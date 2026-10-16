"""Core algorithms. Arrays are (azimuth, range) = (Na, Nr), matching the MATLAB code."""
from __future__ import annotations

import copy
from types import SimpleNamespace

import numpy as np
from numpy.fft import fft, ifft, fft2, ifft2, fftshift, ifftshift
from scipy.ndimage import uniform_filter, label
from scipy.special import i0

C = 299792458.0


def params() -> SimpleNamespace:
    """C-band spaceborne stripmap SAR (Sentinel-1-like, simplified). See matlab/sar_params.m."""
    p = SimpleNamespace()
    p.c = C; p.fc = 5.405e9; p.lam = p.c / p.fc
    p.V = 7000.0; p.R0 = 850e3; p.La = 10.0
    p.Tp = 10e-6; p.Br = 30e6; p.Kr = p.Br / p.Tp; p.fs = 36e6
    p.Ba = 0.886 * 2 * p.V / p.La; p.PRF = 1.25 * p.Ba
    p.Ta = 0.886 * p.lam * p.R0 / (p.La * p.V)
    p.Nr = 1024; p.Na = 2048
    p.win_alpha = 0.75; p.range_filter = "replica"; p.rcmc = True
    p.snr_in_db = 0.0
    p.targets = np.array([[0, 0, 1.0], [-150, -200, 0.8], [200, 250, 0.6], [100, -350, 0.5]])
    p.spec = SimpleNamespace(irw_rg_m=5.5, irw_az_m=6.5, pslr_db=-20.0, islr_db=-15.0,
                             loc_err_m=1.0, model_tol=0.03, nesz_db=-22.0)
    return p


def _sinc(x):
    return np.sinc(x)          # numpy sinc is sin(pi x)/(pi x)


def simulate_raw(p, seed: int = 7):
    """Baseband raw echoes for point targets: hyperbolic range history, LFM chirp,
    two-way sinc^2 azimuth antenna pattern, complex white noise at p.snr_in_db per sample."""
    Rmin = p.R0 - 0.5 * (p.Nr / p.fs) * p.c / 2 + p.Tp * p.c / 4
    t = 2 * Rmin / p.c + np.arange(p.Nr) / p.fs
    eta = (np.arange(p.Na) - p.Na / 2) / p.PRF
    T, ETA = np.meshgrid(t, eta)
    raw = np.zeros((p.Na, p.Nr), complex)
    th_bw = 0.886 * p.lam / p.La
    for xa, dr, A in p.targets:
        R0k = p.R0 + dr
        R = np.sqrt(R0k ** 2 + (p.V * ETA - xa) ** 2)
        tau = T - 2 * R / p.c
        theta = np.arctan((p.V * ETA - xa) / R0k)
        wa = _sinc(0.886 * theta / th_bw) ** 2
        wr = np.abs(tau) <= p.Tp / 2
        raw += A * wr * wa * np.exp(-1j * 4 * np.pi * p.fc * R / p.c) * np.exp(1j * np.pi * p.Kr * tau ** 2)
    rng = np.random.default_rng(seed)
    s = np.sqrt(1 / 10 ** (p.snr_in_db / 10) / 2)
    raw += s * (rng.standard_normal(raw.shape) + 1j * rng.standard_normal(raw.shape))
    return raw, t, eta


def _hamming(u, alpha):
    w = alpha + (1 - alpha) * np.cos(2 * np.pi * u)
    return np.where(np.abs(u) <= 0.5, w, 0.0)


def rda_focus(raw, t, p):
    """Range-Doppler algorithm: range compression, azimuth FFT, RCMC (16-tap Kaiser-windowed sinc),
    azimuth compression with the exact zero-squint filter (Doppler-dependent part only)."""
    Na, Nr = raw.shape
    fr = (np.arange(Nr) - Nr / 2) * p.fs / Nr
    if p.range_filter == "posp":
        Hr = np.exp(1j * np.pi * fr ** 2 / p.Kr) * (np.abs(fr) <= p.Br / 2)
    else:
        tr = (np.arange(Nr) - Nr / 2) / p.fs
        rep = (np.abs(tr) <= p.Tp / 2) * np.exp(1j * np.pi * p.Kr * tr ** 2)
        Rf = fftshift(fft(ifftshift(rep)))
        Hr = np.conj(Rf) / np.maximum(np.abs(Rf), 1e-12) * (np.abs(fr) <= p.Br / 2)
    Hr = Hr * _hamming(fr / p.Br, p.win_alpha)
    src = ifft(ifftshift(fftshift(fft(raw, axis=1), axes=1) * Hr, axes=1), axis=1)
    feta = ((np.arange(Na) - Na / 2) * p.PRF / Na)[:, None]
    SRD = fftshift(fft(src, axis=0), axes=0)
    R = t * p.c / 2
    D = np.sqrt(1 - (p.lam * feta / (2 * p.V)) ** 2)
    dn = R[None, :] * (1 / D - 1) / (p.c / (2 * p.fs))
    if not p.rcmc:
        dn = np.zeros_like(dn)
    pos = np.arange(Nr)[None, :] + dn
    base = np.floor(pos).astype(int); frac = pos - base
    out = np.zeros_like(SRD)
    taps = np.arange(-7, 9); half = len(taps) / 2; beta = 2.5
    rows = np.arange(Na)[:, None]
    for k in taps:
        idx = base + k
        x = frac - k
        u = x / half
        kw = np.where(np.abs(u) <= 1, i0(beta * np.sqrt(np.clip(1 - u ** 2, 0, None))) / i0(beta), 0.0)
        w = _sinc(x) * kw
        valid = (idx >= 0) & (idx < Nr)
        out += np.where(valid, w * SRD[rows, np.clip(idx, 0, Nr - 1)], 0)
    Ha = np.exp(1j * 4 * np.pi * R[None, :] * (D - 1) * p.fc / p.c)
    Ha = Ha * (np.abs(feta) <= p.Ba / 2) * _hamming(feta / p.Ba, p.win_alpha)
    img = ifft(ifftshift(out * Ha, axes=0), axis=0)
    return img, {"max_rcm_cells": float(np.abs(dn).max())}


def cut_metrics(cut, k0, dx):
    """-3 dB width (m), PSLR (dB) and ISLR (dB) of a normalised 1-D power cut peaking at k0."""
    db = 10 * np.log10(cut + 1e-300)
    l = k0
    while l > 0 and db[l] > -3:
        l -= 1
    r = k0
    while r < len(db) - 1 and db[r] > -3:
        r += 1
    xl = l + (-3 - db[l]) / (db[l + 1] - db[l])
    xr = r - 1 + (-3 - db[r - 1]) / (db[r] - db[r - 1])
    a = k0
    while a > 1 and cut[a - 1] < cut[a]:
        a -= 1
    b = k0
    while b < len(cut) - 2 and cut[b + 1] < cut[b]:
        b += 1
    side = np.concatenate([cut[:a], cut[b + 1:]])
    return (xr - xl) * dx, 10 * np.log10(side.max()), 10 * np.log10(side.sum() / cut[a:b + 1].sum())


def measure_irf(img, p, t, eta, xa, dr, n=32, h=16):
    """IRF metrics around a target: 16x FFT upsampling of a 32 x 32 patch, range and azimuth cuts."""
    dx_r = p.c / (2 * p.fs); dx_a = p.V / p.PRF
    R = t * p.c / 2
    ir = int(np.argmin(np.abs(R - (p.R0 + dr)))); ia = int(np.argmin(np.abs(p.V * eta - xa)))
    sub = np.abs(img[ia - 20:ia + 21, ir - 20:ir + 21])
    ka, kr = np.unravel_index(sub.argmax(), sub.shape); ia += ka - 20; ir += kr - 20
    patch = img[ia - n // 2:ia + n // 2, ir - n // 2:ir + n // 2]
    P = fftshift(fft2(patch)); Pz = np.zeros((n * h, n * h), complex); o = (n * h - n) // 2
    Pz[o:o + n, o:o + n] = P
    up = np.abs(ifft2(ifftshift(Pz))) ** 2
    pa, pr = np.unravel_index(up.argmax(), up.shape)
    m = SimpleNamespace()
    m.irw_rg_m, m.pslr_rg_db, m.islr_rg_db = cut_metrics(up[pa, :] / up[pa, pr], pr, dx_r / h)
    m.irw_az_m, m.pslr_az_db, m.islr_az_db = cut_metrics(up[:, pr] / up[pa, pr], pa, dx_a / h)
    m.loc_err_rg_m = R[ir - n // 2] + pr * dx_r / h - (p.R0 + dr)
    m.loc_err_az_m = p.V * eta[ia - n // 2] + pa * dx_a / h - xa
    return m


def irf_model(p, N=2 ** 16):
    """Analytic IRF from the processed spectra (range: rect x Hamming; azimuth: antenna x rect x Hamming)."""
    def from_spec(W, df, scale):
        hh = np.abs(fftshift(ifft(ifftshift(W)))) ** 2; hh /= hh.max()
        return cut_metrics(hh, int(hh.argmax()), scale / (len(W) * df))
    m = SimpleNamespace()
    f = (np.arange(N) - N / 2) / N * 8 * p.Br
    W = (np.abs(f) <= p.Br / 2) * (p.win_alpha + (1 - p.win_alpha) * np.cos(2 * np.pi * f / p.Br))
    m.irw_rg_m, m.pslr_rg_db, m.islr_rg_db = from_spec(W, f[1] - f[0], p.c / 2)
    f = (np.arange(N) - N / 2) / N * 8 * p.Ba
    th = p.lam * f / (2 * p.V); thbw = 0.886 * p.lam / p.La
    W = _sinc(0.886 * th / thbw) ** 2 * (np.abs(f) <= p.Ba / 2) * \
        (p.win_alpha + (1 - p.win_alpha) * np.cos(2 * np.pi * f / p.Ba))
    m.irw_az_m, m.pslr_az_db, m.islr_az_db = from_spec(W, f[1] - f[0], p.V)
    return m


def nesz_model(p, Pt=4000.0, Wa=1.0, eff=0.6, F_db=3.5, L_db=4.0, H=693e3):
    """NESZ across the swath (Curlander and McDonough): returns incidence (deg), slant range (km), NESZ (dB)."""
    k = 1.380649e-23; T0 = 290.0; Re = 6371e3
    G = 4 * np.pi * eff * p.La * Wa / p.lam ** 2
    inc = np.arange(20, 47, 2.0)
    look = np.arcsin(Re / (Re + H) * np.sin(np.radians(inc)))
    R = (Re + H) * np.cos(look) - np.sqrt(Re ** 2 - ((Re + H) * np.sin(look)) ** 2)
    Pavg = Pt * p.Tp * p.PRF
    nesz = 2 * (4 * np.pi) ** 3 * R ** 3 * p.V * k * T0 * 10 ** (F_db / 10) * p.Br * 10 ** (L_db / 10) * \
        np.sin(np.radians(inc)) / (Pavg * G ** 2 * p.lam ** 3 * p.c)
    return inc, R / 1e3, 10 * np.log10(nesz)


def slick_scene(n=512, sigma0_db=-15.0, looks=4, nu=8, contrast_db=6.0, lookalike_db=4.5, seed=1):
    """Image-domain sea clutter (K-distributed: Gamma texture x Gamma speckle) with a slick and a
    broad low-wind look-alike. Returns intensity, slick truth mask, look-alike mask."""
    rng = np.random.default_rng(seed)
    Y, X = np.mgrid[1:n + 1, 1:n + 1].astype(float)
    th = np.radians(35); cx, cy, a, b = 0.45 * n, 0.5 * n, 150, 30
    u = ((X - cx) * np.cos(th) + (Y - cy) * np.sin(th)) / a
    v = (-(X - cx) * np.sin(th) + (Y - cy) * np.cos(th)) / b
    edge = uniform_filter(rng.standard_normal((n, n)), 15) * 6
    truth = (u ** 2 + v ** 2 + 0.35 * edge) < 1
    la = np.exp(-(((X - 0.82 * n) / 60) ** 2 + ((Y - 0.2 * n) / 45) ** 2))
    damp = contrast_db * truth + lookalike_db * la
    sig = 10 ** (sigma0_db / 10) * 10 ** (-damp / 10)
    tex = np.ones((n, n)) if np.isinf(nu) else rng.gamma(nu, 1 / nu, (n, n))
    spk = rng.gamma(looks, 1 / looks, (n, n))
    return sig * tex * spk, truth, la > 0.5


def detect_dark(I, k_db=3.0, win_bg=101, win_sm=7, min_px=200):
    """Adaptive-threshold dark-spot detector: local mean k_db below the background mean, small regions removed."""
    bg = uniform_filter(I, win_bg, mode="nearest"); sm = uniform_filter(I, win_sm, mode="nearest")
    mask = sm < bg * 10 ** (-k_db / 10)
    lab, nlab = label(mask)
    if nlab:
        sizes = np.bincount(lab.ravel()); keep = sizes >= min_px; keep[0] = False
        mask = keep[lab]
    return mask
