"""Make the animated figures in docs/media from the simulator and processor.

    python make_media.py

Output: docs/media/focusing.gif (the processing chain, stage by stage) and
docs/media/slick_detection.gif (detector against simulated truth, by slick contrast).
"""
from pathlib import Path

import numpy as np
from numpy.fft import fft, ifft, fftshift, ifftshift
from PIL import Image, ImageDraw, ImageFont
from scipy.special import i0
import matplotlib
matplotlib.use("Agg")
from matplotlib import cm

from sarsim import params, simulate_raw, slick_scene, detect_dark
from sarsim.core import _hamming, _sinc

OUT = Path(__file__).resolve().parents[1] / "docs" / "media"
OUT.mkdir(parents=True, exist_ok=True)
BG = (12, 18, 30); GOLD = (212, 165, 116); WHITE = (235, 238, 242); GREY = (150, 160, 175)
W, PANEL = 760, 500


def font(size, bold=False):
    for f in ["/usr/share/fonts/truetype/liberation/LiberationSerif-%s.ttf" % ("Bold" if bold else "Regular"),
              "C:/Windows/Fonts/times%s.ttf" % ("bd" if bold else "")]:
        try:
            return ImageFont.truetype(f, size)
        except OSError:
            pass
    return ImageFont.load_default()


def stages(raw, t, p):
    """The same steps as rda_focus, keeping each intermediate product."""
    Na, Nr = raw.shape
    fr = (np.arange(Nr) - Nr / 2) * p.fs / Nr
    tr = (np.arange(Nr) - Nr / 2) / p.fs
    rep = (np.abs(tr) <= p.Tp / 2) * np.exp(1j * np.pi * p.Kr * tr ** 2)
    Rf = fftshift(fft(ifftshift(rep)))
    Hr = np.conj(Rf) / np.maximum(np.abs(Rf), 1e-12) * (np.abs(fr) <= p.Br / 2) * _hamming(fr / p.Br, p.win_alpha)
    src = ifft(ifftshift(fftshift(fft(raw, axis=1), axes=1) * Hr, axes=1), axis=1)
    feta = ((np.arange(Na) - Na / 2) * p.PRF / Na)[:, None]
    SRD = fftshift(fft(src, axis=0), axes=0)
    R = t * p.c / 2
    D = np.sqrt(1 - (p.lam * feta / (2 * p.V)) ** 2)
    dn = R[None, :] * (1 / D - 1) / (p.c / (2 * p.fs))
    pos = np.arange(Nr)[None, :] + dn
    base = np.floor(pos).astype(int); frac = pos - base
    out = np.zeros_like(SRD); rows = np.arange(Na)[:, None]
    for k in np.arange(-7, 9):
        idx = base + k; x = frac - k; u = x / 8
        kw = np.where(np.abs(u) <= 1, i0(2.5 * np.sqrt(np.clip(1 - u ** 2, 0, None))) / i0(2.5), 0.0)
        out += np.where((idx >= 0) & (idx < Nr), _sinc(x) * kw * SRD[rows, np.clip(idx, 0, Nr - 1)], 0)
    Ha = np.exp(1j * 4 * np.pi * R[None, :] * (D - 1) * p.fc / p.c) * (np.abs(feta) <= p.Ba / 2) * _hamming(feta / p.Ba, p.win_alpha)
    img = ifft(ifftshift(out * Ha, axes=0), axis=0)
    return src, SRD, out, img


def to_rgb(a, lo=-45, cmap="inferno"):
    db = 20 * np.log10(np.abs(a) / np.abs(a).max() + 1e-12)
    v = np.clip((db - lo) / -lo, 0, 1)
    rgb = (matplotlib.colormaps[cmap](v)[..., :3] * 255).astype(np.uint8)
    return Image.fromarray(rgb)


def card(panel, step, title, line1, line2=""):
    im = Image.new("RGB", (W, PANEL + 160), BG); d = ImageDraw.Draw(im)
    d.text((30, 18), step, font=font(22), fill=GOLD)
    d.text((30, 46), title, font=font(32, True), fill=WHITE)
    im.paste(panel.resize((W - 60, PANEL - 40)), (30, 100))
    d.text((30, PANEL + 72), line1, font=font(22), fill=WHITE)
    d.text((30, PANEL + 102), line2, font=font(20), fill=GREY)
    d.text((W - 30, PANEL + 136), "Emmanuel Chidi Ibekwe  |  simulated C-band SAR", font=font(17), fill=GOLD, anchor="rs")
    return im


def save_gif(frames, path, hold_ms=2800, fade=2, colors=96):
    seq, dur = [], []
    for i, f in enumerate(frames):
        seq.append(f); dur.append(hold_ms)
        g = frames[(i + 1) % len(frames)]
        for k in range(1, fade + 1):
            seq.append(Image.blend(f, g, k / (fade + 1))); dur.append(70)
    pal = [s.convert("P", palette=Image.ADAPTIVE, colors=colors) for s in seq]
    pal[0].save(path, save_all=True, append_images=pal[1:], duration=dur, loop=0, optimize=True)
    print(path, f"{path.stat().st_size / 1e6:.1f} MB")


def focusing_gif():
    p = params(); p.snr_in_db = 10.0
    raw, t, eta = simulate_raw(p)
    src, SRD, rc, img = stages(raw, t, p)
    A, Rg = slice(520, 1540), slice(150, 560)          # about 4.6 km along track, 1.7 km in range
    Z = (slice(930, 1120), slice(225, 415))            # zoom on the four targets
    irf = Image.open(OUT.parent / "figures" / "fig2_irf_cuts.png").convert("RGB")
    frames = [
        card(to_rgb(raw[A, Rg], -30, "gray"), "What the radar records",
             "Raw echoes from four point targets",
             "Each target is smeared over kilometres in both directions.",
             "Long pulses (10 microseconds) and a 4 km synthetic aperture overlap."),
        card(to_rgb(src[A, Rg], -40, "gray"), "Step 1 of 3: range compression",
             "Matched filter on every pulse",
             "Each target is now sharp across track (range) but still long along track.",
             "Replica chirp filter with Hamming weighting."),
        card(to_rgb(rc[A, Rg], -40, "gray"), "Step 2 of 3: azimuth FFT and migration correction",
             "Range-Doppler domain, migration corrected",
             "A 16-tap Kaiser-windowed sinc interpolator realigns each target's energy.",
             "My verification found that the first 8-tap interpolator failed the sidelobe requirement."),
        card(to_rgb(img[A, Rg], -45), "Step 3 of 3: azimuth compression",
             "Focused image",
             "Four sharp points at the correct positions (error 0.14 m or less).",
             "Same 4.6 km by 1.7 km area as the raw data."),
        card(to_rgb(img[Z], -45), "Close-up",
             "Focused targets, 850 m by 800 m",
             "Each point is about 5 m by 6 m: the resolution the design predicts.",
             "Brightness in dB; faint crosses are sidelobes, kept below -20 dB."),
        card(irf, "Quality check",
             "Measured impulse response against requirements",
             "Width, sidelobes and position are measured for every target. All pass.",
             "Numbers, method and the two defects found: docs/VERIFICATION.md."),
    ]
    save_gif(frames, OUT / "focusing.gif")


def slick_gif():
    frames = []
    for c in [2, 3, 4, 6, 8, 10]:
        I, truth, look = slick_scene(contrast_db=float(c), seed=3)
        m = detect_dark(I)
        pod = (m & truth).sum() / truth.sum(); far = (m & ~truth & ~look).sum() / (~truth & ~look).sum()
        g = 10 * np.log10(I); g = np.clip((g - (-25)) / 17, 0, 1)
        base = np.stack([g] * 3, -1)
        ov = base.copy(); ov[m] = [0.95, 0.65, 0.2]
        edge = truth ^ np.roll(truth, 1, 0) | truth ^ np.roll(truth, 1, 1)
        ov[edge] = [0.3, 0.8, 1.0]
        pair = np.concatenate([base, np.ones((base.shape[0], 8, 3)) * 0.05, ov], 1)
        panel = Image.fromarray((pair * 255).astype(np.uint8))
        frames.append(card(panel, "Oil-slick detection on simulated sea",
                           f"Slick {c} dB darker than the surrounding sea",
                           f"Found {pod * 100:.0f} percent of slick pixels; false alarms on open sea {far:.0e}.",
                           "Left: simulated radar image. Right: detection (orange) and true slick outline (blue)."))
    save_gif(frames, OUT / "slick_detection.gif", hold_ms=1800)


if __name__ == "__main__":
    focusing_gif()
    slick_gif()
