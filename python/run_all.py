"""Run the simulator, processor and verification; write figures to ../docs/figures and a summary."""
import copy
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sarsim import params, simulate_raw, rda_focus, measure_irf, irf_model, nesz_model, slick_scene, detect_dark

FIG = Path(__file__).resolve().parents[1] / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 9})

p = params()
raw, t, eta = simulate_raw(p)
img, info = rda_focus(raw, t, p)
mdl = irf_model(p)
rows = []
for k, (xa, dr, _) in enumerate(p.targets, 1):
    m = measure_irf(img, p, t, eta, xa, dr)
    rows.append([k, xa, dr, m.irw_rg_m, m.irw_az_m, m.pslr_rg_db, m.pslr_az_db, m.islr_rg_db, m.islr_az_db,
                 m.loc_err_rg_m, m.loc_err_az_m])
    print(f"T{k}: IRW rg {m.irw_rg_m:.2f} az {m.irw_az_m:.2f} m | PSLR rg {m.pslr_rg_db:.1f} az {m.pslr_az_db:.1f} | "
          f"ISLR rg {m.islr_rg_db:.1f} az {m.islr_az_db:.1f} | loc {m.loc_err_rg_m:.2f}/{m.loc_err_az_m:.2f} m")
print(f"model: IRW rg {mdl.irw_rg_m:.2f} az {mdl.irw_az_m:.2f}, PSLR rg {mdl.pslr_rg_db:.1f} az {mdl.pslr_az_db:.1f}")
with open(FIG.parent / "python_irf_metrics.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["target", "az_m", "rg_m", "irw_rg_m", "irw_az_m", "pslr_rg_db", "pslr_az_db",
                                   "islr_rg_db", "islr_az_db", "loc_err_rg_m", "loc_err_az_m"])
    w.writerows([[f"{v:.3f}" if isinstance(v, float) else v for v in r] for r in rows])

# Figure 1: raw magnitude, range-compressed, focused image
R = t * p.c / 2
ext = [(R[0] - p.R0) / 1e3, (R[-1] - p.R0) / 1e3, p.V * eta[-1] / 1e3, p.V * eta[0] / 1e3]
fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0), dpi=200)
ax[0].imshow(np.abs(raw), aspect="auto", cmap="gray", extent=ext)
ax[0].set_title("Simulated raw data (magnitude)"); ax[0].set_xlabel("slant range offset (km)"); ax[0].set_ylabel("azimuth (km)")
fdb = 20 * np.log10(np.abs(img) / np.abs(img).max() + 1e-12)
ax[1].imshow(fdb, aspect="auto", cmap="inferno", vmin=-50, vmax=0, extent=ext)
ax[1].set_xlim(-0.5, 0.4); ax[1].set_ylim(0.35, -0.35)
ax[1].set_title("Focused (range-Doppler), dB, zoomed"); ax[1].set_xlabel("slant range offset (km)")
fig.tight_layout(); fig.savefig(FIG / "fig1_raw_vs_focused.png"); plt.close(fig)

# Figure 2: IRF cuts with the model for target 1
n, h = 32, 16
ir = int(np.argmin(np.abs(R - p.R0))); ia = int(np.argmin(np.abs(p.V * eta)))
patch = img[ia - n // 2:ia + n // 2, ir - n // 2:ir + n // 2]
P = np.fft.fftshift(np.fft.fft2(patch)); Pz = np.zeros((n * h, n * h), complex); o = (n * h - n) // 2
Pz[o:o + n, o:o + n] = P; up = np.abs(np.fft.ifft2(np.fft.ifftshift(Pz))) ** 2
pa, pr = np.unravel_index(up.argmax(), up.shape); up /= up.max()
xr = (np.arange(n * h) - pr) * (p.c / (2 * p.fs)) / h; xa = (np.arange(n * h) - pa) * (p.V / p.PRF) / h
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6), dpi=200)
for a_, x, cut, nm, irw, ps in [(ax[0], xr, up[pa, :], "Range", rows[0][3], rows[0][5]),
                                (ax[1], xa, up[:, pr], "Azimuth", rows[0][4], rows[0][6])]:
    a_.plot(x, 10 * np.log10(cut + 1e-12), color="#1f4e79")
    a_.axhline(-3, color="#c55a11", lw=0.8, ls="--"); a_.axhline(p.spec.pslr_db, color="grey", lw=0.8, ls=":")
    a_.set_xlim(-40, 40); a_.set_ylim(-50, 2); a_.set_xlabel(f"{nm.lower()} (m)"); a_.grid(alpha=0.3)
    a_.set_title(f"{nm} IRF: IRW {irw:.2f} m, PSLR {ps:.1f} dB")
ax[0].set_ylabel("power (dB)")
fig.tight_layout(); fig.savefig(FIG / "fig2_irf_cuts.png"); plt.close(fig)

# Figure 3: NESZ across the swath
inc, Rk, nesz = nesz_model(p)
fig, ax = plt.subplots(figsize=(3.4, 2.4), dpi=200)
ax.plot(inc, nesz, "o-", color="#1f4e79"); ax.axhline(p.spec.nesz_db, color="#c55a11", ls="--", lw=0.8)
ax.text(21, p.spec.nesz_db - 0.6, "requirement (-22 dB)", color="#c55a11", fontsize=7)
ax.set_xlabel("incidence angle (deg)"); ax.set_ylabel("NESZ (dB)"); ax.grid(alpha=0.3); ax.set_title("Radar-equation NESZ")
fig.tight_layout(); fig.savefig(FIG / "fig3_nesz.png"); plt.close(fig)

# Figure 4: slick simulator and detector, plus POD vs contrast
I, truth, la = slick_scene(contrast_db=6, seed=3); m = detect_dark(I)
cs = [2, 3, 4, 6, 8, 10]; pod = []; far = []
for c in cs:
    Ic, tr, lac = slick_scene(contrast_db=c, seed=3); mc = detect_dark(Ic); sea = ~tr & ~lac
    pod.append((mc & tr).sum() / tr.sum()); far.append((mc & sea).sum() / sea.sum())
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.5), dpi=200)
ax[0].imshow(10 * np.log10(I), cmap="gray", vmin=-24, vmax=-8); ax[0].set_title("Simulated sea, 6 dB slick")
ax[1].imshow(10 * np.log10(I), cmap="gray", vmin=-24, vmax=-8)
ax[1].contour(truth, levels=[0.5], colors="#ffd400", linewidths=0.7)
ax[1].imshow(np.ma.masked_where(~m, m), cmap="autumn", alpha=0.45); ax[1].set_title("Truth (yellow), detection (red)")
for a_ in ax[:2]: a_.set_xticks([]); a_.set_yticks([])
ax[2].plot(cs, pod, "o-", color="#1f4e79", label="POD"); ax[2].set_xlabel("slick contrast (dB)")
ax[2].set_ylim(0, 1.02); ax[2].grid(alpha=0.3); ax[2].set_title("Detection vs contrast")
ax2 = ax[2].twinx(); ax2.semilogy(cs, np.maximum(far, 1e-6), "s--", color="#c55a11", label="FAR"); ax2.set_ylabel("false-alarm rate (open sea)")
fig.tight_layout(); fig.savefig(FIG / "fig4_slick_simulation.png"); plt.close(fig)
print("POD", np.round(pod, 3).tolist(), "FAR", np.round(far, 5).tolist())
print("figures written to", FIG)
