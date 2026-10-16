"""Assemble docs/media/*.gif into one captioned MP4 reel (docs/media/portfolio_reel.mp4).

    python make_reel.py      (after make_media.py and make_case_media.py)
"""
import imageio_ffmpeg
import numpy as np
from PIL import Image, ImageDraw

from make_media import BG, GOLD, WHITE, GREY, font, OUT

VW, VH, FPS = 1280, 720, 24
NAME = "Emmanuel Chidi Ibekwe"


def fit(im):
    im = im.convert("RGB"); s = min(VW / im.width, VH / im.height)
    im = im.resize((int(im.width * s), int(im.height * s)), Image.LANCZOS)
    c = Image.new("RGB", (VW, VH), BG); c.paste(im, ((VW - im.width) // 2, (VH - im.height) // 2)); return c


def title(kicker, head, sub=""):
    c = Image.new("RGB", (VW, VH), BG); d = ImageDraw.Draw(c)
    d.text((VW // 2, 280), kicker, font=font(30), fill=GOLD, anchor="mm")
    d.text((VW // 2, 345), head, font=font(52, True), fill=WHITE, anchor="mm")
    d.text((VW // 2, 420), sub, font=font(28), fill=GREY, anchor="mm")
    d.text((VW // 2, 660), NAME, font=font(24), fill=GOLD, anchor="mm")
    return c


def gif_frames(name, scale=1.0):
    im = Image.open(OUT / name); out = []
    for i in range(im.n_frames):
        im.seek(i); d = im.info.get("duration", 1000)
        if d >= 1000:
            out.append((fit(im), d / 1000 * scale))
    return out


seq = [(title("SAR processing portfolio", NAME, "From raw radar echoes to monitoring decisions"), 3.5)]
for kick, head, sub, gif, sc in [
    ("Part 1  |  simulation", "Processing raw SAR data", "A range-Doppler processor I wrote and verified", "focusing.gif", 1.3),
    ("Part 2  |  simulation", "Testing a detector against known truth", "Simulated sea with a slick of known contrast", "slick_detection.gif", 1.2),
    ("Part 3  |  real Sentinel-1 data", "Offshore explosion, Nigeria, 2022", "A radar time series before and after the event", "trinity_spirit_2022.gif", 1.3),
    ("Part 4  |  real Sentinel-1 data", "MV Wakashio, Mauritius, 2020", "Radar passes through a grounding and spill", "wakashio_2020.gif", 1.2),
    ("Part 5  |  LillieWatch (personal project)", "Weekly site monitoring", "Radar, thermal, optical and SO2 for one site, week by week", "lilliewatch_bodo_sar.gif", 1.0),
    (None, None, None, "lilliewatch_bodo_all_sensors.gif", 1.0),
    ("Part 6  |  Triton (personal project)", "Offshore watch, Gulf of Guinea", "Slick and vessel candidates for analyst review", "triton_offshore.gif", 1.2),
]:
    if kick:
        seq.append((title(kick, head, sub), 2.5))
    seq += gif_frames(gif, sc)
end = Image.new("RGB", (VW, VH), BG); d = ImageDraw.Draw(end)
d.text((VW // 2, 250), NAME, font=font(54, True), fill=WHITE, anchor="mm")
d.text((VW // 2, 320), "SAR processing  |  satellite remote sensing  |  scientific software", font=font(28), fill=GOLD, anchor="mm")
d.text((VW // 2, 400), "github.com/Ibekwemmanuel7/sar-processing-portfolio", font=font(28), fill=WHITE, anchor="mm")
d.text((VW // 2, 445), "ibekwemmanuel@gmail.com  |  linkedin.com/in/ibekwemmanuel7", font=font(26), fill=GREY, anchor="mm")
d.text((VW // 2, 640), "Personal projects. Simulated data and contains modified Copernicus Sentinel data.", font=font(20), fill=GREY, anchor="mm")
seq.append((end, 5.0))

w = imageio_ffmpeg.write_frames(str(OUT / "portfolio_reel.mp4"), (VW, VH), fps=FPS, codec="libx264",
                                pix_fmt_out="yuv420p", output_params=["-crf", "23", "-movflags", "+faststart"])
w.send(None); fade = int(0.35 * FPS)
for i, (im, sec) in enumerate(seq):
    a = np.asarray(im)
    for _ in range(int(sec * FPS)):
        w.send(a.tobytes())
    if i + 1 < len(seq):
        b = np.asarray(seq[i + 1][0]).astype(float)
        for k in range(1, fade + 1):
            w.send((a * (1 - k / (fade + 1)) + b * k / (fade + 1)).astype(np.uint8).tobytes())
w.close()
print(OUT / "portfolio_reel.mp4", f"{(OUT / 'portfolio_reel.mp4').stat().st_size / 1e6:.1f} MB",
      f"{sum(s for _, s in seq):.0f} s")
