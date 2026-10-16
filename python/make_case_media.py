"""Make the animated case-study figures in docs/media from exported Sentinel-1 frames.

    python make_case_media.py <frames_dir>

<frames_dir> holds the Earth Engine exports (not stored in this repository because of their size):
trinity/, wakashio/, triton/ and the LillieWatch time-lapse GIFs. Contains modified Copernicus
Sentinel data.
"""
import json
import sys
from pathlib import Path

from PIL import Image, ImageDraw

from make_media import BG, GOLD, WHITE, GREY, font, save_gif, OUT

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("media_src")
NAME = "Emmanuel Chidi Ibekwe"
W, BOX_H = 760, 520


def card(img, kicker, title, l1, l2="", footer=""):
    img = img.convert("RGB")
    s = min((W - 60) / img.width, BOX_H / img.height)
    img = img.resize((int(img.width * s), int(img.height * s)), Image.LANCZOS)
    im = Image.new("RGB", (W, BOX_H + 260), BG); d = ImageDraw.Draw(im)
    d.text((30, 18), kicker, font=font(22), fill=GOLD)
    d.text((30, 46), title, font=font(32, True), fill=WHITE)
    im.paste(img, ((W - img.width) // 2, 100 + (BOX_H - img.height) // 2))
    y = 100 + BOX_H + 18
    d.text((30, y), l1, font=font(22), fill=WHITE)
    d.text((30, y + 30), l2, font=font(20), fill=GREY)
    d.text((W - 30, y + 74), f"{NAME}  |  {footer}", font=font(17), fill=GOLD, anchor="rs")
    return im.crop((0, 0, W, y + 86))


def rebrand(gif, out, band, text):
    """Replace the old footer line of a LillieWatch time-lapse with the project name and author."""
    src = Image.open(gif); frames, dur = [], []
    for i in range(src.n_frames):
        src.seek(i); f = src.convert("RGB"); d = ImageDraw.Draw(f)
        x0, y0, x1, y1 = band
        d.rectangle(band, fill=f.getpixel((5, y0 + 2)))
        d.text(((x0 + x1) // 2, (y0 + y1) // 2), text, font=font(17), fill=(150, 160, 175), anchor="mm")
        frames.append(f); dur.append(src.info.get("duration", 1400) * 1.4)
    pal = [f.convert("P", palette=Image.ADAPTIVE, colors=128) for f in frames]
    pal[0].save(out, save_all=True, append_images=pal[1:], duration=dur, loop=0, optimize=True)
    print(out, f"{out.stat().st_size / 1e6:.1f} MB")
    return frames


def trinity():
    m = json.load(open(SRC / "trinity" / "manifest.json"))
    frames = []
    notes = {
        "2022-01-03": ("Before the 2 February 2022 explosion: the vessel is a bright point.", "Radar sees through cloud and at night; it measures how rough the sea surface is."),
        "2022-01-15": ("Before the explosion: the vessel is still a bright point.", "Each pass is 12 days apart on the same orbit."),
        "2022-01-27": ("Six days before the explosion.", "This is the baseline used for comparison."),
        "2022-02-08": ("Six days after: dark streaks consistent with possible surface slicks.", "Smooth (dark) water around the vessel position. Radar alone cannot confirm oil."),
        "2022-02-20": ("Eighteen days after: dark features are still present.", "Repeated passes show how the features persist and drift."),
        "2022-03-04": ("One month after: the dark features have largely dispersed.", "The time series shows the full event, not one snapshot."),
    }
    for s in m["s1"]:
        l1, l2 = notes[s["date"]]
        frames.append(card(Image.open(SRC / "trinity" / s["frame"]),
                           "Real Sentinel-1 radar  |  Offshore Delta State, Nigeria",
                           f"Trinity Spirit FPSO, {s['date']}", l1, l2, "LillieWatch (personal project)"))
    save_gif(frames, OUT / "trinity_spirit_2022.gif", hold_ms=2000, fade=0, colors=48)
    return frames


def wakashio():
    m = json.load(open(SRC / "wakashio" / "manifest.json"))
    frames = []
    for s in sorted(m["s1"], key=lambda r: r["frame"][3:11]):
        f = SRC / "wakashio" / s["frame"]
        if not f.exists() or f.stat().st_size < 500_000:
            continue
        frames.append(card(Image.open(f), "Real Sentinel-1 radar  |  Mauritius",
                           f"MV Wakashio, {s['date']}", "Phase: " + s["event"] + ".",
                           "A time series of radar passes over the reef and lagoon.",
                           "LillieWatch (personal project)"))
    save_gif(frames, OUT / "wakashio_2020.gif", hold_ms=1800, fade=0, colors=48)
    return frames


def triton():
    t = SRC / "triton"
    dash = Image.open(t / "triton_dashboard_evidence.png").convert("RGB"); d = ImageDraw.Draw(dash)
    d.rectangle((25, 116, 520, 156), fill=dash.getpixel((20, 140)))
    d.text((32, 136), f"by {NAME}", font=font(30), fill=(160, 172, 190), anchor="lm")
    frames = [
        card(dash, "Triton (personal project)  |  Gulf of Guinea", "Analyst review dashboard",
             "Radar candidates, vessel AIS context and evidence records on one map.",
             "Demonstration scenario: candidates for review, not confirmed incidents.", "Triton"),
        card(Image.open(t / "slick_20240810.png"), "Triton  |  real Sentinel-1 scene", "Possible slick near vessels",
             "A ResNet34 U-Net outlines a dark patch (cyan); vessels are boxed (orange).",
             "Each output is a candidate that an analyst reviews.", "Triton"),
        card(Image.open(t / "triton_real_S1A_IW_GRDH_1SDV_20240729T051432_20240729T051457_054971_06B24A_224F.png"),
             "Triton  |  real Sentinel-1 scene, 29 July 2024", "Vessel detection",
             "Bright points on dark sea are ships: 73 vessel candidates in this scene.",
             "Each is matched against AIS ship-tracking data.", "Triton"),
    ]
    save_gif(frames, OUT / "triton_offshore.gif", hold_ms=3000, fade=2)
    return frames, dash


if __name__ == "__main__":
    rebrand(SRC / "bodo_sar_timelapse.gif", OUT / "lilliewatch_bodo_sar.gif", (200, 732, 760, 762),
            f"LILLIEWATCH SITE WATCH  ·  {NAME.upper()}")
    rebrand(SRC / "bodo_all_sensors_timelapse.gif", OUT / "lilliewatch_bodo_all_sensors.gif", (300, 910, 1084, 942),
            f"LILLIEWATCH SITE WATCH  ·  {NAME.upper()}")
    trinity(); wakashio(); triton()
