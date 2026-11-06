#!/usr/bin/env python3
"""Synthetic example images for colony-count and invasion-count.

No real lab data: every image is drawn here, with a known ground truth.
  python3 make_examples.py [OUT_DIR]      (default: next to this script)

Writes:
  colony/plates.jpg + colony/design.csv + colony/truth.csv
  invasion/rep{1,2,3}/LineA/{CTRL,TREATED}/field_NN_20X.jpg + invasion/truth.csv
"""
import csv, sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

rng = np.random.default_rng(7)
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
VIOLET = (105, 55, 150)


def colony_plate(n_ctrl=90, n_treated=35):
    w, h, r = 1400, 700, 300
    im = Image.new("RGB", (w, h), (205, 205, 200))
    dr = ImageDraw.Draw(im)
    centres = [(350, 350), (1050, 350)]
    truth = []
    for (cx, cy), n in zip(centres, (n_ctrl, n_treated)):
        dr.ellipse([cx - r, cy - r, cx + r, cy + r], fill=(235, 233, 228), outline=(90, 90, 90), width=6)
        placed = []
        while len(placed) < n:
            ang, rad = rng.uniform(0, 2 * np.pi), r * 0.82 * np.sqrt(rng.uniform())
            x, y, cr = cx + rad * np.cos(ang), cy + rad * np.sin(ang), rng.uniform(4, 10)
            if all((x - px) ** 2 + (y - py) ** 2 > (cr + pr + 4) ** 2 for px, py, pr in placed):
                placed.append((x, y, cr))
                shade = tuple(int(c + rng.integers(-15, 15)) for c in VIOLET)
                dr.ellipse([x - cr, y - cr, x + cr, y + cr], fill=shade)
        truth.append(n)
    im = im.filter(ImageFilter.GaussianBlur(1.2))
    noise = rng.normal(0, 3, (h, w, 3))
    return Image.fromarray(np.clip(np.asarray(im) + noise, 0, 255).astype(np.uint8)), truth


def invasion_field(n):
    """n cells of realistic size (radius ~20-30 px at 1536 px width), some touching, on a pored membrane."""
    w, h = 1536, 1024
    im = Image.new("RGB", (w, h), (228, 216, 192))
    dr = ImageDraw.Draw(im)
    for _ in range(2500):  # membrane pores, drawn under the cells
        x, y = rng.uniform(0, w), rng.uniform(0, h)
        dr.ellipse([x - 4, y - 4, x + 4, y + 4], fill=(205, 195, 172))
    placed = []
    while len(placed) < n:
        x, y = rng.uniform(30, w - 30), rng.uniform(30, h - 30)
        if any((x - px) ** 2 + (y - py) ** 2 < 46 ** 2 for px, py in placed):  # touching is fine, stacked is not
            continue
        placed.append((x, y))
        a, b = rng.uniform(20, 32), rng.uniform(15, 24)
        dr.ellipse([x - a, y - b, x + a, y + b], fill=tuple(int(c + rng.integers(-12, 12)) for c in VIOLET))
    im = im.filter(ImageFilter.GaussianBlur(1.2))
    yy, xx = np.mgrid[0:h, 0:w]
    vignette = 1 - 0.18 * (((xx / w - .5) ** 2 + (yy / h - .5) ** 2) / 0.5)
    return Image.fromarray(np.clip(np.asarray(im) * vignette[..., None], 0, 255).astype(np.uint8))


def main():
    c = OUT / "colony"; c.mkdir(parents=True, exist_ok=True)
    plate, truth = colony_plate()
    plate.save(c / "plates.jpg", quality=90)
    with open(c / "design.csv", "w", newline="") as f:
        csv.writer(f).writerows([["image", "dish", "line", "condition", "dose"],
                                 ["plates.jpg", 1, "LineA", "CTRL", ""],
                                 ["plates.jpg", 2, "LineA", "TREATED", ""]])
    with open(c / "truth.csv", "w", newline="") as f:
        csv.writer(f).writerows([["dish", "condition", "colonies"], [1, "CTRL", truth[0]], [2, "TREATED", truth[1]]])

    rows = [["replicate", "condition", "field", "cells"]]
    for rep, scale in (("rep1", 1.0), ("rep2", 0.85), ("rep3", 1.15)):
        for cond, mean in (("CTRL", 110), ("TREATED", 45)):
            d = OUT / "invasion" / rep / "LineA" / cond; d.mkdir(parents=True, exist_ok=True)
            for i in range(1, 5):
                n = int(round(mean * scale * rng.uniform(0.85, 1.15)))
                invasion_field(n).save(d / f"field_{i:02d}_20X.jpg", quality=88)
                rows.append([rep, cond, i, n])
    with open(OUT / "invasion" / "truth.csv", "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print(f"examples written to {OUT}")


if __name__ == "__main__":
    main()
