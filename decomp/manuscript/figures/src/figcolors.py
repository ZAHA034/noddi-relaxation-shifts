"""Shared colours and colour-vision checks for the manuscript figures (v12; review direction S42 §1.2).

Rules (S42 §1.2): sequential = one hue, light -> dark (blue steps); diverging = red - neutral grey - blue with OKLab
lightness monotone on each arm and similar lightness at the two poles; two categorical colours (agree / oppose);
the theory figures keep the Okabe-Ito colours fixed to meaning (TE blue, TR orange with dark-orange text, tissue
green, relaxation vermillion). Colour is never the only cue: every figure pairs it with outline, glyph, hatch or
position. `cvd_sheet` simulates protanomaly, deuteranomaly and tritanomaly at severity 1.0 (Machado, Oliveira and
Fernandes 2009, IEEE TVCG 15:1291-1298; matrices applied in linear sRGB) for the report.

Import from a figure script in the same folder:
    import importlib.util, pathlib
    spec = importlib.util.spec_from_file_location("figcolors", pathlib.Path(__file__).with_name("figcolors.py"))
    FC = importlib.util.module_from_spec(spec); spec.loader.exec_module(FC)
"""
from __future__ import annotations

import pathlib

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb

# ---------------------------------------------------------------- neutrals
INK, INK2, MUTED, GRID, NEUTRAL = "#0b0b0b", "#52514e", "#8a8984", "#d9d8d4", "#f0efec"

# ---------------------------------------------------------------- Okabe-Ito (theory figures; meaning fixed)
OI = dict(blue="#0072B2", orange="#E69F00", orange_text="#A06E00", green="#009E73", verm="#D55E00",
          sky="#56B4E9", grey="#6E6E6E", ink="#1A1A1A")

# ---------------------------------------------------------------- sequential (magnitude): one hue, light -> dark
SEQ_STEPS = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
# ---------------------------------------------------------------- diverging: red - neutral - blue
RED_ARM = ["#7e1a1c", "#b42f30", "#e34948", "#ef9a92", NEUTRAL]     # dark -> neutral
BLUE_ARM = [NEUTRAL, "#9ec5f4", "#3987e5", "#256abf", "#184f95"]      # neutral -> dark
# ---------------------------------------------------------------- categorical pair
AGREE, OPPOSE = "#2a78d6", "#e34948"


def _lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _enc(c):
    c = np.clip(np.asarray(c, float), 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * np.power(c, 1 / 2.4) - 0.055)


def oklab(hexcol: str) -> tuple[float, float, float]:
    r, g, b = _lin(to_rgb(hexcol))
    l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
    m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
    s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
    l, m, s = np.cbrt([l, m, s])
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def check_palettes() -> dict:
    """The S42 §1.2 assertions; returns the lightness values for the report."""
    out = {}
    Ls = [oklab(c)[0] for c in SEQ_STEPS]
    assert all(np.diff(Ls) < 0), ("sequential steps must darken monotonically", Ls)
    out["seq_L"] = Ls
    for name, arm, rising in (("red_arm", RED_ARM, True), ("blue_arm", BLUE_ARM, False)):
        L = [oklab(c)[0] for c in arm]
        assert all(np.diff(L) > 0) if rising else all(np.diff(L) < 0), (name, L)
        out[name + "_L"] = L
    assert abs(oklab(RED_ARM[0])[0] - oklab(BLUE_ARM[-1])[0]) < 0.06, "diverging poles of similar lightness"
    return out


check_palettes()
SEQ = LinearSegmentedColormap.from_list("seq_blue", SEQ_STEPS)
DIV = LinearSegmentedColormap.from_list("div_red_grey_blue", RED_ARM + BLUE_ARM[1:])

# ---------------------------------------------------------------- colour-vision simulation (Machado et al. 2009, severity 1.0)
MACHADO_1 = {
    "protanomaly (severity 1.0)": np.array([[0.152286, 1.052583, -0.204868],
                                            [0.114503, 0.786281, 0.099216],
                                            [-0.003882, -0.048116, 1.051998]]),
    "deuteranomaly (severity 1.0)": np.array([[0.367322, 0.860646, -0.227968],
                                              [0.280085, 0.672501, 0.047413],
                                              [-0.011820, 0.042940, 0.968881]]),
    "tritanomaly (severity 1.0)": np.array([[1.255528, -0.076749, -0.178779],
                                            [-0.078411, 0.930809, 0.147602],
                                            [0.004733, 0.691367, 0.303900]]),
}


def simulate(rgb: np.ndarray, kind: str) -> np.ndarray:
    """rgb: float array (..., 3) in encoded sRGB 0-1."""
    lin = _lin(rgb[..., :3])
    sim = lin @ MACHADO_1[kind].T
    return _enc(sim)


def cvd_sheet(png_path: str | pathlib.Path, out_path: str | pathlib.Path, title: str | None = None) -> pathlib.Path:
    """Original and the three simulations of one rendered PNG, stacked, labelled; written to out_path (PNG)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    img = mpimg.imread(str(png_path))
    if img.dtype != np.float32 and img.dtype != np.float64:
        img = img.astype(float) / 255.0
    rgb = img[..., :3]
    panels = [("original", rgb)] + [(k, simulate(rgb, k)) for k in MACHADO_1]
    h, w = rgb.shape[:2]
    fig_w = 8.0
    fig_h = fig_w * h / w * 4 + 1.2
    fig, axes = plt.subplots(4, 1, figsize=(fig_w, fig_h))
    for ax, (lab, im) in zip(axes, panels):
        ax.imshow(np.clip(im, 0, 1), interpolation="lanczos")
        ax.set_title(lab, fontsize=9, loc="left")
        ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.45 / fig_h) if title else None)   # keep the title clear of the first panel
    out_path = pathlib.Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps({k: [round(x, 4) for x in v] for k, v in check_palettes().items()}, indent=1))
    for p in sys.argv[1:]:
        print(cvd_sheet(p, pathlib.Path(p).with_name(pathlib.Path(p).stem + "_cvd.png")))
