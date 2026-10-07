"""Print-size typography, export and layout checks shared by all figures."""
from __future__ import annotations
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parents[1]
RESDIR, DATADIR, FIGDIR = ROOT / "results", ROOT / "data", ROOT / "figures"
CONFIG = ROOT / "config" / "figure_style.json"
BLUE, ORANGE, PURPLE, GREEN = "#0072B2", "#D55E00", "#8864A5", "#009E73"
INK, INK2, MUTED = "#202020", "#505050", "#92979D"
COL1, COL2 = 83 / 25.4, 171 / 25.4
EXPORTS, SETTINGS = [], {}


def use_paper_style(*, output_dir=None, dpi=None, formats=None):
    global FIGDIR, SETTINGS
    SETTINGS = json.loads(CONFIG.read_text())
    if output_dir is not None:
        FIGDIR = Path(output_dir).resolve()
    if dpi is not None:
        SETTINGS["raster_dpi"] = dpi
    if formats is not None:
        SETTINGS["formats"] = list(formats)
    EXPORTS.clear()
    font_manager.fontManager.ttflist = [
        font for font in font_manager.fontManager.ttflist
        if font.name != SETTINGS["font_family"]]
    for font_path in sorted((ROOT / "assets/fonts").glob("*.ttf")):
        font_manager.fontManager.addfont(font_path)
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [SETTINGS["font_family"]],
        "font.size": SETTINGS["text_pt"], "axes.labelsize": SETTINGS["axis_label_pt"],
        "axes.labelpad": 4, "axes.titlesize": SETTINGS["text_pt"],
        "xtick.labelsize": SETTINGS["tick_pt"], "ytick.labelsize": SETTINGS["tick_pt"],
        "legend.fontsize": SETTINGS["legend_pt"], "legend.frameon": False,
        "legend.handlelength": 2.1, "legend.handletextpad": .55,
        "legend.columnspacing": 1.25, "axes.edgecolor": INK2, "axes.linewidth": .75,
        "axes.facecolor": "white", "figure.facecolor": "white", "savefig.facecolor": "white",
        "xtick.color": INK2, "ytick.color": INK2, "lines.linewidth": 1.35,
        "mathtext.fontset": "dejavusans", "mathtext.default": "regular",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "svg.hashsalt": "donnan-membrane-partitioning-v0.3.3",
        "axes.unicode_minus": True, "savefig.bbox": None, "savefig.pad_inches": 0,
        "path.simplify": False,
    })


def journal_figure(height_mm, *, single_column=False):
    width = SETTINGS["single_column_mm" if single_column else "double_column_mm"]
    return plt.figure(figsize=(width / 25.4, height_mm / 25.4))


def style_axes(ax, *, grid="y"):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_axisbelow(True)
    ax.tick_params(which="major", direction="out", length=3, width=.75, pad=3)
    ax.tick_params(which="minor", direction="out", length=1.7, width=.5)
    if grid:
        ax.grid(axis=grid, which="major", color="#D8DBDE", linewidth=.45)


def panel_label(ax, letter, *, dx=0, dy=1.025):
    return ax.text(dx, dy, f"({letter})", transform=ax.transAxes,
                   fontsize=SETTINGS["panel_label_pt"], fontweight="bold",
                   ha="left", va="bottom", color=INK)


def legend_axes(fig, bounds, handles, *, ncol=1, loc="center", **kwargs):
    ax = fig.add_axes(bounds)
    ax.set_axis_off()
    ax.legend(handles=handles, ncol=ncol, loc=loc, **kwargs)
    return ax


def _layout_check(fig):
    """Check text at print size and canvas bounds; visual review is also required."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.bbox.width, fig.bbox.height
    sizes, issues = [], []
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text():
            continue
        if text.axes is not None and text in text.axes.get_xticklabels(which="both"):
            x = text.get_position()[0]
            lo, hi = sorted(text.axes.get_xlim())
            if not lo <= x <= hi:
                continue
        if text.axes is not None and text in text.axes.get_yticklabels(which="both"):
            y = text.get_position()[1]
            lo, hi = sorted(text.axes.get_ylim())
            if not lo <= y <= hi:
                continue
        sizes.append(text.get_fontsize())
        if text.get_clip_on():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -.5 or box.y0 < -.5 or box.x1 > width+.5 or box.y1 > height+.5:
            issues.append(text.get_text())
    if issues:
        raise ValueError("Text extends outside figure canvas: " + repr(issues))
    return min(sizes) if sizes else None


def save(fig, stem):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    min_pt = _layout_check(fig)
    if min_pt < SETTINGS["minimum_text_pt"] - .01:
        raise ValueError(f"{stem}: text below configured minimum: {min_pt} pt")
    paths = []
    for fmt in SETTINGS["formats"]:
        suffix = "tif" if fmt in ("tif", "tiff") else fmt
        path = FIGDIR / f"{stem}.{suffix}"
        kwargs = {"dpi": SETTINGS["raster_dpi"]} if suffix in ("png", "tif") else {}
        if suffix == "tif":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        elif suffix == "pdf":
            kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
        elif suffix == "svg":
            kwargs["metadata"] = {"Date": None}
        # Render fully before replacing an existing export.  A failed render
        # must not leave a zero-length or partially written publication file.
        buffer = io.BytesIO()
        fig.savefig(buffer, format=suffix, **kwargs)
        payload = buffer.getvalue()
        with tempfile.NamedTemporaryFile(dir=FIGDIR, prefix=f".{stem}-",
                                         suffix=".tmp", delete=False) as staged:
            staged.write(payload)
            staged.flush()
            os.fsync(staged.fileno())
            staged_path = Path(staged.name)
        try:
            os.replace(staged_path, path)
        finally:
            staged_path.unlink(missing_ok=True)
        paths.append({"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    EXPORTS.append({"stem": stem, "width_mm": fig.get_figwidth()*25.4,
                    "height_mm": fig.get_figheight()*25.4, "minimum_text_pt": min_pt,
                    "layout_bounds_check": "passed", "files": paths})
    plt.close(fig)
    print(f"  wrote {stem}: " + ", ".join(p["path"] for p in paths))


def write_export_manifest():
    manifest = {"style": SETTINGS, "figures": EXPORTS,
                "figure_data": json.loads((FIGDIR / "source_data/manifest.json").read_text())}
    (FIGDIR / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
