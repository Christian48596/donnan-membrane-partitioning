"""Draw manuscript figures with direct legends and brief grey-region labels."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, NullFormatter
from figure_data import export_figure_data
from figstyle import (BLUE, ORANGE, PURPLE, GREEN, INK, INK2, MUTED,
                      use_paper_style, style_axes, panel_label, journal_figure,
                      legend_axes, save, write_export_manifest)

ROOT = Path(__file__).resolve().parents[1]
SALT_LABEL = {"NaCl": "NaCl", "MgCl2": r"MgCl$_2$",
              "Na2SO4": r"Na$_2$SO$_4$", "MgSO4": r"MgSO$_4$"}
MEM_LABEL = {"Nafion117": "Nafion 117", "CR-61": "CR-61"}
SALT_COLOR = dict(zip(SALT_LABEL, (BLUE, ORANGE, GREEN, PURPLE)))
SALT_MARKER = dict(zip(SALT_LABEL, ("o", "s", "^", "D")))
ORDER = [("Nafion117", "NaCl"), ("Nafion117", "MgCl2"),
         ("Nafion117", "Na2SO4"), ("CR-61", "Na2SO4"),
         ("Nafion117", "MgSO4"), ("CR-61", "MgSO4")]
MODELS = [("co_ideal", BLUE, (0, (5, 2.2)), "Ideal"),
          ("co_manning", ORANGE, "-", "Manning-Donnan"),
          ("paired", PURPLE, (0, (5, 1.6, 1.2, 1.6)), "Species equilibrium"),
          ("published_closure_equal_alpha", GREEN, (0, (1.2, 1.5)),
           "Sujanani et al. (2023)\n"+r"$\alpha^{\rm m}=\alpha^{\rm s}$")]
CONCENTRATION = r"(mol L$^{-1}$)"


def series_handles(*, multiline=False, points_only=False):
    handles = []
    for membrane, salt in ORDER:
        color = SALT_COLOR[salt]
        label = MEM_LABEL[membrane] + ("\n" if multiline else " / ") + SALT_LABEL[salt]
        handles.append(Line2D(
            [], [], color=color, marker=SALT_MARKER[salt],
            linestyle="" if points_only else ("--" if membrane == "CR-61" else "-"),
            linewidth=1, markersize=4.5,
            markerfacecolor="white" if membrane == "CR-61" else color,
            label=label))
    return handles


def model_handles():
    return [Line2D([], [], color=color, linestyle=style, label=label)
            for _, color, style, label in MODELS]


def local_legend(ax, handles, *, loc="upper left", ncol=1, **kwargs):
    return ax.legend(handles=handles, loc=loc, ncol=ncol, borderaxespad=.3,
                     labelspacing=.35, handlelength=1.9, **kwargs)


def fig1_screening(source):
    curve = pd.read_csv(source / "fig1_screening_curve.csv")
    fig = journal_figure(85, single_column=True)
    ax = fig.add_axes([.22, .18, .75, .655])
    Istar = curve.I_crossover_M.iloc[0]
    ax.axvspan(Istar, 12, facecolor="#EEF0F2", zorder=0)
    ax.plot(curve.I_M, curve.lambda_D_nm, color=BLUE, linewidth=1.5)
    ax.axhline(.30, color=INK2, linestyle="--", linewidth=1)
    ax.axvline(Istar, color=MUTED, linestyle=":", linewidth=.8)
    ax.plot(Istar, .30, "o", color=INK, markersize=4, markeredgecolor="white")
    ax.text(np.sqrt(Istar*12), 2, "Below\nion size", ha="center", va="center", color=INK2)
    ax.set(xscale="log", yscale="log", xlim=(1e-3, 12), ylim=(.07, 12),
           xlabel=r"$I$ "+CONCENTRATION, ylabel=r"$\lambda_{\rm D}$ (nm)")
    style_axes(ax)
    legend_axes(fig, [.14, .91, .84, .075],
                [Line2D([], [], color=BLUE, label="Debye ref."),
                 Line2D([], [], color=INK2, linestyle="--", label=r"$d=0.30$ nm")],
                ncol=2, columnspacing=.8, handlelength=1.7)
    save(fig, "fig1_screening_regimes")


def fig2_sorption(source):
    table = pd.read_csv(source / "fig2_sorption.csv")
    fig = journal_figure(226)
    legend_axes(fig, [.10, .94, .885, .05], model_handles(), ncol=2,
                columnspacing=1.7, labelspacing=.5)
    gs = fig.add_gridspec(3, 2, left=.13, right=.985, bottom=.067, top=.835,
                         hspace=.78, wspace=.43)
    for k, (membrane, salt) in enumerate(ORDER):
        ax = fig.add_subplot(gs[k//2, k%2])
        group = table[(table.membrane == membrane) & (table.salt == salt)].sort_values("c_ext")
        for key, color, style, _ in MODELS:
            if group[key].notna().any():
                ax.plot(group.c_ext, group[key], color=color, linestyle=style, linewidth=1.3)
        ax.errorbar(group.c_ext, group.co_obs, yerr=group.co_sd, fmt="o", markersize=4.3,
                    markerfacecolor="white", markeredgecolor=INK, markeredgewidth=.9,
                    ecolor=INK, elinewidth=.8, capsize=2, capthick=.8, zorder=5)
        references = group.measurement_reference.unique()
        if len(references) != 1:
            raise ValueError("A panel must identify one measurement source")
        local_legend(ax, [Line2D([], [], color=INK, marker="o", linestyle="",
                                markerfacecolor="white", markersize=4.3,
                                label=MEM_LABEL[membrane]+" / "+SALT_LABEL[salt]+"\n"+references[0])],
                     loc="lower right", bbox_to_anchor=(1, 1.01))
        ax.set(xscale="log", yscale="log", xlim=(.006, 1.6), ylim=(3e-8, 2.5),
               xlabel=r"$C_s^{\rm s}$ "+CONCENTRATION,
               ylabel=(r"Cl$^-$ uptake " if salt in ("NaCl", "MgCl2") else
                       r"Total $\mathrm{SO}_4^{2-}$ ")+CONCENTRATION)
        ax.yaxis.set_major_locator(FixedLocator([1e-6, 1e-4, 1e-2, 1]))
        ax.yaxis.set_minor_formatter(NullFormatter())
        style_axes(ax); panel_label(ax, "abcdef"[k])
    save(fig, "fig2_coion_sorption")


def fig3_residuals(source):
    table = pd.read_csv(source / "fig3_residuals.csv")
    summary = pd.read_csv(source / "fig3_series_means.csv")
    fig = journal_figure(190)
    legend_axes(fig, [.07, .895, .915, .095], series_handles(), ncol=2,
                columnspacing=1.7, labelspacing=.4)
    for left, residual, letter in [(.16, "r_ideal", "a"),
                                   (.655, "r_manning", "b")]:
        ax = fig.add_axes([left, .565, .33, .275])
        ax.axhline(0, color=MUTED, linestyle="--", linewidth=.8)
        for membrane, salt in ORDER:
            group = table[(table.membrane == membrane) & (table.salt == salt)].sort_values("I_ext")
            color = SALT_COLOR[salt]
            ax.plot(group.I_ext, group[residual], color=color, marker=SALT_MARKER[salt],
                    markersize=4.5, markeredgewidth=.85, linewidth=1,
                    linestyle="--" if membrane == "CR-61" else "-",
                    markerfacecolor="white" if membrane == "CR-61" else color)
        ax.set(xscale="log", xlim=(.006, 6), ylim=(-1.15, 8.6), yticks=[0, 2, 4, 6, 8],
               xlabel=r"$I^{\rm s}$ "+CONCENTRATION, ylabel=r"$r$ (-)")
        style_axes(ax); panel_label(ax, letter)
    legend_axes(fig, [.16, .445, .825, .055],
                [Line2D([], [], color=color, marker=marker, linestyle="", label=label)
                 for color, marker, label in [(ORANGE, "o", "Manning-Donnan"),
                 (PURPLE, "D", "Species equilibrium"),
                 (GREEN, "s", "Sujanani et al. (2023)\n"+
                  r"$\alpha^{\rm m}=\alpha^{\rm s}$")]],
                ncol=3, columnspacing=1.2, handlelength=1)
    ax = fig.add_axes([.16, .095, .825, .325])
    ax.axvline(0, color=MUTED, linestyle="--", linewidth=.8)
    labels = []
    for y, (membrane, salt) in enumerate(ORDER[::-1]):
        row = summary[(summary.membrane == membrane) & (summary.salt == salt)].iloc[0]
        for key, offset, color, marker in [
                ("manning", .23, ORANGE, "o"), ("conditional", 0, PURPLE, "D"),
                ("published", -.23, GREEN, "s")]:
            mean = row[key+"_mean"]
            if np.isfinite(mean):
                ax.errorbar(mean, y+offset,
                            xerr=[[mean-row[key+"_min"]], [row[key+"_max"]-mean]],
                            fmt=marker, markersize=4.2, color=color,
                            linewidth=1.1, capsize=2)
        if y < 5:
            ax.axhline(y+.5, color="#E8EAEC", linewidth=.45, zorder=0)
        labels.append(MEM_LABEL[membrane]+"\n"+SALT_LABEL[salt])
    ax.set_yticks(range(6), labels)
    ax.set(ylim=(-.55, 5.6), xlim=(-1.4, 4.4),
           xlabel=r"$\bar r$ (-)", ylabel="Membrane / salt")
    style_axes(ax, grid="x"); panel_label(ax, "c")
    save(fig, "fig3_residuals")


def fig4_identifiability(source):
    data = pd.read_csv(source / "fig4_identifiability.csv")
    fig = journal_figure(90, single_column=True)
    ax = fig.add_axes([.215, .15, .755, .665])
    prospective = data[data.prospective_independent_counterions]
    x, width = np.arange(2), .34
    for side, basis, color in [(-1, "exact", BLUE), (1, "uncertain", ORANGE)]:
        values = prospective[f"power_baseline_{basis}"].to_numpy()
        low = prospective[f"power_{basis}_ci_low"].to_numpy()
        high = prospective[f"power_{basis}_ci_high"].to_numpy()
        xpos = x+side*width/2
        ax.bar(xpos, values, width, color=color)
        ax.errorbar(xpos, values, yerr=np.vstack([values-low, high-values]), fmt="none",
                    ecolor=INK2, elinewidth=.8, capsize=2, capthick=.8)
    ax.set_xticks(x, ["3 salts", "4 salts"])
    ax.set(ylim=(0, 1.05), yticks=np.linspace(0, 1, 6),
           xlabel="Salt panel", ylabel="Power (-)")
    style_axes(ax)
    legend_axes(fig, [.18, .89, .79, .095],
                [Patch(facecolor=BLUE, label="Exact baseline"),
                 Patch(facecolor=ORANGE, label="Uncertain baseline")], ncol=1)
    save(fig, "fig4_identifiability")


def fig5_charge_and_exchange(source):
    charge = pd.read_csv(source / "fig5_apparent_charge.csv")
    transfer = pd.read_csv(source / "fig5_transferability.csv")
    curves = pd.read_csv(source / "fig5_exchange_curves.csv")
    fig = journal_figure(190)
    ax = fig.add_axes([.145, .57, .34, .255])
    ax.axhspan(-1.2, 0, color="#EEF0F2")
    ax.axhline(0, color=INK2, linewidth=.8)
    ax.text(.96, .31, "Negative apparent\ncharge", transform=ax.transAxes,
            ha="right", va="center", color=INK2)
    for (salt, membrane), group in charge.groupby(["salt", "membrane"], sort=False):
        ax.plot(group.C_fixed, group.C_apparent, linestyle="", marker=SALT_MARKER[salt],
                markersize=4.5, color=SALT_COLOR[salt], markeredgewidth=.85,
                markerfacecolor="white" if membrane == "CR-61" else SALT_COLOR[salt])
    ax.plot([3, 4.9], [3, 4.9], color=MUTED, linestyle="--", linewidth=1)
    ax.set_yscale("symlog", linthresh=.05, linscale=.8)
    ax.set(xlim=(3, 4.9), ylim=(-1.2, 10), xticks=[3, 3.5, 4, 4.5],
           xlabel=r"$C_A$ "+CONCENTRATION,
           ylabel=r"$C_{A,\rm app}$ "+CONCENTRATION)
    ax.set_yticks([-1, -.1, 0, .1, 1, 10], ["-1", "-0.1", "0", "0.1", "1", "10"])
    ax.yaxis.set_minor_locator(FixedLocator([]))
    style_axes(ax); panel_label(ax, "a")
    legend_axes(fig, [.06, .875, .46, .105],
                series_handles(points_only=True)+
                [Line2D([], [], color=MUTED, linestyle="--", label="Identity")],
                ncol=2, columnspacing=1.1, labelspacing=.3, handlelength=1.3)
    ax = fig.add_axes([.645, .57, .34, .255])
    ax.plot(transfer.c_ext, transfer.ideal_transferred, color=BLUE, linestyle="--")
    ax.plot(transfer.c_ext, transfer.manning_measured, color=ORANGE)
    ax.errorbar(transfer.c_ext, transfer.obs, yerr=transfer.co_sd, fmt="o", markersize=4.5,
                markerfacecolor="white", markeredgecolor=INK, ecolor=INK,
                elinewidth=.8, capsize=2, zorder=4)
    ax.set(xscale="log", yscale="log", xlim=(.006, 1.5), ylim=(1e-3, 2),
           xlabel=r"$C_s^{\rm s}$ "+CONCENTRATION,
           ylabel=r"$c_{\rm Cl}^{\rm m}$ "+CONCENTRATION)
    style_axes(ax); panel_label(ax, "b")
    legend_axes(fig, [.635, .855, .35, .135],
                [Line2D([], [], color=INK, marker="o", linestyle="", markerfacecolor="white",
                        label="Measured Cl$^-$ uptake\nSujanani et al. (2021)"),
                 Line2D([], [], color=BLUE, linestyle="--", label=r"Ideal ($C_{A,\rm app}$)"),
                 Line2D([], [], color=ORANGE, label=r"Manning-Donnan ($C_A$)")], ncol=1)
    for col in range(2):
        ax = fig.add_axes([.145 if not col else .645, .10, .34, .275])
        handles = []
        cases = (.5, .1, .01) if col == 0 else (0., .5, 1.)
        for value, color in zip(cases, (BLUE, ORANGE, PURPLE)):
            selected = (np.isclose(curves.y, value) & np.isclose(curves.coion_to_fixed_charge, 0)
                        if col == 0 else np.isclose(curves.y, .5) &
                        np.isclose(curves.coion_to_fixed_charge, value))
            group = curves[selected].sort_values("K_ex")
            ax.plot(group.K_ex, group.alpha, color=color, linewidth=1.45)
            ax.axhline(group.alpha_asymptote.iloc[0], color=color, linestyle=(0, (3, 2)), linewidth=.8)
            label = f"$y={value:g}$" if not col else rf"$c_-^{{\rm m}}/C_A^{{\rm m}}={value:g}$"
            handles.append(Line2D([], [], color=color, label=label))
        handles.append(Line2D([], [], color=INK2, linestyle=(0, (3, 2)), label=r"$\alpha_\infty$"))
        ax.set(xscale="log", xlim=(.1, 1e5),
               xlabel=r"$K_{\rm ex}$ (-)", ylabel=r"$\alpha_{\rm Li/Na}$ (-)")
        if not col:
            ax.set(yscale="log", ylim=(.25, 400))
        else:
            ax.set(ylim=(.25, 5.85), yticks=[1, 2, 3, 4, 5])
        style_axes(ax); panel_label(ax, "cd"[col])
        legend_axes(fig, [.125 if not col else .62, .405, .375, .08],
                    handles, ncol=2, columnspacing=1, handlelength=1.8)
    save(fig, "fig5_dle")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"figures")
    parser.add_argument("--formats", nargs="+", choices=["pdf", "svg", "png", "tif"])
    parser.add_argument("--dpi", type=int)
    args = parser.parse_args()
    use_paper_style(output_dir=args.output_dir, dpi=args.dpi, formats=args.formats)
    source = export_figure_data(args.output_dir)
    fig1_screening(source)
    fig2_sorption(source)
    fig3_residuals(source)
    fig4_identifiability(source)
    fig5_charge_and_exchange(source)
    write_export_manifest()


if __name__ == "__main__":
    main()
