"""
Main Analysis Pipeline -- EPIC Array Aging Clocks (Bio-Learn)
=============================================================
Runs the assignment requirements end-to-end:

  (1) Loads two methylation datasets (real GEO data if available, else
      simulated stand-ins).
  (2) Runs >= 8 aging clocks per dataset (we run 14, including
      DunedinPACE which uses pace-of-aging rather than absolute age).
  (3) Saves per-dataset prediction tables to results/tables/.
  (4) Generates the required visualizations:
        - correlation matrix across clocks (per dataset)
        - clock chronological-age deviation heatmap (per dataset)
        - predicted vs chronological age scatter (per dataset)
        - MAE comparison bar chart (both datasets, one plot)
        - predicted-age distribution box plots (both datasets, one plot)
  (5) Writes a results/summary.md report.

Run from the repo root:
    python epic_array/src/analysis.py

Author: Asim Ahmed (BSBI-2023, NUST SINES)
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

# Local imports
SRC_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SRC_DIR))

from simulate_data import make_dataset_a, make_dataset_b  # noqa: E402

from biolearn.data_library import GeoData  # noqa: E402
from biolearn.model_gallery import ModelGallery  # noqa: E402

# --- Patch DunedinPACE: biolearn's quantile_normalize_using_target mutates
#     its input array in-place, but pandas 2.x BlockManager hands out
#     read-only views in some configurations. Wrap the function so it
#     copies the input first. Idempotent and zero-impact on real data.
import biolearn.dunedin_pace as _dp  # noqa: E402

_DP_ORIG = _dp.quantile_normalize_using_target


def _dp_patched(data, target_values):
    return _DP_ORIG(np.array(data, copy=True), target_values)


_dp.quantile_normalize_using_target = _dp_patched


# ------------------------------------------------------------------ Config
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
FIG_DIR = RESULTS_DIR / "figures"
TBL_DIR = RESULTS_DIR / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)

# Fourteen aging clocks. All blood-applicable. The first 13 output age in
# years; DunedinPACE outputs aging rate (years/year, ~1.0 = normal pace).
SELECTED_CLOCKS: list[str] = [
    # --- Classic linear / elastic-net clocks (chronological age) ---
    "Horvathv1",       # Horvath (2013), pan-tissue, 353 CpGs
    "Horvathv2",       # Horvath (2018), skin + blood refinement, 391 CpGs
    "Hannum",          # Hannum (2013), blood-specific, 71 CpGs
    "Lin",             # Lin et al. (2016), 99-CpG blood clock
    "VidalBralo",      # Vidal-Bralo et al. (2018), 8-CpG compact clock
    # --- Biological / mortality-aware ---
    "PhenoAge",        # Levine et al. (2018), 9-biomarker biological age
    "HRSInCHPhenoAge", # Higgins-Chen (2022), PC-imputed PhenoAge
    # --- Causality-decomposed (Ying 2022 trio) ---
    "YingCausAge",     # CpGs causal of aging
    "YingDamAge",      # CpGs reflecting age-related damage
    "YingAdaptAge",    # CpGs reflecting adaptive responses
    # --- Stochastic clocks (Tong 2024) ---
    "StocP",           # Stochastic PhenoAge variant
    "StocH",           # Stochastic Horvath variant
    # --- Deep learning ---
    "AltumAge",        # de Lima Camillo (2022), TabNet deep neural net
    # --- Pace-of-aging (output is years/year, not years) ---
    "DunedinPACE",     # Belsky et al. (2022), Dunedin pace-of-aging clock
]

# Clocks whose predictions are NOT in years.
PACE_CLOCKS: set[str] = {"DunedinPACE", "DunedinPoAm38"}


# --------------------------------------------------------------- Utilities
def load_or_simulate() -> dict[str, GeoData]:
    """
    Try real GEO data first; fall back to simulator if download fails.
    Labels are GSE accession numbers when real data are available;
    simulator labels when not.
    """
    use_real = os.environ.get("USE_REAL_DATA", "0") == "1"
    if use_real:
        try:
            from download_real_data import load_real_datasets  # noqa: WPS433
            real = load_real_datasets()
            return {
                "GSE40279": real["GSE40279"],
                "GSE41169": real["GSE41169"],
            }
        except Exception as exc:  # pragma: no cover
            print(
                f"[analysis] Real-data download failed ({exc}); "
                "falling back to simulator.",
                file=sys.stderr,
            )

    print("[analysis] Generating simulated dataset A (60 samples, 25-75 yr) ...")
    dnam_a, meta_a = make_dataset_a()
    print("[analysis] Generating simulated dataset B (50 samples, 30-90 yr) ...")
    dnam_b, meta_b = make_dataset_b()

    return {
        "DatasetA_sim": GeoData(meta_a, dnam_a),
        "DatasetB_sim": GeoData(meta_b, dnam_b),
    }


def run_all_clocks(
    geo: GeoData, clock_names: Optional[list[str]] = None
) -> pd.DataFrame:
    """Run each clock; return DataFrame indexed by sample, one column per clock."""
    if clock_names is None:
        clock_names = SELECTED_CLOCKS
    gallery = ModelGallery()
    out = pd.DataFrame(index=geo.metadata.index)
    for name in clock_names:
        try:
            t0 = time.time()
            model = gallery.get(name)
            pred = model.predict(geo)
            out[name] = pred.iloc[:, 0].reindex(out.index)
            print(f"  - {name:<18} done in {time.time() - t0:5.1f}s")
        except Exception as exc:
            print(f"  - {name:<18} FAILED: {exc}")
            out[name] = np.nan
    return out


# ----------------------------------------------------------- Visualizations
def plot_correlation_matrix(
    preds: pd.DataFrame, ages: pd.Series, sex: pd.Series,
    out_path: Path, title: str
) -> None:
    """Pearson correlation between all clock predictions, age, and sex."""
    df = preds.copy()
    df["ChronologicalAge"] = ages.reindex(df.index)
    df["Sex"] = sex.reindex(df.index)
    corr = df.corr(method="pearson")
    plt.figure(figsize=(11, 9))
    sns.heatmap(
        corr, annot=True, fmt=".2f", cmap="RdBu_r",
        vmin=-1, vmax=1, square=True,
        cbar_kws={"label": "Pearson r"},
        linewidths=0.4, linecolor="white",
        annot_kws={"fontsize": 8},
    )
    plt.title(f"Clock-vs-Clock Correlation Matrix\n{title}", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_deviation_heatmap(
    preds: pd.DataFrame, ages: pd.Series, out_path: Path, title: str
) -> None:
    """
    Heatmap of (predicted_age - chronological_age) per sample x clock.
    Pace clocks (DunedinPACE) are excluded since their unit is years/year.
    """
    age_clocks = [c for c in preds.columns if c not in PACE_CLOCKS]
    deviations = preds[age_clocks].sub(ages, axis=0)
    order = ages.sort_values().index
    deviations = deviations.reindex(order)
    vmax = float(np.nanpercentile(np.abs(deviations.values), 98))
    plt.figure(figsize=(11, max(6, 0.20 * len(deviations))))
    sns.heatmap(
        deviations, cmap="RdBu_r", center=0,
        vmin=-vmax, vmax=vmax,
        cbar_kws={"label": "Predicted - Chronological Age (yr)"},
        yticklabels=False,
    )
    plt.title(
        f"Per-sample Age Deviation by Clock\n{title} "
        f"(rows = samples, sorted ascending by chronological age; "
        f"pace clocks excluded)",
        fontsize=12,
    )
    plt.xlabel("Clock")
    plt.ylabel("Samples (youngest at top -> oldest at bottom)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_age_scatter_grid(
    preds: pd.DataFrame, ages: pd.Series, out_path: Path, title: str
) -> None:
    """Predicted vs chronological age scatter for every clock (excl. pace)."""
    age_clocks = [c for c in preds.columns if c not in PACE_CLOCKS]
    n_clocks = len(age_clocks)
    n_cols = 3
    n_rows = int(np.ceil(n_clocks / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.4 * n_cols, 4.0 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    age_min = float(ages.min()) - 5
    age_max = float(ages.max()) + 5
    pred_min = float(np.nanmin(preds[age_clocks].values)) - 5
    pred_max = float(np.nanmax(preds[age_clocks].values)) + 5
    lim_min = min(age_min, pred_min)
    lim_max = max(age_max, pred_max)

    for i, clock in enumerate(age_clocks):
        ax = axes[i]
        x = ages.values
        y = preds[clock].values
        valid = ~np.isnan(y)
        if valid.sum() >= 3:
            r, _ = pearsonr(x[valid], y[valid])
            mae = float(np.mean(np.abs(y[valid] - x[valid])))
        else:
            r = float("nan")
            mae = float("nan")

        ax.scatter(x, y, alpha=0.65, s=22, edgecolor="white", linewidth=0.4)
        ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=1, alpha=0.6,
                label="y = x")
        ax.set_title(f"{clock}\n(r = {r:.3f}, MAE = {mae:.1f} yr)", fontsize=10)
        ax.set_xlabel("Chronological age (yr)")
        ax.set_ylabel("Predicted age (yr)")
        ax.set_xlim(lim_min, lim_max)
        ax.set_ylim(lim_min, lim_max)
        ax.grid(alpha=0.25)

    for j in range(n_clocks, len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"Aging Clock Predictions vs Chronological Age -- {title}",
                 fontsize=13, y=1.00)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_mae_comparison(
    summaries: dict[str, pd.DataFrame], out_path: Path
) -> None:
    """Grouped bar chart of MAE per clock across both datasets (excl. pace)."""
    plt.figure(figsize=(13, 5.5))
    labels = list(summaries.keys())
    palette = sns.color_palette("Set2", n_colors=len(labels))

    common: Optional[list[str]] = None
    for s in summaries.values():
        clocks = [c for c in s["Clock"].tolist() if c not in PACE_CLOCKS]
        common = clocks if common is None else [c for c in common if c in clocks]
    common = common or []

    x = np.arange(len(common))
    width = 0.8 / max(1, len(labels))

    for i, (label, s) in enumerate(summaries.items()):
        s_idx = s.set_index("Clock")
        mae_vals = [s_idx.loc[c, "MAE_years"] for c in common]
        offset = (i - (len(labels) - 1) / 2) * width
        bars = plt.bar(x + offset, mae_vals, width=width, label=label,
                       color=palette[i], edgecolor="black", linewidth=0.5)
        for bar, v in zip(bars, mae_vals):
            plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                     f"{v:.1f}", ha="center", va="bottom", fontsize=7)

    plt.xticks(x, common, rotation=40, ha="right")
    plt.ylabel("Mean Absolute Error (years)")
    plt.title("MAE per Aging Clock Across Both Datasets", fontsize=12)
    plt.legend(loc="upper left", fontsize=9, framealpha=0.9)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_predicted_age_distribution(
    predictions: dict[str, pd.DataFrame],
    metadata: dict[str, pd.Series],
    out_path: Path,
) -> None:
    """Box plots of predicted-age distributions per clock per dataset."""
    rows = []
    for ds_label, preds in predictions.items():
        for clock in preds.columns:
            if clock in PACE_CLOCKS:
                continue
            for sample, val in preds[clock].dropna().items():
                rows.append({"Dataset": ds_label, "Clock": clock,
                             "PredictedAge": float(val)})
    df = pd.DataFrame(rows)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(13, 6))
    sns.boxplot(data=df, x="Clock", y="PredictedAge", hue="Dataset",
                ax=ax, palette="Set2", linewidth=0.8, fliersize=2)

    cmap = sns.color_palette("Set2", n_colors=len(metadata))
    for i, (ds_label, ages) in enumerate(metadata.items()):
        ax.axhspan(ages.min(), ages.max(), alpha=0.10, color=cmap[i], zorder=0)
        ax.axhline(ages.median(), ls="--", color=cmap[i], alpha=0.55,
                   linewidth=1, label=f"{ds_label} median true age")

    ax.set_ylabel("Predicted age (yr)")
    ax.set_xlabel("")
    ax.set_title("Predicted-age Distribution per Clock\n"
                 "(shaded bands = chronological-age range; "
                 "dashed lines = chronological-age median)",
                 fontsize=11)
    plt.xticks(rotation=40, ha="right")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9, ncol=2)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------------------- Driver
def analyse_dataset(
    label: str, geo: GeoData
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Run clocks + plots for a single dataset."""
    print(f"\n{'=' * 70}\nDataset: {label}\n{'=' * 70}")
    print(f"  shape: {geo.dnam.shape[1]} samples x {geo.dnam.shape[0]} CpGs")
    if "age" in geo.metadata.columns:
        print(f"  age range: {geo.metadata['age'].min():.1f} - "
              f"{geo.metadata['age'].max():.1f} yr "
              f"(median {geo.metadata['age'].median():.1f})")

    preds = run_all_clocks(geo)
    preds.to_csv(TBL_DIR / f"predictions_{label}.csv")

    ages = geo.metadata["age"].astype(float)
    sex = geo.metadata.get("sex", pd.Series(np.nan, index=ages.index))
    rows = []
    for clock in preds.columns:
        y = preds[clock].values
        x = ages.values
        valid = ~np.isnan(y)
        if valid.sum() < 3:
            rows.append([clock, np.nan, np.nan, np.nan, np.nan])
            continue
        r, _ = pearsonr(x[valid], y[valid])
        if clock in PACE_CLOCKS:
            rows.append([clock, r, np.nan, np.nan, np.nan])
            continue
        mae = float(np.mean(np.abs(y[valid] - x[valid])))
        bias = float(np.mean(y[valid] - x[valid]))
        rmse = float(np.sqrt(np.mean((y[valid] - x[valid]) ** 2)))
        rows.append([clock, r, mae, bias, rmse])
    summary = pd.DataFrame(
        rows,
        columns=["Clock", "Pearson_r", "MAE_years",
                 "Mean_bias_years", "RMSE_years"],
    )
    summary.to_csv(TBL_DIR / f"summary_{label}.csv", index=False)
    print("\n  Per-clock summary:")
    print(summary.round(3).to_string(index=False))

    plot_correlation_matrix(preds, ages, sex,
                            FIG_DIR / f"correlation_matrix_{label}.png", label)
    plot_deviation_heatmap(preds, ages,
                           FIG_DIR / f"age_deviation_heatmap_{label}.png", label)
    plot_age_scatter_grid(preds, ages,
                          FIG_DIR / f"age_prediction_{label}.png", label)
    print(f"  -> Per-dataset figures saved under {FIG_DIR}/")
    return preds, summary, ages


def write_summary_report(per_dataset: dict[str, dict]) -> None:
    out = RESULTS_DIR / "summary.md"
    lines: list[str] = ["# EPIC Array Aging Clocks -- Analysis Summary", ""]
    lines.append(
        "Generated by `epic_array/src/analysis.py`. "
        "Each section reports clock-by-clock performance against the "
        "chronological age provided in dataset metadata. DunedinPACE "
        "outputs aging *rate* (years/year) so MAE/bias/RMSE are not "
        "meaningful for it -- only Pearson correlation is reported."
    )
    lines.append("")
    for label, payload in per_dataset.items():
        lines.append(f"## {label}")
        lines.append("")
        lines.append("| Clock | Pearson r | MAE (yr) | Mean bias (yr) | RMSE (yr) |")
        lines.append("|---|---:|---:|---:|---:|")
        for _, row in payload["summary"].iterrows():
            mae = "—" if pd.isna(row["MAE_years"]) else f"{row['MAE_years']:.2f}"
            bias = "—" if pd.isna(row["Mean_bias_years"]) else f"{row['Mean_bias_years']:.2f}"
            rmse = "—" if pd.isna(row["RMSE_years"]) else f"{row['RMSE_years']:.2f}"
            lines.append(
                f"| {row['Clock']} | {row['Pearson_r']:.3f} | "
                f"{mae} | {bias} | {rmse} |"
            )
        lines.append("")
    lines.append("## Cross-dataset comparison")
    lines.append("")
    lines.append("- `figures/mae_comparison.png` -- MAE per clock across both datasets")
    lines.append("- `figures/predicted_age_distribution.png` -- box plots per clock per dataset")
    lines.append("")
    out.write_text("\n".join(lines))
    print(f"\n[analysis] Wrote {out}")


def main() -> None:
    print("[analysis] Loading datasets ...")
    datasets = load_or_simulate()
    per_dataset: dict[str, dict] = {}
    predictions: dict[str, pd.DataFrame] = {}
    summaries: dict[str, pd.DataFrame] = {}
    metadata_ages: dict[str, pd.Series] = {}
    for label, geo in datasets.items():
        preds, summary, ages = analyse_dataset(label, geo)
        per_dataset[label] = {"preds": preds, "summary": summary}
        predictions[label] = preds
        summaries[label] = summary
        metadata_ages[label] = ages

    print("\n[analysis] Generating cross-dataset comparison plots ...")
    plot_mae_comparison(summaries, FIG_DIR / "mae_comparison.png")
    plot_predicted_age_distribution(
        predictions, metadata_ages, FIG_DIR / "predicted_age_distribution.png"
    )
    print("  -> mae_comparison.png + predicted_age_distribution.png saved")

    write_summary_report(per_dataset)
    print("\n[analysis] Done.")


if __name__ == "__main__":
    main()
