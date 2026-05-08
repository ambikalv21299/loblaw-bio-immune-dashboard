"""
analysis.py
-----------
Parts 2-4: frequency summary, statistical comparison, and subset queries.
Writes CSV tables to outputs/ and plots to plots/.
"""

import os
import sqlite3
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

DB_PATH = "immune_trial.db"
OUT_DIR = "outputs"
PLT_DIR = "plots"
CELL_POPS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_LABELS = {
    "b_cell": "B Cell",
    "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell",
    "nk_cell": "NK Cell",
    "monocyte": "Monocyte",
}

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PLT_DIR, exist_ok=True)


# ── helpers ─────────────────────────────────────────────────────────────────

def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def query_df(sql: str, params=()) -> pd.DataFrame:
    with get_conn() as conn:
        return pd.read_sql_query(sql, conn, params=params)


# ── Part 2: relative frequency table ────────────────────────────────────────

def part2_frequency_table() -> pd.DataFrame:
    print("\n── Part 2: Cell-population frequency table ──")

    sql = """
        SELECT
            cc.sample_id                                                            AS sample,
            SUM(cc.count) OVER (PARTITION BY cc.sample_id)                         AS total_count,
            cc.population                                                           AS population,
            cc.count                                                                AS count,
            ROUND(cc.count * 100.0 / SUM(cc.count) OVER (PARTITION BY cc.sample_id), 4)
                                                                                    AS percentage
        FROM cell_counts cc
        ORDER BY cc.sample_id, cc.population
    """
    df = query_df(sql)
    out = os.path.join(OUT_DIR, "part2_frequency_table.csv")
    df.to_csv(out, index=False)
    print(f"  Rows: {len(df):,}  →  {out}")
    print(df.head(10).to_string(index=False))
    return df


# ── Part 3: responder vs non-responder comparison ───────────────────────────

def part3_statistics(freq_df: pd.DataFrame) -> None:
    print("\n── Part 3: Responder vs Non-Responder (melanoma, miraclib, PBMC) ──")

    # Filter: melanoma, miraclib, PBMC only
    sql = """
        SELECT s.sample_id, s.response
        FROM   samples s
        JOIN   subjects sub ON sub.subject_id = s.subject_id
        WHERE  sub.condition = 'melanoma'
          AND  s.treatment   = 'miraclib'
          AND  s.sample_type = 'PBMC'
          AND  s.response    IN ('yes','no')
    """
    meta = query_df(sql)
    merged = freq_df.merge(meta, left_on="sample", right_on="sample_id")

    # Mann-Whitney U per population (non-parametric, appropriate for this data)
    results = []
    for pop in CELL_POPS:
        sub = merged[merged["population"] == pop]
        resp    = sub[sub["response"] == "yes"]["percentage"].values
        nonresp = sub[sub["response"] == "no"]["percentage"].values
        stat, pval = stats.mannwhitneyu(resp, nonresp, alternative="two-sided")
        results.append({
            "population":       pop,
            "n_responders":     len(resp),
            "n_nonresponders":  len(nonresp),
            "median_resp":      round(float(np.median(resp)), 3),
            "median_nonresp":   round(float(np.median(nonresp)), 3),
            "mw_statistic":     round(float(stat), 3),
            "p_value":          round(float(pval), 6),
            "significant":      pval < 0.05,
        })

    stats_df = pd.DataFrame(results)

    # Bonferroni correction
    stats_df["p_bonferroni"] = (stats_df["p_value"] * len(CELL_POPS)).clip(upper=1.0).round(6)
    stats_df["significant_corrected"] = stats_df["p_bonferroni"] < 0.05

    out = os.path.join(OUT_DIR, "part3_statistics.csv")
    stats_df.to_csv(out, index=False)
    print(stats_df.to_string(index=False))
    print(f"\n  → {out}")

    sig = stats_df[stats_df["significant"]]["population"].tolist()
    sig_corr = stats_df[stats_df["significant_corrected"]]["population"].tolist()
    print(f"\n  Significant (p<0.05):              {sig if sig else 'none'}")
    print(f"  Significant (Bonferroni-corrected): {sig_corr if sig_corr else 'none'}")

    # Boxplot
    _plot_boxplot(merged)

    return stats_df


def _plot_boxplot(merged: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, len(CELL_POPS), figsize=(18, 6), sharey=False)
    fig.patch.set_facecolor("#0f1117")

    colors = {"yes": "#4ade80", "no": "#f87171"}
    labels = {"yes": "Responder", "no": "Non-Responder"}

    for ax, pop in zip(axes, CELL_POPS):
        sub = merged[merged["population"] == pop]
        data_resp    = sub[sub["response"] == "yes"]["percentage"].values
        data_nonresp = sub[sub["response"] == "no"]["percentage"].values

        bp = ax.boxplot(
            [data_resp, data_nonresp],
            patch_artist=True,
            widths=0.5,
            medianprops=dict(color="white", linewidth=2),
            whiskerprops=dict(color="#94a3b8"),
            capprops=dict(color="#94a3b8"),
            flierprops=dict(marker="o", color="#94a3b8", markersize=4, alpha=0.6),
        )
        bp["boxes"][0].set_facecolor(colors["yes"])
        bp["boxes"][0].set_alpha(0.75)
        bp["boxes"][1].set_facecolor(colors["no"])
        bp["boxes"][1].set_alpha(0.75)

        ax.set_facecolor("#1e2130")
        ax.set_title(POP_LABELS[pop], color="white", fontsize=11, pad=8)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["R", "NR"], color="#94a3b8")
        ax.tick_params(colors="#94a3b8")
        ax.spines[:].set_color("#334155")
        ax.yaxis.label.set_color("#94a3b8")
        ax.set_ylabel("Relative Frequency (%)", color="#94a3b8", fontsize=9)

    legend_handles = [
        mpatches.Patch(color=colors["yes"], alpha=0.75, label="Responder"),
        mpatches.Patch(color=colors["no"],  alpha=0.75, label="Non-Responder"),
    ]
    fig.legend(handles=legend_handles, loc="upper right",
               facecolor="#1e2130", edgecolor="#334155", labelcolor="white")
    fig.suptitle(
        "Cell Population Frequencies: Responders vs Non-Responders\n"
        "(Melanoma · Miraclib · PBMC)",
        color="white", fontsize=13, y=1.01,
    )
    plt.tight_layout()
    path = os.path.join(PLT_DIR, "part3_boxplot.png")
    plt.savefig(path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"  Boxplot saved → {path}")


# ── Part 4: subset analysis ──────────────────────────────────────────────────

def part4_subset() -> None:
    print("\n── Part 4: Subset – melanoma PBMC baseline miraclib ──")

    # 4a: Identify samples
    sql_samples = """
        SELECT s.sample_id, sub.project, sub.condition, sub.sex,
               s.response, s.time_from_treatment_start, s.sample_type, s.treatment
        FROM   samples s
        JOIN   subjects sub ON sub.subject_id = s.subject_id
        WHERE  sub.condition              = 'melanoma'
          AND  s.sample_type              = 'PBMC'
          AND  s.time_from_treatment_start = 0
          AND  s.treatment                = 'miraclib'
    """
    base = query_df(sql_samples)
    out_base = os.path.join(OUT_DIR, "part4a_baseline_samples.csv")
    base.to_csv(out_base, index=False)
    print(f"\n  Total baseline samples: {len(base)}")
    print(f"  → {out_base}")

    # 4b-i: Samples per project
    by_project = base.groupby("project").size().reset_index(name="n_samples")
    print("\n  Samples per project:")
    print(by_project.to_string(index=False))
    by_project.to_csv(os.path.join(OUT_DIR, "part4b_samples_per_project.csv"), index=False)

    # 4b-ii: Responders / non-responders (unique subjects)
    # Need unique subjects per sample (baseline = 1 sample per subject here)
    by_response = base[base["response"].isin(["yes","no"])].groupby("response").size().reset_index(name="n_subjects")
    print("\n  Subjects by response:")
    print(by_response.to_string(index=False))
    by_response.to_csv(os.path.join(OUT_DIR, "part4b_response_counts.csv"), index=False)

    # 4b-iii: Males / females
    by_sex = base.groupby("sex").size().reset_index(name="n_subjects")
    print("\n  Subjects by sex:")
    print(by_sex.to_string(index=False))
    by_sex.to_csv(os.path.join(OUT_DIR, "part4b_sex_counts.csv"), index=False)

    # 4c: Average B cells – melanoma, male, responders, time=0
    sql_bcell = """
        SELECT ROUND(AVG(cc.count), 2) AS avg_b_cells
        FROM   cell_counts cc
        JOIN   samples  s   ON s.sample_id  = cc.sample_id
        JOIN   subjects sub ON sub.subject_id = s.subject_id
        WHERE  sub.condition              = 'melanoma'
          AND  sub.sex                   = 'M'
          AND  s.treatment               = 'miraclib'
          AND  s.response                = 'yes'
          AND  s.time_from_treatment_start = 0
          AND  s.sample_type             = 'PBMC'
          AND  cc.population             = 'b_cell'
    """
    avg_df = query_df(sql_bcell)
    avg_bcell = avg_df["avg_b_cells"].iloc[0]
    print(f"\n  Avg B cells (melanoma, male, responder, t=0): {avg_bcell:.2f}")
    avg_df.to_csv(os.path.join(OUT_DIR, "part4c_avg_bcell.csv"), index=False)

    return base, by_project, by_response, by_sex, avg_bcell


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    freq_df   = part2_frequency_table()
    stats_df  = part3_statistics(freq_df)
    part4_subset()
    print("\n✓ All outputs written to outputs/ and plots/")
