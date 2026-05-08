# Loblaw Bio · Immune Cell Dashboard

**Analyzes how miraclib affects immune cell populations in a clinical trial.**

---

## Quick Start (GitHub Codespaces)

```bash
make setup     # install dependencies
make pipeline  # run full data pipeline (Parts 1–4)
make dashboard # start dashboard at http://localhost:8050
```

---

## Project Structure

```
.
├── cell-count.csv        # Input data (converted from xlsx)
├── load_data.py          # Part 1: SQLite DB initialisation & loading
├── analysis.py           # Parts 2–4: frequency table, stats, subset analysis
├── dashboard.py          # Interactive Dash dashboard
├── requirements.txt
├── Makefile
├── outputs/              # CSV output tables
│   ├── part2_frequency_table.csv
│   ├── part3_statistics.csv
│   ├── part4a_baseline_samples.csv
│   ├── part4b_samples_per_project.csv
│   ├── part4b_response_counts.csv
│   ├── part4b_sex_counts.csv
│   └── part4c_avg_bcell.csv
└── plots/
    └── part3_boxplot.png
```

---

## Database Schema

Three normalized tables in `immune_trial.db`:

### `subjects`
| Column | Type | Notes |
|---|---|---|
| `subject_id` | TEXT PK | Unique patient identifier |
| `project` | TEXT | Clinical project (prj1/prj2/prj3) |
| `condition` | TEXT | melanoma / carcinoma / healthy |
| `age` | REAL | Age at enrollment |
| `sex` | TEXT | M / F |

### `samples`
| Column | Type | Notes |
|---|---|---|
| `sample_id` | TEXT PK | Unique biological sample |
| `subject_id` | TEXT FK→subjects | Parent subject |
| `sample_type` | TEXT | PBMC / WB |
| `treatment` | TEXT | miraclib / phauximab / none |
| `response` | TEXT | yes / no / NULL |
| `time_from_treatment_start` | INTEGER | Days since treatment start |

### `cell_counts`
| Column | Type | Notes |
|---|---|---|
| `id` | INTEGER PK | Auto-increment |
| `sample_id` | TEXT FK→samples | Parent sample |
| `population` | TEXT | b_cell / cd8_t_cell / cd4_t_cell / nk_cell / monocyte |
| `count` | REAL | Raw cell count |

### Design Rationale

**Normalization** — Subject demographics (age, sex, condition, project) live in `subjects` and are not repeated on every sample row. This eliminates update anomalies and reduces storage proportionally to the number of samples per subject.

**Long format for cell counts** — Each (sample, population) pair is a separate row in `cell_counts`. Compared to wide format (one column per population), this means:
- Adding new populations (e.g. dendritic cells, regulatory T cells) requires no schema changes — just new rows.
- `GROUP BY population` aggregations are natural SQL.
- The table scales linearly: 10 500 samples × 5 populations = 52 500 rows today; 100 000 samples × 20 populations = 2M rows — trivially handled by SQLite with proper indexes.

**Indexes** — `idx_cc_sample` (cell_counts.sample_id), `idx_cc_pop` (cell_counts.population), `idx_smp_filter` (treatment, response, time), `idx_sub_project` (project, condition) keep analytic queries fast.

**Scaling to hundreds of projects** — A `projects` table would be added (project_id PK, name, PI, start_date, protocol), with `subjects.project` becoming a FK. `condition` and `treatment` strings could be normalised into lookup tables to prevent free-text inconsistencies. For millions of samples, migrating to DuckDB or a columnar store (e.g. Parquet + Polars) would maintain the same schema design with far higher analytical throughput.

---

## Analysis Summary

### Part 2 — Frequency Table
Relative frequency of each cell population per sample, computed as:
`percentage = (population_count / total_sample_count) × 100`

### Part 3 — Statistical Comparison
Filters: melanoma · miraclib · PBMC only.  
Test: Mann-Whitney U (two-sided, non-parametric — appropriate for count-derived frequencies that are not guaranteed normal).  
Correction: Bonferroni across 5 populations.

**Result:** CD4 T Cell shows a nominally significant difference (p ≈ 0.013) between responders and non-responders, but this does not survive Bonferroni correction. No population reaches corrected significance at α = 0.05.

### Part 4 — Baseline Subset
Melanoma + PBMC + miraclib + time = 0:
- **prj1:** 384 samples | **prj3:** 272 samples
- Responders: 331 | Non-responders: 325
- Male: 344 | Female: 312
- **Average B cells (melanoma, male, responders, t=0): 10401.28**

---

## Dashboard

Run `make dashboard` and open **http://localhost:8050**.

The dashboard includes:
- **KPI strip** — total samples, subjects, projects, melanoma count
- **Part 2 tab** — filterable bar chart, pie chart, violin plot, and data table
- **Part 3 tab** — per-population boxplot, strip plot, all-populations comparison, stats table
- **Part 4 tab** — baseline subset donuts (project, response, sex) and key metrics
