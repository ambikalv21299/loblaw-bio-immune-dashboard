"""
dashboard.py
------------
Interactive Dash dashboard for Bob Loblaw's immune-cell clinical trial data.
Run: python dashboard.py
"""

import sqlite3

import dash
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, dash_table, dcc, html
from scipy import stats

DB_PATH = "immune_trial.db"

# Auto-build DB if missing
import subprocess, os as _os
if not _os.path.exists(DB_PATH):
    subprocess.run(["python", "load_data.py"], check=True)
CELL_POPS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_LABELS = {
    "b_cell": "B Cell",
    "cd8_t_cell": "CD8 T Cell",
    "cd4_t_cell": "CD4 T Cell",
    "nk_cell": "NK Cell",
    "monocyte": "Monocyte",
}

DARK_BG   = "#0f1117"
CARD_BG   = "#1a1d2e"
BORDER    = "#2d3154"
TEXT      = "#e2e8f0"
MUTED     = "#94a3b8"
ACCENT    = "#6366f1"
GREEN     = "#4ade80"
RED       = "#f87171"
YELLOW    = "#fbbf24"

PLOTLY_LAYOUT = dict(
    paper_bgcolor=CARD_BG,
    plot_bgcolor=CARD_BG,
    font=dict(color=TEXT, family="Inter, system-ui, sans-serif"),
    margin=dict(l=50, r=20, t=40, b=50),
    legend=dict(bgcolor=CARD_BG, bordercolor=BORDER),
    xaxis=dict(gridcolor=BORDER, linecolor=BORDER),
    yaxis=dict(gridcolor=BORDER, linecolor=BORDER),
)

# ── data loaders ────────────────────────────────────────────────────────────

def qdf(sql, params=()):
    with sqlite3.connect(DB_PATH) as conn:
        return pd.read_sql_query(sql, conn, params=params)


def load_freq_table():
    return qdf("""
        SELECT cc.sample_id AS sample,
               SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
               cc.population, cc.count,
               ROUND(cc.count*100.0/SUM(cc.count) OVER (PARTITION BY cc.sample_id),4) AS percentage,
               s.treatment, s.response, s.sample_type, s.time_from_treatment_start,
               sub.condition, sub.project, sub.sex, sub.age
        FROM cell_counts cc
        JOIN samples  s   ON s.sample_id   = cc.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
    """)


def load_part4():
    base = qdf("""
        SELECT s.sample_id, sub.project, sub.condition, sub.sex,
               s.response, s.time_from_treatment_start, s.sample_type, s.treatment
        FROM   samples s JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE  sub.condition='melanoma' AND s.sample_type='PBMC'
          AND  s.time_from_treatment_start=0 AND s.treatment='miraclib'
    """)
    avg_bc = qdf("""
        SELECT ROUND(AVG(cc.count),2) AS avg_b_cells
        FROM cell_counts cc
        JOIN samples  s   ON s.sample_id   = cc.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition='melanoma' AND sub.sex='M'
          AND s.treatment='miraclib' AND s.response='yes'
          AND s.time_from_treatment_start=0 AND s.sample_type='PBMC'
          AND cc.population='b_cell'
    """)["avg_b_cells"].iloc[0]
    return base, avg_bc


FREQ   = load_freq_table()
BASE4, AVG_BCELL = load_part4()

# ── app layout ───────────────────────────────────────────────────────────────

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.CYBORG,
                          "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"],
    title="Loblaw Bio · Immune Dashboard",
)
app.config.suppress_callback_exceptions = True


def stat_card(title, value, subtitle="", color=ACCENT):
    return dbc.Card(
        dbc.CardBody([
            html.P(title, className="mb-1", style={"color": MUTED, "fontSize": "0.75rem", "textTransform": "uppercase", "letterSpacing": "0.05em"}),
            html.H3(value, style={"color": color, "fontWeight": "700", "margin": "0"}),
            html.P(subtitle, style={"color": MUTED, "fontSize": "0.75rem", "margin": "0"}),
        ]),
        style={"background": CARD_BG, "border": f"1px solid {BORDER}", "borderRadius": "12px"},
    )


def section_header(title, subtitle=""):
    return html.Div([
        html.H4(title, style={"color": TEXT, "fontWeight": "600", "marginBottom": "2px"}),
        html.P(subtitle, style={"color": MUTED, "fontSize": "0.85rem", "marginBottom": "16px"}),
    ])


nav = dbc.Navbar(
    dbc.Container([
        html.Span("🧬", style={"fontSize": "1.4rem", "marginRight": "8px"}),
        dbc.NavbarBrand("Loblaw Bio · Immune Cell Dashboard", style={"color": TEXT, "fontWeight": "700"}),
        dbc.Nav([
            dbc.NavItem(dbc.NavLink("Overview", href="#overview", style={"color": MUTED})),
            dbc.NavItem(dbc.NavLink("Statistics", href="#stats", style={"color": MUTED})),
            dbc.NavItem(dbc.NavLink("Subset", href="#subset", style={"color": MUTED})),
        ], className="ms-auto", navbar=True),
    ], fluid=True),
    style={"background": CARD_BG, "borderBottom": f"1px solid {BORDER}"},
    dark=True,
)

# Top KPI cards
n_samples   = FREQ["sample"].nunique()
n_subjects  = qdf("SELECT COUNT(*) AS n FROM subjects")["n"].iloc[0]
n_projects  = qdf("SELECT COUNT(DISTINCT project) AS n FROM subjects")["n"].iloc[0]
melanoma_n  = qdf("SELECT COUNT(*) AS n FROM subjects WHERE condition='melanoma'")["n"].iloc[0]

kpi_row = dbc.Row([
    dbc.Col(stat_card("Total Samples",  f"{n_samples:,}",  "across all projects"), md=3),
    dbc.Col(stat_card("Total Subjects", f"{n_subjects:,}", "unique patients"), md=3),
    dbc.Col(stat_card("Projects",       str(n_projects),   "clinical projects"), md=3),
    dbc.Col(stat_card("Melanoma Subjects", str(melanoma_n), "target indication"), md=3),
], className="g-3 mb-4")

# Part 2 tab controls
part2_controls = dbc.Row([
    dbc.Col([
        html.Label("Condition", style={"color": MUTED, "fontSize": "0.8rem"}),
        dcc.Dropdown(
            id="p2-condition",
            options=[{"label": "All", "value": "All"}] +
                    [{"label": c.title(), "value": c} for c in sorted(FREQ["condition"].dropna().unique())],
            value="All", clearable=False,
            style={"background": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}"},
        ),
    ], md=3),
    dbc.Col([
        html.Label("Treatment", style={"color": MUTED, "fontSize": "0.8rem"}),
        dcc.Dropdown(
            id="p2-treatment",
            options=[{"label": "All", "value": "All"}] +
                    [{"label": t, "value": t} for t in sorted(FREQ["treatment"].dropna().unique())],
            value="All", clearable=False,
            style={"background": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}"},
        ),
    ], md=3),
    dbc.Col([
        html.Label("Sample Type", style={"color": MUTED, "fontSize": "0.8rem"}),
        dcc.Dropdown(
            id="p2-stype",
            options=[{"label": "All", "value": "All"}] +
                    [{"label": s, "value": s} for s in sorted(FREQ["sample_type"].dropna().unique())],
            value="All", clearable=False,
            style={"background": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}"},
        ),
    ], md=3),
    dbc.Col([
        html.Label("Timepoint (days)", style={"color": MUTED, "fontSize": "0.8rem"}),
        dcc.Dropdown(
            id="p2-time",
            options=[{"label": "All", "value": "All"}] +
                    [{"label": f"Day {int(t)}", "value": t}
                     for t in sorted(FREQ["time_from_treatment_start"].dropna().unique())],
            value="All", clearable=False,
            style={"background": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}"},
        ),
    ], md=3),
], className="g-2 mb-3")

# Part 3 controls
part3_controls = dbc.Row([
    dbc.Col([
        html.Label("Cell Population", style={"color": MUTED, "fontSize": "0.8rem"}),
        dcc.Dropdown(
            id="p3-pop",
            options=[{"label": POP_LABELS[p], "value": p} for p in CELL_POPS],
            value="b_cell", clearable=False,
            style={"background": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}"},
        ),
    ], md=4),
], className="g-2 mb-3")

layout = dbc.Container([
    nav,
    html.Div(style={"height": "24px"}),

    # KPI strip
    kpi_row,

    # ── Part 2 ─────────────────────────────────────────────────────
    html.Div(id="overview"),
    dbc.Card([
        dbc.CardBody([
            section_header("Part 2 · Cell Population Frequencies",
                           "Relative frequency (%) of each immune cell population per sample"),
            part2_controls,
            dbc.Row([
                dbc.Col(dcc.Graph(id="p2-bar",    config={"displayModeBar": False}), md=8),
                dbc.Col(dcc.Graph(id="p2-pie",    config={"displayModeBar": False}), md=4),
            ], className="mb-3"),
            dbc.Row([
                dbc.Col(dcc.Graph(id="p2-violin", config={"displayModeBar": False}), md=12),
            ], className="mb-3"),
            html.H6("Frequency Table (first 200 rows)", style={"color": MUTED, "fontSize": "0.8rem"}),
            html.Div(id="p2-table"),
        ])
    ], style={"background": CARD_BG, "border": f"1px solid {BORDER}", "borderRadius": "12px", "marginBottom": "24px"}),

    # ── Part 3 ─────────────────────────────────────────────────────
    html.Div(id="stats"),
    dbc.Card([
        dbc.CardBody([
            section_header("Part 3 · Responders vs Non-Responders",
                           "Melanoma · Miraclib · PBMC · Mann-Whitney U test with Bonferroni correction"),
            part3_controls,
            dbc.Row([
                dbc.Col(dcc.Graph(id="p3-box",  config={"displayModeBar": False}), md=6),
                dbc.Col(dcc.Graph(id="p3-strip", config={"displayModeBar": False}), md=6),
            ], className="mb-3"),
            dbc.Row([
                dbc.Col(dcc.Graph(id="p3-allbox", config={"displayModeBar": False}), md=12),
            ], className="mb-3"),
            html.H6("Statistical Summary (all populations)", style={"color": MUTED, "fontSize": "0.8rem"}),
            html.Div(id="p3-stats-table"),
        ])
    ], style={"background": CARD_BG, "border": f"1px solid {BORDER}", "borderRadius": "12px", "marginBottom": "24px"}),

    # ── Part 4 ─────────────────────────────────────────────────────
    html.Div(id="subset"),
    dbc.Card([
        dbc.CardBody([
            section_header("Part 4 · Baseline Subset Analysis",
                           "Melanoma · PBMC · Miraclib · Time = 0"),
            dbc.Row([
                dbc.Col(stat_card("Total Baseline Samples", str(len(BASE4)),   "melanoma + PBMC + miraclib + t=0"), md=3),
                dbc.Col(stat_card("Responders",   str((BASE4["response"]=="yes").sum()), "of baseline samples"), md=3),
                dbc.Col(stat_card("Non-Responders", str((BASE4["response"]=="no").sum()), "of baseline samples"), md=3),
                dbc.Col(stat_card("Avg B Cells (M, R, t=0)", f"{AVG_BCELL:,.2f}", "melanoma males, responders"), md=3),
            ], className="g-3 mb-4"),
            dbc.Row([
                dbc.Col(dcc.Graph(id="p4-project", config={"displayModeBar": False}), md=4),
                dbc.Col(dcc.Graph(id="p4-response", config={"displayModeBar": False}), md=4),
                dbc.Col(dcc.Graph(id="p4-sex",     config={"displayModeBar": False}), md=4),
            ]),
        ])
    ], style={"background": CARD_BG, "border": f"1px solid {BORDER}", "borderRadius": "12px", "marginBottom": "24px"}),

    html.Footer(
        "Loblaw Bio · Internal Use Only",
        style={"textAlign": "center", "color": MUTED, "padding": "24px 0", "fontSize": "0.75rem"},
    ),
], fluid=True, style={"background": DARK_BG, "minHeight": "100vh", "padding": "0"})

app.layout = layout


# ── callbacks: Part 2 ────────────────────────────────────────────────────────

def filter_freq(condition, treatment, stype, time):
    df = FREQ.copy()
    if condition != "All": df = df[df["condition"] == condition]
    if treatment != "All": df = df[df["treatment"] == treatment]
    if stype     != "All": df = df[df["sample_type"] == stype]
    if time      != "All": df = df[df["time_from_treatment_start"] == float(time)]
    return df


@app.callback(
    Output("p2-bar", "figure"),
    Output("p2-pie", "figure"),
    Output("p2-violin", "figure"),
    Output("p2-table", "children"),
    Input("p2-condition", "value"),
    Input("p2-treatment", "value"),
    Input("p2-stype", "value"),
    Input("p2-time", "value"),
)
def update_part2(cond, trt, stype, time):
    df = filter_freq(cond, trt, stype, time)
    df["pop_label"] = df["population"].map(POP_LABELS)

    # Mean % per population
    mean_pct = df.groupby("pop_label")["percentage"].mean().reset_index()
    mean_pct.columns = ["population", "mean_pct"]

    bar = go.Figure(go.Bar(
        x=mean_pct["population"], y=mean_pct["mean_pct"],
        marker=dict(color=[ACCENT, GREEN, "#38bdf8", YELLOW, RED],
                    line=dict(color=DARK_BG, width=1)),
        text=mean_pct["mean_pct"].round(2).astype(str) + "%",
        textposition="outside", textfont=dict(color=TEXT, size=11),
    ))
    bar.update_layout(**PLOTLY_LAYOUT,
                      title=dict(text="Mean Relative Frequency per Population", font=dict(color=TEXT)),
                      yaxis_title="Mean %", xaxis_title="")

    pie = go.Figure(go.Pie(
        labels=mean_pct["population"], values=mean_pct["mean_pct"],
        hole=0.45,
        marker=dict(colors=[ACCENT, GREEN, "#38bdf8", YELLOW, RED],
                    line=dict(color=DARK_BG, width=2)),
        textfont=dict(color=TEXT),
    ))
    pie.update_layout(**{**PLOTLY_LAYOUT, "margin": dict(l=10, r=10, t=40, b=10)},
                      title=dict(text="Composition", font=dict(color=TEXT)),
                      showlegend=False)

    violin = px.violin(
        df, x="pop_label", y="percentage", color="pop_label",
        box=True, points=False,
        color_discrete_sequence=[ACCENT, GREEN, "#38bdf8", YELLOW, RED],
        labels={"pop_label": "", "percentage": "Relative Frequency (%)"},
        title="Distribution of Relative Frequencies",
    )
    violin.update_layout(**PLOTLY_LAYOUT)
    violin.update_traces(meanline_visible=True)

    # Table
    tbl_df = df[["sample","population","count","total_count","percentage"]].head(200).copy()
    tbl_df["percentage"] = tbl_df["percentage"].round(3)
    tbl = dash_table.DataTable(
        data=tbl_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in tbl_df.columns],
        page_size=10,
        style_table={"overflowX": "auto"},
        style_header={"backgroundColor": DARK_BG, "color": TEXT, "fontWeight": "600", "border": f"1px solid {BORDER}"},
        style_cell={"backgroundColor": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}", "fontSize": "12px"},
        style_data_conditional=[{"if": {"row_index": "odd"}, "backgroundColor": DARK_BG}],
        sort_action="native", filter_action="native",
    )
    return bar, pie, violin, tbl


# ── callbacks: Part 3 ────────────────────────────────────────────────────────

def get_part3_data():
    return FREQ[
        (FREQ["condition"] == "melanoma") &
        (FREQ["treatment"] == "miraclib") &
        (FREQ["sample_type"] == "PBMC") &
        (FREQ["response"].isin(["yes", "no"]))
    ].copy()


@app.callback(
    Output("p3-box",   "figure"),
    Output("p3-strip", "figure"),
    Output("p3-allbox", "figure"),
    Output("p3-stats-table", "children"),
    Input("p3-pop", "value"),
)
def update_part3(pop):
    df3 = get_part3_data()
    df3["Response"] = df3["response"].map({"yes": "Responder", "no": "Non-Responder"})
    df3["pop_label"] = df3["population"].map(POP_LABELS)

    sub = df3[df3["population"] == pop].copy()

    # Box
    box = px.box(sub, x="Response", y="percentage", color="Response",
                 color_discrete_map={"Responder": GREEN, "Non-Responder": RED},
                 points="outliers",
                 labels={"percentage": "Relative Frequency (%)"},
                 title=f"{POP_LABELS[pop]} – Responders vs Non-Responders")
    box.update_layout(**PLOTLY_LAYOUT, showlegend=False)

    # Strip / jitter
    strip = px.strip(sub, x="Response", y="percentage", color="Response",
                     color_discrete_map={"Responder": GREEN, "Non-Responder": RED},
                     labels={"percentage": "Relative Frequency (%)"},
                     title=f"{POP_LABELS[pop]} – Individual Samples")
    strip.update_traces(jitter=0.35, marker_size=4, opacity=0.6)
    strip.update_layout(**PLOTLY_LAYOUT, showlegend=False)

    # All populations box
    allbox = px.box(df3, x="pop_label", y="percentage", color="Response",
                    color_discrete_map={"Responder": GREEN, "Non-Responder": RED},
                    labels={"pop_label": "", "percentage": "Relative Frequency (%)"},
                    title="All Populations – Responders vs Non-Responders")
    allbox.update_layout(**PLOTLY_LAYOUT)

    # Stats table
    rows = []
    for p in CELL_POPS:
        s = df3[df3["population"] == p]
        r  = s[s["response"] == "yes"]["percentage"].values
        nr = s[s["response"] == "no"]["percentage"].values
        if len(r) > 1 and len(nr) > 1:
            stat, pval = stats.mannwhitneyu(r, nr, alternative="two-sided")
            p_bon = min(pval * len(CELL_POPS), 1.0)
        else:
            stat, pval, p_bon = np.nan, np.nan, np.nan
        rows.append({
            "Population":      POP_LABELS[p],
            "N Responders":    len(r),
            "N Non-Responders":len(nr),
            "Median R (%)":    round(float(np.median(r)), 3) if len(r) else np.nan,
            "Median NR (%)":   round(float(np.median(nr)), 3) if len(nr) else np.nan,
            "MW Statistic":    round(float(stat), 1) if not np.isnan(stat) else "—",
            "p-value":         f"{pval:.4f}" if not np.isnan(pval) else "—",
            "p (Bonferroni)":  f"{p_bon:.4f}" if not np.isnan(p_bon) else "—",
            "Significant*":    "✓" if (not np.isnan(pval) and pval < 0.05) else "",
        })
    stats_df = pd.DataFrame(rows)

    style_cond = [
        {"if": {"filter_query": '{Significant*} = "✓"', "column_id": "Significant*"},
         "color": GREEN, "fontWeight": "bold"},
        {"if": {"row_index": "odd"}, "backgroundColor": DARK_BG},
    ]
    tbl = dash_table.DataTable(
        data=stats_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in stats_df.columns],
        style_table={"overflowX": "auto"},
        style_header={"backgroundColor": DARK_BG, "color": TEXT, "fontWeight": "600", "border": f"1px solid {BORDER}"},
        style_cell={"backgroundColor": CARD_BG, "color": TEXT, "border": f"1px solid {BORDER}", "fontSize": "12px", "textAlign": "center"},
        style_data_conditional=style_cond,
    )
    note = html.P("* p < 0.05, uncorrected Mann-Whitney U (two-sided). Bonferroni correction applied across 5 populations.",
                  style={"color": MUTED, "fontSize": "0.75rem", "marginTop": "6px"})
    return box, strip, allbox, html.Div([tbl, note])


# ── callbacks: Part 4 ────────────────────────────────────────────────────────

@app.callback(
    Output("p4-project",  "figure"),
    Output("p4-response", "figure"),
    Output("p4-sex",      "figure"),
    Input("p4-project", "id"),  # dummy trigger on load
)
def update_part4(_):
    by_proj = BASE4.groupby("project").size().reset_index(name="n")
    by_resp = BASE4[BASE4["response"].isin(["yes","no"])].groupby("response").size().reset_index(name="n")
    by_sex  = BASE4.groupby("sex").size().reset_index(name="n")

    by_resp["label"] = by_resp["response"].map({"yes": "Responder", "no": "Non-Responder"})
    by_sex["label"]  = by_sex["sex"].map({"M": "Male", "F": "Female"})

    def donut(labels, values, title, colors):
        fig = go.Figure(go.Pie(
            labels=labels, values=values, hole=0.5,
            marker=dict(colors=colors, line=dict(color=DARK_BG, width=2)),
            textfont=dict(color=TEXT),
        ))
        fig.update_layout(**{**PLOTLY_LAYOUT, "margin": dict(l=10, r=10, t=50, b=10)},
                          title=dict(text=title, font=dict(color=TEXT, size=13)))
        return fig

    proj_colors = [ACCENT, "#38bdf8", YELLOW]
    fig_proj = donut(by_proj["project"].tolist(), by_proj["n"].tolist(),
                     "Samples by Project", proj_colors[:len(by_proj)])

    fig_resp = donut(by_resp["label"].tolist(), by_resp["n"].tolist(),
                     "Responders vs Non-Responders", [GREEN, RED])

    fig_sex  = donut(by_sex["label"].tolist(), by_sex["n"].tolist(),
                     "Sex Distribution", ["#38bdf8", RED])

    return fig_proj, fig_resp, fig_sex


# ── run ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8050))
app.run(debug=False, host="0.0.0.0", port=port)
