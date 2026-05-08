"""
load_data.py
------------
Initialises the SQLite database and loads cell-count.csv.

Run:
    python load_data.py

Schema design
─────────────
Three tables:

  subjects(subject_id PK, project, condition, age, sex)
    – One row per biological subject.  Demographic attributes belong here,
      not repeated on every sample row.  At scale (thousands of subjects /
      hundreds of projects) this table can carry FK→projects if projects grow
      into their own entity with metadata.

  samples(sample_id PK, subject_id FK, sample_type, treatment, response,
          time_from_treatment_start)
    – One row per biological sample.  Treatment and response live here because
      in a real trial a subject may receive different treatments across arms or
      have per-sample response calls.  sample_type (PBMC / WB) is a sample
      attribute, not a subject attribute.

  cell_counts(id PK, sample_id FK, population, count)
    – Long / tidy format: one row per (sample, cell population).  This makes
      GROUP BY population, aggregate across populations, and adding new
      populations trivial without schema changes.

Scaling rationale
──────────────────
• Normalisation avoids repeating subject demographics on every sample row
  (saves storage, prevents update anomalies).
• Long format for cell_counts means adding a new population (e.g. dendritic
  cells) requires no DDL change – just new rows.
• Indexes on sample_id (cell_counts), subject_id (samples), and the most
  common filter columns (condition, treatment, response, time) keep analytic
  queries fast at millions of rows.
• A future projects table would add project-level metadata (PI, protocol,
  start date) with a FK from subjects.project.
"""

import csv
import os
import sqlite3

DB_PATH = "immune_trial.db"
CSV_PATH = "cell-count.csv"
CELL_POPS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


# ── schema ──────────────────────────────────────────────────────────────────

DDL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project     TEXT NOT NULL,
    condition   TEXT,
    age         REAL,
    sex         TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                   TEXT PRIMARY KEY,
    subject_id                  TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                 TEXT,
    treatment                   TEXT,
    response                    TEXT,
    time_from_treatment_start   INTEGER
);

CREATE TABLE IF NOT EXISTS cell_counts (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id   TEXT    NOT NULL REFERENCES samples(sample_id),
    population  TEXT    NOT NULL,
    count       REAL    NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_cc_sample   ON cell_counts(sample_id);
CREATE INDEX IF NOT EXISTS idx_cc_pop      ON cell_counts(population);
CREATE INDEX IF NOT EXISTS idx_smp_subject ON samples(subject_id);
CREATE INDEX IF NOT EXISTS idx_smp_filter  ON samples(treatment, response, time_from_treatment_start);
CREATE INDEX IF NOT EXISTS idx_sub_project ON subjects(project, condition);
"""


# ── helpers ─────────────────────────────────────────────────────────────────

def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(DDL)
    conn.commit()


def load_csv(conn: sqlite3.Connection, csv_path: str) -> None:
    subjects_seen: set = set()
    samples_seen: set = set()

    sub_rows, smp_rows, cc_rows = [], [], []

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            subj = row["subject"]
            samp = row["sample"]

            # subjects (dedup)
            if subj not in subjects_seen:
                subjects_seen.add(subj)
                sub_rows.append((
                    subj,
                    row["project"],
                    row["condition"],
                    float(row["age"]) if row["age"] else None,
                    row["sex"],
                ))

            # samples (dedup)
            if samp not in samples_seen:
                samples_seen.add(samp)
                smp_rows.append((
                    samp,
                    subj,
                    row["sample_type"],
                    row["treatment"],
                    row["response"] if row["response"] else None,
                    int(float(row["time_from_treatment_start"])) if row["time_from_treatment_start"] else 0,
                ))

            # cell counts (one row per population)
            for pop in CELL_POPS:
                cc_rows.append((samp, pop, float(row[pop]) if row[pop] else 0.0))

    cur = conn.cursor()
    cur.executemany(
        "INSERT OR IGNORE INTO subjects VALUES (?,?,?,?,?)", sub_rows
    )
    cur.executemany(
        "INSERT OR IGNORE INTO samples VALUES (?,?,?,?,?,?)", smp_rows
    )
    cur.executemany(
        "INSERT INTO cell_counts(sample_id, population, count) VALUES (?,?,?)", cc_rows
    )
    conn.commit()
    print(f"  Loaded {len(sub_rows)} subjects, {len(smp_rows)} samples, "
          f"{len(cc_rows)} cell-count rows.")


# ── entry point ─────────────────────────────────────────────────────────────

def main() -> None:
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        print(f"Removed existing {DB_PATH}")

    conn = sqlite3.connect(DB_PATH)
    try:
        print("Initialising schema …")
        init_db(conn)
        print(f"Loading {CSV_PATH} …")
        load_csv(conn, CSV_PATH)
        print(f"Database ready: {DB_PATH}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
