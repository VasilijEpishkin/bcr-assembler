#!/usr/bin/env python3
"""Build AIRR QC statistics, empirical quality distributions, and combined igbrowser HTML."""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from collections import Counter
from pathlib import Path

import polars as pl

DATASET = "PRJEB30386"
RUNS = {
    "ERR3004229": {"library": "IgM", "expected_locus": "IGH"},
    "ERR3004230": {"library": "IgG", "expected_locus": "IGH"},
    "ERR3004231": {"library": "IgK", "expected_locus": "IGK"},
    "ERR3004232": {"library": "IgL", "expected_locus": "IGL"},
}

# Альтернативные конфигурации наборов данных. Выбираются по имени каталога
# dataset-root, чтобы один и тот же скрипт обслуживал и человека, и мышь.
DATASET_CONFIGS = {
    "PRJEB30386": {
        "dataset": "PRJEB30386",
        "runs": RUNS,
        "lib_colors": {
            "IgM": "#3267d6",
            "IgG": "#d65f5f",
            "IgK": "#2b8a6e",
            "IgL": "#8e63ce",
        },
    },
    "ERP003950_fastp_q30_u40": {
        "dataset": "ERP003950_fastp_q30_u40",
        "runs": {
            "ERR346596": {"library": "Mouse IGH 596", "expected_locus": "IGH"},
            "ERR346597": {"library": "Mouse IGH 597", "expected_locus": "IGH"},
            "ERR346598": {"library": "Mouse IGH 598", "expected_locus": "IGH"},
            "ERR346599": {"library": "Mouse IGH 599", "expected_locus": "IGH"},
            "ERR346600": {"library": "Mouse IGH 600", "expected_locus": "IGH"},
            "ERR346601": {"library": "Mouse IGH 601", "expected_locus": "IGH"},
        },
        "lib_colors": {
            "Mouse IGH 596": "#3267d6",
            "Mouse IGH 597": "#d65f5f",
            "Mouse IGH 598": "#2b8a6e",
            "Mouse IGH 599": "#8e63ce",
            "Mouse IGH 600": "#c7862a",
            "Mouse IGH 601": "#3aa0b8",
        },
    },
}
BOOL_TRUE = ["t", "true", "1", "yes"]

# Conservative, dataset-specific post-annotation QC cutoffs.
# Current report: minimum P01 across libraries is 83.391% for V identity
# and 88.372% for J identity. Therefore 80/85% sit below the observed
# 1% tails and avoid a strong bias against somatically hypermutated IgG.
MIN_V_IDENTITY = 80.0
MIN_J_IDENTITY = 85.0
MAX_V_SUPPORT = 1e-5
MAX_J_SUPPORT = 1e-5
MIN_VJ_SPAN = 300

QUALITY_SPECS = {
    "v_identity": {
        "label": "V identity",
        "x_label": "identity, %",
        "x_min": 0.0,
        "x_max": 100.0,
        "bin_width": 1.0,
        "threshold": MIN_V_IDENTITY,
        "threshold_label": f"cutoff {MIN_V_IDENTITY:.0f}%",
    },
    "j_identity": {
        "label": "J identity",
        "x_label": "identity, %",
        "x_min": 0.0,
        "x_max": 100.0,
        "bin_width": 1.0,
        "threshold": MIN_J_IDENTITY,
        "threshold_label": f"cutoff {MIN_J_IDENTITY:.0f}%",
    },
    "v_support_score": {
        "label": "V support",
        "x_label": "-log10(E-value); higher is better",
        "x_min": 0.0,
        "x_max": 200.0,
        "bin_width": 2.0,
        "threshold": -math.log10(MAX_V_SUPPORT),
        "threshold_label": f"E={MAX_V_SUPPORT:.0e}",
    },
    "j_support_score": {
        "label": "J support",
        "x_label": "-log10(E-value); higher is better",
        "x_min": 0.0,
        "x_max": 60.0,
        "bin_width": 1.0,
        "threshold": -math.log10(MAX_J_SUPPORT),
        "threshold_label": f"E={MAX_J_SUPPORT:.0e}",
    },
}

def repo_root_from_script() -> Path:
    here = Path(__file__).resolve().parent
    for p in (here, *here.parents):
        if (p / "scripts").is_dir():
            return p
    return Path.cwd().resolve()


def default_dataset_root() -> Path:
    return repo_root_from_script() / "results" / DATASET




def resolve_dataset_config(dataset_root: Path) -> dict:
    name = dataset_root.name
    if name not in DATASET_CONFIGS:
        raise ValueError(
            f"Unknown dataset root name {name!r}; expected one of {sorted(DATASET_CONFIGS)}"
        )
    return DATASET_CONFIGS[name]

def wilson(k: int, n: int) -> tuple[float, float]:
    if not n:
        return 0.0, 0.0
    z = 1.959963984540054
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return 100 * (centre - half), 100 * (centre + half)


def scan_airr(path: Path) -> pl.LazyFrame:
    with path.open() as handle:
        cols = handle.readline().rstrip("\n").split("\t")
    schema = {c: pl.String for c in cols}
    return pl.scan_csv(path, separator="\t", schema=schema, has_header=True, null_values=None)


def is_true(name: str) -> pl.Expr:
    return pl.col(name).fill_null("").str.to_lowercase().is_in(BOOL_TRUE)


def nonempty(name: str) -> pl.Expr:
    return pl.col(name).fill_null("").str.len_chars().gt(0)


def number(name: str) -> pl.Expr:
    return pl.col(name).cast(pl.Float64, strict=False)


def first_gene(name: str) -> pl.Expr:
    return (
        pl.col(name).fill_null("")
        .str.split(",").list.first()
        .str.split("*").list.first()
        .str.strip_chars()
    )


def quality_histograms(lf: pl.LazyFrame, run: str) -> list[dict]:
    """Build normalized empirical histograms from all annotated rows of one library."""
    df = lf.select(
        number("v_identity").alias("v_identity"),
        number("j_identity").alias("j_identity"),
        number("v_support").alias("v_support"),
        number("j_support").alias("j_support"),
    ).collect(engine="streaming")

    counts: dict[str, Counter] = {k: Counter() for k in QUALITY_SPECS}
    valid = Counter()

    def add(metric: str, value: float | None) -> None:
        if value is None or not math.isfinite(value):
            return
        spec = QUALITY_SPECS[metric]
        x_min, x_max, width = spec["x_min"], spec["x_max"], spec["bin_width"]
        # First/last bins absorb values outside the plotting range.
        value = min(max(float(value), x_min), math.nextafter(x_max, x_min))
        low = x_min + math.floor((value - x_min) / width) * width
        low = round(low, 10)
        counts[metric][low] += 1
        valid[metric] += 1

    for vi, ji, vs, js in df.iter_rows():
        add("v_identity", vi)
        add("j_identity", ji)
        if vs is not None and math.isfinite(vs) and vs >= 0:
            add("v_support_score", 300.0 if vs == 0 else -math.log10(vs))
        if js is not None and math.isfinite(js) and js >= 0:
            add("j_support_score", 300.0 if js == 0 else -math.log10(js))

    rows: list[dict] = []
    for metric, counter in counts.items():
        width = QUALITY_SPECS[metric]["bin_width"]
        denom = valid[metric]
        for low in sorted(counter):
            n = int(counter[low])
            rows.append({
                "run": run,
                "library": RUNS[run]["library"],
                "metric": metric,
                "bin_low": low,
                "bin_high": min(low + width, QUALITY_SPECS[metric]["x_max"]),
                "count": n,
                "pct_of_valid": 100.0 * n / denom if denom else 0.0,
                "valid_n": int(denom),
            })
    return rows


def collect_metrics(path: Path, run: str) -> tuple[dict, list[dict], list[dict], list[dict], list[dict]]:
    lf = scan_airr(path)
    expected = RUNS[run]["expected_locus"]

    vj_span_ok = (number("j_sequence_end") - number("v_sequence_start") + 1).ge(MIN_VJ_SPAN).fill_null(False)
    v_start_ok = number("v_germline_start").eq(1).fill_null(False)
    complete_ok = is_true("complete_vdj")
    productive_ok = is_true("productive")
    no_stop = ~is_true("stop_codon")
    locus_ok = pl.col("locus").fill_null("").eq(expected)
    v_id_ok = number("v_identity").ge(MIN_V_IDENTITY).fill_null(False)
    j_id_ok = number("j_identity").ge(MIN_J_IDENTITY).fill_null(False)
    v_e_ok = number("v_support").le(MAX_V_SUPPORT).fill_null(False)
    j_e_ok = number("j_support").le(MAX_J_SUPPORT).fill_null(False)
    proposed_ok = vj_span_ok & v_start_ok & complete_ok & productive_ok & no_stop & locus_ok & v_id_ok & j_id_ok & v_e_ok & j_e_ok

    row = lf.select(
        pl.len().alias("total"),
        productive_ok.sum().alias("productive"),
        complete_ok.sum().alias("complete_vdj"),
        is_true("stop_codon").sum().alias("stop_codon"),
        no_stop.sum().alias("no_stop"),
        is_true("vj_in_frame").sum().alias("vj_in_frame"),
        nonempty("v_call").sum().alias("v_call_present"),
        nonempty("d_call").sum().alias("d_call_present"),
        nonempty("j_call").sum().alias("j_call_present"),
        number("v_identity").is_not_null().sum().alias("v_identity_present"),
        number("j_identity").is_not_null().sum().alias("j_identity_present"),
        number("v_support").is_not_null().sum().alias("v_support_present"),
        number("j_support").is_not_null().sum().alias("j_support_present"),
        v_e_ok.sum().alias("v_support_le_1e5"),
        j_e_ok.sum().alias("j_support_le_1e5"),
        v_id_ok.sum().alias("v_identity_ge_80"),
        j_id_ok.sum().alias("j_identity_ge_85"),
        vj_span_ok.sum().alias("vj_span_ge_300"),
        v_start_ok.sum().alias("v_germline_start_eq1"),
        locus_ok.sum().alias("expected_locus_count"),
        proposed_ok.sum().alias("proposed_filter_pass"),
        nonempty("cdr3_aa").sum().alias("cdr3_aa_present"),
        nonempty("fwr1").sum().alias("fwr1_present"),
        nonempty("cdr1").sum().alias("cdr1_present"),
        number("v_germline_start").gt(1).sum().alias("v_germline_start_gt1"),
        number("v_germline_start").gt(15).sum().alias("v_germline_start_gt15"),
        pl.col("sequence_id").str.contains(r"(?:^|\|)BARCODE=$").sum().alias("empty_barcode"),
        number("v_identity").median().alias("v_identity_median"),
        number("v_identity").quantile(0.01).alias("v_identity_p01"),
        number("v_identity").quantile(0.05).alias("v_identity_p05"),
        number("v_identity").quantile(0.95).alias("v_identity_p95"),
        number("v_identity").quantile(0.99).alias("v_identity_p99"),
        number("j_identity").median().alias("j_identity_median"),
        number("j_identity").quantile(0.01).alias("j_identity_p01"),
        number("j_identity").quantile(0.05).alias("j_identity_p05"),
        number("j_identity").quantile(0.95).alias("j_identity_p95"),
        number("j_identity").quantile(0.99).alias("j_identity_p99"),
        number("v_support").median().alias("v_support_median"),
        number("v_support").quantile(0.95).alias("v_support_p95"),
        number("v_support").quantile(0.99).alias("v_support_p99"),
        number("j_support").median().alias("j_support_median"),
        number("j_support").quantile(0.95).alias("j_support_p95"),
        number("j_support").quantile(0.99).alias("j_support_p99"),
        pl.col("cdr3_aa").filter(nonempty("cdr3_aa")).n_unique().alias("unique_cdr3_aa"),
        pl.col("locus").filter(nonempty("locus")).n_unique().alias("n_loci"),
        (pl.col("locus").fill_null("") != expected).sum().alias("unexpected_locus"),
    ).collect(engine="streaming").row(0, named=True)

    row.update(run=run, library=RUNS[run]["library"], expected_locus=expected)
    total = int(row["total"])
    count_fields = [
        "productive", "complete_vdj", "stop_codon", "no_stop", "vj_in_frame",
        "v_call_present", "d_call_present", "j_call_present", "v_identity_present",
        "j_identity_present", "v_support_present", "j_support_present",
        "v_support_le_1e5", "j_support_le_1e5", "v_identity_ge_80",
        "j_identity_ge_85", "vj_span_ge_300", "v_germline_start_eq1",
        "expected_locus_count", "proposed_filter_pass", "cdr3_aa_present",
        "fwr1_present", "cdr1_present", "v_germline_start_gt1",
        "v_germline_start_gt15", "empty_barcode", "unexpected_locus",
    ]
    for field in count_fields:
        k = int(row[field])
        lo, hi = wilson(k, total)
        row[field + "_pct"] = 100 * k / total if total else 0.0
        row[field + "_ci95_low"] = lo
        row[field + "_ci95_high"] = hi

    def top_calls(field: str, n: int = 15) -> list[dict]:
        x = (
            lf.with_columns(first_gene(field).alias("gene"))
            .filter(pl.col("gene") != "")
            .group_by("gene").len()
            .sort("len", descending=True).head(n)
            .collect(engine="streaming")
        )
        return [
            {
                "run": run,
                "library": RUNS[run]["library"],
                "gene": r["gene"],
                "count": int(r["len"]),
                "pct": 100 * int(r["len"]) / total,
            }
            for r in x.iter_rows(named=True)
        ]

    hist_rows = quality_histograms(lf, run)
    return row, top_calls("v_call"), top_calls("d_call"), top_calls("j_call"), hist_rows


def select_representatives(path: Path, run: str, n: int = 5) -> tuple[list[str], list[dict]]:
    required = [
        "sequence", "sequence_id", "v_sequence_start", "v_sequence_end",
        "j_sequence_start", "j_sequence_end", "cdr3_start", "cdr3_end",
        "fwr4", "fwr4_start", "fwr4_end",
    ]
    selected = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        for row in reader:
            if all(row.get(c, "") not in ("", "NA", "NaN") for c in required) and row.get("productive", "").lower() in BOOL_TRUE:
                row["report_source_run"] = run
                row["report_library"] = RUNS[run]["library"]
                selected.append(row)
                if len(selected) == n:
                    break
    if len(selected) < n:
        raise RuntimeError(f"{run}: only {len(selected)} rows valid for igbrowser")
    return fields + ["report_source_run", "report_library"], selected


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise RuntimeError(f"No rows for {path}")
    fields = fields or list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def stats(dataset_root: Path) -> None:
    global DATASET, RUNS, LIB_COLORS
    DATASET, RUNS, LIB_COLORS = (lambda c: (c["dataset"], c["runs"], c["lib_colors"]))(resolve_dataset_config(dataset_root))
    annotation = dataset_root / "annotation" / "igblast"
    out = dataset_root / "annotation" / "report"
    if not annotation.is_dir():
        raise FileNotFoundError(f"Annotation directory not found: {annotation}")
    out.mkdir(parents=True, exist_ok=True)

    summaries, top_v, top_d, top_j, hist, reps = [], [], [], [], [], []
    rep_fields = None
    for run in RUNS:
        path = annotation / f"{run}.airr.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"[stats] {run}: {path}", flush=True)
        s, v, d, j, h = collect_metrics(path, run)
        summaries.append(s)
        top_v.extend(v)
        top_d.extend(d)
        top_j.extend(j)
        hist.extend(h)
        fields, selected = select_representatives(path, run)
        rep_fields = rep_fields or fields
        reps.extend(selected)
        print(
            f"  rows={s['total']:,} productive={s['productive_pct']:.2f}% "
            f"complete={s['complete_vdj_pct']:.2f}% proposed-pass={s['proposed_filter_pass_pct']:.2f}%",
            flush=True,
        )

    write_tsv(out / "annotation_summary_by_library.tsv", summaries)
    write_tsv(out / "top_v_calls.tsv", top_v)
    write_tsv(out / "top_d_calls.tsv", top_d)
    write_tsv(out / "top_j_calls.tsv", top_j)
    write_tsv(out / "annotation_quality_distributions.tsv", hist)
    write_tsv(out / "igbrowser_representatives.airr.tsv", reps, rep_fields)
    (out / "annotation_summary.json").write_text(
        json.dumps({"dataset": DATASET, "libraries": summaries}, indent=2) + "\n"
    )
    print(f"[stats] wrote {out}", flush=True)


def fmt_int(x) -> str:
    return f"{int(x):,}".replace(",", " ")


def fmt_pct(x) -> str:
    return f"{float(x):.2f}%"


def sci(x) -> str:
    try:
        return f"{float(x):.2e}"
    except (TypeError, ValueError):
        return "NA"


def metric_bar(label: str, rows: list[dict], field: str, color: str) -> str:
    parts = [f"<h3>{html.escape(label)}</h3><div class='bars'>"]
    for r in rows:
        v = float(r[field])
        parts.append(
            f"<div class='barrow'><span>{r['library']}</span><div class='track'>"
            f"<div class='fill' style='width:{min(v,100):.3f}%;background:{color}'></div>"
            f"</div><b>{v:.2f}%</b></div>"
        )
    return "".join(parts) + "</div>"


def top_table(rows: list[dict], run: str) -> str:
    x = [r for r in rows if r["run"] == run][:10]
    return (
        "<table><tr><th>Gene</th><th>Reads</th><th>%</th></tr>"
        + "".join(
            f"<tr><td>{html.escape(r['gene'])}</td><td>{fmt_int(r['count'])}</td>"
            f"<td>{float(r['pct']):.2f}</td></tr>" for r in x
        )
        + "</table>"
    )


def quality_table(rows: list[dict]) -> str:
    body = []
    for r in rows:
        body.append(
            f"<tr><td>{r['library']}</td><td>{r['run']}</td>"
            f"<td>{float(r['v_identity_p01'] or 0):.3f}</td><td>{float(r['v_identity_p05'] or 0):.3f}</td>"
            f"<td>{float(r['v_identity_median'] or 0):.3f}</td><td>{float(r['v_identity_p95'] or 0):.3f}</td>"
            f"<td>{float(r['v_identity_p99'] or 0):.3f}</td>"
            f"<td>{float(r['j_identity_p01'] or 0):.3f}</td><td>{float(r['j_identity_p05'] or 0):.3f}</td>"
            f"<td>{float(r['j_identity_median'] or 0):.3f}</td><td>{float(r['j_identity_p95'] or 0):.3f}</td>"
            f"<td>{float(r['j_identity_p99'] or 0):.3f}</td>"
            f"<td>{sci(r['v_support_median'])}</td><td>{sci(r['v_support_p95'])}</td><td>{sci(r['v_support_p99'])}</td>"
            f"<td>{sci(r['j_support_median'])}</td><td>{sci(r['j_support_p95'])}</td><td>{sci(r['j_support_p99'])}</td></tr>"
        )
    return (
        "<div class='table-scroll'><table><tr><th>Library</th><th>Run</th>"
        "<th>V id P01</th><th>V id P05</th><th>V id median</th><th>V id P95</th><th>V id P99</th>"
        "<th>J id P01</th><th>J id P05</th><th>J id median</th><th>J id P95</th><th>J id P99</th>"
        "<th>V E median</th><th>V E P95</th><th>V E P99</th>"
        "<th>J E median</th><th>J E P95</th><th>J E P99</th></tr>"
        + "".join(body)
        + "</table></div>"
    )


def filter_preview_table(rows: list[dict]) -> str:
    body = []
    for r in rows:
        body.append(
            f"<tr><td>{r['library']}</td>"
            f"<td>{fmt_pct(r['vj_span_ge_300_pct'])}</td>"
            f"<td>{fmt_pct(r['v_germline_start_eq1_pct'])}</td>"
            f"<td>{fmt_pct(r['complete_vdj_pct'])}</td>"
            f"<td>{fmt_pct(r['productive_pct'])}</td>"
            f"<td>{fmt_pct(r['no_stop_pct'])}</td>"
            f"<td>{fmt_pct(r['expected_locus_count_pct'])}</td>"
            f"<td>{fmt_pct(r['v_identity_ge_80_pct'])}</td>"
            f"<td>{fmt_pct(r['j_identity_ge_85_pct'])}</td>"
            f"<td>{fmt_pct(r['v_support_le_1e5_pct'])}</td>"
            f"<td>{fmt_pct(r['j_support_le_1e5_pct'])}</td>"
            f"<td><b>{fmt_pct(r['proposed_filter_pass_pct'])}</b></td></tr>"
        )
    return (
        "<div class='table-scroll'><table><tr><th>Library</th><th>V–J span ≥300</th>"
        "<th>V starts at germline 1</th><th>Complete VDJ / J end</th><th>Productive</th>"
        "<th>No stop</th><th>Expected locus</th><th>V id ≥80%</th><th>J id ≥85%</th>"
        "<th>V E≤1e−5</th><th>J E≤1e−5</th><th>All criteria</th></tr>"
        + "".join(body)
        + "</table></div>"
    )


def distribution_svg(hist_rows: list[dict], metric: str) -> str:
    spec = QUALITY_SPECS[metric]
    width, height = 760, 300
    left, right, top, bottom = 62, 18, 24, 58
    plot_w, plot_h = width - left - right, height - top - bottom
    x_min, x_max, bw = spec["x_min"], spec["x_max"], spec["bin_width"]
    y_min_log, y_max_log = -4.0, 2.0  # 0.0001% .. 100%, logarithmic frequency axis

    def xpix(x: float) -> float:
        return left + (x - x_min) / (x_max - x_min) * plot_w

    def ypix(pct: float) -> float:
        lv = math.log10(max(pct, 1e-4))
        lv = min(max(lv, y_min_log), y_max_log)
        return top + (y_max_log - lv) / (y_max_log - y_min_log) * plot_h

    # Grid and y ticks.
    parts = [f"<svg viewBox='0 0 {width} {height}' role='img' aria-label='{html.escape(spec['label'])} distribution'>"]
    for val in (1e-4, 1e-2, 1.0, 100.0):
        y = ypix(val)
        parts.append(f"<line x1='{left}' y1='{y:.1f}' x2='{width-right}' y2='{y:.1f}' stroke='#e1e6ef'/>")
        parts.append(f"<text x='{left-8}' y='{y+4:.1f}' text-anchor='end' font-size='10'>{val:g}%</text>")

    # x ticks.
    nticks = 5
    for i in range(nticks + 1):
        xval = x_min + (x_max - x_min) * i / nticks
        x = xpix(xval)
        parts.append(f"<line x1='{x:.1f}' y1='{top+plot_h}' x2='{x:.1f}' y2='{top+plot_h+5}' stroke='#8993a5'/>")
        parts.append(f"<text x='{x:.1f}' y='{height-35}' text-anchor='middle' font-size='10'>{xval:g}</text>")

    # Threshold line.
    thr = float(spec["threshold"])
    tx = xpix(thr)
    parts.append(f"<line x1='{tx:.1f}' y1='{top}' x2='{tx:.1f}' y2='{top+plot_h}' stroke='#111' stroke-dasharray='5 4' stroke-width='1.5'/>")
    parts.append(f"<text x='{min(tx+5,width-150):.1f}' y='{top+12}' font-size='10'>{html.escape(spec['threshold_label'])}</text>")

    # One normalized empirical series per library.
    legend_x = left + 8
    for idx, (run, meta) in enumerate(RUNS.items()):
        library = meta["library"]
        color = LIB_COLORS[library]
        row_map = {float(r["bin_low"]): float(r["pct_of_valid"]) for r in hist_rows if r["run"] == run and r["metric"] == metric}
        pts = []
        x = x_min
        while x < x_max - 1e-9:
            pct = row_map.get(round(x, 10), 0.0)
            pts.append(f"{xpix(x + bw/2):.1f},{ypix(pct):.1f}")
            x += bw
        parts.append(f"<polyline points='{' '.join(pts)}' fill='none' stroke='{color}' stroke-width='2'/>")
        ly = top + 18 + idx * 16
        parts.append(f"<line x1='{legend_x}' y1='{ly}' x2='{legend_x+18}' y2='{ly}' stroke='{color}' stroke-width='3'/>")
        parts.append(f"<text x='{legend_x+24}' y='{ly+4}' font-size='10'>{library}</text>")

    parts.append(f"<line x1='{left}' y1='{top+plot_h}' x2='{width-right}' y2='{top+plot_h}' stroke='#8993a5'/>")
    parts.append(f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top+plot_h}' stroke='#8993a5'/>")
    parts.append(f"<text x='{left+plot_w/2:.1f}' y='{height-8}' text-anchor='middle' font-size='11'>{html.escape(spec['x_label'])}</text>")
    parts.append(f"<text x='14' y='{top+plot_h/2:.1f}' transform='rotate(-90 14 {top+plot_h/2:.1f})' text-anchor='middle' font-size='11'>reads per bin, % (log scale)</text>")
    parts.append("</svg>")
    return "".join(parts)


def extract_body(text: str) -> str:
    m = re.search(r"<body[^>]*>(.*)</body>", text, flags=re.I | re.S)
    return m.group(1) if m else text


def assemble(dataset_root: Path) -> Path:
    global DATASET, RUNS, LIB_COLORS
    DATASET, RUNS, LIB_COLORS = (lambda c: (c["dataset"], c["runs"], c["lib_colors"]))(resolve_dataset_config(dataset_root))
    out = dataset_root / "annotation" / "report"
    summaries = json.loads((out / "annotation_summary.json").read_text())["libraries"]
    with (out / "top_v_calls.tsv").open() as h:
        top_v = list(csv.DictReader(h, delimiter="\t"))
    with (out / "top_d_calls.tsv").open() as h:
        top_d = list(csv.DictReader(h, delimiter="\t"))
    with (out / "top_j_calls.tsv").open() as h:
        top_j = list(csv.DictReader(h, delimiter="\t"))
    with (out / "annotation_quality_distributions.tsv").open() as h:
        hist_rows = list(csv.DictReader(h, delimiter="\t"))

    native_path = out / "igbrowser_representatives_native.html"
    if not native_path.exists():
        raise FileNotFoundError(f"Native igbrowser HTML missing: {native_path}")

    total = sum(int(r["total"]) for r in summaries)
    cards = "".join(
        f"<div class='card'><h3>{r['library']} · {r['expected_locus']}</h3>"
        f"<div class='big'>{fmt_int(r['total'])}</div><p>productive {fmt_pct(r['productive_pct'])}"
        f"<br>complete VDJ {fmt_pct(r['complete_vdj_pct'])}"
        f"<br>proposed filter pass {fmt_pct(r['proposed_filter_pass_pct'])}</p></div>"
        for r in summaries
    )
    overview_rows = "".join(
        f"<tr><td>{r['run']}</td><td>{r['library']}</td><td>{r['expected_locus']}</td>"
        f"<td>{fmt_int(r['total'])}</td><td>{fmt_pct(r['productive_pct'])}</td>"
        f"<td>{fmt_pct(r['complete_vdj_pct'])}</td><td>{fmt_pct(r['stop_codon_pct'])}</td>"
        f"<td>{float(r['v_identity_median'] or 0):.3f}</td>"
        f"<td>{fmt_pct(r['v_support_le_1e5_pct'])}</td><td>{fmt_pct(r['j_support_le_1e5_pct'])}</td>"
        f"<td>{fmt_pct(r['v_germline_start_gt1_pct'])}</td><td>{fmt_int(r['unexpected_locus'])}</td>"
        f"<td>{fmt_pct(r['empty_barcode_pct'])}</td></tr>" for r in summaries
    )
    calls = "".join(
        f"<section><h3>{r['library']} ({r['run']})</h3><div class='threecol'>"
        f"<div><h4>Top V</h4>{top_table(top_v,r['run'])}</div>"
        f"<div><h4>Top D</h4>{top_table(top_d,r['run']) if r['expected_locus']=='IGH' else '<p>Not applicable for light chain</p>'}</div>"
        f"<div><h4>Top J</h4>{top_table(top_j,r['run'])}</div></div></section>" for r in summaries
    )
    dist = "".join(
        f"<div class='dist-card'><h3>{QUALITY_SPECS[m]['label']}</h3>{distribution_svg(hist_rows,m)}</div>"
        for m in ("v_identity", "j_identity", "v_support_score", "j_support_score")
    )
    native = extract_body(native_path.read_text())

    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>{DATASET} AIRR annotation report</title><style>
body{{font-family:Inter,system-ui,sans-serif;max-width:1500px;margin:0 auto;padding:28px;color:#172033;background:#f5f7fb}}h1,h2,h3{{color:#13213c}}section{{background:white;padding:20px;margin:18px 0;border-radius:12px;box-shadow:0 2px 10px #18304a12}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}}.card{{background:white;padding:16px;border-radius:12px;border-top:5px solid #3267d6}}.big{{font-size:28px;font-weight:750}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{border-bottom:1px solid #dce2ec;padding:7px;text-align:right;white-space:nowrap}}th:first-child,td:first-child{{text-align:left}}.table-scroll{{overflow-x:auto}}.twocol{{display:grid;grid-template-columns:1fr 1fr;gap:24px}}.threecol{{display:grid;grid-template-columns:1fr 1fr 1fr;gap:18px}}.dist-grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}.dist-card{{border:1px solid #e1e6ef;border-radius:10px;padding:12px}}svg{{width:100%;height:auto}}.barrow{{display:grid;grid-template-columns:55px 1fr 75px;align-items:center;gap:10px;margin:8px 0}}.track{{height:18px;background:#e7ebf3;border-radius:9px;overflow:hidden}}.fill{{height:100%}}.note{{background:#fff4cf;border-left:5px solid #e5ae20;padding:12px}}.native{{overflow-x:auto;background:#fff}}code{{background:#eef2f8;padding:2px 5px}}@media(max-width:900px){{.cards{{grid-template-columns:1fr 1fr}}.twocol,.threecol,.dist-grid{{grid-template-columns:1fr}}}}
</style></head><body><h1>{DATASET} · IgBLAST AIRR annotation</h1><p>Единый отчёт по всем библиотекам набора. Всего <b>{fmt_int(total)}</b> аннотированных merged sequences.</p><div class='cards'>{cards}</div>
<section><h2>Основные показатели</h2><div class='table-scroll'><table><tr><th>Run</th><th>Library</th><th>Locus</th><th>Reads</th><th>Productive</th><th>Complete VDJ</th><th>Stop codon</th><th>Median V identity</th><th>V e≤1e−5</th><th>J e≤1e−5</th><th>V germline start &gt;1</th><th>Unexpected locus</th><th>Empty UMI</th></tr>{overview_rows}</table></div></section>
<section><h2>Качество V/J-аннотации</h2><p>Identity показан в процентах. Для support используется E-value IgBLAST; на графиках он преобразован в <code>-log10(E-value)</code>, поэтому правее = более уверенное назначение. Частота по Y логарифмическая, чтобы был виден слабый хвост, а не только основной пик.</p>{quality_table(summaries)}<div class='dist-grid'>{dist}</div><p class='note'>Пороговые линии — консервативные dataset-specific QC-кандидаты: V identity ≥ {MIN_V_IDENTITY:.0f}%, J identity ≥ {MIN_J_IDENTITY:.0f}%, V/J E-value ≤ {MAX_V_SUPPORT:.0e}. Они не являются универсальными биологическими порогами.</p></section>
<section><h2>Предпросмотр post-annotation filter</h2><p>Это расчёт по AIRR-таблицам до записи нового FASTQ. V–J span = <code>j_sequence_end - v_sequence_start + 1</code>. <code>complete_vdj=True</code> используется для требования полного конца J; дополнительное <code>v_germline_start=1</code> явно требует начало germline V.</p>{filter_preview_table(summaries)}</section>
<section><h2>Annotation completeness</h2>{metric_bar('Productive',summaries,'productive_pct','#2b8a6e')}{metric_bar('Complete VDJ',summaries,'complete_vdj_pct','#3267d6')}{metric_bar('FWR1 present',summaries,'fwr1_present_pct','#8e63ce')}{metric_bar('V germline start > 1',summaries,'v_germline_start_gt1_pct','#d67b32')}</section>
<section><h2>V/J usage</h2>{calls}</section>
<section><h2>Native igblastr::igbrowser() — representative sequences</h2><p>По 5 reproducibly selected productive rows с полными координатами из каждой библиотеки. Это per-sequence browser; статистика выше рассчитана по всем {fmt_int(total)} AIRR rows.</p><div class='native'>{native}</div></section>
</body></html>"""

    final = out / f"{DATASET}_annotation_igbrowser_report.html"
    final.write_text(report)
    manifest = {
        "dataset": DATASET,
        "airr_rows": total,
        "libraries": {r["run"]: int(r["total"]) for r in summaries},
        "quality_cutoffs": {
            "min_v_identity": MIN_V_IDENTITY,
            "min_j_identity": MIN_J_IDENTITY,
            "max_v_support": MAX_V_SUPPORT,
            "max_j_support": MAX_J_SUPPORT,
            "min_vj_span": MIN_VJ_SPAN,
        },
        "native_igbrowser": str(native_path),
        "report": str(final),
    }
    (out / "report_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(final)
    return final


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="all", choices=["stats", "assemble", "all"])
    ap.add_argument("--dataset-root", type=Path, default=default_dataset_root())
    args = ap.parse_args()
    dataset_root = args.dataset_root.expanduser().resolve()
    print(f"dataset_root={dataset_root}", flush=True)
    if args.mode in {"stats", "all"}:
        stats(dataset_root)
    if args.mode in {"assemble", "all"}:
        assemble(dataset_root)


if __name__ == "__main__":
    main()
