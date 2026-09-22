"""Completion and summary checks for pRESTO AssemblePairs outputs."""

from __future__ import annotations

import gzip
from pathlib import Path
import re


def parse_presto_summary(text: str) -> dict[str, int]:
    values: dict[str, int] = {}
    for key in ("PAIRS", "PASS", "FAIL"):
        matches = re.findall(rf"^\s*{key}>\s*([0-9,]+)\s*$", text, flags=re.MULTILINE)
        if not matches:
            raise ValueError(f"Missing terminal {key}> summary")
        values[key.lower()] = int(matches[-1].replace(",", ""))
    if values["pass"] + values["fail"] != values["pairs"]:
        raise ValueError(f"Inconsistent pRESTO summary: {values}")
    return values


def count_fastq_records(path: str | Path) -> int:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    lines = 0
    with opener(path, "rt") as handle:
        for _ in handle:
            lines += 1
    if lines % 4:
        raise ValueError(f"FASTQ line count is not divisible by four: {path}")
    return lines // 4


def is_complete_presto_output(pass_fastq: str | Path, stdout_log: str | Path) -> bool:
    pass_fastq = Path(pass_fastq)
    stdout_log = Path(stdout_log)
    if not pass_fastq.is_file() or not stdout_log.is_file():
        return False
    try:
        text = stdout_log.read_text(errors="replace")
        if "END> AssemblePairs" not in text:
            return False
        summary = parse_presto_summary(text)
        return count_fastq_records(pass_fastq) == summary["pass"]
    except (OSError, EOFError, ValueError):
        return False
