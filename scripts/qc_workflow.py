"""Input resolution and execution helpers for the canonical QC notebook."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import Iterable


STAGES = {"raw", "trimmed", "pr_trimmed", "merged"}
VIEWS = {"primary", "branch_diagnostics"}


@dataclass(frozen=True)
class QCInput:
    path: Path
    sample: str
    mate: str | None
    technical_branch: str | None
    role: str


@dataclass(frozen=True)
class QCOutputPaths:
    stage: Path
    fastqc: Path
    multiqc: Path
    report: Path
    input_manifest: Path
    run_manifest: Path


def _stage_dir(volume: Path, dataset: str, stage: str, branch: str | None) -> Path:
    results_root = volume / "results" / dataset
    is_branch_dataset = branch is not None or (results_root / "branches").is_dir()
    if stage == "raw":
        return results_root / "raw"
    if is_branch_dataset:
        if not branch:
            raise ValueError(f"branch is required for {dataset} stage {stage}")
        return results_root / "branches" / branch / stage
    if branch:
        raise ValueError(f"branch is not supported for conventional dataset {dataset}")
    return results_root / stage


def _mate(name: str) -> str | None:
    if re.search(r"(?:_|-)(?:1|R1)(?:[_.-]|$)", name, re.IGNORECASE):
        return "R1"
    if re.search(r"(?:_|-)(?:2|R2)(?:[_.-]|$)", name, re.IGNORECASE):
        return "R2"
    return None


def _sample(name: str) -> str:
    return re.split(r"(?:_|-)(?:R?[12])(?:[_.-]|$)", name, maxsplit=1, flags=re.IGNORECASE)[0]


def _items(paths: Iterable[Path], role: str, technical_branch: str | None = None) -> list[QCInput]:
    return [
        QCInput(path=p, sample=_sample(p.name), mate=_mate(p.name), technical_branch=technical_branch, role=role)
        for p in sorted(paths)
    ]


def resolve_qc_inputs(
    volume: str | Path,
    dataset: str,
    stage: str,
    branch: str | None = None,
    view: str = "primary",
) -> list[QCInput]:
    """Resolve file-level QC inputs without mixing aggregates and components."""
    volume = Path(volume)
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {sorted(STAGES)}")
    if view not in VIEWS:
        raise ValueError(f"view must be one of {sorted(VIEWS)}")

    results_root = volume / "results" / dataset
    is_branch_dataset = (results_root / "branches").is_dir()
    stage_dir = _stage_dir(volume, dataset, stage, branch)

    if stage == "raw":
        if view != "primary":
            raise ValueError("branch_diagnostics view is only valid for branch-aware pr_trimmed/merged stages")
        source = volume / "raw" / dataset
        paths = sorted(source.glob("*.fastq.gz")) + sorted(source.glob("*.fastq"))
        items = _items(paths, "raw_mate")
    elif is_branch_dataset and stage == "trimmed":
        if view != "primary":
            raise ValueError("branch_diagnostics view is not valid for trimmed")
        source = stage_dir / "fastq"
        paths = sorted(source.glob("*.fastq.gz")) + sorted(source.glob("*.fastq"))
        items = _items(paths, "trimmed_mate")
    elif is_branch_dataset and stage == "pr_trimmed":
        if view not in {"primary", "branch_diagnostics"}:
            raise ValueError("invalid view")
        source = stage_dir / "sync"
        items = []
        for branch_dir in sorted(path for path in source.iterdir() if path.is_dir()):
            branch_inputs = list(branch_dir.glob("*pair-pass.fastq.gz"))
            if branch_inputs:
                items.extend(_items(branch_inputs, "synchronised_primer_branch_mate", branch_dir.name))
    elif is_branch_dataset and stage == "merged":
        if view == "primary":
            items = _items(stage_dir.glob("*_all_assemble-pass.fastq.gz"), "combined_merged")
        else:
            items = []
            for branch_dir in sorted(path for path in stage_dir.iterdir() if path.is_dir()):
                branch_inputs = list(branch_dir.glob("*_assemble-pass.fastq.gz"))
                if branch_inputs:
                    items.extend(_items(branch_inputs, "merged_branch", branch_dir.name))
    else:
        if view != "primary":
            raise ValueError("branch_diagnostics view is only valid for branch-aware pr_trimmed/merged stages")
        source = stage_dir / "fastq"
        if stage == "merged":
            # AssemblePairs may leave mate-specific *assemble-fail* files next
            # to the consensus. Primary merged QC must contain pass consensus
            # only, otherwise MultiQC mixes failed mates with assembled reads.
            paths = sorted(source.glob("*_assemble-pass.fastq.gz")) + sorted(
                source.glob("*_assemble-pass.fastq")
            )
        else:
            paths = sorted(source.glob("*.fastq.gz")) + sorted(source.glob("*.fastq"))
        role = "merged" if stage == "merged" else f"{stage}_mate"
        items = _items(paths, role)

    if not items:
        raise FileNotFoundError(
            f"No QC inputs for dataset={dataset}, branch={branch}, stage={stage}, view={view}"
        )
    return sorted(items, key=lambda item: str(item.path))


def qc_output_paths(
    volume: str | Path,
    dataset: str,
    stage: str,
    branch: str | None = None,
    view: str = "primary",
) -> QCOutputPaths:
    volume = Path(volume)
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {sorted(STAGES)}")
    if view not in VIEWS:
        raise ValueError(f"view must be one of {sorted(VIEWS)}")
    stage_dir = _stage_dir(volume, dataset, stage, branch)
    suffix = "" if view == "primary" else "_branch_diagnostics"
    fastqc = stage_dir / f"fastqc{suffix}"
    multiqc = stage_dir / f"multiqc{suffix}"
    stem = "_".join(x for x in (dataset, branch, stage, None if view == "primary" else view) if x)
    report = multiqc / f"{stem}_multiqc.html"
    return QCOutputPaths(
        stage=stage_dir,
        fastqc=fastqc,
        multiqc=multiqc,
        report=report,
        input_manifest=multiqc / "qc_input_manifest.tsv",
        run_manifest=multiqc / "qc_run_manifest.json",
    )


def _run_with_heartbeat(command: list[str], log_path: Path, heartbeat: int = 30) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with log_path.open("w") as handle:
        process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
        print(f"PID={process.pid} log={log_path}", flush=True)
        while process.poll() is None:
            print(f"PID={process.pid} elapsed={(time.monotonic() - started) / 60:.1f} min", flush=True)
            time.sleep(heartbeat)
    if process.returncode:
        raise RuntimeError(f"Command failed with rc={process.returncode}; inspect {log_path}")


def _tool(name: str, env_dir: str | Path | None = None) -> str:
    if env_dir:
        candidate = Path(env_dir) / "bin" / name
        if candidate.is_file():
            return str(candidate)
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(name)
    return found


def run_qc(
    volume: str | Path,
    dataset: str,
    stage: str,
    branch: str | None = None,
    view: str = "primary",
    threads: int = 4,
    expected_count: int | None = None,
    env_dir: str | Path | None = None,
    heartbeat: int = 30,
) -> QCOutputPaths:
    """Run FastQC/MultiQC after exact input resolution and write provenance manifests."""
    inputs = resolve_qc_inputs(volume, dataset, stage, branch=branch, view=view)
    if expected_count is not None and len(inputs) != expected_count:
        raise RuntimeError(f"Expected {expected_count} inputs, found {len(inputs)}")
    outputs = qc_output_paths(volume, dataset, stage, branch=branch, view=view)
    shutil.rmtree(outputs.fastqc, ignore_errors=True)
    shutil.rmtree(outputs.multiqc, ignore_errors=True)
    outputs.fastqc.mkdir(parents=True)
    outputs.multiqc.mkdir(parents=True)

    with outputs.input_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "sample", "mate", "technical_branch", "role"], delimiter="\t")
        writer.writeheader()
        for item in inputs:
            row = asdict(item)
            row["path"] = str(row["path"])
            writer.writerow(row)

    _run_with_heartbeat(
        [_tool("fastqc", env_dir), "--threads", str(threads), "--quiet", "--noextract", "--outdir", str(outputs.fastqc), *[str(x.path) for x in inputs]],
        outputs.fastqc / "fastqc.log",
        heartbeat,
    )
    _run_with_heartbeat(
        [_tool("multiqc", env_dir), "-f", "-o", str(outputs.multiqc), "-n", outputs.report.name, str(outputs.fastqc)],
        outputs.multiqc / "multiqc.log",
        heartbeat,
    )
    if not outputs.report.is_file():
        raise RuntimeError(f"MultiQC report missing: {outputs.report}")

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": dataset,
        "branch": branch,
        "stage": stage,
        "view": view,
        "input_count": len(inputs),
        "inputs": [str(item.path) for item in inputs],
        "report": str(outputs.report),
        "python": sys.version,
        "fastqc": _tool("fastqc", env_dir),
        "multiqc": _tool("multiqc", env_dir),
    }
    outputs.run_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(outputs.report, flush=True)
    return outputs
