import ast
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANNOT = ROOT / "results/ERP003950/notebooks/annotate_mouse.ipynb"
FILTER = ROOT / "results/ERP003950/notebooks/filter_mouse_post_annotation.ipynb"


def code_source(path):
    nb = json.loads(path.read_text())
    return "\n\n".join(
        "".join(cell.get("source", []))
        for cell in nb["cells"]
        if cell.get("cell_type") == "code"
    )


def function_source(source, name):
    tree = ast.parse(source)
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    segment = ast.get_source_segment(source, node)
    assert segment is not None
    return segment


def test_orientation_helper_matches_igblast_airr_semantics():
    ns = {
        "_COMPLEMENT": str.maketrans(
            "ACGTRYKMSWBDHVN",
            "TGCAYRMKSWVHDBN",
        )
    }
    source = code_source(ANNOT)
    exec(function_source(source, "reverse_complement"), ns)
    exec(function_source(source, "is_reverse_complemented"), ns)
    exec(function_source(source, "expected_airr_sequence"), ns)

    query = "AACCGT"
    expected = ns["expected_airr_sequence"]
    assert callable(expected)
    assert expected(query, {"rev_comp": "F"}) == query
    assert expected(query, {"rev_comp": "T"}) == "ACGGTT"


def test_projection_preserves_oriented_airr_sequence_and_indexes_original_query():
    source = code_source(ANNOT)
    build_index = function_source(source, "build_sqlite_index")
    project = function_source(source, "project_to_per_read_airr")
    assert "query_sequence" in build_index
    assert "INSERT INTO ann(sequence,row_json)" in build_index
    assert '(query_sequence, json.dumps(row' in build_index
    assert project.count('row["sequence"] = sequence') == 1


def test_projection_keeps_ambiguous_reads_as_explicit_unannotated_rows():
    ns = {}
    source = code_source(ANNOT)
    exec(function_source(source, "has_ambiguous_bases"), ns)
    checker = ns["has_ambiguous_bases"]
    assert callable(checker)
    assert checker("ACGTN") is True
    assert checker("ACGTR") is True
    assert checker("ACGT") is False
    project = function_source(source, "project_to_per_read_airr")
    assert "unannotated_ambiguous_reads" in project
    assert "has_ambiguous_bases(sequence)" in project


def test_igblast_subprocess_receives_bcr_env_library_path():
    source = code_source(ANNOT)
    runner = function_source(source, "run_unique_igblast")
    assert 'env["LD_LIBRARY_PATH"]' in runner
    assert 'BCR_ENV / "lib"' in runner


def test_filter_compares_fastq_to_orientation_aware_airr_sequence():
    source = code_source(FILTER)
    assert "expected_airr_sequence" in source
    body = function_source(source, "filter_sample")
    assert "expected_airr_sequence(fq_seq, row)" in body
    evaluate = function_source(source, "evaluate")
    assert 'reasons.append("ambiguous_sequence")' in evaluate
