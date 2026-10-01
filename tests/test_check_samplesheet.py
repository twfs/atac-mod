"""Unit tests for bin/check_samplesheet.py (FastQ and BAM entry points).

Run with:  python -m pytest tests/test_check_samplesheet.py
"""
import csv
import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "check_samplesheet.py"
spec = importlib.util.spec_from_file_location("check_samplesheet", SCRIPT)
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)

PE1, PE2, SE = "s3://b/a_1.fastq.gz", "s3://b/a_2.fastq.gz", "s3://b/c.fq.gz"


def run(tmp_path, text, with_control=False):
    fin, fout = tmp_path / "in.csv", tmp_path / "out.csv"
    fin.write_text(text)
    cs.check_samplesheet(str(fin), str(fout), with_control)
    with open(fout) as fh:
        return list(csv.DictReader(fh))


def fails(tmp_path, text, message, capsys, with_control=False):
    with pytest.raises(SystemExit):
        run(tmp_path, text, with_control)
    assert message in capsys.readouterr().out


def test_fastq_only_matches_legacy_columns(tmp_path):
    rows = run(tmp_path, f"sample,fastq_1,fastq_2,replicate\nA,{PE1},{PE2},1\nB,{SE},,1\n")
    assert rows[0] == {"sample": "A_REP1_T1", "fastq_1": PE1, "fastq_2": PE2, "replicate": "1", "single_end": "0", "control": "", "bam": ""}
    assert rows[1]["single_end"] == "1"


def test_bam_only_sheet_without_fastq_columns(tmp_path):
    rows = run(tmp_path, "sample,replicate,bam\nA,1,s3://b/a.bam\nA,2,s3://b/a2.bam\n")
    assert [(r["sample"], r["bam"], r["single_end"], r["fastq_1"]) for r in rows] == [
        ("A_REP1_T1", "s3://b/a.bam", "", ""),
        ("A_REP2_T1", "s3://b/a2.bam", "", ""),
    ]


def test_bam_column_present_but_empty_is_fastq_sheet(tmp_path):
    rows = run(tmp_path, f"sample,fastq_1,fastq_2,replicate,bam\nA,{PE1},{PE2},1,\nB,{PE1},{PE2},1,\n")
    assert [r["bam"] for r in rows] == ["", ""]


def test_multiple_bam_runs_per_replicate(tmp_path):
    rows = run(tmp_path, "sample,replicate,bam\nA,1,s3://b/a_run1.bam\nA,1,s3://b/a_run2.bam\n")
    assert [r["sample"] for r in rows] == ["A_REP1_T1", "A_REP1_T2"]


def test_error_mixed_fastq_and_bam_same_replicate(tmp_path, capsys):
    fails(tmp_path, f"sample,fastq_1,fastq_2,replicate,bam\nA,{PE1},{PE2},1,\nA,,,1,s3://b/a_run2.bam\n", "cannot be mixed", capsys)


def test_error_mixed_fastq_and_bam_different_samples(tmp_path, capsys):
    fails(tmp_path, f"sample,fastq_1,replicate,bam\nA,,1,s3://b/a.bam\nB,{SE},1,\n", "cannot be mixed", capsys)


def test_column_order_does_not_matter(tmp_path):
    rows = run(tmp_path, "replicate,bam,sample\n1,s3://b/a.bam,A\n1,s3://b/b.bam,B\n")
    assert rows[0]["sample"] == "A_REP1_T1" and rows[1]["bam"] == "s3://b/b.bam"


def test_bam_with_controls(tmp_path):
    rows = run(
        tmp_path,
        f"sample,fastq_1,fastq_2,replicate,control,control_replicate,bam\nT,,,1,IN,1,s3://b/t.bam\nIN,,,1,,,s3://b/in.bam\n",
        with_control=True,
    )
    assert {r["sample"]: r["control"] for r in rows} == {"IN_REP1_T1": "", "T_REP1_T1": "IN_REP1"}


def test_extra_columns_passed_through(tmp_path):
    rows = run(tmp_path, "sample,replicate,bam,batch\nA,1,s3://b/a.bam,x\n")
    assert rows[0]["batch"] == "x"


def test_error_fastq_and_bam_on_same_row(tmp_path, capsys):
    fails(tmp_path, f"sample,fastq_1,replicate,bam\nA,{SE},1,s3://b/a.bam\n", "not both", capsys)


def test_error_bad_bam_extension(tmp_path, capsys):
    fails(tmp_path, "sample,replicate,bam\nA,1,s3://b/a.cram\n", "extension '.bam'", capsys)


def test_error_no_input_column(tmp_path, capsys):
    fails(tmp_path, "sample,replicate\nA,1\n", "fastq_1 and/or bam", capsys)


def test_error_empty_row_inputs(tmp_path, capsys):
    fails(tmp_path, "sample,fastq_1,replicate,bam\nA,,1,\n", "Invalid combination", capsys)


def test_error_mixed_fastq_datatypes(tmp_path, capsys):
    fails(tmp_path, f"sample,fastq_1,fastq_2,replicate\nA,{PE1},{PE2},1\nA,{SE},,2\n", "same datatype", capsys)


def test_error_replicates_not_contiguous(tmp_path, capsys):
    fails(tmp_path, "sample,replicate,bam\nA,2,s3://b/a.bam\n", "Replicate ids must start", capsys)


def test_error_duplicate_rows(tmp_path, capsys):
    fails(tmp_path, "sample,replicate,bam\nA,1,s3://b/a.bam\nA,1,s3://b/a.bam\n", "duplicate rows", capsys)


def test_error_unknown_control(tmp_path, capsys):
    fails(tmp_path, "sample,replicate,control,control_replicate,bam\nA,1,NOPE,1,s3://b/a.bam\n", "Control identifier", capsys, with_control=True)
