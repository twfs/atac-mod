"""Unit tests for bin/bam_dedup_status.sh (has duplicate removal already been applied to a BAM?).

Run with:  python -m pytest tests/test_bam_dedup_status.py
"""
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "bam_dedup_status.sh"

HD = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n"
BWA = "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem -M -t 2 ./bwa/genome r1.fq.gz r2.fq.gz\n"
SORT = "@PG\tID:samtools\tPN:samtools\tPP:bwa\tVN:1.16.1\tCL:samtools sort -@ 2 -o x.sorted.bam x.bam\n"


def pg(cl, pid="p"):
    return f"@PG\tID:{pid}\tPN:{pid}\tCL:{cl}\n"


def status(header):
    out = subprocess.run(["bash", str(SCRIPT)], input=header, capture_output=True, text=True, check=True)
    return out.stdout.rstrip("\n").split("\t")


# Header of OSMOTIC_STRESS_T15_PE_REP1.mLb.clN.sorted.bam from a previous nf-core/atacseq 2.1.2 run (trimmed)
ATACSEQ_CLN = (
    HD
    + BWA
    + SORT
    + pg("MarkDuplicates --INPUT S.mLb.sorted.bam --OUTPUT S.mLb.mkD.sorted.bam --METRICS_FILE S.metrics.txt "
         "--ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT", "MarkDuplicates")
    + pg("samtools view -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -q 1 -L genome.include_regions.bed -b S.mLb.mkD.sorted.bam", "samtools.2")
    + pg("samtools sort -@ 2 -o S.mLb.clN.sorted.bam -T S.mLb.clN.sorted S.mLb.clN.bam", "samtools.7")
)


def test_nfcore_atacseq_cln_bam_is_deduplicated():
    result = status(ATACSEQ_CLN)
    assert result[0] == "true"
    assert "-F 0x0400" in result[1]


def test_raw_aligner_output_is_not_deduplicated():
    assert status(HD + BWA + SORT) == ["false"]


def test_duplicate_marked_only_is_not_deduplicated():
    mkd = pg("MarkDuplicates --INPUT a.bam --OUTPUT b.bam --REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES false")
    assert status(HD + BWA + mkd) == ["false"]


def test_stripped_header_cannot_be_classified():
    assert status(HD) == ["false"]


@pytest.mark.parametrize(
    "cl",
    [
        "MarkDuplicates --INPUT a.bam --OUTPUT b.bam --REMOVE_DUPLICATES true",
        "picard MarkDuplicates INPUT=a.bam OUTPUT=b.bam REMOVE_DUPLICATES=true",
        "MarkDuplicates -I a.bam -O b.bam --REMOVE_SEQUENCING_DUPLICATES true",
        "gatk MarkDuplicatesSpark -I a.bam -O b.bam --remove-all-duplicates true",
        "gatk MarkDuplicatesSpark -I a.bam -O b.bam --remove-sequencing-duplicates",
        "samtools markdup -r -@ 4 in.bam out.bam",
        "sambamba markdup --remove-duplicates -t 4 in.bam out.bam",
        "umi_tools dedup -I in.bam -S out.bam",
        "samtools view -b -F 1024 in.bam",
        "samtools view -b -F0x400 in.bam",
        "samtools view -b -F 3844 in.bam",
        "samtools view -b --exclude-flags DUP,UNMAP in.bam",
        "samtools view -b --exclude-flags=0xF04 in.bam",
    ],
)
def test_duplicate_removal_detected(cl):
    assert status(HD + BWA + pg(cl))[0] == "true"


@pytest.mark.parametrize(
    "cl",
    [
        "samtools markdup -@ 4 in.bam out.bam",
        "sambamba markdup -t 4 in.bam out.bam",
        "gatk MarkDuplicatesSpark -I a.bam -O b.bam --remove-all-duplicates false",
        "samtools view -b -F 4 in.bam",
        "samtools view -b -F 0x904 in.bam",
        "samtools view -b -f 1024 in.bam",
        "samtools view -b -q 30 -F UNMAP,SECONDARY in.bam",
    ],
)
def test_no_duplicate_removal(cl):
    assert status(HD + BWA + pg(cl)) == ["false"]
