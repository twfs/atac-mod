#!/usr/bin/env python3

import os
import sys
import csv
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/atacseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument("--with_control", action="store_true", help="shows output")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# Columns understood by the pipeline. Any other column is passed through unchanged.
KNOWN_COLS = ["sample", "fastq_1", "fastq_2", "replicate", "bam", "control", "control_replicate"]


def check_samplesheet(file_in, file_out, with_control=False):
    """
    Check the samplesheet and write a normalised copy with one row per sequencing run.

    Columns are matched by name, so their order does not matter.
    Each row starts from EITHER FastQ files OR an existing BAM file:

    sample,fastq_1,fastq_2,replicate,bam
    WT,s3://bucket/WT_R1_1.fastq.gz,s3://bucket/WT_R1_2.fastq.gz,1,
    WT,,,2,s3://bucket/WT_REP2.bam
    KO,s3://bucket/KO_R1.fastq.gz,,1,

    - `fastq_2` and `bam` columns are optional. The `fastq_1` column may be omitted if every row uses `bam`.
    - For BAM rows single-end / paired-end is detected from the BAM itself later in the pipeline,
      so it is written as an empty `single_end` value here.
    - With --with_control, `control` and `control_replicate` columns are required.

    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/atacseq/samplesheet/v2.1/samplesheet_test.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig", newline="") as fin:
        reader = csv.reader(fin)
        try:
            header = [x.strip().strip('"') for x in next(reader)]
        except StopIteration:
            print_error("Samplesheet is empty!")

        ## Check header
        required = ["sample", "replicate"] + (["control", "control_replicate"] if with_control else [])
        missing = [x for x in required if x not in header]
        if missing or not ("fastq_1" in header or "bam" in header):
            if not ("fastq_1" in header or "bam" in header):
                missing.append("fastq_1 and/or bam")
            print(
                "ERROR: Please check samplesheet header -> {}\nMissing column(s): {}".format(
                    ",".join(header), ", ".join(missing)
                )
            )
            sys.exit(1)
        if len(set(header)) != len(header):
            print(f"ERROR: Please check samplesheet header -> duplicate column names: {','.join(header)}")
            sys.exit(1)
        extra_cols = [x for x in header if x not in KNOWN_COLS]

        ## Check sample entries
        for lspl in reader:
            line = ",".join(lspl)
            if not line.strip(", \t"):
                continue
            lspl = [x.strip().strip('"') for x in lspl]
            if len(lspl) < len(header):
                lspl += [""] * (len(header) - len(lspl))
            elif len(lspl) > len(header):
                print_error("Invalid number of columns (expected {})!".format(len(header)), "Line", line)
            row = dict(zip(header, lspl))

            sample = row.get("sample", "")
            fastq_1 = row.get("fastq_1", "")
            fastq_2 = row.get("fastq_2", "")
            replicate = row.get("replicate", "")
            bam = row.get("bam", "")
            control = row.get("control", "") if with_control else ""
            control_replicate = row.get("control_replicate", "") if with_control else ""

            ## Check sample name entries
            if sample.find(" ") != -1:
                print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Check BAM file extension
            if bam:
                if bam.find(" ") != -1:
                    print_error("BAM file contains spaces!", "Line", line)
                if not bam.endswith(".bam"):
                    print_error("BAM file does not have extension '.bam'!", "Line", line)

            ## Check replicate column is integer
            if not replicate.isdecimal():
                print_error("Replicate id not an integer!", "Line", line)

            if with_control and control:
                if control.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for control: {control}")
                    control = control.replace(" ", "_")
                if not control_replicate.isdecimal():
                    print_error("Control replicate id not an integer!", "Line", line)
                control = "{}_REP{}".format(control, control_replicate)

            ## Auto-detect input type: paired-end FastQ / single-end FastQ / BAM
            ## sample_info = [ fastq_1, fastq_2, replicate, single_end, control, bam ]
            if bam and (fastq_1 or fastq_2):
                print_error("Provide either FastQ file(s) or a BAM file for a row, not both!", "Line", line)
            elif bam:
                sample_info = ["", "", replicate, "", control, bam]
            elif fastq_1 and fastq_2:
                sample_info = [fastq_1, fastq_2, replicate, "0", control, ""]
            elif fastq_1 and not fastq_2:
                sample_info = [fastq_1, fastq_2, replicate, "1", control, ""]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = {sample: {replicate: [[ fastq_1, fastq_2, replicate, single_end, control, bam, extras... ]]}}
            replicate = int(replicate)
            sample_info = sample_info + [row[x] for x in extra_cols]
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = {}
            if replicate not in sample_mapping_dict[sample]:
                sample_mapping_dict[sample][replicate] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample][replicate]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample][replicate].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(["sample", "fastq_1", "fastq_2", "replicate", "single_end", "control", "bam"] + extra_cols)
                + "\n"
            )

            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that replicate ids are in format 1..<num_replicates>
                uniq_rep_ids = sorted(list(set(sample_mapping_dict[sample].keys())))
                if len(uniq_rep_ids) != max(uniq_rep_ids) or 1 != min(uniq_rep_ids):
                    print_error(
                        "Replicate ids must start with 1..<num_replicates>!",
                        "Sample",
                        "{}, replicate ids: {}".format(sample, ",".join([str(x) for x in uniq_rep_ids])),
                    )

                ## Check that FastQ runs of a sample are of the same datatype i.e. single-end / paired-end.
                ## BAM rows (single_end == "") are checked once the BAM has been inspected in the pipeline.
                datatypes = set(
                    x[3] for runs in sample_mapping_dict[sample].values() for x in runs if x[3] != ""
                )
                if len(datatypes) > 1:
                    print_error(
                        f"Multiple replicates and runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )

                for replicate in sorted(sample_mapping_dict[sample].keys()):
                    for val in sample_mapping_dict[sample][replicate]:
                        control = "_REP".join(val[4].split("_REP")[:-1])
                        control_replicate = val[4].split("_REP")[-1]
                        if control and (
                            control not in sample_mapping_dict.keys()
                            or int(control_replicate) not in sample_mapping_dict[control].keys()
                        ):
                            print_error(
                                f"Control identifier and replicate has to match a provided sample identifier and replicate!",
                                "Control",
                                val[4],
                            )

                    ## Write to file
                    for idx, val in enumerate(sample_mapping_dict[sample][replicate]):
                        sample_id = "{}_REP{}_T{}".format(sample, replicate, idx + 1)
                        fout.write(",".join([sample_id] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet", file_in)


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT, args.with_control)


if __name__ == "__main__":
    sys.exit(main())
