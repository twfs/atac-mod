#!/usr/bin/env bash
#
# Decide from a SAM/BAM header whether duplicate reads have already been REMOVED from the BAM.
#
# Usage:  samtools view -H in.bam | bam_dedup_status.sh
# Prints: "true<TAB><program command line that removed duplicates>"  or  "false"
#
# Duplicate removal is recognised from the @PG command lines (CL:) of:
#   - Picard / GATK MarkDuplicates with REMOVE_DUPLICATES or REMOVE_SEQUENCING_DUPLICATES true
#   - GATK MarkDuplicatesSpark with --remove-all-duplicates or --remove-sequencing-duplicates
#   - samtools markdup -r / sambamba markdup -r|--remove-duplicates
#   - umi_tools dedup
#   - samtools view -F / --exclude-flags whose value includes 0x400 (1024, "DUP")
# Duplicate MARKING alone (flags set, reads kept) is not removal and returns "false".
# BAMs whose @PG history has been stripped cannot be classified and return "false".

awk -F'\t' '
function hex2dec(h,    i, c, v) {
    h = tolower(h); sub(/^0x/, "", h); v = 0
    for (i = 1; i <= length(h); i++) {
        c = index("0123456789abcdef", substr(h, i, 1))
        if (c == 0) return -1
        v = v * 16 + c - 1
    }
    return v
}
function has_dup_bit(val,    v) {
    if (toupper(val) ~ /DUP/) return 1
    if (val ~ /^0[xX][0-9a-fA-F]+$/) v = hex2dec(val)
    else if (val ~ /^[0-9]+$/)       v = val + 0
    else return 0
    return (int(v / 1024) % 2 == 1)
}
function removes_dups(cl,    n, t, i, v, lc) {
    lc = tolower(cl)
    if (lc ~ /markduplicates/ && lc ~ /remove_(sequencing_)?duplicates[ =:]+true/) return 1
    if (lc ~ /markduplicatesspark/ && lc ~ /remove-(all|sequencing)-duplicates( +true|=true| +-|$)/) return 1
    if (lc ~ /umi_tools(\.py)? +dedup/) return 1
    n = split(cl, t, / +/)
    if (lc ~ /(samtools|sambamba)[^ ]* +markdup/) {
        for (i = 1; i <= n; i++) if (t[i] == "-r" || t[i] == "--remove-duplicates") return 1
    }
    if (lc ~ /samtools[^ ]* +view/) {
        for (i = 1; i <= n; i++) {
            v = ""
            if (t[i] == "-F" || t[i] == "--exclude-flags") v = t[i + 1]
            else if (t[i] ~ /^-F./)                        v = substr(t[i], 3)
            else if (t[i] ~ /^--exclude-flags=/)            v = substr(t[i], 17)
            if (v != "" && has_dup_bit(v)) return 1
        }
    }
    return 0
}
# Read the whole header (no early exit) so the upstream "samtools view -H" never gets SIGPIPE under pipefail
$1 == "@PG" && found == "" {
    cl = ""
    for (i = 2; i <= NF; i++) if ($i ~ /^CL:/) cl = substr($i, 4)
    if (cl != "" && removes_dups(cl)) found = cl
}
END {
    if (found != "") print "true\t" found
    else             print "false"
}
'
