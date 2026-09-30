process BAM_INPUT_CHECK {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam)
    path  fai

    output:
    tuple val(meta), path(bam), env('SINGLE_END'), emit: bam
    path  "*.reference_check.txt"              , emit: report
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def n_reads  = task.ext.n_reads ?: 100000
    """
    ## 1. Reference check: every @SQ in the BAM header must be in the genome .fai with the same length
    samtools view -H $bam | awk -F'\\t' '\$1 == "@SQ" { sn = ""; ln = ""; for (i = 2; i <= NF; i++) { if (\$i ~ /^SN:/) sn = substr(\$i, 4); if (\$i ~ /^LN:/) ln = substr(\$i, 4) } print sn "\\t" ln }' | sort > bam_contigs.txt
    cut -f1,2 $fai | sort > genome_contigs.txt

    n_bam=\$(wc -l < bam_contigs.txt)
    if [ "\$n_bam" -eq 0 ]; then
        echo "ERROR: BAM file '$bam' (sample ${meta.id}) has no @SQ lines in its header - is it an unaligned BAM?" >&2
        exit 1
    fi

    comm -23 bam_contigs.txt genome_contigs.txt > mismatched_contigs.txt
    n_mismatch=\$(wc -l < mismatched_contigs.txt)
    n_genome_only=\$(cut -f1 bam_contigs.txt | comm -13 - <(cut -f1 genome_contigs.txt) | wc -l)

    {
        echo -e "sample\\t${meta.id}"
        echo -e "bam\\t$bam"
        echo -e "bam_contigs\\t\$n_bam"
        echo -e "bam_contigs_not_matching_genome\\t\$n_mismatch"
        echo -e "genome_contigs_absent_from_bam\\t\$n_genome_only"
    } > ${prefix}.reference_check.txt

    if [ "\$n_mismatch" -gt 0 ]; then
        echo "ERROR: BAM file '$bam' (sample ${meta.id}) was not aligned to the genome supplied to this run." >&2
        echo "\$n_mismatch of \$n_bam contigs in the BAM header are missing from the genome .fai or differ in length. First mismatches (name<TAB>length):" >&2
        head -n 10 mismatched_contigs.txt >&2
        echo "Re-run with the matching --genome / --fasta, or start this sample from FastQ instead." >&2
        exit 1
    fi
    if [ "\$n_genome_only" -gt 0 ]; then
        echo "WARNING: \$n_genome_only genome contigs are not in the BAM header for sample ${meta.id}; continuing." >&2
    fi

    ## 2. Detect single-end / paired-end from the first $n_reads primary alignments
    set +o pipefail
    n_total=\$(samtools view -F 0x900 $bam | head -n $n_reads | wc -l)
    n_paired=\$(samtools view -F 0x900 -f 0x1 $bam | head -n $n_reads | wc -l)
    set -o pipefail

    if [ "\$n_total" -eq 0 ]; then
        echo "ERROR: BAM file '$bam' (sample ${meta.id}) contains no primary alignments." >&2
        exit 1
    elif [ "\$n_paired" -eq 0 ]; then
        SINGLE_END=true
    elif [ "\$n_paired" -eq "\$n_total" ]; then
        SINGLE_END=false
    else
        echo "ERROR: BAM file '$bam' (sample ${meta.id}) mixes paired and unpaired reads (\$n_paired of \$n_total checked reads are paired)." >&2
        exit 1
    fi
    echo -e "single_end\\t\$SINGLE_END" >> ${prefix}.reference_check.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    SINGLE_END=false
    touch ${prefix}.reference_check.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.17
    END_VERSIONS
    """
}
