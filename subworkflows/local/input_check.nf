//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet  // file: /path/to/samplesheet.csv
    seq_center   // string: sequencing center for read group
    with_control // boolean: samplesheet contains controls


    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .branch { row ->
            bam  : row.bam
            fastq: true
        }
        .set { ch_rows }

    ch_rows
        .fastq
        .map { row -> create_fastq_channel(row, seq_center) }
        .set { reads }

    ch_rows
        .bam
        .map { row -> create_bam_channel(row) }
        .set { bams }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    bams                                      // channel: [ val(meta), bam ]  (meta.single_end not yet known)
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row, String seq_center) {
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.control    = row.control

    def read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id - ~/_T\d+$/}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = read_group

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

// Function to get list of [ meta, bam ] for samples that start from an existing alignment.
// meta.single_end is filled in by BAM_INPUT_CHECK, which inspects the reads in the BAM.
def create_bam_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = null
    meta.control    = row.control

    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
    }
    return [ meta, file(row.bam) ]
}
