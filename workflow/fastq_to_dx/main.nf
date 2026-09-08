#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Publication entry point. The analysis implementation is preserved unchanged.
include { CTDNA_FASTQ_TO_DX } from './workflows/ctdna_fastq_to_dx_dedup'

workflow {
    if (!params.input) {
        error 'Missing --input: provide the tab-separated FASTQ samplesheet.'
    }
    def required_files = [
        'ref_fasta', 'ref_fasta_fai', 'ref_dict',
        'panel_intervals', 'bait_intervals', 'target_intervals'
    ]
    required_files.each { key ->
        if (!params[key]) {
            error "Missing required parameter: --${key}"
        }
        file(params[key], checkIfExists: true)
    }
    if (!params.bwamem2_index) {
        error 'Missing --bwamem2_index: provide the BWA-MEM2 index-file glob.'
    }
    if (workflow.profile.tokenize(',').contains('conda') && !params.dedup_conda_env) {
        error 'The conda profile requires --dedup_conda_env (an environment prefix or YAML containing gatk and samtools).'
    }
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
    CTDNA_FASTQ_TO_DX(ch_samplesheet)
}
