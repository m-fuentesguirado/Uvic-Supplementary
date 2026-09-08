/*
========================================================================================
    Subworkflow: CTDNA_FASTQ_TO_DX

    Description:

    Converts FASTQ files into analysis-ready duplex consensus BAMs.

    Steps in this subworkflow:
    1-FASTQ QC with FastQC
    2- MULTIQC -Collecting fastqc reports into a single multiqc report 
    3-FASTQ to unmapped BAM using fgbio FastqToBam
    4-Convert unmapped BAM back to interleaved FASTQ using samtools fastq
    5-Align interleaved FASTQ using bwa-mem2
    6-Combine BWA alignment with original uBAM metadata using fgbio ZipperBams
       
========================================================================================
*/

include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_FASTQC } from '../modules/nf-core/multiqc/main'
include { FGBIO_FASTQTOBAM } from '../modules/nf-core/fgbio/fastqtobam/main'
include { SAMTOOLS_FASTQ } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_UBAM_QNAME } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_BWA_QNAME } from '../modules/nf-core/samtools/sort/main'
include { BWAMEM2_MEM } from '../modules/nf-core/bwamem2/mem/main'
include { FGBIO_ZIPPERBAMS } from '../modules/nf-core/fgbio/zipperbams/main'
include { FGBIO_GROUPREADSBYUMI } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS } from '../modules/nf-core/fgbio/collectduplexseqmetrics/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS } from '../modules/nf-core/fgbio/callduplexconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS } from '../modules/nf-core/fgbio/filterconsensusreads/main'
//To remap the filtered consensus bams:
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_DX } from '../modules/nf-core/samtools/fastq/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_DX } from '../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DX_UBAM_QNAME } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DX_BWA_QNAME } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_DX_COORD } from '../modules/nf-core/samtools/sort/main'
include { FGBIO_ZIPPERBAMS as FGBIO_ZIPPERBAMS_DX } from '../modules/nf-core/fgbio/zipperbams/main'
// add multiqc
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DX_FINAL } from '../modules/nf-core/samtools/index/main'
include { MULTIQC as MULTIQC_DX_FINAL_MAPPING } from '../modules/nf-core/multiqc/main'
include { CTDNA_FASTQ_TO_DX_QC } from '../subworkflows/local/ctdna_fastq_to_dx_qc/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MAPPED_COORD } from '../modules/nf-core/samtools/sort/main'
include { MAPPED_TO_DEDUP } from '../subworkflows/local/mapped_to_dedup'


workflow CTDNA_FASTQ_TO_DX {

    take:
    ch_samplesheet

    main:

    /*
     * Read FASTQ samplesheet.
     */
    ch_fastqs = ch_samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row ->

            /*
             * Extract run/lane information from the FASTQ filename.
             *This is what is being done here:
             * 
             * 
             */
            def r1_name = file(row.fastq_1).getName()
            def name_no_r1 = r1_name.replaceFirst(/_R1\.fastq\.gz$/, '')
            def parts = name_no_r1.tokenize('_')

            def run_id = parts[4]
            def lane_full = parts[6]
            def lane = lane_full.replaceFirst(/^L00/, '')

            tuple(
                [
                    id: row.sample,
                    patient: row.patient,
                    endpoint: row.endpoint,
                    run_id: run_id,
                    lane: lane,
                    lane_full: lane_full
                ],
                [
                    file(row.fastq_1),
                    file(row.fastq_2)
                ]
            )
        }




    /*
     * Reference channels.
     *
     * BWAMEM2_MEM requires:
     * tuple val(meta2), path(index)
     * tuple val(meta3), path(fasta)
     *
     * FGBIO_ZIPPERBAMS requires:
     * tuple val(meta2), path(fasta), path(index), path(dict)
     *
     * These reference files come from params in nextflow.config.
     */
    ch_bwamem2_index = Channel
        .fromPath(params.bwamem2_index, checkIfExists: true)
        .collect()
        .map { index_files ->
            tuple(
                [ id: 'hg38_bwamem2_index' ],
                index_files
            )
        }

    ch_ref_fasta_bwa = Channel.value(
        tuple(
            [ id: 'hg38_fasta' ],
            file(params.ref_fasta)
        )
    )

    ch_ref_fasta_sort = Channel.value(
        tuple(
            [ id: 'hg38_fasta_sort' ],
            file(params.ref_fasta),
            file(params.ref_fasta_fai)
        )
    )

    ch_ref_zipperbams = Channel.value(
        tuple(
            [ id: 'hg38_reference' ],
            file(params.ref_fasta),
            file(params.ref_fasta_fai),
            file(params.ref_dict)
        )
    )

        /*
     * Step 1:
     * Run FastQC on raw FASTQ files.
     */
    FASTQC(
        ch_fastqs
    )

    /*
     * MultiQC report for FastQC only.
     */
ch_multiqc_fastqc_files = FASTQC.out.zip
    .map { meta, zip -> zip }
    .collect()

ch_multiqc_fastqc_input = ch_multiqc_fastqc_files.map { files ->
    tuple(
        [ id: 'fastqc' ],
        files,
        [],
        [],
        [],
        []
    )
}

MULTIQC_FASTQC(ch_multiqc_fastqc_input)

    /*
     * Step 3:
     * Convert FASTQ files to unmapped BAM using fgbio FastqToBam.
     */
    FGBIO_FASTQTOBAM(
        ch_fastqs
    )
    /*
     * Queryname-sort the unmapped BAM.
     *
     * ZipperBams later requires the unmapped BAM and mapped BAM
     * to have the same queryname/read-name order.
     */
    SAMTOOLS_SORT_UBAM_QNAME(
        FGBIO_FASTQTOBAM.out.bam,
        ch_ref_fasta_sort,
        ''
    )
    /*
     * Step 4:
     * Convert unmapped BAM back to interleaved FASTQ.
     *
     * The SAMTOOLS_FASTQ module expects:
     * tuple val(meta), path(input)
     * val interleave
     *
     * FGBIO_FASTQTOBAM.out.bam gives:
     * tuple val(meta), path(ubam)
     *
     * So the first input already has the correct format.
     *
     * interleave = true because bwa-mem2 will later use -p.
     *
     * Interleaved means R1 and R2 are in a single FASTQ file:
     * R1, R2, R1, R2...
     */
   
SAMTOOLS_FASTQ(
    SAMTOOLS_SORT_UBAM_QNAME.out.bam,
    true
)
    /*
     * Step 5:
     * Align interleaved FASTQ using bwa-mem2.
     *
     * BWAMEM2_MEM expects:
     * tuple val(meta), path(reads)
     * tuple val(meta2), path(index)
     * tuple val(meta3), path(fasta)
     * val sort_bam
     *
     * SAMTOOLS_FASTQ.out.interleaved gives:
     * tuple val(meta), path(interleaved_fastq)
     *
     * sort_bam = false because the original bash workflow pipes raw BWA output
     * directly into fgbio ZipperBams.
     */
    BWAMEM2_MEM(
        SAMTOOLS_FASTQ.out.interleaved,
        ch_bwamem2_index,
        ch_ref_fasta_bwa,
        false
    )
    /*
     * Queryname-sort the BWA BAM.
     *
     * This makes the mapped BAM read-name order match the unmapped BAM.
     */
    SAMTOOLS_SORT_BWA_QNAME(
        BWAMEM2_MEM.out.bam,
        ch_ref_fasta_sort,
        ''
    )
    /*
     * Step 6:
     * Combine BWA alignment with original uBAM metadata using fgbio ZipperBams.
     */
   ch_zipperbams_input = SAMTOOLS_SORT_BWA_QNAME.out.bam.join(
    SAMTOOLS_SORT_UBAM_QNAME.out.bam
    )

    FGBIO_ZIPPERBAMS(
        ch_zipperbams_input,
        ch_ref_zipperbams
    )

  /*
 * Create a separate coordinate-sorted and indexed BAM
 * for pre-consensus mapping and capture QC.
 *
 * The original queryname-sorted ZipperBams output still
 * continues into GroupReadsByUmi.
 */
ch_mapped_bam_for_qc = FGBIO_ZIPPERBAMS.out.bam.map { meta, bam ->

    def qc_meta = meta + [
        id         : "${meta.id}_raw_mapped",
        original_id: meta.id,
        qc_stage   : 'raw_mapped'
    ]

    tuple(
        qc_meta,
        bam
    )
}

/*
 * Coordinate-sort the pre-consensus mapped BAM
 * and create its BAI index.
 */
SAMTOOLS_SORT_MAPPED_COORD(
    ch_mapped_bam_for_qc,
    ch_ref_fasta_sort,
    'bai'
)

/*
 * Deduplicate pre-consensus mapped BAMs for fragmentomics.
 */
MAPPED_TO_DEDUP(
    SAMTOOLS_SORT_MAPPED_COORD.out.bam
)

/*
 * Join the coordinate-sorted BAM and BAI.
 *
 * Structure:
 * tuple(meta, bam, bai)
 */
ch_mapped_bam_bai = SAMTOOLS_SORT_MAPPED_COORD.out.bam.join(
    SAMTOOLS_SORT_MAPPED_COORD.out.index
)
     /*
     * Step 7:
     * Group reads by UMI.
     */
    FGBIO_GROUPREADSBYUMI(
        FGBIO_ZIPPERBAMS.out.bam,
        'paired'
    )

    /*
     * MultiQC report after UMI grouping.
     */

    /*
     * Prepare input for CollectDuplexSeqMetrics.
     *
     * Input required:
     * tuple val(meta), path(grouped_bam), path(interval_list)
     */
    ch_grouped_bam_for_metrics = FGBIO_GROUPREADSBYUMI.out.bam
        .map { meta, grouped_bam ->
            tuple(
                meta,
                grouped_bam,
                file(params.panel_intervals)
            )
        }

    /*
     * Step 8:
     * Collect duplex sequencing metrics.
     */
    FGBIO_COLLECTDUPLEXSEQMETRICS(
        ch_grouped_bam_for_metrics
    )

 
 /*
 * Combine the text-based fgbio QC outputs into one channel.
 *
 * strucuture
 * tuple(meta, metric_file)
 */
ch_fgbio_qc_files = Channel.empty()
    .mix(FGBIO_GROUPREADSBYUMI.out.histogram)
    .mix(FGBIO_GROUPREADSBYUMI.out.read_metrics)
    .mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.family_sizes)
    .mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_family_sizes)
    .mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_yield_metrics)
    .mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.umi_counts)
    .mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_umi_counts)
     


/*
 * Build duplex consensus reads.
 *
 * Input 
 * - grouped BAM from GroupReadsByUmi
 *tuple val(meta), path(grouped_bam)
    val min_reads
    val min_baseq
 * Output:
 * - unmapped duplex consensus BAM
 tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | tr -d "[:cntrl:]" | sed -e "s/^.*Version: //;s/\\[.*$//"'), topic: versions, emit: versions_fgbio
 */
 FGBIO_CALLDUPLEXCONSENSUSREADS(
   FGBIO_GROUPREADSBYUMI.out.bam,
    params.duplex_min_reads,
    params.duplex_min_input_baseq
)
/*
 * Filter duplex consensus reads.
 *
 * input:
    tuple val(meta), path(bam) #this is the sample channel, taken from the output from above
    tuple val(meta2), path(fasta), path(index), path(dict) # this is the ref channel!
    #the ref channel is the same as the one used for the bwa-mem2 alignment, 
    'ch_ref_zipperbams' so let's us this
    val min_reads
    val min_baseq
    val max_base_error_rate

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | tr -d "[:cntrl:]" | sed -e "s/^.*Version: //;s/\\[.*$//"'), topic: versions, emit: versions_fgbio

 */
FGBIO_FILTERCONSENSUSREADS(
    FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam,
    ch_ref_zipperbams,
    params.duplex_min_reads,
    params.duplex_min_base_quality,
    params.duplex_max_base_error_rate
)



/*****************************************************************************
* BLOCK FOR REMAPPING THE FILTERED CONSENSUS BAMS  using alisases!
******************************************************************************
*/

/*
 * Convert filtered unmapped duplex consensus BAM to interleaved FASTQ.
 */
SAMTOOLS_FASTQ_DX(
    FGBIO_FILTERCONSENSUSREADS.out.bam,
    true
)

/*
 * Remap duplex consensus reads with bwa-mem2.
 */
BWAMEM2_MEM_DX(
    SAMTOOLS_FASTQ_DX.out.interleaved,
    ch_bwamem2_index,
    ch_ref_fasta_bwa,
    false
)

/*
 * Queryname-sort the filtered consensus uBAM.
 */
SAMTOOLS_SORT_DX_UBAM_QNAME(
    FGBIO_FILTERCONSENSUSREADS.out.bam,
    ch_ref_fasta_sort,
    ''
)

/*
 * Queryname-sort the BWA-mapped consensus BAM.
 */
SAMTOOLS_SORT_DX_BWA_QNAME(
    BWAMEM2_MEM_DX.out.bam,
    ch_ref_fasta_sort,
    ''
)

/*
 * Match mapped consensus reads with original consensus metadata.
 */
ch_dx_zipperbams_input = SAMTOOLS_SORT_DX_BWA_QNAME.out.bam.join(
    SAMTOOLS_SORT_DX_UBAM_QNAME.out.bam
)

FGBIO_ZIPPERBAMS_DX(
    ch_dx_zipperbams_input,
    ch_ref_zipperbams
)

/*
 * Final coordinate sort.
 */
SAMTOOLS_SORT_DX_COORD(
    FGBIO_ZIPPERBAMS_DX.out.bam,
    ch_ref_fasta_sort,
    ''
)
/*********************************************************************************
 * multiqc report for the final mapping of the filtered consensus bams
 **********************************************************************************
 */
 /*
 * Index final mapped duplex consensus BAMs.
 */
SAMTOOLS_INDEX_DX_FINAL(
    SAMTOOLS_SORT_DX_COORD.out.bam
)

/*
 * Join final BAMs with their BAI indexes.
 *
 * This creates:
 * tuple(meta, bam, bai)
 */
ch_dx_final_bam_bai = SAMTOOLS_SORT_DX_COORD.out.bam.join(
    SAMTOOLS_INDEX_DX_FINAL.out.index
)

/*
 * Run final FASTQ-to-duplex-BAM QC USING the subworkflow:
 *(raw fasq, fgbio Umi and duplex metrics, mapped reads before consenesus and dinal duplex consensus mapped reads)
 */
CTDNA_FASTQ_TO_DX_QC(
    FASTQC.out.zip,
    ch_fgbio_qc_files,
    ch_mapped_bam_bai,
    ch_dx_final_bam_bai,
    ch_ref_zipperbams
)

       
        emit:

    /*
     * Initial FASTQ-to-uBAM and alignment outputs.
     */
    ubam = FGBIO_FASTQTOBAM.out.bam
    ubam_interleaved_fastq = SAMTOOLS_FASTQ.out.interleaved
    bwa_bam = BWAMEM2_MEM.out.bam

    /*
     *  ONCE ZIPPERBAMS : THERE ARE TWO BRAHCES : QC AND GROUPREADSBYUMI
     Queryname-sorted mapped BAM from the first ZipperBams step.
     *
     * This BAM continues into GroupReadsByUmi.
     */
    mapped_bam = FGBIO_ZIPPERBAMS.out.bam


    /*
     * Coordinate-sorted and indexed copy of the mapped BAM.
     *
     * This branch is used for pre-consensus QC.
     */
    mapped_coord_bam = SAMTOOLS_SORT_MAPPED_COORD.out.bam
    mapped_coord_bai = SAMTOOLS_SORT_MAPPED_COORD.out.index

    /*
     * Deduplicated pre-consensus BAMs for fragmentomics.
     */
    dedup_mapped_bam = MAPPED_TO_DEDUP.out.bam
    dedup_mapped_bai = MAPPED_TO_DEDUP.out.bai
    dedup_metrics = MAPPED_TO_DEDUP.out.metrics


    /*
     * UMI grouping outputs.
     */
    umi_grouped_bam = FGBIO_GROUPREADSBYUMI.out.bam
    umi_family_histogram = FGBIO_GROUPREADSBYUMI.out.histogram
    umi_grouping_metrics = FGBIO_GROUPREADSBYUMI.out.read_metrics


    /*
     * Duplex sequencing metrics.
     */
    umi_family_sizes = FGBIO_COLLECTDUPLEXSEQMETRICS.out.family_sizes
    duplex_family_sizes = FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_family_sizes
    duplex_yield_metrics = FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_yield_metrics
    umi_counts = FGBIO_COLLECTDUPLEXSEQMETRICS.out.umi_counts
    duplex_qc = FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_qc
    duplex_umi_counts = FGBIO_COLLECTDUPLEXSEQMETRICS.out.duplex_umi_counts
    versions_fgbio = FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions_fgbio
    versions_ggplot2 = FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions_ggplot2

    /*
     * Duplex-consensus generation outputs.
     */
    consensus_unmapped_bam = FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam
    filtered_consensus_bam = FGBIO_FILTERCONSENSUSREADS.out.bam


    /*
     * Duplex-consensus remapping outputs.
     */
    consensus_fastq = SAMTOOLS_FASTQ_DX.out.interleaved
    consensus_bwa_bam = BWAMEM2_MEM_DX.out.bam
    consensus_mapped_unsorted_bam = FGBIO_ZIPPERBAMS_DX.out.bam
    consensus_mapped_bam = SAMTOOLS_SORT_DX_COORD.out.bam
    consensus_mapped_bai = SAMTOOLS_INDEX_DX_FINAL.out.index


    /*
     * QC for the mapped BAM before UMI grouping
     * and duplex-consensus generation.
     */
    mapped_flagstat = CTDNA_FASTQ_TO_DX_QC.out.mapped_flagstat
    mapped_idxstats = CTDNA_FASTQ_TO_DX_QC.out.mapped_idxstats
    mapped_stats = CTDNA_FASTQ_TO_DX_QC.out.mapped_stats


    /*
     * QC for the final mapped duplex-consensus BAM.
     */
    dx_final_flagstat = CTDNA_FASTQ_TO_DX_QC.out.dx_flagstat
    dx_final_idxstats = CTDNA_FASTQ_TO_DX_QC.out.dx_idxstats
    dx_final_stats = CTDNA_FASTQ_TO_DX_QC.out.dx_stats


    /*
     * Raw FASTQ QC.
     */
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
    multiqc_fastqc_report = MULTIQC_FASTQC.out.report


}
