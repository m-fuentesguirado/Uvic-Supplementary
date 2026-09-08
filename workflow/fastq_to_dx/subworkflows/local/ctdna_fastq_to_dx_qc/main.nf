/*
========================================================================================
    Subworkflow: CTDNA_FASTQ_TO_DX_QC

    Runs quality control for the FASTQ-to-duplex-consensus-BAM workflow.

    Things inclused in this workflow:
    - samtools flagstat
    - samtools idxstats
    - samtools stats

    Ideas to include ??
    - Mosdepth panel coverage
    - Picard CollectHsMetrics -- needs to be chechcked
    - fgbio duplex-yield summary
    - Consensus survival summary
    - multi MultiQC report
========================================================================================
*/


include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_MAPPED } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_MAPPED } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_MAPPED } from '../../../modules/nf-core/samtools/stats/main'

include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DX } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_DX } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_DX } from '../../../modules/nf-core/samtools/stats/main'


include { MULTIQC as MULTIQC_FASTQ_TO_DX_QC } from '../../../modules/nf-core/multiqc/main'
include { BAM_QC_PICARD as BAM_QC_PICARD_MAPPED } from '../../nf-core/bam_qc_picard/main'
include { BAM_QC_PICARD as BAM_QC_PICARD_DX } from '../../nf-core/bam_qc_picard/main'
workflow CTDNA_FASTQ_TO_DX_QC {

       take:

    /*
     * Raw FastQC ZIP files.
     * tuple(meta, fastqc.zip)
     */
    ch_fastqc_zip

    /*
     * fgbio QC, which includes histogram and other metrics
     * Structure:
     * tuple(meta, histogram.txt)
     */
    ch_fgbio_qc_files

    /*
     * coordiante-sorted mapped BAM before UMI grouping and duplex-consensus generation.
     *
     * Structure:
     * tuple(meta, bam, bai)
     *  
     */
    ch_mapped_bam_bai

    /*
     * Final coordinate-sorted duplex-consensus BAM and index.
     * tuple(meta, bam, bai)
     *
     * Used for samtools QC and the second CollectHsMetrics run.
     */
    ch_dx_final_bam_bai

    /*
     * Reference FASTA, FASTA index and sequence dictionary.
     *
     * Structure:
     * tuple(meta, fasta, fai, dict)
     */
    ch_ref_zipperbams
    
            main:

    /*
     * Rename the final duplex-consensus samples for QC.
     *
     * This avoids MultiQC confusing the pre-consensus BAM
     * with the final duplex-consensus BAM.
     */
    ch_dx_qc_bam_bai = ch_dx_final_bam_bai.map { meta, bam, bai ->

        def qc_meta = meta + [
            id         : "${meta.id}_dx_mapped",
            original_id: meta.id,
            qc_stage   : 'dx_mapped'
        ]

        tuple(
            qc_meta,
            bam,
            bai
        )
    }


    /*
     * Prepare the reference channel for samtools stats.
     *
     * Input:
     * tuple(meta, fasta, fai, dict)
     *
     * Output:
     * tuple(meta, fasta, fai)
     */
    ch_ref_fasta_stats = ch_ref_zipperbams.map { meta, fasta, fai, dict ->
        tuple(
            meta,
            fasta,
            fai
        )
    }


    /*
    ====================================================================================
        PRE-CONSENSUS MAPPED BAM QC
    ====================================================================================
    */

    /*
     * Overall mapping statistics before UMI grouping
     * and duplex-consensus generation.
     */
    SAMTOOLS_FLAGSTAT_MAPPED(
        ch_mapped_bam_bai
    )

    /*
     * Mapped and unmapped reads per chromosome
     * before consensus generation.
     */
    SAMTOOLS_IDXSTATS_MAPPED(
        ch_mapped_bam_bai
    )

    /*
     * Detailed alignment statistics before
     * consensus generation.
     */
    SAMTOOLS_STATS_MAPPED(
        ch_mapped_bam_bai,
        ch_ref_fasta_stats
    )


    /*
    ====================================================================================
        FINAL DUPLEX-CONSENSUS BAM QC
    ====================================================================================
    */

    /*
     * Overall mapping statistics for the final
     * mapped duplex-consensus BAM.
     */
    SAMTOOLS_FLAGSTAT_DX(
        ch_dx_qc_bam_bai
    )

    /*
     * Mapped and unmapped reads per chromosome
     * in the final duplex-consensus BAM.
     */
    SAMTOOLS_IDXSTATS_DX(
        ch_dx_qc_bam_bai
    )

    /*
     * Detailed alignment statistics for the
     * final duplex-consensus BAM.
     */
    SAMTOOLS_STATS_DX(
        ch_dx_qc_bam_bai,
        ch_ref_fasta_stats
    )
    /*
    ====================================================================================
        PREPARE INPUTS FOR NF-CORE BAM_QC_PICARD ---> MULTIQC
    ====================================================================================
    */

    /*
     * Add bait and target intervals to the pre-consensus BAM.
     *
     * Output:
     * tuple(meta, bam, bai, bait_intervals, target_intervals)
     */
    ch_picard_mapped_input = ch_mapped_bam_bai.map { meta, bam, bai ->
        tuple(
            meta,
            bam,
            bai,
            file(params.bait_intervals, checkIfExists: true),
            file(params.target_intervals, checkIfExists: true)
        )
    }

    /*
     * Add bait and target intervals to the final duplex BAM.
     */
    ch_picard_dx_input = ch_dx_qc_bam_bai.map { meta, bam, bai ->
        tuple(
            meta,
            bam,
            bai,
            file(params.bait_intervals, checkIfExists: true),
            file(params.target_intervals, checkIfExists: true)
        )
    }

    /*
     * Split the bundled reference channel into the separate
     * channels required by BAM_QC_PICARD.
     */
    ch_ref_fasta_picard = ch_ref_zipperbams.map { meta, fasta, fai, dict ->
        tuple(meta, fasta)
    }

    ch_ref_fai_picard = ch_ref_zipperbams.map { meta, fasta, fai, dict ->
        tuple(meta, fai)
    }

    ch_ref_dict_picard = ch_ref_zipperbams.map { meta, fasta, fai, dict ->
        tuple(meta, dict)
    }

    /*
     * No GZI is required because the reference is an
     * uncompressed .fasta file.
     */
    ch_ref_gzi_picard = ch_ref_zipperbams.map { meta, fasta, fai, dict ->
        tuple(meta, [])
    }


    /*
    ====================================================================================
        RUN NF-CORE PICARD QC
    ====================================================================================
    */

    
    /*
     * Picard QC before UMI grouping and consensus generation.
     */
    BAM_QC_PICARD_MAPPED(
        ch_picard_mapped_input,
        ch_ref_fasta_picard,
        ch_ref_fai_picard,
        ch_ref_dict_picard,
        ch_ref_gzi_picard
    )

    /*
     * Picard QC on the final mapped duplex-consensus BAM.
     */
    BAM_QC_PICARD_DX(
        ch_picard_dx_input,
        ch_ref_fasta_picard,
        ch_ref_fai_picard,
        ch_ref_dict_picard,
        ch_ref_gzi_picard
    )
        

    /*
    ====================================================================================
        PREPARE MULTIQC INPUT
    ====================================================================================
    */

    /*
     * Using nf-Core BAMQC_PICARD outputs, we need to extract the individual files for MultiQC.
     * flatMap converts either form into individual files for MultiQC.
     */
    ch_picard_mapped_coverage_files =
    BAM_QC_PICARD_MAPPED.out.coverage_metrics.flatMap {
        meta, metric_files -> metric_files
    }

ch_picard_mapped_multiple_files =
    BAM_QC_PICARD_MAPPED.out.multiple_metrics.flatMap {
        meta, metric_files -> metric_files
    }

ch_picard_dx_coverage_files =
    BAM_QC_PICARD_DX.out.coverage_metrics.flatMap {
        meta, metric_files -> metric_files
    }

ch_picard_dx_multiple_files =
    BAM_QC_PICARD_DX.out.multiple_metrics.flatMap {
        meta, metric_files -> metric_files
    }


    /*
     * Combine all QC files into one channel.
     */
    ch_multiqc_files = Channel.empty()
        .mix(
            ch_fastqc_zip.map { meta, file ->
                file
            }
        )
        .mix(
            ch_fgbio_qc_files.map { meta, file ->
                file
            }
        )
        .mix(
            SAMTOOLS_FLAGSTAT_MAPPED.out.flagstat.map { meta, file ->
                file
            }
        )
        .mix(
            SAMTOOLS_IDXSTATS_MAPPED.out.idxstats.map { meta, file ->
                file
            }
        )
        .mix(
            SAMTOOLS_STATS_MAPPED.out.stats.map { meta, file ->
                file
            }
        )
        .mix(
            SAMTOOLS_FLAGSTAT_DX.out.flagstat.map { meta, file ->
                file
            }
        )
        .mix(
            SAMTOOLS_IDXSTATS_DX.out.idxstats.map { meta, file ->
                file
            }
        )
        .mix(
            SAMTOOLS_STATS_DX.out.stats.map { meta, file ->
                file
            }
        )
        .mix(ch_picard_mapped_coverage_files)
        .mix(ch_picard_mapped_multiple_files)
        .mix(ch_picard_dx_coverage_files)
        .mix(ch_picard_dx_multiple_files)
        .collect()


    /*
     * Build the input tuple required by the nf-core MultiQC module.
     */
    ch_multiqc_input = ch_multiqc_files.map { files ->
        tuple(
            [ id: 'ctdna_fastq_to_dx_qc' ],
            files,
            [],
            [],
            [],
            []
        )
    }


    /*
     * Generate the combined QC report.
     */
    MULTIQC_FASTQ_TO_DX_QC(
        ch_multiqc_input
    )

    /*
     * Temporary emit block.
     *
     * Picard and MultiQC outputs will be added after
     * their processes are connected.
     */
    emit:

    mapped_flagstat = SAMTOOLS_FLAGSTAT_MAPPED.out.flagstat
    mapped_idxstats = SAMTOOLS_IDXSTATS_MAPPED.out.idxstats
    mapped_stats    = SAMTOOLS_STATS_MAPPED.out.stats

    dx_flagstat = SAMTOOLS_FLAGSTAT_DX.out.flagstat
    dx_idxstats = SAMTOOLS_IDXSTATS_DX.out.idxstats
    dx_stats    = SAMTOOLS_STATS_DX.out.stats


    /*
     * Pre-consensus Picard QC.
     */
    mapped_picard_coverage = BAM_QC_PICARD_MAPPED.out.coverage_metrics
    mapped_picard_multiple_metrics = BAM_QC_PICARD_MAPPED.out.multiple_metrics

    /*
     * Final duplex-consensus Picard QC.
     */
    dx_picard_coverage = BAM_QC_PICARD_DX.out.coverage_metrics
    dx_picard_multiple_metrics = BAM_QC_PICARD_DX.out.multiple_metrics

    /*
     * MultiQC outputs.
     */
    multiqc_report = MULTIQC_FASTQ_TO_DX_QC.out.report
    multiqc_data   = MULTIQC_FASTQ_TO_DX_QC.out.data
    multiqc_plots  = MULTIQC_FASTQ_TO_DX_QC.out.plots
    
}

