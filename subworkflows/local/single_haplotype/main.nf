include { GFASTATS as GFASTATS_SH                               } from '../../../modules/nf-core/gfastats/main'
include { GFASTATS2 as GFASTATS2_SH                             } from '../../../modules/local/gfastats2/main'
include { BUSCO_BUSCO as BUSCO_BUSCO_SH                         } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as BUSCO_GENERATEPLOT_SH           } from '../../../modules/nf-core/busco/generateplot/main'
include { BUSCO_BUSCO as BUSCO_BUSCO_SH_FINAL                   } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as BUSCO_GENERATEPLOT_SH_FINAL     } from '../../../modules/nf-core/busco/generateplot/main'
include { MERQURY as MERQURY_SH                                 } from '../../../modules/nf-core/merqury/main'
include { OMNIC as OMNIC_SH                                     } from '../../../modules/local/omnic/main'
include { OMNIC as OMNIC_SH_FINAL                               } from '../../../modules/local/omnic/main'
include { YAHS as YAHS_SH                                       } from '../../../modules/nf-core/yahs/main'
include { FCS_FCSGX as FCS_FCSGX_SH                             } from '../../../modules/nf-core/fcs/fcsgx/main'
include { FCSGX_CLEANGENOME as FCSGX_CLEANGENOME_SH             } from '../../../modules/nf-core/fcsgx/cleangenome/main'
include { TIARA_TIARA as TIARA_TIARA_SH                         } from '../../../modules/nf-core/tiara/tiara/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_SH           } from '../../../modules/local/bbmap/filterbyname/main'
include { RENAME_SCAFFOLDS                                       } from '../../../modules/local/rename_scaffolds/main'
include { COVERAGE_TRACKS                                        } from '../coverage_tracks/main'
include { TELO_FINDER                                            } from '../telo_finder/main'
include { PRETEXTMAP as PRETEXTMAP_SH                            } from '../../../modules/nf-core/pretextmap/main'
include { PRETEXTMAP2 as PRETEXTMAP_SH_HIGH_RES                  } from '../../../modules/nf-core/pretextmap2/main'
include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_SH                  } from '../../../modules/nf-core/pretextsnapshot/main'
include { PRETEXTGRAPH as PRETEXTGRAPH_SH_GAPS                   } from '../../../modules/local/pretextgraph/main'
include { PRETEXTGRAPH as PRETEXTGRAPH_SH_COVERAGE               } from '../../../modules/local/pretextgraph/main'
include { PRETEXTGRAPH as PRETEXTGRAPH_SH_TELOMERE               } from '../../../modules/local/pretextgraph/main'

workflow SINGLE_HAPLOTYPE {

    take:
    ch_fasta               // channel: [val(meta), path(fasta)] - precomputed assembly from PREPARE_PRECOMPUTED_ASSEMBLY
    ch_hifi_reads          // channel: [val(meta), path(reads)] - from HIFIADAPTERFILT
    ch_hic_reads           // channel: [val(meta), path(reads)] - from FASTP_HIC
    ch_meryl_db            // channel: [val(meta), path(db)]    - from MERYL_COUNT
    ch_genomescope_summary // channel: [val(meta), path(summary)] - from GENOMESCOPE2
    gxdb                   // val: path to GX database
    tempdir                // val: temp directory
    ch_teloseq             // channel.value: telomere sequence string
    buscomode              // val: busco mode
    buscodb                // val: busco lineage path

    main:
    ch_versions = Channel.empty()

    //
    // STEP 1: GFASTATS - contig assembly QC
    //
    ch_gfastats_sh_in = ch_fasta.join(ch_genomescope_summary)

    GFASTATS_SH (
        ch_gfastats_sh_in,
        "fasta",
        "",
        "hap1",
        "1.hifiasm",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_SH.out.versions.first())

    //
    // STEP 2: BUSCO - contig assembly completeness
    //
    ch_busco_sh_contig_in = GFASTATS_SH.out.assembly
        .map { meta, fasta -> [ meta, [ fasta ] ] }

    BUSCO_BUSCO_SH (
        ch_busco_sh_contig_in,
        buscomode,
        buscodb,
        "1.hifiasm",
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO_SH.out.versions.first())

    BUSCO_GENERATEPLOT_SH (
        BUSCO_BUSCO_SH.out.short_summaries_txt,
        "1.hifiasm"
    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT_SH.out.versions.first())

    //
    // STEP 3: MERQURY - contig assembly quality value
    //
    ch_merqury_sh_in = ch_meryl_db.join(GFASTATS_SH.out.assembly)

    MERQURY_SH (
        ch_merqury_sh_in,
        "contigs",
        "1.hifiasm"
    )
    ch_versions = ch_versions.mix(MERQURY_SH.out.versions.first())

    //
    // STEP 4: OMNIC - map HiC reads to contig assembly (pre-scaffolding)
    //
    ch_omnic_sh_in = ch_hic_reads.join(GFASTATS_SH.out.assembly)
        .map { meta, reads, fasta -> [ meta, reads, fasta ] }

    OMNIC_SH (
        ch_omnic_sh_in,
        "hap1",
        tempdir
    )
    ch_versions = ch_versions.mix(OMNIC_SH.out.versions.first())

    //
    // STEP 5: YAHS - scaffold with HiC (single haplotype)
    //
    ch_yahs_sh_in = OMNIC_SH.out.omnic_bam
        .join(GFASTATS_SH.out.assembly)
        .join(OMNIC_SH.out.omnic_fai)

    YAHS_SH (
        ch_yahs_sh_in,
        ".2.yahs.hap1"
    )
    ch_versions = ch_versions.mix(YAHS_SH.out.versions.first())

    //
    // STEP 6: DECONTAMINATION - FCS-GX, clean genome, Tiara, BBmap
    //
    FCS_FCSGX_SH (
        YAHS_SH.out.scaffolds_fasta,
        gxdb,
        ".2.yahs.hap1.NCBI"
    )
    ch_versions = ch_versions.mix(FCS_FCSGX_SH.out.versions.first())

    FCSGX_CLEANGENOME_SH (
        YAHS_SH.out.scaffolds_fasta.join(FCS_FCSGX_SH.out.fcs_gx_report)
    )
    ch_versions = ch_versions.mix(FCSGX_CLEANGENOME_SH.out.versions.first())

    TIARA_TIARA_SH (
        FCSGX_CLEANGENOME_SH.out.cleaned,
        ".2.yahs.hap1.tiara"
    )
    ch_versions = ch_versions.mix(TIARA_TIARA_SH.out.versions.first())

    ch_bbmap_sh_in = FCSGX_CLEANGENOME_SH.out.cleaned
        .join(TIARA_TIARA_SH.out.contig_removal)

    BBMAP_FILTERBYNAME_SH (
        ch_bbmap_sh_in,
        "2.tiara.hap1"
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_SH.out.versions.first())

    //
    // STEP 7: GFASTATS2 - scaffold assembly QC (post-decontamination)
    //
    ch_gfastats2_sh_in = BBMAP_FILTERBYNAME_SH.out.scaffolds
        .join(ch_genomescope_summary)

    GFASTATS2_SH (
        ch_gfastats2_sh_in,
        "fa",
        "",
        "hap1_scaffolds",
        "2.tiara",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS2_SH.out.versions.first())

    //
    // STEP 8: BUSCO - scaffold assembly completeness (post-decontamination)
    //
    ch_busco_sh_final_in = BBMAP_FILTERBYNAME_SH.out.scaffolds
        .map { meta, fasta -> [ meta, [ fasta ] ] }

    BUSCO_BUSCO_SH_FINAL (
        ch_busco_sh_final_in,
        buscomode,
        buscodb,
        "2.tiara",
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO_SH_FINAL.out.versions.first())

    BUSCO_GENERATEPLOT_SH_FINAL (
        BUSCO_BUSCO_SH_FINAL.out.short_summaries_txt,
        "2.tiara."
    )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT_SH_FINAL.out.versions.first())

    //
    // STEP 9: RENAME_SCAFFOLDS - rename sequences with sequential SCAFFOLD_N names
    //
    // Pass gfastats summary as a hard dependency so that if BBMAP re-runs on resume,
    // gfastats re-runs too and invalidates the RENAME_SCAFFOLDS cache. Prevents stale output.
    RENAME_SCAFFOLDS (
        BBMAP_FILTERBYNAME_SH.out.scaffolds,
        GFASTATS2_SH.out.assembly_summary
    )
    ch_versions = ch_versions.mix(RENAME_SCAFFOLDS.out.versions.first())

    //
    // STEP 10: COVERAGE_TRACKS - map HiFi reads to renamed scaffolds
    //
    ch_coverage_sh_in = ch_hifi_reads.join(RENAME_SCAFFOLDS.out.renamed_fasta)
        .map { meta, reads, fasta -> [ meta, reads, fasta ] }

    COVERAGE_TRACKS (
        ch_coverage_sh_in
    )
    ch_versions = ch_versions.mix(COVERAGE_TRACKS.out.versions)

    //
    // STEP 11: TELO_FINDER - find telomeric sequences
    //
    TELO_FINDER (
        RENAME_SCAFFOLDS.out.renamed_fasta,
        ch_teloseq,
        false,
        false
    )
    ch_versions = ch_versions.mix(TELO_FINDER.out.versions)

    ch_telomere_sh_for_pretext = TELO_FINDER.out.bedgraph_file
        .map { meta, bedgraphs -> [ meta, bedgraphs[0] ] }

    //
    // STEP 12: OMNIC_SH_FINAL - map HiC reads to final renamed scaffolds (for Pretext)
    //
    ch_omnic_sh_final_in = ch_hic_reads.join(RENAME_SCAFFOLDS.out.renamed_fasta)
        .map { meta, reads, fasta -> [ meta, reads, fasta ] }

    OMNIC_SH_FINAL (
        ch_omnic_sh_final_in,
        "hap1",
        tempdir
    )
    ch_versions = ch_versions.mix(OMNIC_SH_FINAL.out.versions.first())

    //
    // STEP 13: PRETEXTMAP - generate contact map
    //
    PRETEXTMAP_SH (
        OMNIC_SH_FINAL.out.omnic_bam,
        "hap1",
        "3.final"
    )
    ch_versions = ch_versions.mix(PRETEXTMAP_SH.out.versions.first())

    PRETEXTMAP_SH_HIGH_RES (
        OMNIC_SH_FINAL.out.omnic_bam,
        "hap1-hi-res",
        "3.final"
    )
    ch_versions = ch_versions.mix(PRETEXTMAP_SH_HIGH_RES.out.versions.first())

    //
    // STEP 14: PRETEXTSNAPSHOT - generate contact map image
    //
    PRETEXTSNAPSHOT_SH (
        PRETEXTMAP_SH.out.pretext_map,
        "hap1",
        "3.final"
    )
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT_SH.out.versions.first())

    //
    // STEP 15: PRETEXTGRAPH - inject coverage, gaps and telomere tracks
    //
    PRETEXTGRAPH_SH_GAPS (
        PRETEXTMAP_SH.out.pretext_map.join(COVERAGE_TRACKS.out.genomecov),
        "coverage"
    )
    ch_versions = ch_versions.mix(PRETEXTGRAPH_SH_GAPS.out.versions.first())

    PRETEXTGRAPH_SH_COVERAGE (
        PRETEXTGRAPH_SH_GAPS.out.pretext_graph.join(COVERAGE_TRACKS.out.gaps),
        "gaps"
    )
    ch_versions = ch_versions.mix(PRETEXTGRAPH_SH_COVERAGE.out.versions.first())

    PRETEXTGRAPH_SH_TELOMERE (
        PRETEXTGRAPH_SH_COVERAGE.out.pretext_graph.join(ch_telomere_sh_for_pretext),
        "telomere"
    )
    ch_versions = ch_versions.mix(PRETEXTGRAPH_SH_TELOMERE.out.versions.first())

    emit:
    versions = ch_versions
}
