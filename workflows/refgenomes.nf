/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIFIADAPTERFILT                                } from '../modules/local/hifiadapterfilt/main'
include { FASTQC as FASTQC_HIFI                          } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_HIC                           } from '../modules/nf-core/fastqc/main'
include { FASTP as FASTP_HIC                             } from '../modules/nf-core/fastp/main'
include { PREPARE_PRECOMPUTED_ASSEMBLY                   } from '../modules/local/prepare_precomputed_assembly/main'
include { PREPARE_DUAL_ASSEMBLY                          } from '../modules/local/prepare_dual_assembly/main'
include { SINGLE_HAPLOTYPE                               } from '../subworkflows/local/single_haplotype/main'
include { MERYL_COUNT                                    } from '../modules/nf-core/meryl/count/main'
include { MERYL_HISTOGRAM                                } from '../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2                                   } from '../modules/nf-core/genomescope2/main'
include { HIFIASM1 as HIFIASM_SOLO                       } from '../modules/local/hifiasm1/main'
include { GFASTATS as GFASTATS_HIFI_PRIMARY              } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_HIFI_ALT                  } from '../modules/nf-core/gfastats/main'
include { CAT_HIC                                        } from '../modules/local/cat_hic/main'
include { HIFIASM                                        } from '../modules/nf-core/hifiasm/main'
include { GFASTATS as GFASTATS_HAP1                      } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_HAP2                      } from '../modules/nf-core/gfastats/main'
include { BUSCO_BUSCO                                    } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT                             } from '../modules/nf-core/busco/generateplot/main'
include { MERQURY                                        } from '../modules/nf-core/merqury/main'
include { OMNIC as OMNIC_HAP1                            } from '../modules/local/omnic/main'
include { OMNIC as OMNIC_HAP2                            } from '../modules/local/omnic/main'
include { SCAFFOLDING                                    } from '../subworkflows/local/scaffolding/main'
include { DECONTAMINATION                                } from '../subworkflows/local/decontamination/main'
include { GFASTATS2 as GFASTATS_HAP1_FINAL                } from '../modules/local/gfastats2/main'
include { GFASTATS2 as GFASTATS_HAP2_FINAL                } from '../modules/local/gfastats2/main'
include { BUSCO_BUSCO as BUSCO_BUSCO_FINAL               } from '../modules/nf-core/busco/busco/main'
include { BUSCO_GENERATEPLOT as BUSCO_GENERATEPLOT_FINAL } from '../modules/nf-core/busco/generateplot/main'
include { CAT_SCAFFOLDS                                  } from '../modules/local/cat_scaffolds/main'
include { VALIDATE_SCAFFOLD_COUNTS                       } from '../modules/local/validate_scaffold_counts/main'
// PREPARE_SINGLE_HAPLOTYPE_OUTPUTS replaced by SINGLE_HAPLOTYPE subworkflow
include { TAR                                            } from '../modules/local/tar/main'
include { COVERAGE_TRACKS                                } from '../subworkflows/local/coverage_tracks/main'
include { TELO_FINDER                                   } from '../subworkflows/local/telo_finder/main'
include { OMNIC as OMNIC_HAP1_FINAL                      } from '../modules/local/omnic/main'
include { OMNIC as OMNIC_HAP2_FINAL                      } from '../modules/local/omnic/main'
include { OMNIC as OMNIC_DUAL_HAP                        } from '../modules/local/omnic_dual/main'
include { PRETEXTMAP as PRETEXTMAP_HAP_1                 } from '../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_HAP_2                 } from '../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_DUAL_HAP              } from '../modules/nf-core/pretextmap/main'
include { PRETEXTMAP2 as PRETEXTMAP_HIGH_RES            } from '../modules/nf-core/pretextmap2/main'
include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_HAP1        } from '../modules/nf-core/pretextsnapshot/main' 
include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_HAP2        } from '../modules/nf-core/pretextsnapshot/main' 
include { PRETEXTSNAPSHOT as PRETEXTSNAPSHOT_DUAL_HAP    } from '../modules/nf-core/pretextsnapshot/main'
include { PRETEXTGRAPH as PRETEXTGRAPH_GAPS              } from '../modules/local/pretextgraph/main'
include { PRETEXTGRAPH as PRETEXTGRAPH_COVERAGE          } from '../modules/local/pretextgraph/main'
include { PRETEXTGRAPH as PRETEXTGRAPH_TELOMERE          } from '../modules/local/pretextgraph/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                               } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                         } from '../subworkflows/local/utils_oceangenomesrefgenomes_pipeline'
include { databaseVersionsToYAML                         } from '../subworkflows/local/utils_oceangenomesrefgenomes_pipeline'
include { CAT_HIFI                                      } from '../modules/local/cat_hifi/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REFGENOMES {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    
    // Validate scaffolder parameter
    if (params.scaffolder != 'yahs' && params.scaffolder != 'salsa') {
        error "Invalid scaffolder parameter: ${params.scaffolder}. Must be 'yahs' or 'salsa'"
    }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_rclone_in = Channel.empty()

    ch_hifi = ch_samplesheet
        .map {
            meta ->
                meta = meta[0]
                meta.id = meta.sample + "_" + meta.date + "." + meta.version
                return [ meta, meta.hifi_dir ]
        }

    ch_hic = ch_samplesheet
        .map {
            tuple ->
                def meta = tuple[0]  // ← Extract meta from the tuple FIRST
                if (meta.hic_dir != null) {
                    meta.id = meta.sample + "_" + meta.date + "." + meta.version
                    return [ meta, meta.hic_dir ]
                }
    }

    ch_species = ch_samplesheet
        .map {
            meta ->
                meta = meta[0]
                meta.id = meta.sample + "_" + meta.date + "." + meta.version
                return [ meta, meta.species ]
        }
    //
    // Prepare pre-computed assembly channel if provided
    // NOTE: We do NOT add any extra flags to meta here — ch_samplesheet is broadcast
    // to multiple subscribers that share the same underlying Groovy map object.
    // Mutating meta here would corrupt the meta seen by ch_hifi and downstream
    // processes (MERYL, GENOMESCOPE2), breaking channel joins.
    // Precomputed vs. hifiasm branching is done via meta.primary_assembly at the
    // point of branching (line ~579).
    //
    // Helper closure: resolve an assembly path (file or directory) to a single FASTA/GFA file
    def resolve_assembly_file = { String path_str, String label ->
        def p = file(path_str)
        if (!p.exists()) { error "Assembly path not found for ${label}: ${path_str}" }
        if (p.isDirectory()) {
            def files = p.listFiles().findAll { f -> f.name =~ /\.(fa|fasta|gfa)(\.gz)?$/ }
            if (files.isEmpty()) { error "No assembly file (.fa/.fasta/.gfa) found in ${label} directory: ${path_str}" }
            if (files.size() > 1) { log.warn "Multiple assembly files in ${label} directory ${path_str}, using: ${files[0].name}" }
            return files[0]
        }
        return p
    }

    // Dual-precomputed assembly channel (hap1 + hap2 paths provided → full dual-hap pipeline, skip Hifiasm)
    ch_dual_precomputed = ch_samplesheet
        .map { meta ->
            meta = meta[0]
            meta.id = meta.sample + "_" + meta.date + "." + meta.version
            def has_hap1 = meta.hap1_assembly != null && meta.hap1_assembly != ''
            def has_hap2 = meta.hap2_assembly != null && meta.hap2_assembly != ''
            if (has_hap1 && has_hap2) {
                def hap1_file = resolve_assembly_file(meta.hap1_assembly, "${meta.sample} hap1")
                def hap2_file = resolve_assembly_file(meta.hap2_assembly, "${meta.sample} hap2")
                log.info "Sample ${meta.sample}: Using precomputed hap1=${hap1_file.name}, hap2=${hap2_file.name}"
                return [ meta, hap1_file, hap2_file ]
            } else if (has_hap1 ^ has_hap2) {
                error "Sample ${meta.sample}: hap1_assembly and hap2_assembly must both be provided, or both left empty"
            } else {
                return null
            }
        }
        .filter { it != null }

    ch_precomputed_assembly = ch_samplesheet
        .map { meta ->
            meta = meta[0]
            meta.id = meta.sample + "_" + meta.date + "." + meta.version
            def has_assembly = meta.primary_assembly != null && meta.primary_assembly != ''
            if (has_assembly) {
                def assembly_path = file(meta.primary_assembly)

                // Check if path exists
                if (!assembly_path.exists()) {
                    error "Assembly path not found: ${meta.primary_assembly}"
                }

                // If it's a directory, find the assembly file inside
                def assembly_file
                if (assembly_path.isDirectory()) {
                    // Find assembly files (FASTA or GFA, optionally compressed)
                    def assembly_files = assembly_path.listFiles().findAll { f ->
                        f.name =~ /\.(fa|fasta|gfa)(\.gz)?$/
                    }

                    if (assembly_files.isEmpty()) {
                        error "No assembly files (.fa/.fasta/.gfa) found in directory: ${meta.primary_assembly}"
                    }

                    if (assembly_files.size() > 1) {
                        log.warn "Multiple assembly files found in ${meta.primary_assembly}, using: ${assembly_files[0].name}"
                    }

                    assembly_file = assembly_files[0]
                    log.info "Sample ${meta.sample}: Auto-detected assembly ${assembly_file.name} from ${meta.primary_assembly}"
                } else {
                    // It's a file, use it directly
                    assembly_file = assembly_path
                    log.info "Sample ${meta.sample}: Using assembly file ${assembly_file.name}"
                }

                return [ meta, assembly_file ]
            } else {
                return null
            }
        }
        .filter { it != null }

    //
    // MODULE: Run HiFiAdapterFilt
    //
    HIFIADAPTERFILT (
        ch_hifi
    )
    ch_versions = ch_versions.mix(HIFIADAPTERFILT.out.versions.first())


    CAT_HIFI (
        HIFIADAPTERFILT.out.reads)
    ch_versions = ch_versions.mix(CAT_HIFI.out.versions.first())


    //
    // MODULE: Run FastQC on HiFi fastqc files
    //
    FASTQC_HIFI (
        HIFIADAPTERFILT.out.reads,
        "hifi"
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HIFI.out.zip)
    ch_versions = ch_versions.mix(FASTQC_HIFI.out.versions.first())

    //
    // MODULE: Run Meryl
    //
    MERYL_COUNT (
        HIFIADAPTERFILT.out.reads,
        params.kvalue
    )
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())

    //
    // MODULE: Run Meryl histogram
    //
    MERYL_HISTOGRAM (
        MERYL_COUNT.out.meryl_db,
        params.kvalue
    )
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

    //
    // MODULE: Run Genomescope2
    //
    GENOMESCOPE2 (
        MERYL_HISTOGRAM.out.hist
    )
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())


    //
    // MODULE: Run hifi only assembly (conditional on assembly_mode)
    //
    if (params.assembly_mode == 'hifi_only') {
        HIFIASM_SOLO(
            HIFIADAPTERFILT.out.reads,
            "0.hifiasm"
        )
        ch_versions = ch_versions.mix(HIFIASM_SOLO.out.versions.first())

        //
        // MODULE: gfa stats primary and alternate
        //
        ch_gfastats_hifi_only_primary = HIFIASM_SOLO.out.primary_contigs.join(GENOMESCOPE2.out.summary)
        ch_gfastats_hifi_only_alternate = HIFIASM_SOLO.out.alternate_contigs.join(GENOMESCOPE2.out.summary)

        GFASTATS_HIFI_PRIMARY(
            ch_gfastats_hifi_only_primary,
            "fasta",
            "",
            "p_ctg",
            "0.hifiasm",
            [],
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(GFASTATS_HIFI_PRIMARY.out.versions.first())

        GFASTATS_HIFI_ALT(
            ch_gfastats_hifi_only_alternate,
            "fasta",
            "",
            "a_ctg",
            "0.hifiasm",
            [],
            [],
            [],
            []
        )
        ch_versions = ch_versions.mix(GFASTATS_HIFI_ALT.out.versions.first())
        
        // Create ch_contig_assemblies for hifi_only mode
        ch_contig_assemblies = GFASTATS_HIFI_PRIMARY.out.assembly
            .join(GFASTATS_HIFI_ALT.out.assembly)
            .map { meta, primary, alternate ->
                return [ meta, [ primary, alternate ] ]
            }
    }

    //
    // HiFi+HiC assembly mode: Run HiC-assisted assembly (conditional on assembly_mode)
    //
    if (params.assembly_mode == 'hifi_hic') {
        //
        // MODULE: Concatenate Hi-C files together for cases when there is multiple R1 and multiple R2 files
        //
        CAT_HIC (
            ch_hic
        )
        ch_versions = ch_versions.mix(CAT_HIC.out.versions.first())


    //
    // MODULE: Run FastQC on Hi-C fastqc files
    //
    FASTQC_HIC (
        CAT_HIC.out.cat_files,
        "hic"
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HIC.out.zip)
    ch_versions = ch_versions.mix(FASTQC_HIC.out.versions.first())

    FASTP_HIC (
        CAT_HIC.out.cat_files,
        [],
        [],
        [],
        [],
        "hic"
    )
    ch_versions = ch_versions.mix(FASTP_HIC.out.versions.first())

    //
    // CONDITIONAL: Use pre-computed assemblies OR run hifiasm
    //
    
    // Branch: exclude any sample that already has a precomputed assembly from Hifiasm
    ch_assembly_source = ch_samplesheet
        .map { meta ->
            meta = meta[0]
            meta.id = meta.sample + "_" + meta.date + "." + meta.version
            def has_single  = meta.primary_assembly != null && meta.primary_assembly != ''
            def has_dual    = meta.hap1_assembly != null && meta.hap1_assembly != '' &&
                              meta.hap2_assembly != null && meta.hap2_assembly != ''
            return [ meta, has_single || has_dual ]
        }
        .branch {
            precomputed: it[1] == true
            hifiasm:     it[1] == false
        }

    // Path 1a: Single precomputed assembly → single-haplotype subworkflow
    PREPARE_PRECOMPUTED_ASSEMBLY (
        ch_precomputed_assembly
    )
    ch_versions = ch_versions.mix(PREPARE_PRECOMPUTED_ASSEMBLY.out.versions.first())

    // Path 1b: Dual precomputed assemblies → full dual-haplotype pipeline
    PREPARE_DUAL_ASSEMBLY (
        ch_dual_precomputed
    )
    ch_versions = ch_versions.mix(PREPARE_DUAL_ASSEMBLY.out.versions.first())
    
    // Path 2: Run hifiasm for samples without pre-computed assemblies
    ch_hifiasm_samples = ch_assembly_source.hifiasm
        .map { meta, flag -> meta }
        .join(HIFIADAPTERFILT.out.reads)
        .join(FASTP_HIC.out.fastp_hic)
        .map { meta, hifi, hic ->
            return [ meta, hifi, hic[0], hic[1] ]
        }
    
    HIFIASM (
        ch_hifiasm_samples,
        "1.hifiasm",
        [],
        []
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions.first())
    
    // Define teloseq channel once — used by both SINGLE_HAPLOTYPE and TELO_FINDER below
    ch_teloseq = channel.value(params.teloseq)

    //
    // SUBWORKFLOW: Single-haplotype pipeline for precomputed assemblies
    // Runs all QC, scaffolding, decontamination and visualisation for samples
    // that already have an assembly — bypassing HIFIASM and the dual-hap pipeline.
    //
    SINGLE_HAPLOTYPE (
        PREPARE_PRECOMPUTED_ASSEMBLY.out.hap1_fasta,
        HIFIADAPTERFILT.out.reads.filter { meta, reads ->
            meta.primary_assembly != null && meta.primary_assembly != '' &&
            !(meta.hap1_assembly != null && meta.hap1_assembly != '')
        },
        FASTP_HIC.out.fastp_hic.filter { meta, reads ->
            meta.primary_assembly != null && meta.primary_assembly != '' &&
            !(meta.hap1_assembly != null && meta.hap1_assembly != '')
        },
        MERYL_COUNT.out.meryl_db.filter { meta, db ->
            meta.primary_assembly != null && meta.primary_assembly != '' &&
            !(meta.hap1_assembly != null && meta.hap1_assembly != '')
        },
        GENOMESCOPE2.out.summary.filter { meta, summary ->
            meta.primary_assembly != null && meta.primary_assembly != '' &&
            !(meta.hap1_assembly != null && meta.hap1_assembly != '')
        },
        params.gxdb,
        params.tempdir,
        ch_teloseq,
        params.buscomode,
        params.buscodb
    )
    ch_versions = ch_versions.mix(SINGLE_HAPLOTYPE.out.versions)

    //
    // Dual-haplotype pipeline: HIFIASM samples + dual-precomputed assembly samples
    // (Single-precomputed assembly samples are fully handled by SINGLE_HAPLOTYPE above)
    //
    ch_hap1_assemblies = HIFIASM.out.hap1_contigs.mix(PREPARE_DUAL_ASSEMBLY.out.hap1_fasta)
    ch_hap2_assemblies = HIFIASM.out.hap2_contigs.mix(PREPARE_DUAL_ASSEMBLY.out.hap2_fasta)

    ch_gfastats_hap1_in = ch_hap1_assemblies.join(GENOMESCOPE2.out.summary)
    ch_gfastats_hap2_in = ch_hap2_assemblies.join(GENOMESCOPE2.out.summary)

    GFASTATS_HAP1 (
        ch_gfastats_hap1_in,
        "fasta",
        "",
        "hap1",
        "1.hifiasm",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_HAP1.out.versions.first())

    GFASTATS_HAP2 (
        ch_gfastats_hap2_in,
        "fasta",
        "",
        "hap2",
        "1.hifiasm",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_HAP2.out.versions.first())


    // Create ch_contig_assemblies for hifi_hic mode
    ch_contig_assemblies = GFASTATS_HAP1.out.assembly
        .join(GFASTATS_HAP2.out.assembly)
        .map { meta, hap1, hap2 ->
            return [ meta, [ hap1, hap2 ] ]
        }

    //
    // MODULE: Run Busco
    //
    BUSCO_BUSCO (
        ch_contig_assemblies,
        params.buscomode,
        params.buscodb,
        "1.hifiasm",
        []
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run Busco generate_plot
    //
    BUSCO_GENERATEPLOT (
        BUSCO_BUSCO.out.short_summaries_txt,
        "1.hifiasm"
    )
   ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT.out.versions.first())


    //
    // MODULE: Run Merqury
    //
    ch_merqury_in = MERYL_COUNT.out.meryl_db.join(ch_contig_assemblies)

    MERQURY (
        ch_merqury_in,
        "contigs",
        "1.hifiasm"
    )
    ch_versions = ch_versions.mix(MERQURY.out.versions.first())

    //
    // SUBWORKFLOW: Run omnic workflow
    //
    ch_omnic_hap1_in = FASTP_HIC.out.fastp_hic.join(ch_contig_assemblies)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies[0] ]
        }
    
    
    ch_omnic_hap2_in = FASTP_HIC.out.fastp_hic.join(ch_contig_assemblies)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies[1] ]
        }


    OMNIC_HAP1 (
        ch_omnic_hap1_in,
        "hap1",
        params.tempdir
    )
    ch_versions = ch_versions.mix(OMNIC_HAP1.out.versions.first())

    OMNIC_HAP2 (
        ch_omnic_hap2_in,
        "hap2",
        params.tempdir
    )
    ch_versions = ch_versions.mix(OMNIC_HAP2.out.versions.first())


    // Prepare input channels for YAHS
    ch_yahs_hap1_in = OMNIC_HAP1.out.omnic_bam.join(GFASTATS_HAP1.out.assembly).join(OMNIC_HAP1.out.omnic_fai)
    ch_yahs_hap2_in = OMNIC_HAP2.out.omnic_bam.join(GFASTATS_HAP2.out.assembly).join(OMNIC_HAP2.out.omnic_fai)
    
    // Prepare input channels for SALSA
    ch_bamtobed_hap1_in = OMNIC_HAP1.out.omnic_bam
    ch_bamtobed_hap2_in = OMNIC_HAP2.out.omnic_bam
    ch_fasta_index_hap1 = GFASTATS_HAP1.out.assembly.join(OMNIC_HAP1.out.omnic_fai)
    ch_fasta_index_hap2 = GFASTATS_HAP2.out.assembly.join(OMNIC_HAP2.out.omnic_fai)
    ch_hap1_contigs = HIFIASM.out.hap1_contigs
    ch_hap2_contigs = HIFIASM.out.hap2_contigs
    
    // Set up scaffolder suffix for file naming
    scaffolder_suffix = params.scaffolder == 'yahs' ? '.2.yahs' : '.2.salsa'
    
    //
    // SUBWORKFLOW: Run scaffolding (YAHS or SALSA based on parameter)
    //
    SCAFFOLDING (
        ch_yahs_hap1_in,
        ch_yahs_hap2_in,
        ch_bamtobed_hap1_in,
        ch_bamtobed_hap2_in,
        ch_fasta_index_hap1,
        ch_fasta_index_hap2,
        ch_hap1_contigs,  // Add HIFIASM contigs for SALSA
        ch_hap2_contigs,  // Add HIFIASM contigs for SALSA
        params.scaffolder  // 'yahs' or 'salsa'
    )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)

    //
    // SUBWORKFLOW: Run decontamination pipeline
    //
    DECONTAMINATION (
        SCAFFOLDING.out.hap1_scaffolds,
        SCAFFOLDING.out.hap2_scaffolds,
        params.gxdb,
        scaffolder_suffix
    )
    ch_versions = ch_versions.mix(DECONTAMINATION.out.versions)

    
    //
    // MODULE: Run Gfastats again
    //


    ch_gfastats2_hap1_in = DECONTAMINATION.out.hap1_clean_scaffolds.join(GENOMESCOPE2.out.summary)
    ch_gfastats2_hap2_in = DECONTAMINATION.out.hap2_clean_scaffolds.join(GENOMESCOPE2.out.summary)


    GFASTATS_HAP1_FINAL (
        ch_gfastats2_hap1_in,
        "fa",
        "",
        "hap1_scaffolds",
        "2.tiara",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_HAP1_FINAL.out.versions.first())


    GFASTATS_HAP2_FINAL (
        ch_gfastats2_hap2_in,
        "fa",
        "",
        "hap2_scaffolds",
        "2.tiara",
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GFASTATS_HAP2_FINAL.out.versions.first())

    
    
    //
    // MODULE: Run Busco again
    //

    busco_final_assemblies_ch = DECONTAMINATION.out.hap1_clean_scaffolds.join(DECONTAMINATION.out.hap2_clean_scaffolds)
    .map {
        meta, hap1_scaffolds, hap2_scaffolds ->
            return [ meta, [ hap1_scaffolds, hap2_scaffolds ] ]
    }

    BUSCO_BUSCO_FINAL (
        busco_final_assemblies_ch,
        params.buscomode,
        params.buscodb,
        "2.tiara",
        []
    )
   ch_versions = ch_versions.mix(BUSCO_BUSCO_FINAL.out.versions.first())

    //
    // MODULE: Run Busco generate_plot again
    
    BUSCO_GENERATEPLOT_FINAL (
        BUSCO_BUSCO_FINAL.out.short_summaries_txt,
        "2.tiara."
        )
    ch_versions = ch_versions.mix(BUSCO_GENERATEPLOT_FINAL.out.versions.first())

    //
    // MODULE: Rename and concatenate scaffolds (dual-haplotype HIFIASM samples only)
    // Precomputed single-haplotype samples are handled by SINGLE_HAPLOTYPE subworkflow.
    //
    ch_filtered_scaffolds = DECONTAMINATION.out.hap1_clean_scaffolds
        .join(DECONTAMINATION.out.hap2_clean_scaffolds)
        .map { meta, hap1_scaffolds, hap2_scaffolds ->
            return [ meta, [ hap1_scaffolds, hap2_scaffolds ] ]
        }

    // Join gfastats outputs so CAT_SCAFFOLDS has a hard dependency on them.
    // This prevents stale-cache mismatches: if BBMAP re-runs (new scaffold count),
    // gfastats re-runs too, which invalidates the CAT_SCAFFOLDS cache and forces it
    // to re-run with the new input rather than returning a cached result.
    ch_cat_gfastats = GFASTATS_HAP1_FINAL.out.assembly_summary
        .join(GFASTATS_HAP2_FINAL.out.assembly_summary)

    CAT_SCAFFOLDS (
        ch_filtered_scaffolds,
        ".2.tiara.hap1.hap2",
        ch_cat_gfastats
    )
    ch_versions = ch_versions.mix(CAT_SCAFFOLDS.out.versions.first())

    ch_combined_scaffolds = CAT_SCAFFOLDS.out.cat_file
    ch_hap1_scaffolds     = CAT_SCAFFOLDS.out.hap1_scaffold
    ch_hap2_scaffolds     = CAT_SCAFFOLDS.out.hap2_scaffold

    // Cross-validate gfastats scaffold counts vs scaffold_counts.txt.
    // Fails loudly if they disagree (indicates a resume cache inconsistency).
    ch_validate = CAT_SCAFFOLDS.out.count_report
        .join(GFASTATS_HAP1_FINAL.out.assembly_summary)
        .join(GFASTATS_HAP2_FINAL.out.assembly_summary)

    VALIDATE_SCAFFOLD_COUNTS ( ch_validate )
    ch_versions = ch_versions.mix(VALIDATE_SCAFFOLD_COUNTS.out.versions.first())

    //
    // SUBWORKFLOW: Run TELO_FINDER
    //

    TELO_FINDER (
        ch_combined_scaffolds,
        ch_teloseq,
        false,
        false
    )
    ch_versions = ch_versions.mix(TELO_FINDER.out.versions)

    // Prepare telomere channel for PretextGraph
    ch_telomere_for_pretext = TELO_FINDER.out.bedgraph_file
        .map { meta, bedgraphs ->
            [ meta, bedgraphs[0] ]
        }

    //
    // SUBWORKFLOW: Run coverage tracks
    //

    ch_coverage_tracks_in = HIFIADAPTERFILT.out.reads.join(ch_combined_scaffolds)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies ]
            }
             
    
    COVERAGE_TRACKS (
        ch_coverage_tracks_in
    )
    ch_versions = ch_versions.mix(COVERAGE_TRACKS.out.versions)


    
    //
    // SUBWORFLOW: Run omnic again
    //

    ch_omnic_hap1_in = FASTP_HIC.out.fastp_hic.join(ch_hap1_scaffolds)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies ]
        }

    ch_omnic_hap2_in = FASTP_HIC.out.fastp_hic.join(ch_hap2_scaffolds)
        .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies ]
        }

    ch_omic_dual_in = FASTP_HIC.out.fastp_hic.join(ch_combined_scaffolds)
            .map {
            meta, reads, assemblies ->
                return [ meta, reads, assemblies ]
            }

   OMNIC_HAP1_FINAL (
        ch_omnic_hap1_in,
        "hap1",
        params.tempdir
         )
  ch_versions = ch_versions.mix(OMNIC_HAP1_FINAL.out.versions.first())

      OMNIC_HAP2_FINAL (
       ch_omnic_hap2_in,
        "hap2",
        params.tempdir
    )
    ch_versions = ch_versions.mix(OMNIC_HAP2_FINAL.out.versions.first())

    OMNIC_DUAL_HAP (
        ch_omic_dual_in,
        "dual",
        params.tempdir
    )

    ch_versions = ch_versions.mix(OMNIC_DUAL_HAP.out.versions.first())

    //
    // MODULE: Run Pretext Map
    ///

    PRETEXTMAP_HAP_1 (OMNIC_HAP1_FINAL.out.omnic_bam,
                        "hap1",
                        "2.tiara")
    
    ch_versions = ch_versions.mix(PRETEXTMAP_HAP_1.out.versions.first())
    
    PRETEXTMAP_HAP_2 (OMNIC_HAP2_FINAL.out.omnic_bam,
                    "hap2",
                    "2.tiara")
    
    ch_versions = ch_versions.mix(PRETEXTMAP_HAP_2.out.versions.first())
    
    PRETEXTMAP_DUAL_HAP (OMNIC_DUAL_HAP.out.omnic_bam, 
                    "dual",
                    "2.tiara")
    
    ch_versions = ch_versions.mix(PRETEXTMAP_DUAL_HAP.out.versions.first())


    PRETEXTMAP_HIGH_RES (OMNIC_DUAL_HAP.out.omnic_bam, 
                    "dual-hi-res",
                    "2.tiara")
    
    ch_versions = ch_versions.mix(PRETEXTMAP_HIGH_RES.out.versions.first())
    
    //
    // MODULE: Run Pretext snapshot
    //

    PRETEXTSNAPSHOT_HAP1 (PRETEXTMAP_HAP_1.out.pretext_map,
                        "hap1",
                        "2.tiara")  
    
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT_HAP1.out.versions.first())

    PRETEXTSNAPSHOT_HAP2 (PRETEXTMAP_HAP_2.out.pretext_map,
                    "hap2",
                    "2.tiara")
    
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT_HAP2.out.versions.first())

    PRETEXTSNAPSHOT_DUAL_HAP (PRETEXTMAP_DUAL_HAP.out.pretext_map, 
                    "dual",
                    "2.tiara")
    
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT_DUAL_HAP.out.versions.first())


    //Modual run pretext graph

    PRETEXTGRAPH_GAPS(PRETEXTMAP_DUAL_HAP.out.pretext_map.join(COVERAGE_TRACKS.out.genomecov),
            "coverage")

    ch_versions = ch_versions.mix(PRETEXTGRAPH_GAPS.out.versions.first())

    PRETEXTGRAPH_COVERAGE(PRETEXTMAP_DUAL_HAP.out.pretext_map.join(COVERAGE_TRACKS.out.gaps),
                "gaps")

    ch_versions = ch_versions.mix(PRETEXTGRAPH_COVERAGE.out.versions.first())

// Inject telomere track
    PRETEXTGRAPH_TELOMERE(
        PRETEXTGRAPH_COVERAGE.out.pretext_graph.join(ch_telomere_for_pretext),
        "telomere"
    )
    ch_versions = ch_versions.mix(PRETEXTGRAPH_TELOMERE.out.versions.first())

    } // End of HiFi+HiC assembly mode conditional block

 //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .mix(Channel.of(databaseVersionsToYAML()).filter { yaml -> yaml })
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

   // ch_multiqc_files                      = ch_multiqc_files.mix(BUSCO_GENERATEPLOT.out.png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(BUSCO_GENERATEPLOT_FINAL.out.png)
    ch_multiqc_files                      = ch_multiqc_files.mix(GENOMESCOPE2.out.linear_plot_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_fl_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_ln_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_cn_st_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_fl_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_ln_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(MERQURY.out.spectra_asm_st_png)
    //ch_multiqc_files                      = ch_multiqc_files.mix(GFASTATS_HAP1.out.assembly_summary)
    //ch_multiqc_files                      = ch_multiqc_files.mix(GFASTATS_HAP2.out.assembly_summary)
    //ch_multiqc_files                      = ch_multiqc_files.mix(GFASTATS_HAP1_FINAL.out.assembly_summary)
   // ch_multiqc_files                      = ch_multiqc_files.mix(GFASTATS_HAP2_FINAL.out.assembly_summary)
    ch_multiqc_wf_summary                 = ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    ch_multiqc_versions                   = ch_collated_versions
    ch_multiqc_method_desc                = ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false)

    ch_multiqc_files = ch_multiqc_files
        .groupTuple()
        .map {
            meta, files ->
                return [ meta, files[0] ]
        }



  //  ch_multiqc_files = ch_multiqc_files
  //      .groupTuple()
  //      .map {
  //          meta, files ->
  //              return [ meta, files[0], files[1], files[2], files[3], files[4], files[5], files[6], files[7], files[8], files[9], files[10], files[11], files[12], files[13], files[14] ]
   //     }

    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_multiqc_wf_summary.first(),
        ch_multiqc_versions.first(),
        ch_multiqc_method_desc.first()
   )


    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
