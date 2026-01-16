/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { PYTHON3_CIRCULARIZE_MTDNA } from '../modules/local/python3/circularize_mtdna/main'
include { CIRCULARMAPPER_REALIGNSAMFILE } from '../modules/nf-core/circularmapper/realignsamfile/main'
include { ADAPTERREMOVAL           } from '../modules/local/adapterremoval/main'
include { SEQKIT_SPLIT2          } from '../modules/nf-core/seqkit/split2/main'
include { DEDUP                  } from '../modules/nf-core/dedup/main'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FLAGSTAT as RAWFLAGSTAT } from '../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_FLAGSTAT as FILTEREDFLAGSTAT } from '../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_FLAGSTAT as DEDUPFLAGSTAT } from '../modules/nf-core/samtools/flagstat'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { ANGSD_DOCOUNTS         } from '../modules/local/angsd/docounts/main'
include { BWA_ALN                } from '../modules/local/bwa/aln/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BAMUTIL_TRIMBAM        } from '../modules/nf-core/bamutil/trimbam/main'
include { SAMTOOLS_MERGE as MERGELIB_SAMTOOLS         } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MERGE as MERGESAMPLE_SAMTOOLS         } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_DEPTH         } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_SORT as SORTMERGEDLIB_SAMTOOLS          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SORTMERGEDSAMPLE_SAMTOOLS          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SORTDEDUP_SAMTOOLS          } from '../modules/nf-core/samtools/sort/main'
include { ENDORSPY               } from "../modules/nf-core/endorspy/main"
include { QUALIMAP_BAMQC         } from '../modules/nf-core/qualimap/bamqc/main'
include { AWK_REPORT_DOFASTA     } from '../modules/local/awk/report_dofasta/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_assemble_mtdna_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLE_MTDNA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    
    //ch_samplesheet.view()
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    reference_fasta = Channel.fromPath(params.fasta)
    prech_reference_fasta = reference_fasta.map{fasta->tuple([id:"reference"],fasta)}
    
    if(!params.skip_adapterremoval){
            ADAPTERREMOVAL(
                ch_samplesheet
            )
            reads = ADAPTERREMOVAL.out.read1
        }
    else{
            reads = ch_samplesheet
        }


    SEQKIT_SPLIT2(
        reads
    )
    

    split_reads = SEQKIT_SPLIT2.out.reads
    .flatMap { meta, fileList ->
        // The spread operator * is implied in the channel definition above,
        // grouping all files into 'fileList'.

        // Ensure fileList is always a List (Nextflow does this automatically for multiple files,
        // but it's good practice to ensure for single-file cases, too, though often not needed here)
        def files = fileList instanceof List ? fileList : [fileList]

        // Transform [meta, file1, file2, ...] into [[meta, file1], [meta, file2], ...]
        return files.collect { file ->
            [meta, file]
        }
    }
   if(params.create_circularmtdna && params.mtdna_consense == true){
		//
		// MODULE: PYTHON3_CIRCULARIZE_MTDNA
		//
		PYTHON3_CIRCULARIZE_MTDNA(
		    prech_reference_fasta
		)
		fa = PYTHON3_CIRCULARIZE_MTDNA.out.fa
	}else{
		fa = prech_reference_fasta
	}

   if(!params.bwa_idx){
		//
		//MODULE: BWA_INDEX
		//
		BWA_INDEX(
		    fa
		)
		m1_fa_m2_idx = fa.combine(BWA_INDEX.out.index)
	}else{
		index = Channel.fromPath(params.bwa_idx)
        meta_idx = index.map{idx_path -> tuple([id:"reference"],idx_path)}
		m1_fa_m2_idx = fa.combine(meta_idx)
	}
    prech_bwa_aln = split_reads.combine(m1_fa_m2_idx)
    //
    //MODULE:BWA_ALN
    //
    BWA_ALN(
        prech_bwa_aln.map{meta, fastq, meta_f, fa, meta_i, idx -> tuple(meta, fastq)},
        prech_bwa_aln.map{meta, fastq, meta_f, fa, meta_i, idx -> tuple(meta_i, idx)}
    )
    //
    //MODULE: SAMTOOLS_MERGE
    //
    grp_bwa_aln_out = BWA_ALN.out.bam.groupTuple()

    br_grp_bwa_aln_out = grp_bwa_aln_out.branch { id, bams ->
        to_merge: bams.size() > 1 
        singles:  bams.size() == 1 
        } // only merged the fastq files from the same library at this stage

    MERGELIB_SAMTOOLS(
        br_grp_bwa_aln_out.to_merge.map{meta, bam -> 
            def new_meta = meta + [id:"${meta.id}.${meta.lib}.merged"]
            tuple(new_meta, bam)
            },
        [[],[]],
        [[],[]],
        [[],[]]
    )
    //
    //MODULE: SAMTOOLS_SORT
    //
    SORTMERGEDLIB_SAMTOOLS(
        MERGELIB_SAMTOOLS.out.bam.map{meta, bam -> 
            def new_meta = meta + [id:"${meta.id}.sorted"]
            tuple(new_meta, bam)}.groupTuple(),
        [[],[]],
        Channel.value("bai")
    )


    ch_combined_raw_bam =  SORTMERGEDLIB_SAMTOOLS.out.bam.map { meta, bam ->
                def cleaned_id = meta.id.replaceFirst("\\.${meta.lib}\\.merged\\.sorted\$", "")
                def new_meta =  meta + [ id: cleaned_id ]
                tuple(new_meta, bam)
            }


    //br_grp_bwa_aln_out.singles.map{meta, bam -> tuple([id:meta.id],bam[0])}.view()
        
    ///// this is exclusively to get the statistics of the raw bam files
    /////
    /////
        ch_merge_raw_bams = ch_combined_raw_bam
            .map { meta, bam -> 
                tuple([id:meta.id], bam) 
            }
            .join(br_grp_bwa_aln_out.singles.map{ meta, bam ->
                tuple([id:meta.id], bam[0])
            }, remainder: true)
            .map { it.findAll { it != null } } // This line removes the nulls

        ch_merge_raw_bams = ch_merge_raw_bams.map { it ->
            // it[0] is the meta map
            // it[1..-1] is a list of all remaining elements (the BAMs)
            [ it[0], it[1..-1] ]
        }


        ch_merge_raw_bams = ch_merge_raw_bams.branch { id, bams ->
            to_merge: bams.size() > 1 
            singles:  bams.size() == 1 
            }

            //
            // MODULE: MERGE_RAW_BAMS_SAMPLE
            //

            MERGESAMPLE_SAMTOOLS(
                ch_merge_raw_bams.to_merge.map{meta, bam -> 
                    def new_meta = meta + [id:"${meta.id}.raw"]
                    tuple(new_meta, bam)
                    },
                [[],[]],
                [[],[]],
                [[],[]]
            )

            //
            // SORT_RAW_BAMS_SAMPLE
            //

            SORTMERGEDSAMPLE_SAMTOOLS(
                MERGESAMPLE_SAMTOOLS.out.bam.map{meta, bam -> 
                    def new_meta = meta + [id:"${meta.id}.sorted"]
                    tuple(new_meta, bam)}.groupTuple(),
                [[],[]],
                Channel.value("bai")
            )

            ch_raw_sample_flagstat = SORTMERGEDSAMPLE_SAMTOOLS.out.bam.mix(ch_merge_raw_bams.singles)
            
            // 
            // RAW_SAMPLES_FLAGSTAT
            //

            RAWFLAGSTAT(
              ch_raw_sample_flagstat.map{meta, bam -> tuple(meta, bam, [])} 
            )
            //
            //MODULE: SAMTOOLS_VIEW
            //
            SAMTOOLS_VIEW(
                ch_raw_sample_flagstat.map{meta, bam -> tuple(meta,bam,[])},
                [[],[]],
                [],
                Channel.value("bai")
            )

            //
            // MODULE: FILTERED_FLAGSTAT
            //
            FILTEREDFLAGSTAT(
                SAMTOOLS_VIEW.out.bam.map{meta, bam -> tuple(meta, bam, [])}
            )
            SAMTOOLS_VIEW.out.bam.view()

            if (params.dedup_tool == "picard" ){
                    PICARD_MARKDUPLICATES(
                        SAMTOOLS_VIEW.out.bam,
                        [[],[]],
                        [[],[]]
                    )

                ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})

                }
                ch_dedup_flagstat = PICARD_MARKDUPLICATES.out.bam
            if (params.dedup_tool == "dedup" ){
                    DEDUP(
                        SAMTOOLS_VIEW.out.bam
                    )
                    ///
                    /// MODULE: SORT_DEDUP
                    ///
                    SORTDEDUP_SAMTOOLS(
                        DEDUP.out.bam.map{meta, bam -> 
                            def new_meta = meta + [id:"${meta.id}.sorted"]
                            tuple(new_meta, bam)}.groupTuple(),
                        [[],[]],
                        Channel.value("bai")
                    )
                ch_dedup_flagstat = SORTDEDUP_SAMTOOLS.out.bam
                }
            ///
            /// MODULE: DEDUPFLAGSTAT
            ///
            DEDUPFLAGSTAT(
                ch_dedup_flagstat.map{meta, bam -> tuple(meta, bam, [])}
            )
            
            //
            // MODULE: ENDORSPY
            //
            ENDORSPY(
                RAWFLAGSTAT.out.flagstat.combine(FILTEREDFLAGSTAT.out.flagstat,by:0).combine(DEDUPFLAGSTAT.out.flagstat,by:0)
            )
            ch_multiqc_files = ch_multiqc_files.mix(ENDORSPY.out.json.collect{it[1]})

            //
            // MODULE: SAMTOOLS_DEPTH
            //
            SAMTOOLS_DEPTH(
                ch_dedup_flagstat,
                [[],[]]
            )

    /*
    /////
    /////
    /////
    
    prepare_ch_dup = ch_combined_raw_bam.map{meta,bam -> tuple(meta, [bam])}.mix(br_grp_bwa_aln_out.singles)

    prepare_ch_dup.collect().view()

    //prepare_ch_dup.view()

    prepare_ch_dup = prepare_ch_dup.branch{ meta, bam -> 
                picard_b: meta.single_end == true
                dedup_b: meta.single_end == false
        }

    //
    // MODULE: PICARD_MARKDUPLICATES
    //
    PICARD_MARKDUPLICATES(
        prepare_ch_dup.picard_b,
        [[],[]],
        [[],[]]
    )
    //
    // MODULE: DEDUP
    //
    DEDUP(
        prepare_ch_dup.dedup_b
    )
    
        ///
        /// MODULE: SORT_DEDUP
        ///
            SORT_DEDUP(
                DEDUP.out.bam.map{meta, bam -> 
                    def new_meta = meta + [id:"${meta.id}.sorted"]
                    tuple(new_meta, bam)}.groupTuple(),
                [[],[]],
                Channel.value("bai")
            )


    // MODULE: BAMUTIL_TRIMBAM
    //
    //
    ch_pre_bamtrim = PICARD_MARKDUPLICATES.out.bam.mix(SORT_DEDUP.out.bam).branch{meta, bam -> 
                udg: meta.udg == true
                non_udg: meta.udg == false
        }

    ch_bamtrim = ch_pre_bamtrim.non_udg.combine([params.trim_left]).combine([params.trim_right])

    BAMUTIL_TRIMBAM(
        ch_bamtrim
    )

    //
    // MODULE: MERGE_SAMPLE
    //
    ch_merge_sample = BAMUTIL_TRIMBAM.out.bam.map{meta, bam -> tuple([id:meta.id], bam)}.mix(ch_pre_bamtrim.udg.map{meta, bam -> tuple([id:meta.id], bam)}).groupTuple()

    ch_merge_sample = ch_merge_sample.branch{ meta, bams ->
            to_merge: bams.size() > 1 
            singles:  bams.size() == 1 
            }

    MERGE_SAMPLE(
        ch_merge_sample.to_merge.map{meta, bam -> 
            def new_meta = meta + [id:"${meta.id}.dedup"]
            tuple(new_meta, bam)
            },
        [[],[]],
        [[],[]],
        [[],[]]
    )

    //
    // MODULE: SORT_DEDUP_MERGED
    //

    SORT_DEDUP_MERGED(
        MERGE_SAMPLE.out.bam.map{meta, bam -> 
            def new_meta = meta + [id:"${meta.id}.sorted"]
            tuple(new_meta, bam)}.groupTuple(),
        [[],[]],
        Channel.value("bai")
    )

    //
    // DEDUP_FLAGSTAT
    //
    
    ch_dedup_flagstat = SORT_DEDUP_MERGED.out.bam.mix(ch_merge_sample.singles)

    SAMTOOLS_FLAGSTAT(
        ch_dedup_flagstat.map{meta, bam-> tuple(meta, bam, [])}
    )
    
    //
    //MODULE: SAMTOOLS_VIEW
    //
    SAMTOOLS_VIEW(
        ch_dedup_flagstat.map{meta, bam -> tuple(meta,bam,[])},
        [[],[]],
        [],
        Channel.value("bai")
    )

    //
    // MODULE: FILTERED_FLAGSTAT
    //
    FILTERED_FLAGSTAT(
        SAMTOOLS_VIEW.out.bam.map{meta, bam -> tuple(meta, bam, [])}
    )

    //
    // MODULE: ENDORSPY
    //
    RAW_SAMPLES_FLAGSTAT.out.flagstat.view()
    SAMTOOLS_FLAGSTAT.out.flagstat.view()
    FILTERED_FLAGSTAT.out.flagstat.view()


    if(params.mtdna_consense == true ){
        prech_realignsamfile = SAMTOOLS_SORT.out.bam.combine(fa)
        //
        //
        //MODULE:CIRCULARMAPPER_REALIGNSAMFILE
        //
        CIRCULARMAPPER_REALIGNSAMFILE(
            prech_realignsamfile.map{meta, bam, meta_r, fa -> tuple(meta, bam)},
            prech_realignsamfile.map{meta, bam, meta_r, fa -> tuple(meta_r, fa)},
            [[],params.circle_nbp],
            [[],[]]
        )
        bam_to_dup = CIRCULARMAPPER_REALIGNSAMFILE.out.bam
    }
    else{
            bam_to_dup = SAMTOOLS_SORT.out.bam
        }
    

    


    */


    if(params.mtdna_consense == true ){
        ch_angsd_docounts = SAMTOOLS_MERGE.out.bam.mix(singles)
        //
        // MODULE: ANGSD_DOCOUNTS
        //
        ANGSD_DOCOUNTS(
            ch_angsd_docounts.map{meta, bam ->tuple(meta, bam, [], [])}
        )
        //
        //MODULE: AWK_REPORT_DOFASTA
        //
        AWK_REPORT_DOFASTA(
            ANGSD_DOCOUNTS.out.log.map{meta, fasta -> fasta}.collect()
        )
        //
        //QUALIMAP_BAMQC
        //
        QUALIMAP_BAMQC(
            ch_angsd_docounts,
            []
        )
    }
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]})
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'assemble_mtdna_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )


    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
