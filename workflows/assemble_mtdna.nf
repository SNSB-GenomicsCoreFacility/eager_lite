/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { PYTHON3_CIRCULARIZE_MTDNA } from '../modules/local/python3/circularize_mtdna/main'
include { CIRCULARMAPPER_REALIGNSAMFILE } from '../modules/nf-core/circularmapper/realignsamfile/main'
include { ADAPTERREMOVAL           } from '../modules/local/adapterremoval/main'
include { ADAPTERREMOVALFIXPREFIX  } from '../modules/nf-core/adapterremovalfixprefix/main'
include { DEDUP                  } from '../modules/nf-core/dedup/main'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view/main'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { ANGSD_DOCOUNTS         } from '../modules/local/angsd/docounts/main'
include { BWA_ALN                } from '../modules/local/bwa/aln/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_MERGE         } from '../modules/nf-core/samtools/merge/main'
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

   if(params.create_circularmtdna){
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

   if(!params.bwamem_idx){
		//
		//MODULE: BWA_INDEX
		//
		BWA_INDEX(
		fa
		)
		m1_fa_m2_idx = fa.combine(BWA_INDEX.out.index)
	}else{
		index = Channel.fromPath(params.bwamem_idx)
		m1_fa_m2_idx = fa.combine(index)
	}
    prech_bwa_aln = reads.combine(m1_fa_m2_idx)
    //
    //MODULE:BWA_ALN
    //
    BWA_ALN(
        prech_bwa_aln.map{meta, fastq, meta_f, fa, meta_i, idx -> tuple(meta, fastq)},
        prech_bwa_aln.map{meta, fastq, meta_f, fa, meta_i, idx -> tuple(meta_i, idx)}
    )
    //
    //MODULE: SAMTOOLS_VIEW
    //
    SAMTOOLS_VIEW(
        BWA_ALN.out.bam.map{meta, bam -> tuple(meta,bam,[])},
        [[],[]],
        [],
        Channel.value("bai")
    )
    prech_realignsamfile = SAMTOOLS_VIEW.out.bam.combine(fa)
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
    // MODULE: PICARD_MARKDUPLICATES
    //
    PICARD_MARKDUPLICATES(
        SAMTOOLS_VIEW.out.bam,
        [[],[]],
        [[],[]]
    )

    bam_groups = PICARD_MARKDUPLICATES.out.bam.map{meta,bam -> tuple([id:meta.id],bam)}.groupTuple()

    // Split into two named branches
    bam_split = bam_groups.branch {
        to_merge: { id, bams -> bams.size() > 1 }
        singles:  { id, bams -> bams.size() == 1 }
    }

    // Now you can access them as separate channels:
    to_merge = bam_split.to_merge
    singles  = bam_split.singles.map { id, bams -> tuple(id, bams[0]) }

    //
    // MODULE: SAMTOOLS_MERGE
    //
    SAMTOOLS_MERGE(
        to_merge,
        [[],[]],
        [[],[]],
        [[],[]]
    )

    ch_angsd_docounts = SAMTOOLS_MERGE.out.bam.mix(singles)


    //
    // MODULE: ANGSD_DOCOUNTS
    //
    ANGSD_DOCOUNTS(
        ch_angsd_docounts.map{meta, bam ->tuple(meta, bam, [], [])}
    )
    //
    //QUALIMAP_BAMQC
    //
    QUALIMAP_BAMQC(
        ch_angsd_docounts,
        []
    )
    //
    //MODULE: AWK_REPORT_DOFASTA
    //
    AWK_REPORT_DOFASTA(
        ANGSD_DOCOUNTS.out.log.map{meta, fasta -> fasta}.collect()
    )
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]})
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
