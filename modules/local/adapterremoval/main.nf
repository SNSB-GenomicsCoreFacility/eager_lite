process ADAPTERREMOVAL {
    label 'process_medium'
    tag "${sample}_L${lib}"
    conda "${moduleDir}/environment.yml"
    //publishDir "${params.outdir}/adapterremoval", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("output/*{combined.fq,.se.truncated,pair1.truncated}.gz"),    emit: read1
    tuple val(meta), path("output/*pair2.truncated.gz"), emit: read2, optional: true
    tuple val(meta), path("output/*.settings"), emit: logs
    
    when: 
    !params.skip_adapterremoval

    script:
    sample = "${meta.id}"
    lib = "${meta.lib}"
    def base = "${reads[0].baseName}_${meta.lib}"
    def single_end = "${meta.single_end}"
    def adapters_to_remove = !params.clip_adapters_list ? "--adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor}" : "--adapter-list ${adapterlist}"
    //This checks whether we skip trimming and defines a variable respectively
    def preserve5p = params.preserve5p ? '--preserve5p' : '' // applies to any AR command - doesn't affect output file combination
    r1 = reads[0]
    if(single_end == false){
        r2=reads[1]
    }
    
    if ( single_end == false   && !params.skip_collapse && !params.skip_trim  && !params.mergedonly && !params.preserve5p ) {
    """
    mkdir -p output

    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}

    cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
    
    mv *.settings output/

    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix -Xmx${task.memory.toGiga()}g output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus - 1} > output/${base}.pe.combined.fq.gz
    
    """
    //PE mode, collapse and trim, outputting all reads, preserving 5p
    } else if ( single_end == false && !params.skip_collapse && !params.skip_trim  && !params.mergedonly && params.preserve5p) {
    """
    mkdir -p output

    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}


    cat *.collapsed.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz

    mv *.settings output/

    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix -Xmx${task.memory.toGiga()}g output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus - 1} > output/${base}.pe.combined.fq.gz

    """
    // PE mode, collapse and trim but only output collapsed reads
    } else if ( single_end == false && !params.skip_collapse && !params.skip_trim && params.mergedonly && !params.preserve5p ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe  --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}

    
    cat *.collapsed.gz *.collapsed.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
        
    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix -Xmx${task.memory.toGiga()}g output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus - 1} > output/${base}.pe.combined.fq.gz

    mv *.settings output/
    """
    // PE mode, collapse and trim but only output collapsed reads, preserving 5p
    } else if ( single_end == false && !params.skip_collapse && !params.skip_trim && params.mergedonly && params.preserve5p ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe  --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}

    
    cat *.collapsed.gz > output/${base}.pe.combined.tmp.fq.gz
    
    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix -Xmx${task.memory.toGiga()}g output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus - 1} > output/${base}.pe.combined.fq.gz

    mv *.settings output/
    """
    // PE mode, collapsing but skip trim, (output all reads). Note: seems to still generate `truncated` files for some reason, so merging for safety.
    // Will still do default AR length filtering I guess
    } else if ( single_end == false   && !params.skip_collapse && params.skip_trim && !params.mergedonly ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} --collapse ${preserve5p} --adapter1 "" --adapter2 ""

    
    cat *.collapsed.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
        
    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix -Xmx${task.memory.toGiga()}g output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus - 1} > output/${base}.pe.combined.fq.gz

    mv *.settings output/
    """
    // PE mode, collapsing but skip trim, and only output collapsed reads. Note: seems to still generate `truncated` files for some reason, so merging for safety.
    // Will still do default AR length filtering I guess
    } else if ( single_end == false   && !params.skip_collapse && params.skip_trim && params.mergedonly ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} --collapse ${preserve5p}  --adapter1 "" --adapter2 ""
    
    cat *.collapsed.gz > output/${base}.pe.combined.tmp.fq.gz
    
    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix -Xmx${task.memory.toGiga()}g output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus - 1} > output/${base}.pe.combined.fq.gz

    mv *.settings output/
    """
    // PE mode, skip collapsing but trim (output all reads, as merging not possible) - activates paired-end mapping!
    } else if ( single_end == false   && params.skip_collapse && !params.skip_trim ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}
    
    mv ${base}.pe.pair*.truncated.gz *.settings output/
    """
    } else if ( single_end == true  && !params.skip_trim ) {
    //SE, collapse not possible, trim reads only
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --basename ${base}.se --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} ${preserve5p} --trimns --trimqualities ${adapters_to_remove} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}
    mv *.settings *.se.truncated.gz output/
    """
    } else if ( single_end == true  && params.skip_trim ) {
    //SE, collapse not possible, trim reads only
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --basename ${base}.se --gzip --threads ${task.cpus} --qualitymax ${params.qualitymax} ${preserve5p} --adapter1 "" --adapter2 ""
    mv *.settings *.se.truncated.gz output/
    """
    }
}
