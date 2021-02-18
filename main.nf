nextflow.enable.dsl = 2

include { create_index } from './modules/index'

process build_design_species {

    container 'amancevice/pandas:1.2.1'
    publishDir path: "${params.outDir}/designs", pattern: '{single,paired}.csv', mode: 'copy'

    input:
        path script

    output:
	    path 'single.csv', optional: true
	    path 'paired.csv', optional: true
        path 'single_*.csv', optional: true, emit: single
        path 'paired_*.csv', optional: true, emit: paired
    
    script:
        """
        python3 ${params.script}/chipseq_input_match.py '${params.species}'
        """
}

process run_chipseq_genotoul {
    publishDir path: "${params.outDir}", mode: 'copy'
       
    input:
        tuple file(design), val(single)
    
    output:
        path 'results_*'
    
    when:
        params.indexExists
        
    script:
        """
        srun nextflow run nf-core/chipseq -profile genotoul \
             --input $design \
             $single \
             --bwa_index ${params.index}/${fasta.baseName}.fa \
             --skip_diff_analysis \
             --fasta ${params.fasta} \
             --gtf ${params.gtf} \
             --macs2_gsize ${params.gsize}
             --save_reference \
             --outdir 'results_${design.baseName}'\
             &
        """
}

index = file("${params.index}/${params.species}", type:"dir")

workflow {
    params.indexExists = index.exists()
    if (params.indexExists) {
        create_index()
    }
    
    build_design_species()
    all_chipseq = Channel.empty()
    all_chipseq = all_chipseq.concat(build_design_species.out.single.flatten().combine(Channel.from('--single_end')))
    all_chipseq = all_chipseq.concat(build_design_species.out.paired.flatten().combine(Channel.from('')))
    
    
    run_chipseq_genotoul(all_chipseq)
}

