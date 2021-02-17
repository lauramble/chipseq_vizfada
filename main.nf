nextflow.enable.dsl = 2

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
	pip3 install requests
        python3 $script '${params.species}'
        """
}

process run_chipseq_genotoul {
    publishDir path: "${params.outDir}", mode: 'copy'
       
    input:
        tuple file(design), val(single)
        path fasta
        path gtf
        path bwa_index
    
    output:
        path 'results_*'
        
    script:
        """
        srun nextflow run nf-core/chipseq -profile genotoul \
             --input $design \
             $single \
             --bwa_index ${bwa_index}/${fasta.baseName}.fa \
             --skip_diff_analysis \
             --fasta $fasta \
             --gtf $gtf \
             --save_reference \
             --outdir 'results_${design.baseName}'\
             &
        """
}

workflow {
    build_design_species("${baseDir}/scripts/chipseq_input_match.py")
    all_chipseq = Channel.empty()
    all_chipseq = all_chipseq.concat(build_design_species.out.single.flatten().combine(Channel.from('--single_end')))
    all_chipseq = all_chipseq.concat(build_design_species.out.paired.flatten().combine(Channel.from('')))
    run_chipseq_genotoul(all_chipseq, "${baseDir}/${params.fasta}", "${baseDir}/${params.gtf}", "${baseDir}/${params.bwa_index}")
}

