nextflow.enable.dsl = 2

/*
process index_bwa_mem {
	publishDir "${params.index}"
        
    label 'process_high'

	container 'nfcore/chipseq:1.2.1'
	
	stageInMode 'copy'

	input:
	path fasta
	path gtf

	output:
	path "$params.species", emit: index_chipseq
	path "${params.species}/${fasta.baseName}", emit: fasta_chipseq
	path "${params.species}/${gtf.baseName}", emit: gtf_chipseq

	script:
	"""
	gunzip *.gz
	bwa index -a bwtsw *.fa
	mkdir ${params.species} && mv ${fasta.baseName}* $params.species && mv *.gtf ${params.species}
	"""
}
*/

process unzip {

    input:
    path fasta
    path gtf
    
    output:
    path "*.fa", emits: fasta_chipseq
    path "*.gtf", emits: gtf_chipseq
    
    script:
    """
    gunzip *.gz
    """
    
}

process build_design_species {

    container 'amancevice/pandas:1.2.1'
    publishDir path: "${params.outDir}/designs", pattern: '{single,paired}.csv', mode: 'copy'

    label 'process_low'

    output:
	    path 'single.csv', optional: true
	    path 'paired.csv', optional: true
        path 'single_*.csv', optional: true, emit: single
        path 'paired_*.csv', optional: true, emit: paired
    
    script:
        """
        python3 ${params.scripts}/chipseq_input_match.py '${params.species}'
        """
}

process run_chipseq_genotoul {
    tag "$design"
    publishDir pattern: "results_*/genome", path: "${params.index}/${params.species}", mode: 'copy', enabled: params.makeIndex
    publishDir path: "${params.outDir}"
    
    label "process_high"
native
    input:
        tuple file(design), val(single)
        path index
        path fasta
        path gtf
    
    output:
        path 'results_*'
        
    script:
        def index = ${params.makeIndex} ? '--save-reference' : "--bwa_index $index/${file(fasta).baseName}"
        """
        srun nextflow run nf-core/chipseq -profile genotoul \
             --input $design \
             $single \
             $index \
             --fasta $fasta \native
             --gtf $gtf \
             --macs_gsize ${params.gsize} \
             --narrow_peaks \
             --outdir 'results_${design.baseName}'\
             --skip_consensus_peaks \
             --skip_diff_analysis \
             --skip_igv \
             &
        """
}

index = file("${params.index}/${params.species}", type:"dir")

workflow {

    //Parsing reference file
    awkCMD = ['awk', 'BEGIN{FS=\",\"; OFS=\"\\n\"} \$1 ~ species {print}', "species=$params.species", "$params.reference"]
    data = awkCMD.execute()
    data.waitFor()
    (species, fasta, gtf, gsize) = data.text.split(',')
    //params.fastaBase = file(fasta).baseName
    //params.gtfBase = file(gtf).baseName
    params.gsize = gsize
   
    if (!index.exists()) {
        params.makeIndex = true
    } else {
        params.makeIndex = false
    }
    
        
    unzip(fasta, gtf)
    build_design_species()
    
    all_chipseq = Channel.empty()
    all_chipseq = all_chipseq.concat(build_design_species.out.single.flatten().combine(Channel.from('--single_end')))
    all_chipseq = all_chipseq.concat(build_design_species.out.paired.flatten().combine(Channel.from('')))
    
    run_chipseq_genotoul(all_chipseq, "$index/${file(fasta).baseName}", unzip.out.fasta_chipseq, unzip.out.gtf_chipseq)
}

