nextflow.enable.dsl = 2

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
	def species = params.species.replace(" ", '_')
	"""
	gunzip *.gz
	bwa index -a bwtsw *.fa
	mkdir $species && mv ${fasta.baseName}* $species && mv *.gtf $species
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
    publishDir path: "${params.outDir}", mode: 'copy'

    
    label "process_high"

    input:
        tuple file(design), val(single)
        path index
        path fasta
        path gtf
    
    output:
        path 'results_*'
        
    script:
        """
        srun nextflow run nf-core/chipseq -profile genotoul \
             --input $design \
             $single \
             --bwa_index ${index}/${fasta.baseName} \
             --skip_diff_analysis \
             --fasta $fasta \
             --gtf $gtf \
             --macs_gsize ${params.gsize} \
             --outdir 'results_${design.baseName}'\
             --narrow_peaks \
             --skip_consensus_peaks \
             --skip_igv \
             &
        """
}

workflow {

    //Parsing reference file
    awkCMD = ['awk', 'BEGIN{FS=\",\"; OFS=\"\\n\"} \$1 ~ species {print}', "species=$params.species", "$params.reference"]
    data = awkCMD.execute()
    data.waitFor()
    (species, fasta, gtf, gsize) = data.text.split(',')
    //params.fastaBase = file(fasta).baseName
    //params.gtfBase = file(gtf).baseName
    params.gsize = gsize

    index = file("${params.index}/${params.species.replace(" ", '_')}", type:"dir")
   
    if (!index.exists()) {
        index_bwa_mem(fasta, gtf)
        index_chipseq = index_bwa_mem.out.index_chipseq
        fasta_chipseq = index_bwa_mem.out.fasta_chipseq
        gtf_chipseq = index_bwa_mem.out.gtf_chipseq
    } else {
        index_chipseq = index
        fasta_chipseq = "$index/${file(fasta).baseName}"
        gtf_chipseq = "$index/${file(gtf).baseName}"
    }
    
    build_design_species()
    all_chipseq = Channel.empty()
    all_chipseq = all_chipseq.concat(build_design_species.out.single.flatten().combine(Channel.from('--single_end')))
    all_chipseq = all_chipseq.concat(build_design_species.out.paired.flatten().combine(Channel.from('')))
    
    
    run_chipseq_genotoul(all_chipseq, index_chipseq, fasta_chipseq, gtf_chipseq)
}

