process INDEX_BWA_MEM {
	publishDir "${params.index}", mode: "copy"
        
  label 'process_high'

	container 'nfcore/chipseq:1.2.1'
	
	stageInMode 'copy'

	input:
	path fasta
	path gtf

	output:
	path "${params.species.replace(" ", '_')}", emit: index_chipseq
	path "${params.species.replace(" ", '_')}/${fasta.baseName}", emit: fasta_chipseq
	path "${params.species.replace(" ", '_')}/${gtf.baseName}", emit: gtf_chipseq

	script:
	def species = params.species.replace(" ", '_')
	"""
	gunzip *.gz
	bwa index -a bwtsw *.fa
	mkdir $species && mv ${fasta.baseName}* $species && mv *.gtf $species
	"""
}
