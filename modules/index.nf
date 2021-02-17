process get_reference {
	input:
	val species
	path reference

	output:
	env fasta, emit: fasta
	env gsize, emit: gsize

	shell:
	'''
	fasta=$(awk 'BEGIN{FS=","} $1 ~ species {print $2}' species="!{species}" !{reference})
	gsize=$(awk 'BEGIN{FS=","} $1 ~ species {print $3}' species="!{species}" !{reference})
	'''

}


process download_reference {
	input:
	val reference

	output:
	path '*.fa'

	script:
	"""
	wget $reference
	gunzip *.gz
	"""
}

process index_bwa_mem {
	publishDir "${params.index}"

	input:
	path fasta

	output:
	path "${species}"

	script:
	"""
	bwa index -a bwtsw $fasta
	mkdir $species && mv ${fasta}* $species	
	"""
}

workflow create_index {
	take: species, reference
	main:
	get_reference(species, reference)
	download_reference(get_reference.out.fasta)
	index_bwa_mem(download_reference.out)
}

