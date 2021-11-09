process download_reference {
    publishDir "${params.index}/${params.species}", pattern: "*.gtf"

	input:
	val reference

	output:
	path '*.fa', emit: fasta
	path '*.gtf'

	script:
	"""
	wget ${params.fasta}
	wget ${params.gtf}
	gunzip *.gz
	"""
}

process index_bwa_mem {
	publishDir "${params.index}/${params.species}"

	input:
	path fasta

	output:
	path "$fasta*"
	
	afterScript:
	params.indexExists = true

	script:
	"""
	bwa index -a bwtsw $fasta
	mkdir $params.species && mv ${params.fasta}* $params.species	
	"""
}

workflow create_index {
    main:
    
        //Parsing reference file
        awkCMD = ['awk', 'BEGIN{FS=\",\"; OFS=\"\\n\"} \$1 ~ species {print}', "species=$params.species", "$params.reference"]
        data = awkCMD.execute()
        data.waitFor()
        (species, params.fasta, params.gtf, params.gsize) = data.text.split(',')
        
	    download_reference()
	    index_bwa_mem(download_reference.out.fasta)

}

