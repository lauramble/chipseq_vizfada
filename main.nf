nextflow.enable.dsl = 2

include { CHIPSEQ_VIZFADA } from './workflows/chipseq_vizfada'

workflow {
	CHIPSEQ_VIZFADA()
}
