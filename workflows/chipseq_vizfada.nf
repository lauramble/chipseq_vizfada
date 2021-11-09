nextflow.enable.dsl = 2

/*
process build_design_species {

    container 'amancevice/pandas:1.2.1'
    publishDir path: "${params.outDir}/designs", pattern: '*.csv', mode: 'copy'
    publishDir path: "${params.outDir}/metadata", pattern: '*.tsv', mode: 'copy'

    label 'process_low'

    output:
		path 'single.csv', optional: true
		path 'paired.csv', optional: true
    path 'single_*.csv', optional: true, emit: single
    path 'paired_*.csv', optional: true, emit: paired
		path 'metadata.tsv'

    script:
    """
    python3 ${params.scripts}/chipseq_input_match.py '${params.species}'
    """
}
*/
process run_chipseq_genotoul {
    tag "$expID"
    publishDir path: "${params.outDir}/${expID}", mode: 'move', saveAs: {
      filename -> {
        if (filename.endsWith("bigwig")) "bigwig"
        else if (filename.endsWith("macs")) "macs"
        else filename
      }
    }
    errorStrategy 'retry'
    maxRetries 3
    maxForks 3

    label "process_high"
    label "process_long"

    input:
        tuple file(design), val(single)
        path index
        path fasta
        path gtf
    
    output:
        path 'results_*/bwa/mergedLibrary/bigwig'
        path 'results_*/bwa/mergedLibrary/macs'
        path 'results_*/pipeline_info'
	      path 'results_*/multiqc'
  
    afterScript "rm -rf work"
        
    script:
      matches = "${design.baseName}" =~ /.*_(.*)_samplesheet/
      expID = matches[0][1]
        """
        srun --mem=8G nextflow run nf-core/chipseq -profile genotoul \
             --input $design \
             $single \
             -c conf/chipseq.config \
             --bwa_index ${index}/${fasta} \
             --skip_diff_analysis \
             --skip_trimming \
             --fasta $fasta \
             --gtf $gtf \
             --macs_gsize ${params.macs_gsize} \
             --outdir 'results_${expID}'\
             --narrow_peak \
             --skip_consensus_peaks \
             --skip_igv \
             --skip_spp
        """
}

include { FETCHNGS as FETCHNGS_SINGLE; FETCHNGS as FETCHNGS_PAIRED } from '../subworkflows/local/fetchngs' addParams(outdir: "${params.outdir}/fetchngs")
include { INDEX_BWA_MEM } from "../modules/local/index_bwa_mem"
include { GET_FAANG } from '../modules/local/get_faang'

workflow CHIPSEQ_VIZFADA {

    /*
    //Parsing reference file
    awkCMD = ['awk', 'BEGIN{FS=\",\"; OFS=\"\\n\"} \$1 ~ species {print}', "species=$params.species", "$params.reference"]
    data = awkCMD.execute()
    data.waitFor()
    (species, fasta, gtf, gsize) = data.text.split(',')
    //params.fastaBase = file(fasta).baseName
    //params.gtfBase = file(gtf).baseName
    params.gsize = gsize.replace("\n", "")
    */
    
    index = file("${params.index}/${params.species.replace(" ", '_')}", type:"dir")
    
    if (!index.exists()) {
        INDEX_BWA_MEM(file(params.fasta), file(params.gtf))
        index_chipseq = INDEX_BWA_MEM.out.index_chipseq
        fasta_chipseq = INDEX_BWA_MEM.out.fasta_chipseq
        gtf_chipseq = INDEX_BWA_MEM.out.gtf_chipseq
    } else {
        index_chipseq = index
        fasta_chipseq = "$index/${file(params.fasta).baseName}"
        gtf_chipseq = "$index/${file(params.gtf).baseName}"
    }
    
		GET_FAANG()
    GET_FAANG.out.single.view()
    GET_FAANG.out.paired.view()
    FETCHNGS_SINGLE( GET_FAANG.out.single.flatten() )
		FETCHNGS_PAIRED( GET_FAANG.out.paired.flatten() )
		
    Channel
        .empty()
        .concat(FETCHNGS_SINGLE.out.samplesheet.flatten().combine(Channel.from('--single_end')))
        .concat(FETCHNGS_PAIRED.out.samplesheet.flatten().combine(Channel.from('')))
        .multiMap{it -> run: print: it}
        .set{ all_chipseq }
    all_chipseq.print.view{"$it"}
    
    run_chipseq_genotoul(all_chipseq.run, index_chipseq, fasta_chipseq, gtf_chipseq)
}
