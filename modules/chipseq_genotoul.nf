process CHIPSEQ_GENOTOUL {
    publishDir path: 'results'
    
    scratch true
    
    input:
        path design
        path fasta
        path gtf
        val singleOrPaired
    
    output:
        path 'results_*'
        
    script:
        """
        srun nextflow run nf-core/chipseq -profile genotoul \
             --input design \
             ${singleOrPaired} \
             --skip_diff_analysis \
             --fasta fasta \
             --gtf gtf \
             --save_reference \
             --outdir 'results_${design.baseName =~ /.*_(.*)\.txt/}'
             &
        """
}

