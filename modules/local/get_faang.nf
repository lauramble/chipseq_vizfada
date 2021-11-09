process GET_FAANG {
  publishDir "${params.outdir}/metadata",
      mode: params.publish_dir_mode,
      pattern: "*.tsv"

  container "lauramble/python-vizfada"
    
  output:
  path "paired_*.csv", emit: paired, optional: true
  path "single_*.csv", emit: single, optional: true
  path "*.tsv", emit: metadata
  
  script:
  def species = "${params.species}".capitalize().replace("_", " ")
  if (params.ids) {
    def idsFile = file(params.ids)
    """
    chipseq_input_match.py "$species" $idsFile
    """
  } else {
    """
    chipseq_input_match.py "$species"
    """
  }
}
