process FastTree {

    tag "${id}:${orthogroup_id}"

    publishDir( 
      "${params.outputs}/alignments/trees", 
      mode: 'copy',
      saveAs: { "${id}_${orthogroup_id}_${it}"},
    )

    input:
    tuple val( id ), val( orthogroup_id ), path( codon_aln )

    output:
    tuple val( id ), val( orthogroup_id ), path( "tree.nwk" )

    script:
    """
    FastTree -nt -gtr -quiet "${codon_aln}" > tree.nwk
    """
}
