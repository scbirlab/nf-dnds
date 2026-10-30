process MAFFT {
    tag "${id}:${orthogroup_id}"
    label "big_cpu"

    publishDir( 
      "${params.outputs}/alignments/protein", 
      mode: 'copy',
      saveAs: { "${id}_${orthogroup_id}_${it}"},
    )

    input:
    tuple val( id ), val( orthogroup_id ), path( prot_faa )

    output:
    tuple val( id ), val( orthogroup_id ), path( "aln.faa" )

    script:
    """
    mafft --auto --thread "${task.cpus}" "${prot_faa}" > aln.faa
    """
}
