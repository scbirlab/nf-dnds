process MAFFT {
    tag "${id}:${orthogroup_id}"
    cpus 2

    input:
    tuple val( id ), val( orthogroup_id ), path( prot_faa ), path( cds_fna )

    output:
    tuple val( id ), val( orthogroup_id ), path( "aln.faa" ), path( cds_fna )

    script:
    """
    mafft --auto --thread "${task.cpus}" "${prot_faa}" > aln.faa
    """
}
