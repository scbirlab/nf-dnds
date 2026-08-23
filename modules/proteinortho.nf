process ProteinOrtho {
    tag "$id"

    input:
    tuple val( id ), path( proteins_dir )

    output:
    tuple val( id ), path( "ortho.proteinortho.tsv" )

    script:
    """
    proteinortho -project=ortho -cpus=${task.cpus} ${proteins_dir}/*.faa
    """
}