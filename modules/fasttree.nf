process FastTree {
    // FastTree is substantially faster than IQ-TREE for genome-wide sweeps.
    // Tree accuracy doesn't need to be high — it's used only to parameterise
    // HyPhy branch lengths, not for topology inference.

    tag "${id}:${orthogroup_id}"

    input:
    tuple val( id ), val( orthogroup_id ), path( codon_aln )

    output:
    tuple val( id ), val( orthogroup_id ), path("tree.nwk")

    script:
    """
    FastTree -nt -gtr -quiet ${codon_aln} > tree.nwk
    """
}
