process Pal2Nal {
    tag "${id}:${orthogroup_id}"

    publishDir( 
      "${params.outputs}/alignments/cds", 
      mode: 'copy',
      saveAs: { "${id}_${orthogroup_id}_${it}"},
    )

    input:
    tuple val( id ), val( orthogroup_id ), path( prot_aln ), path( cds_fna )

    output:
    tuple val( id ), val( orthogroup_id ), path( "codon.fna" )

    script:
    """
    pal2nal.pl "${prot_aln}" "${cds_fna}" \\
        -output fasta \\
        -nogap \\
        -nomismatch \\
    > codon.fna
    """
}