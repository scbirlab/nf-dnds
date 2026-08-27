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
    > codon.fna \\
    2> pal2nal.err

    cat pal2nal.err >&2

    if [ ! -s codon.fna ]
    then
        >&2 echo "[ERROR] PAL2NAL produced an empty alignment"
        exit 1
    fi

    if ! grep -q '^>' codon.fna
    then
        >&2 echo "[ERROR] PAL2NAL output is not FASTA"
        exit 1
    fi
    """
}