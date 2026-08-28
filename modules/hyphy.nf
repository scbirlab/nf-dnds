process HyPhy {
    // SLAC: fast, tree-aware, per-site dN/dS with FDR.
    // Outputs JSON parsed downstream.
    // Alternative: swap --method slac for fel or fubar if you want more power.

    tag "${id}:${orthogroup_id}"

    // errorStrategy "ignore"

    publishDir( 
      "${params.outputs}/hyphy/slac", 
      mode: 'copy',
      saveAs: { "${id}_${orthogroup_id}_${it}"},
    )

    input:
    tuple val( id ), val( orthogroup_id ), path( codon_aln ), path( tree )
    val pvalue

    output:
    tuple val( id ), val( orthogroup_id ), path( "hyphy.json" )

    script:
    """
    hyphy slac \\
        --alignment "${codon_aln}" \\
        --tree "${tree}" \\
        --branches All \\
        --pvalue "${pvalue}" \\
        --output hyphy.json \\
        --ci Yes
    """
}

process ParseHyPy {
    // bin/parse_slac.py:
    //   extracts from SLAC JSON:
    //     og_id, n_sites, dN, dS, omega (gene-level averages),
    //     n_positively_selected_sites (p < hyphy_pvalue),
    //     n_negatively_selected_sites

    tag "$taxon_id"
    publishDir "${params.outdir}/dnds", mode: 'copy'

    input:
    tuple val( id ), val( orthogroup_id ), path(slac_jsons)

    output:
    tuple val( id ), val( orthogroup_id ), path ( "dnds.tsv" )

    script:
    """
    python3 ${projectDir}/bin/parse_slac.py \\
        --jsons  ${slac_jsons} \\
        --taxon  ${taxon_id} \\
        --pvalue ${params.hyphy_pvalue} \\
        --out    ${taxon_id}.dnds.tsv
    """
}
