process ProteinOrtho {
    tag "$id"
    cpus 32

    publishDir( 
      "${params.outputs}/alignment", 
      mode: 'copy',
      saveAs: { "${id}_${it}"},
   )

    input:
    tuple val( id ), path( proteins, stageAs: 'protein-????.faa' )

    output:
    tuple val( id ), path( "ortho.proteinortho.tsv" ), emit: tsv
    tuple val( id ), path( "*.html" ), emit: html

    script:
    """
    set -euox pipefail
    FILES=(${proteins})

    for f in "\${FILES[@]}"
    do
      sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^\$/d' \$f
    done
    proteinortho -project=ortho -cpus=${task.cpus} \${FILES[@]}
    """
}