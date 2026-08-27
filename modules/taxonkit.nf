process FetchTaxonkitDB {

    tag "${url}"

    publishDir( 
        "${params.outputs}/taxonkit", 
        mode: 'copy',
    )

    input:
    val url

    output:
    path '.taxonkit'

    script:
    """
    set -euxo pipefail
    
    curl -v "${url}" -o taxdump.tar.gz
    tar -zxvf taxdump.tar.gz

    mkdir .taxonkit
    mv names.dmp nodes.dmp delnodes.dmp merged.dmp .taxonkit
    rm *.{dmp,prt,txt}

    """
}


process FetchTaxonomicRanks {

    tag "${table}"

    publishDir( 
        "${params.outputs}/taxonkit", 
        mode: 'copy',
    )

    input:
    tuple path( table ), path( taxonkit_db )
    val column

    output:
    path "taxonomy.csv"

    script:
    """
    set -euox pipefail

    tax_id_col=\$(head -n1 "${table}" | tr , \$'\\n' | grep -nxF "${column}" | cut -d: -f1)

    head -n1 "${table}" \
    | awk -F, -v OFS=, -v col="\$tax_id_col" '
        { 
            print \$col, "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "strain"
        }
        ' \
    > "taxonomy.csv"

    tail -n+2 "${table}" \
    | cut -f"\$tax_id_col" -d, \
    | sort -u -n \
    | taxonkit reformat \
        --data-dir "${taxonkit_db}" \
        --threads "${task.cpus}" \
        -I 1 \
        --pseudo-strain \
        -f '{d},{K},{p},{c},{o},{f},{g},{s},{S},{t}' \
    | tr \$'\\t' , \
    >> "taxonomy.csv"

    """
}

process ResolveSpecies {
    tag "$id"

    input:
    tuple val( id ), path( taxonkit_db )

    output:
    tuple val( id ), path( "result.csv" ), emit: file
    tuple val( id ), stdout, emit: stdout

    script:
    """
    set -euox pipefail

    echo "${id}" \
    | taxonkit lineage -t \
        --data-dir "${taxonkit_db}" \
    > lineage.tsv

    awk -F'\\t' '{
        n=split(\$3, a, ";"); 
        for(i=1; i<=n; i++) print a[i]
    }' lineage.tsv \
    > lineage_ids.txt

    tail -n+2 lineage_ids.txt \
    | taxonkit filter \
        --equal-to species \
        --data-dir "${taxonkit_db}" \
    > species.txt

    species_id=\$(cat species.txt | cut -f1)

    if [ -z "\$species_id" ]; then
        >&2 echo "[WARN] No species ancestor for ${id}, using as-is"
        species_id="${id}"
    elif [ "\$species_id" != "${id}" ]; 
    then
        >&2 echo "Resolved ${id} = \$species_id" 
    fi

    echo "taxon_id,species_id" > result.csv
    echo "${id},\$species_id" >> result.csv
    printf "\$species_id"
    """
}
