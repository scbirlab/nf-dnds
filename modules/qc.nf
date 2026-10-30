process AssessCodonAlignment {

    tag "${id}:${orthogroup_id}"

    publishDir(
        "${params.outputs}/qc/codon",
        mode: 'copy',
        saveAs: { "${id}_${orthogroup_id}_${it}" },
    )

    input:
    tuple val(id), val(orthogroup_id), path(codon_aln)

    val min_sequences
    val min_unique_sequences
    val min_codons
    val min_variable_codons
    val max_ambiguous_fraction

    output:
    tuple val(id),
          val(orthogroup_id),
          path("analysable.fna"),
          emit: analysable,
          optional: true

    tuple val(id),
          val(orthogroup_id),
          path("assessment.tsv"),
          emit: assessment

    tuple val(id),
          val(orthogroup_id),
          path("sequence_groups.tsv"),
          emit: groups

    script:
    """
    python3 ${projectDir}/bin/assess_codon_alignment.py \
        --alignment "${codon_aln}" \
        --taxon-id "${id}" \
        --orthogroup-id "${orthogroup_id}" \
        --assessment assessment.tsv \
        --groups sequence_groups.tsv \
        --analysable analysable.fna \
        --min-sequences "${min_sequences}" \
        --min-unique-sequences "${min_unique_sequences}" \
        --min-codons "${min_codons}" \
        --min-variable-codons "${min_variable_codons}" \
        --max-ambiguous-fraction "${max_ambiguous_fraction}"
    """
}