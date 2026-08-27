process ProteinOrtho {
    tag "${id}"

    publishDir( 
      "${params.outputs}/orthogroups", 
      mode: 'copy',
      saveAs: { "${id}_${it}"},
   )

    input:
    tuple val( id ), path( proteins )

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

process ExtractOrthogroups {
    tag "${id}"

    publishDir( 
      "${params.outputs}/orthogroups/proteins", 
      mode: 'copy',
      saveAs: { "${id}_${it}"},
   )

    input:
    tuple val( id ), path( orthogroups ), path( proteins ), path( cds )

    output:
    tuple val( id ), path( "orthogroup-*.faa" ), emit: protein
    tuple val( id ), path( "orthogroup-*.fna" ), emit: cds

    script:
    """
    #!/usr/bin/env python

    from itertools import chain

    from carabiner import print_err
    import pandas as pd

    df = pd.read_csv(
      "${orthogroups}",
      sep="\\t",
    )
    print_err(f"{df=}")

    input_fastas = [
      col.removesuffix("-protein.faa") for col in df.columns.tolist()[3:]
    ]
    print_err(f"{input_fastas=}")


    for i, row in enumerate(df.itertuples(index=False, name=None)):
      new_cds_ids, new_protein_ids = set(), set()
      print_err(">>>", f"{row=}")
      n_species, n_genes,	alg_conn,	*protein_ids = row
      protein_ids = [
          [] if pd.isna(ids)
          else [
            x.strip() 
            for x in str(ids).split(",")
          ] for ids in protein_ids
      ]
      
      print_err(f"{protein_ids=}")
      assert isinstance(protein_ids[0], (list, tuple))

      orthogroup_name = f"OG{i:04d}-" + str(next(id[0] for id in protein_ids if len(id) > 0 and id[0] != "*"))
      print_err(f"{orthogroup_name=}")
      these_orthologs = dict(zip(input_fastas, protein_ids))

      _line = ">{genome_id}-{seqname}\\n"

      with open(f"orthogroup-{orthogroup_name}.faa", "w") as out:
        for genome_id, _proteins in these_orthologs.items():
          seq_names_to_keep = {name: f">{name} " for name in _proteins}
          print_err(f"{seq_names_to_keep=}")
          with open(f"{genome_id}-protein.faa", "r") as f:
            save = False
            for line in f:
              if line.startswith(tuple(seq_names_to_keep.values())):
                save = True
                seqname = next(
                  name for name, string in seq_names_to_keep.items() 
                  if line.startswith(string)
                )
                line = _line.format(genome_id=genome_id, seqname=seqname)
                new_protein_ids.add(line)
              elif line.startswith((">", "\\n", "*")):
                save = False

              if save:
                print(line, file=out, end="")
              else:
                continue

      with open(f"orthogroup-{orthogroup_name}.fna", "w") as out:
        for genome_id, _proteins in these_orthologs.items():
          seq_names_to_keep = {name: f"[protein_id={name}]" for name in _proteins}
          print_err(f"{seq_names_to_keep=}")
          with open(f"{genome_id}-cds_from_genomic.fna", "r") as f:
            save = False
            for line in f:
              if line.startswith(">") and any(name in line for name in seq_names_to_keep.values()):
                save = True
                seqname = next(
                  name for name, string in seq_names_to_keep.items() 
                  if string in line
                )
                line = _line.format(genome_id=genome_id, seqname=seqname)
                new_cds_ids.add(line)
              elif line.startswith((">", "\\n", "*")):
                save = False

              if save:
                print(line, file=out, end="")
              else:
                continue
    
      assert new_cds_ids == new_protein_ids, f"CDS and protein IDs dont all match:\\n{(new_cds_ids - new_protein_ids)=}\\n{(new_protein_ids - new_cds_ids)=}"
    
    """
}
