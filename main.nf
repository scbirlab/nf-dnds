#!/usr/bin/env nextflow

/*
========================================================================================
   Variant calling Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-dnds
   Contact  : Eachan Johnson <user@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/

def pipeline_name = """\
         S C B I R   P I P E L I N E
         ===========================
         """.stripIndent()

if ( params.help ) {
   println """${pipeline_name}
         Nextflow pipeline to ....

         Usage:
            nextflow run scbirlab/nf-dnds --sample_sheet <csv> --inputs <dir>
            nextflow run scbirlab/nf-dnds -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 to be processed

         Optional parameters (with defaults):  
            inputs             Directory containing inputs. Default: "./inputs".

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}

working_dir = params.outputs

log.info """${pipeline_name}
         inputs
            input dir.     : ${params.inputs}
            sample sheet   : ${params.sample_sheet}
            max. genomes   : ${params.max_strains}
         output            : ${params.outputs}
         """
         .stripIndent()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

include {
   FastTree;
} from './modules/fasttree.nf'
include {
   HyPhy;
   ParseHyPy;
} from './modules/hyphy.nf'
include {
   MAFFT;
} from './modules/mafft.nf'
include {
   Fetch_genome_from_NCBI2 as Fetch_genome_from_NCBI;
} from './modules/ncbi.nf'
include {
   Pal2Nal;
} from './modules/pal2nal.nf'
include {
   ProteinOrtho;
   ExtractOrthogroups;
} from './modules/proteinortho.nf'
include {
   FetchTaxonkitDB;
   FetchTaxonomicRanks;
   ResolveSpecies;
} from './modules/taxonkit.nf'
include {
   MultiQC;
} from './modules/multiqc.nf'

workflow {

   Channel.fromPath( 
      params.sample_sheet, 
      checkIfExists: true 
   )
      .tap { sample_sheet_ch }
      .splitCsv( header: true )
      .set { csv_ch }

   csv_ch
      .map { v -> v.taxon_id }
      .set { taxon_ch }

   FetchTaxonkitDB(
      Channel.of( params.taxonkit_db_url ),
   )

   FetchTaxonomicRanks(
      sample_sheet_ch
         .combine( FetchTaxonkitDB.out ),
      Channel.value( "taxon_id" ),
   )

   taxon_ch
      .combine( FetchTaxonkitDB.out ) 
   | ResolveSpecies

   Fetch_genome_from_NCBI( 
      ResolveSpecies.out.stdout, 
      Channel.value( params.max_strains ),
   )

   ProteinOrtho(
      Fetch_genome_from_NCBI.out.protein,
   )

   ExtractOrthogroups(
      ProteinOrtho.out.tsv
         .combine( Fetch_genome_from_NCBI.out.protein, by: 0 )
         .combine( Fetch_genome_from_NCBI.out.cds, by: 0 )
   )

   MAFFT(
      ExtractOrthogroups.out.protein
         .transpose()
         .map { v -> [
            v[0], v[1].simpleName.split("orthogroup-")[-1], v[1]] 
         },
   )

   Pal2Nal(
      MAFFT.out
         .combine(
            ExtractOrthogroups.out.cds
               .transpose()
               .map { v -> [
                  v[0], v[1].simpleName.split("orthogroup-")[-1], v[1]
               ] },
            by: [0, 1],
         ),
   )
   | FastTree
   
   HyPhy(
      Pal2Nal.out.combine( FastTree.out, by: [0, 1] ),
      Channel.value( params.pvalue ),
   )

   // HyPhy.out
   //    .groupTuple()
   //    .set { grouped_ch }

   // ParseHyPhy( grouped_ch )

   // outputs
   //    .concat( other_outputs )
   //    .flatten()
   //    .unique()
   //    .collect()
   //    | multiQC

}


/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/