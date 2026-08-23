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
} from './modules/proteinortho.nf'
include {
   Fetch_taxonkit_db;
   Fetch_taxonomic_ranks;
} from './modules/taxonkit.nf'
include {
   MultiQC;
} from './modules/multiqc.nf'

workflow {

   Channel.fromPath( 
      params.sample_sheet, 
      checkIfExists: true 
   )
      .splitCsv( header: true )
      .set { csv_ch }

   csv_ch
      .map { v -> v.taxon_id }
      .set { taxon_ch }

   Fetch_genome_from_NCBI( taxon_ch )
   EXTRACT_CDS( DOWNLOAD_GENOMES.out )
   ProteinOrtho( EXTRACT_CDS.out )
   EXTRACT_ORTHOGROUPS( ProteinOrtho.out )

   // Fan-out: one tuple per orthogroup
   og_ch = EXTRACT_ORTHOGROUPS.out
      .flatMap { taxon_id, manifest, og_dir ->
         manifest.readLines().drop(1)
               .collect { og_id ->
                  [
                     taxon_id,
                     og_id,
                     og_dir.resolve("${og_id}/prot.faa"),
                     og_dir.resolve("${og_id}/cds.fna")
                  ]
               }
      }

   MAFFT(og_ch)
   pal2nal(MAFFT.out)
   FastTree(PAL2NAL.out)
   HyPhy(FASTTREE.out)

   HyPhy.out
      .groupTuple()
      .set { grouped_ch }

   PARSE_HYPHY( grouped_ch )

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