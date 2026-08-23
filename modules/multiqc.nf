process multiQC {

   errorStrategy 'retry'
   maxRetries 1

   publishDir( 
      "${params.outputs}/multiqc", 
      mode: 'copy',
   )

   input:
   path( '*', stageAs: '?/*' )

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}