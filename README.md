
## Sub Folders
+ alignment_pipeline_scripts

  primarily shell scripts for processing fastq files, aligning reads
  and estimating transcript abundance.
  
+ GEO_Split_Files

  Processed data files (read counts / TPM) for each sample. These files
  are submitted to GEO
  
+ iSyte_Data

  Lens enrichment data manually retrieved from the iSyte 2.0 website
  
+ manually_reviewed_results

  Analysis output that has been manually reviewed and incorporated into
  a report / paper / presentation (figures and tables)
 
+ transcriptomic_analysis_scripts

  primarily R scripts, used in the bulk of this analysis. 
  
## Transcriptomic Analysis Scripts:
===================================

1. BuildDataMatrix.R

   This script includes a set of functions to assemble the raw output from
   htseq-count and stringtie into data matrices that are suitable for
   down-stream analysese (clustering, differential expression, ANOVA)

2. PreprocessingFunctions.R

   Apply various pre-processing functions to data matrices; these include
   filtering, scaling, noramlization and transformation.


  