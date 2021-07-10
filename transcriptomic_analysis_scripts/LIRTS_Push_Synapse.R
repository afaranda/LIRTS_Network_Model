library('synapser')
library('dplyr')


## Synapse Login 
synapser::synLogin()
syn_project <- ""

## Specify Count Data and scripts Location
count_data_dir <- "~/Documents/LEC_Time_Series/LIRTS_Raw_Data/"
scripts_dir <- "~/Documents/LEC_Time_Series/transcriptomic_analysis_scripts/"
count_files <- list.files(pattern="_rf_GeneCount.txt", count_data_dir)
annotation_path <- paste0(
  "~/Documents/LEC_Time_Series/",
  "LIRTS_Raw_Data/Mouse_Gene_Annotations.csv"
)

# Find Synapse Project and push data to it
syn_project <- Project(
  name="LIRTS_Network_Model"
)
syn_project <- synStore(
  syn_project
)

# Add a count directory to the new project
syn_count_dir <- Folder(name="Counts", parent=syn_project)
syn_count_dir <- synStore(
  syn_count_dir
)


# Add a gene metadata directory to the new project
syn_gene_dir <- Folder(name="Gene_Metadata", parent=syn_project)
syn_gene_dir <- synStore(
  syn_gene_dir
)

# Add a scripts directory to the new folder
syn_code_dir <- Folder(name="Code", parent=syn_project)
syn_code_dir <- synStore(
  syn_code_dir
)

# Add this script to the code dir
synapse_push <- File(
  path=paste0(scripts_dir, "/LIRTS_Push_Synapse.R"),
  parent=syn_code_dir
)
synapse_push <- synStore(
  synapse_push
)

# Add the mouse annotation file to the gene metadata folder
synapse_annot <- File(
  path=annotation_path,
  parent=syn_gene_dir
)
synapse_annot <- synStore(
  synapse_annot
)

syn_act <- Activity(
  name="upload_counts",
  description="upload HTSeq count files to Synapse Project Counts folder"
)

syn_act$executed(synapse_push)

synSetProvenance(
  synapse_annot,
  syn_act
)

# Upload files in the counts folder
cf <- c()
for(f in list.files(count_data_dir, pattern="_rf_GeneCount")){
  print(f)
  count_file <- File(
    path=paste(count_data_dir, f, sep="/"),
    parent=synGet(syn_count_dir)
  )
  count_file <- synStore(count_file)
  synSetProvenance(
    count_file,
    syn_act
  )
}
