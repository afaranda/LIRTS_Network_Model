library('synapser')
library('dplyr')


## Synapse Login 
synapser::synLogin()

## Specify Count Data and scripts Location
count_data_dir <- "~/Documents/LEC_Time_Series/LIRTS_Raw_Data/"
scripts_dir <- "~/Documents/LEC_Time_Series/transcriptomic_analysis_scripts/"
count_files <- list.files(pattern="_rf_GeneCount.txt", count_data_dir)
annotation_path <- paste0(
  "~/Documents/LEC_Time_Series/",
  "LIRTS_Raw_Data/Mouse_Gene_Annotations.csv"
)

## Construct Sample Covariate Table
ft <- data.frame(
  sample=gsub("_rf_GeneCount.txt", "", count_files),
  files=count_files,
  group=factor(
      ifelse(
          grepl("Cre_N", count_files),
	  ifelse(
	      grepl("^0H", count_files),
	      "WT0H", "WT72H"
	  ),
	  ifelse(
	      grepl("^0H", count_files),
	      "R10H", "R172H"
	  )
      ),
      levels=c("WT0H", "WT72H", "R10H", "R172H")
  ),
  genotype=factor(
      ifelse(
          grepl("Cre_N", count_files),
	  "WT", "R1"
      ),
      levels=c("WT", "R1")
  ),
  interval=factor(
      ifelse(
          grepl("^0H", count_files),
	      "0H", "72H"
	  ),
      levels=c("0H", "72H")
  )
) %>%
group_by(group) %>%
  mutate(label = paste(group, row_number(), sep="_")) %>%
  as.data.frame()


# Create New Project and push data to it
syn_project <- Project(
  name="Runx1_Study"
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
  path=paste0(scripts_dir, "/Runx1_Push_Synapse.R"),
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
for(f in list.files(count_data_dir)){
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

## Define schema for synapse sample table
syn_cols <- list(                                 
  Column(name="sample", columnType='STRING', maximumSize =50), 
  Column(name="files", columnType='STRING', maximumSize =100),
  Column(name="group", columnType='STRING', maximumSize =20),
  Column(name="genotype", columnType='STRING', maximumSize =10),
  Column(name="interval", columnType='STRING', maximumSize =10),
  Column(name="label", columnType='STRING', maximumSize =10)
)
syn_schema <- Schema(
  name="Runx1_Study_Sample_Metadata",
  columns=syn_cols,
  parent=syn_project
)

## Create synapse sample table
syn_sample_table <- Table(
  schema = syn_schema,
  values = ft
)

syn_sample_table <- synStore(syn_sample_table)
synSetProvenance(
    syn_sample_table,
    syn_act
)
