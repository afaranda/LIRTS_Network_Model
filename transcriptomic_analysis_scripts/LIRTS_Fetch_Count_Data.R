##############################################################################
#                                                                            #
# File: LIRTS_Fetch_Count_Data.R                                             #
# Puropose: Fetch gene level htseq-count data from synapse project and       #
#           assemble an annotated master DGEList.                            #
#                                                                            #
#                                                                            #
#                                                                            #
#                                                                            #
#                                                                            #
#                                                                            #
# Created: July 9, 2021                                                      #
# Author: Adam Faranda                                                       #
#                                                                            #
##############################################################################

############ Load Libraries and Source Dependencies ##########################
wd <- "/Users/adam/Documents/LEC_Time_Series"
library('synapser')
library('edgeR')

synapser::synLogin()
syn_project <- "syn25579691"                 ## Synapse ID for this project
syn_count_dir <- "syn25979754"               ## Synapse folder with counts
syn_code_dir <- "syn25976329"                ## Synapse folder with code
syn_sample_table <- "syn25582010"            ## Synapse sample table
syn_gene_meta <- "syn25976328"               ## Gene Annotations
local_data_dir <- paste0(wd,"/LIRTS_Raw_Data")        ## Local data directory

# Add this script to the code dir
synapse_push <- File(
  path="transcriptomic_analysis_scripts/LIRTS_Fetch_Count_Data.R",
  parent=syn_code_dir
)
synapse_push <- synStore(
  synapse_push
)

# Retrieve files from a synapse folder and store locally
count_files <- as.list(synGetChildren(syn_count_dir))
if(!dir.exists(local_data_dir))
{
  dir.create(local_data_dir)
}
ft <- data.frame()
for(f in count_files) {
  if(grepl("_rf_GeneCount", f$name)){
    ft <- rbind(
      ft, 
      data.frame(
        sample=gsub("_rf_GeneCount.txt", "", f$name),
        files=f$name,
        count_file_id=f$id
      )
    )
    synapser::synGet(f$id, downloadLocation=local_data_dir, downloadFile=T)
  }
}

# Retrieve sample information table from synapse and join 
# to table of count files
ft <- merge(
  ft, 
  synTableQuery(
      sprintf(
        'select sample, genotype, batch, hours_pcs from %s',
      synFindEntityId(name='Samples', parent=syn_project))
  )$asDataFrame(),
  by="sample"
)

# Read files into a dgelist object and store
master <- readDGE(
  files=subset(ft, select=-c(ROW_ID, ROW_VERSION)),
  header=F, path=local_data_dir
)

# Factorize covariates
master$samples$hours_pcs<-factor(
  paste(master$samples$hours_pcs,"H",sep=""),
  levels=c("0H", "6H", "24H", "48H", "72H", "120H")
)

master$samples$genotype<-factor(
  master$samples$genotype,
  levels=c("WT","FN", "B8")
)

master$samples$group<-factor(
  paste(
    master$samples$genotype,
    master$samples$hours_pcs,
    master$samples$batch,
    sep="_"
  ),
  levels=c(
    "WT_0H_DBI","WT_0H_DNA1","WT_0H_DNA2","WT_0H_DNA3",
    "WT_6H_DNA1",
    "WT_24H_DBI","WT_24H_DNA1",
    "WT_48H_DBI",
    "WT_72H_DNA3",
    "WT_120H_DNA2",
    "FN_0H_DBI", "FN_48H_DBI",
    "B8_0H_DNA1", "B8_24H_DNA1"
  )
)

master$samples %>%
  group_by(group) %>%
  mutate(
    label=paste(group, row_number(),sep="_")
  ) %>% as.data.frame() -> master$samples

# Strip QC Rows from Count Matrix
master <- master[!grepl("__", row.names(master$counts)),,keep.lib.sizes=F]

## If the annotation table was already downloaded, add it to master DGEList
fn<-'LIRTS_Raw_Data/Mouse_Gene_Annotations.csv'

## If the annotation table has not yet been generated, generate it
if(file.exists(fn)){
  lt<-read.csv(fn, row.names = 1)
  row.names(lt) <- lt$gene_id
  master$genes <- lt[row.names(master$counts),]
  print(TRUE)
} else {
  lt <- synGet(
    syn_gene_meta,
    downloadLocation=local_data_dir
  )
  
  lt<-read.csv(fn, row.names = 1)
  row.names(lt) <- lt$gene_id
  master$genes <- lt[row.names(master$counts),]
}

save(master, file="LIRTS_Raw_Data/LIRTS_master_dgelist.Rdata")
syn_dge <- File(
  path="LIRTS_Raw_Data/LIRTS_master_dgelist.Rdata",
  parent=syn_count_dir
)
syn_dge <- synStore(
  syn_dge
)

syn_act <- Activity(
  name="upload_dgelist",
  description="retrieved counts from synapse and assembled a DGEList"
)

syn_act$executed(synapse_push)

synSetProvenance(
  syn_dge,
  syn_act
)


