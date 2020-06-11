################################################################################
# File: Overlap_Comparison_Functions.R                                         #
# Author: Adam Faranda                                                         #
# Created: June 27, 2019                                                       #
# Purpose: Compare differential expression results from various ocular tissues #
#          in pairwise contrasts between wildtype and pax6 heterozygous mice   #
#                                                                              #
################################################################################

# Setup Workspace
library('openxlsx')
library('dplyr')
wd<-getwd()

# Given a data frame and a vector of column names; return a vector
# of each column's corresponding number in the data frame. 
colNum<-function(df, n){
  if(is.character(n)){
    positions<-c()
    for(i in n){
      positions<-c(
        positions,
        grep(i, colnames(df))
      )
    }
  } else {
    positions<-n
  }
  positions
}

# Function recieves a data frame and a list of name pairs, returns a data frame
# with names updated based on the list
dfRename<-function(df, name_pairs){
  if (!is.null(names(name_pairs))){
    for(np in names(name_pairs)){
      names(df)[grep(name_pairs[np], names(df))]<-np
    }
  }
  df
}

dfSubname<-function(df, name_pairs, pf="^", po="$"){
  if (!is.null(names(name_pairs))){
    for(np in names(name_pairs)){
      r<-paste(pf,name_pairs[np], po,sep="")
      names(df)<-gsub(r, np, names(df))
    }
  }
  df
}

# Query for Pairwise overlap between two Deglists -- assumes both tables
# have the same column headers, and that "MGI.symbol" is a unique id in both
# lists
query<-function(
  dg1=ss_master, dg2=ss_master, id_col="MGI.symbol",
  cols=c("logFC", "PValue", "FDR", "Group_1", "Group_2", "Avg1", "Avg2")
){
  print(names(dg2))
  print(paste(id_col, cols))
 dg1<-dg1[,c(id_col, cols)]
 dg2<-dg2[,c(id_col, cols)]
 
 for(c in cols){
   names(dg1)[grep(c, names(dg1))]<-paste("dg1.",c, sep="")
   names(dg2)[grep(c, names(dg1))]<-paste("dg2.",c, sep="")
 }
 out<-inner_join(dg1, dg2, by=id_col)
 out
}


# Function takes a set of DEGs tablutes Total, Up and Down Genes
# For statistical and biological significance
degSummary<-function(
  df, lfc_min=1, fdr_max=0.05, 
  Avg1="Avg1", Avg2="Avg2",
  lfc="logFC", fdr="FDR", 
  bioFun = bioSigRNASeq, minExp = 2
){
  names(df)[grep(lfc, names(df))]<-"lfc"
  names(df)[grep(fdr, names(df))]<-"fdr"
  names(df)[grep(Avg1, names(df))]<-"Avg1"
  names(df)[grep(Avg2, names(df))]<-"Avg2"
  
  dg<-data.frame(
    criteria = c("Statistically Significant", "Biologically Significant"),
    Total = c(
      nrow(df %>% filter(abs(lfc) > lfc_min, fdr < fdr_max)), 
      nrow(
      	bioFun(df, lfc ="lfc", stat="fdr", a1="Avg1", a2="Avg2", 
  		minExp= minExp, maxStat=fdr_max, minLfc=lfc_min)
      )
    ),
    Up = c(
      nrow(df %>% filter(lfc > lfc_min, fdr < fdr_max)), 
      nrow(
      	bioFun(df, lfc ="lfc", stat="fdr", a1="Avg1", a2="Avg2", 
  		minExp= minExp, maxStat=fdr_max, minLfc=lfc_min) %>%
        filter(lfc > lfc_min)
      )
    ),
    Down = c(
      nrow(df %>% filter(lfc < -lfc_min, fdr < fdr_max)),
      nrow(
      	bioFun(df, lfc ="lfc", stat="fdr", a1="Avg1", a2="Avg2", 
  		minExp= minExp, maxStat=fdr_max, minLfc=lfc_min) %>%
        filter(lfc < lfc_min)
      ) 
    )
  )
  dg
}

# Extract directional subsets of statistically significant genes
subsetTables<-function(
  df,                       # Data frame with a joined pair of results
  id_col="MGI.symbol",      # Unique Identifier for this gene
  Contrast_1 = "LE",        # Name of the first contrast in df
  Contrast_2 = "PCO",       # Name of the second contrast in df
  dg1="dg1",                # Prefix for contrast 1
  dg2="dg2",                # Prefix for contrast 2
  lfc="logFC",              # Column with log 2 fold change values
  pvl="p_value",            # Column with p value for pairwise test
  fdr="FDR",                # Column with FDR values
  g1 = "Group_1",           # Column with Group_1 label
  g2 = "Group_2",           # Column with Group_2 label
  a1 = "Avg1",              # Column with average values for Group_1
  a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
  unlog = T,                # Whether to report absolute or log2 fold changes
  descname = F,             # Use original, or descriptive attribute names
  annot = NULL,             # Optionally provide table (keyed on ID)
  dropGroup = T             # Drop or keep group label columns
){
  # Standardize column headers
  cols<-c(lfc=lfc, pvl=pvl, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
  df<-dfSubname(df, cols, pf="")
  
  prf<-c(dg1=dg1, dg2=dg2)
  df<-dfSubname(df, prf, po="")
  
  if(unlog){
    df$dg1.lfc<-ifelse(df$dg1.lfc >= 0, 2^df$dg1.lfc, -1/2^df$dg1.lfc)
    df$dg2.lfc<-ifelse(df$dg2.lfc >= 0, 2^df$dg2.lfc, -1/2^df$dg2.lfc)
  }
  if(!is.null(annot)){
    if(any(grepl(id_col, names(annot)))){
       if(length(unique(annot[,id_col])) == nrow(annot)){
        acols<-setdiff(names(annot), id_col)
        dcols<-setdiff(names(df), id_col)
        df<-left_join(df, annot, by=id_col)
        df<-data.frame(df, stringsAsFactors = F)
        df<-df[, c(id_col, acols, dcols)]
      }
    }
  }
  # Subset results
  if (stat){
    tables<-list(
      `Stat Sig Intersection`=df,
      `SS Dn C1 Up C2` = df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      `SS Dn C1 Dn C2` = df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      `SS Up C1 Up C2` = df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      `SS Up C1 Dn C2` = df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
    )
    
  # Use if the data tables submitted via df are biologically significant
  } else {
    tables<-list(
      `Bio Sig Intersection`=df,
      `BS Dn C1 Up C2`= df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      `BS Dn C1 Dn C2`= df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      `BS Up C1 Up C2`= df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      `BS Up C1 Dn C2`= df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
    )
  }
  
  # Rename Contrasts
  names(tables)<-gsub("C1", Contrast_1, names(tables))
  names(tables)<-gsub("C2", Contrast_2, names(tables))
  if(descname){
    cols<-c(
      dg1.g1 = paste(Contrast_1, df$dg1.g1[1], sep="_"),
      dg1.g2 = paste(Contrast_1, df$dg1.g2[1], sep="_"),
      dg2.g1 = paste(Contrast_2, df$dg2.g1[1], sep="_"),
      dg2.g2 = paste(Contrast_2, df$dg2.g2[1], sep="_"),
      dg1.a1 = paste(Contrast_1, df$dg1.g1[1], "Avg", sep="_"),
      dg1.a2 = paste(Contrast_1, df$dg1.g2[1], "Avg", sep="_"),
      dg2.a1 = paste(Contrast_2, df$dg2.g1[1], "Avg", sep="_"),
      dg2.a2 = paste(Contrast_2, df$dg2.g2[1], "Avg", sep="_"),
      dg1.lfc = ifelse(
        unlog, paste(Contrast_1, "Fold_Change",sep="_"),
        paste(Contrast_1, lfc,sep="_")
      ),
      dg2.lfc=ifelse(
        unlog, paste(Contrast_2, "Fold_Change",sep="_"),
        paste(Contrast_2, lfc,sep="_")
      ),
      dg1.pvl = paste(Contrast_1, pvl,sep="_"),
      dg2.pvl = paste(Contrast_2, pvl,sep="_"),
      dg1.fdr = paste(Contrast_1, fdr,sep="_"),
      dg2.fdr = paste(Contrast_2, fdr,sep="_")
    )
    
    for (t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        names(sloc)<-cols
        keep<-setdiff(
          names(tables[[t]]), 
          sloc[c(grep("\\.g1$", sloc), grep("\\.g2$", sloc))]
        )
        sloc<-sloc[-c(grep("\\.g1$", sloc), grep("\\.g2$", sloc)) ]
        print(keep)
        tables[[t]]<-tables[[t]][, keep]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        
      } else {
        sloc<-names(cols)
        names(sloc)<-cols
        tables[[t]]<-dfSubname(tables[[t]], sloc)
      }
    }
    
  } else {
    for(t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        
        drop<-c(
          grep("\\.g1$", names(tables[[t]])), 
          grep("\\.g2$", names(tables[[t]]))
        )
        keep<-names(tables[[t]])[-drop]
        print(keep)
        sloc<-sloc[-c(grep("g1$", sloc), grep("g2$", sloc))]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
        
      } else {
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
      }
    }
  }
  tables<-append(
    tables, list(Contrasts=c(Contrast_1=Contrast_1, Contrast_2=Contrast_2))
  )
  tables
}

################################################################

# Extract directional subsets of statistically significant genes
compareHits<-function(
  df,                       # Data frame with a joined pair of results
  id_col="MGI.symbol",      # Unique Identifier for this gene
  Contrast_1 = "LE",        # Name of the first contrast in df
  Contrast_2 = "PCO",       # Name of the second contrast in df
  dg1="dg1",                # Prefix for contrast 1
  dg2="dg2",                # Prefix for contrast 2
  lfc="logFC",              # Column with log 2 fold change values
  pvl="p_value",            # Column with p value for pairwise test
  fdr="FDR",                # Column with FDR values
  g1 = "Group_1",           # Column with Group_1 label
  g2 = "Group_2",           # Column with Group_2 label
  a1 = "Avg1",              # Column with average values for Group_1
  a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
  unlog = T,                # Whether to report absolute or log2 fold changes
  descname = F,             # Use original, or descriptive attribute names
  annot = NULL,             # Optionally provide table (keyed on ID)
  dropGroup = T,            # Drop or keep group label columns
  lfcmin = 1,               # Fold Change threshold
  minDiff = 2,              # Minimum difference for biological significance
  minAvg  = 2               # Minimum Abundance for biologuical significance
){
  print(length(stat))

  # Standardize column headers
  cols<-c(lfc=lfc, pvl=pvl, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
  df<-dfSubname(df, cols, pf="")
  
  prf<-c(dg1=dg1, dg2=dg2)
  df<-dfSubname(df, prf, po="")
  
  if(unlog){
    df$dg1.lfc<-ifelse(df$dg1.lfc >= 0, 2^df$dg1.lfc, -1/2^df$dg1.lfc)
    df$dg2.lfc<-ifelse(df$dg2.lfc >= 0, 2^df$dg2.lfc, -1/2^df$dg2.lfc)
  }
  if(!is.null(annot)){
    print(head(annot))
    if(any(grepl(id_col, names(annot)))){
      print("found ID")
      if(length(unique(annot[,id_col])) == nrow(annot)){
        acols<-setdiff(names(annot), id_col)
        dcols<-setdiff(names(df), id_col)
        df<-left_join(df, annot, by=id_col)
        df<-data.frame(df, stringsAsFactors = F)
        df<-df[, c(id_col, acols, dcols)]
      } else {
        print("Length Mismatch")
      }
    }
  }
  # Subset results
  print(length(stat))
  if (stat){
    tables<-list(
      `Stat Sig Intersection`= df %>%
        filter(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05)%>%
        filter(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05),
      `SS Dn C1 Up C2` = df %>% 
        filter(
          dg1.lfc < 0-lfcmin & dg1.fdr < 0.05, 
          dg2.lfc > lfcmin & dg2.fdr < 0.05),
      
      `SS Dn C1 Dn C2` = df %>% 
        filter(
          dg1.lfc < 0-lfcmin & dg1.fdr < 0.05, 
          dg2.lfc < 0-lfcmin & dg2.fdr < 0.05),
      
      `SS Up C1 Up C2` =df %>% 
        filter(dg1.lfc > lfcmin & dg1.fdr < 0.05, 
               dg2.lfc > lfcmin & dg2.fdr < 0.05),
      
      `SS Up C1 Dn C2` = df %>% 
        filter(dg1.lfc > lfcmin & dg1.fdr < 0.05, 
               dg2.lfc < 0-lfcmin & dg2.fdr < 0.05),
      
      `SS Up C1 Only`  = df %>% 
        filter(dg1.lfc > lfcmin & dg1.fdr < 0.05, 
               !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05)),
      
      `SS Dn C1 Only`  = df %>% 
        filter(dg1.lfc < 0-lfcmin & dg1.fdr < 0.05, 
               !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05)),
      
      `SS Up C2 Only`  = df %>% 
        filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05), 
               dg2.lfc > lfcmin & dg2.fdr < 0.05),
      
      `SS Dn C2 Only`  = df %>% 
        filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05), 
               dg2.lfc < 0-lfcmin & dg2.fdr < 0.05)
    )
    
    # Use if the data tables submitted via df are biologically significant
  } else {
    dg <- df
    df <- df %>% 
      filter(abs(dg1.a1 - dg1.a2) > minDiff) %>%
      filter(dg1.a1 > minAvg | dg1.a2 > minAvg) %>%
      filter(abs(dg2.a1 - dg2.a2) > minDiff) %>%
      filter(dg2.a1 > minAvg | dg2.a2 > minAvg)
    
    tables<-list(
      # `Bio Sig Intersection`=df,
      # `BS Dn C1 Up C2`= df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      # `BS Dn C1 Dn C2`= df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      # `BS Up C1 Up C2`= df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      # `BS Up C1 Dn C2`= df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
      
      
      `Bio Sig Intersection`=df %>%
        filter(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05)%>%
        filter(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05),
      `BS Dn C1 Up C2` = df %>% 
        filter(
          dg1.lfc < 0-lfcmin & dg1.fdr < 0.05, 
          dg2.lfc > lfcmin & dg2.fdr < 0.05),
      
      `BS Dn C1 Dn C2` = df %>% 
        filter(
          dg1.lfc < 0-lfcmin & dg1.fdr < 0.05, 
          dg2.lfc < 0-lfcmin & dg2.fdr < 0.05),
      
      `BS Up C1 Up C2` =df %>% 
        filter(dg1.lfc > lfcmin & dg1.fdr < 0.05, 
               dg2.lfc > lfcmin & dg2.fdr < 0.05),
      
      `BS Up C1 Dn C2` = df %>% 
        filter(dg1.lfc > lfcmin & dg1.fdr < 0.05, 
               dg2.lfc < 0-lfcmin & dg2.fdr < 0.05),
      
      `BS Up C1 Only`  = df %>% 
        filter(dg1.lfc > lfcmin & dg1.fdr < 0.05, 
               !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05)),
      
      `BS Dn C1 Only`  = df %>% 
        filter(dg1.lfc < 0-lfcmin & dg1.fdr < 0.05, 
               !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05)),
      
      `BS Up C2 Only`  = df %>% 
        filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05), 
               dg2.lfc > lfcmin & dg2.fdr < 0.05),
      
      `BS Dn C2 Only`  = df %>% 
        filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05), 
               dg2.lfc < 0-lfcmin & dg2.fdr < 0.05),
      `All Observations` = dg
    )
  }
  
  # Rename Contrasts
  names(tables)<-gsub("C1", Contrast_1, names(tables))
  names(tables)<-gsub("C2", Contrast_2, names(tables))
  if(descname){
    cols<-c(
      dg1.g1 = paste(Contrast_1, df$dg1.g1[1], sep="_"),
      dg1.g2 = paste(Contrast_1, df$dg1.g2[1], sep="_"),
      dg2.g1 = paste(Contrast_2, df$dg2.g1[1], sep="_"),
      dg2.g2 = paste(Contrast_2, df$dg2.g2[1], sep="_"),
      dg1.a1 = paste(Contrast_1, df$dg1.g1[1], "Avg", sep="_"),
      dg1.a2 = paste(Contrast_1, df$dg1.g2[1], "Avg", sep="_"),
      dg2.a1 = paste(Contrast_2, df$dg2.g1[1], "Avg", sep="_"),
      dg2.a2 = paste(Contrast_2, df$dg2.g2[1], "Avg", sep="_"),
      dg1.lfc = ifelse(
        unlog, paste(Contrast_1, "Fold_Change",sep="_"),
        paste(Contrast_1, lfc,sep="_")
      ),
      dg2.lfc=ifelse(
        unlog, paste(Contrast_2, "Fold_Change",sep="_"),
        paste(Contrast_2, lfc,sep="_")
      ),
      dg1.pvl = paste(Contrast_1, pvl,sep="_"),
      dg2.pvl = paste(Contrast_2, pvl,sep="_"),
      dg1.fdr = paste(Contrast_1, fdr,sep="_"),
      dg2.fdr = paste(Contrast_2, fdr,sep="_")
    )
    
    for (t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        names(sloc)<-cols
        keep<-setdiff(
          names(tables[[t]]), 
          sloc[c(grep("\\.g1$", sloc), grep("\\.g2$", sloc))]
        )
        sloc<-sloc[-c(grep("\\.g1$", sloc), grep("\\.g2$", sloc)) ]
        print(keep)
        tables[[t]]<-tables[[t]][, keep]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        
      } else {
        sloc<-names(cols)
        names(sloc)<-cols
        tables[[t]]<-dfSubname(tables[[t]], sloc)
      }
    }
    
  } else {
    for(t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        
        drop<-c(
          grep("\\.g1$", names(tables[[t]])), 
          grep("\\.g2$", names(tables[[t]]))
        )
        keep<-names(tables[[t]])[-drop]
        print(keep)
        sloc<-sloc[-c(grep("g1$", sloc), grep("g2$", sloc))]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
        
      } else {
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
      }
    }
  }
  tables<-append(
    tables, list(Contrasts=c(Contrast_1=Contrast_1, Contrast_2=Contrast_2))
  )
  tables
}
################################################################

# Tabulate Directional Intersections between the two data sets
# Recieves set of tables generated by "subsetTables()" returns a
# data frame with row counts for each directional intersect
tabulateOverlap<-function(tables, rename=F){
  stat.Intersect<-data.frame(
    C1_UP=c(nrow(tables[[4]]), nrow(tables[[5]])),
    C1_Down=c(nrow(tables[[2]]), nrow(tables[[3]]))
  )
  row.names(stat.Intersect)<-c("C2_Up", "C2_Down")
  if(rename){
    names(stat.Intersect)<-gsub(
      "C1", tables[['Contrasts']][1], names(stat.Intersect)
    )
    row.names(stat.Intersect)<-gsub(
      "C2", tables[['Contrasts']][2], row.names(stat.Intersect)
    )
  }
  stat.Intersect
}

# Biosig Filter -- filter a pair of deg-lists for biologically significant
# genes in both contrasts (only  applicable for a pair of RNASeq contrasts)
bioSig<-function(ar){
  ar %>%
    filter(dg1.Avg1 > 2 | dg1.Avg2 > 2) %>%
    filter(abs(dg1.Avg1 - dg1.Avg2) > 2) %>%
    filter(dg2.Avg1 > 2 | dg2.Avg2 > 2) %>%
    filter(abs(dg2.Avg1 - dg2.Avg2) > 2)
}


bioSigRNASeq<-function(
  df, lfc = "logFC", stat="FDR", a1="Avg1", a2="Avg2", 
  minExp=2, maxStat=0.05, minLfc=1
  ){
    cols<-c(lfc=lfc, stat=stat, a1=a1, a2=a2)
    df<-dfSubname(df, cols)
    df<-df %>% 
      filter(abs(lfc) >= minLfc & stat <= maxStat) %>% 
      filter(a1 >= minExp | a2 >= minExp) %>%
      filter(abs( a1 - a2) > minExp)
    
    sloc<-names(cols)
    names(sloc)<-cols
    df<-dfSubname(df, sloc)
    df
}

bioSigArray<-function(
  df, lfc = "logFC", stat="FDR", a1="Avg1", a2="Avg2", 
  minExp=0, maxStat=0.05, minLfc=1
){
  cols<-c(lfc=lfc, stat=stat, a1=a1, a2=a2)
  df<-dfSubname(df, cols)
  df<-df %>% 
    filter(abs(lfc) >= minLfc & stat <= maxStat) %>% 
    filter(a1 >= minExp | a2 >= minExp)
  sloc<-names(cols)
  names(sloc)<-cols
  df<-dfSubname(df, sloc)
  df
}


bioSigNone<-function(
  df, lfc = "logFC", stat="FDR", a1="Avg1", a2="Avg2", 
  minExp=0, maxStat=0.05, minLfc=1
){
  # Pass-through function that can handle the same parameters
  # as the other sig filters; returns the same data that it is
  # passed
  df
}



# Returns true if all values in a vector have the same sign, false otherwise
dxn<-function(x){
  x<-x[!is.na(x)]
  for(i in sign(x)){
    if(any(i != sign(x))){
      return(FALSE)
    }
  }
  return(TRUE)
}


# Given a grouped data set, returns the value of f1 that corresponds
# to the row in f2 matching the "key" value -- this function can be used
# to pivot data in a "summarize" query

pivot<-function(f1, f2, key){ 
  if(length(grep(key, f2)) == 1){
    return(
      nth(f1, grep(key, f2))
    )
  } else {
    return (NA)
  }
}

# General Purpose "Unlog" function. Given a data frame and a range or list of
# columns, convert log2 fold changes in the given columns to fold changes. 
unlog<-function(x, cols){
  for(c in cols){
    x[,c]<-ifelse(x[,c] >= 0, 2^x[,c], -1/2^x[,c])
  }
  x
}


# Function to select one row from a set of duplicates.  For any
# duplicate symbol, the record with the greatest absolute logFC is
# retained and all others are discarded
uniqueMaxLfc<-function(
  df, idc="MGI.symbol", lfc="logFC", fdr="FDR", fdr_min=0.05){
  names(df)[grep(idc, names(df))]<-"idc"
  names(df)[grep(lfc, names(df))]<-"lfc"
  names(df)[grep(fdr, names(df))]<-"fdr"
  
  df <- bind_rows(  
    df %>% 
      group_by(idc) %>% filter(min(fdr) > fdr_min) %>%
      group_by(idc) %>% filter(abs(lfc) == max(abs(lfc))) %>% 
      filter(row_number() == 1),
    df %>% 
      group_by(idc) %>% filter(min(fdr) <= fdr_min ) %>%
      group_by(idc) %>% filter(fdr <= fdr_min) %>%
      group_by(idc) %>% filter(abs(lfc) ==  max(abs(lfc)))
  )
  
  names(df)[grep("idc", names(df))]<-idc
  names(df)[grep("lfc", names(df))]<-lfc
  names(df)[grep("fdr", names(df))]<-fdr
  return (data.frame(df, stringsAsFactors=F))
}

# Helper function convert fold change to Log2 fold change
logify<-function(x, base=2){
  if(x >= 0){
    return (log(x, base))
  }
  else if( x < 0){
    return( log( 1/abs(x), base ))
  }
}

# Function to detect duplicates in a data frame and report them based
# on one or more columns given as unique ID's
printDups<-function(df, idc='MGI.symbol'){
	df$pk<-apply(data.frame(df[,idc]), 1, paste, collapse = "")
	dups<-unique(df[duplicated(df$pk), 'pk'])
	dg<-df[df$pk %in% dups,]
	dg<-dg[order(dg$pk),setdiff(names(dg), 'pk')]
	dg
}


# Function to select one row from a set of duplicates.  For any
# duplicate symbol, the record with the greatest total abundance based
# on the sum of group averages is selected for downstream analysis
uniqueTotalExp<-function(
  df, idc="MGI.symbol", av1="Avg1", av2="Avg2"){
  names(df)[grep(paste("^",idc,"$", sep=""), names(df))]<-"idc"
  names(df)[grep(paste("^",av1,"$", sep=""), names(df))]<-"av1"
  names(df)[grep(paste("^",av2,"$", sep=""), names(df))]<-"av2"
    

  df <- bind_rows(  
    df %>% 
      group_by(idc) %>% filter(n() == 1),
    df %>% 
      group_by(idc) %>% filter(n() > 1 ) %>%
      group_by(idc) %>% filter(av1 + av2 == max(av1 + av2)) %>%
      group_by(idc) %>% filter(row_number() == 1)
  )
  
  names(df)[grep("idc", names(df))]<-idc
  names(df)[grep("av1", names(df))]<-av1
  names(df)[grep("av2", names(df))]<-av2
  return (data.frame(df, stringsAsFactors=F))
}

# vennIntersections: Given a joined pair of DEG tables, this function
# tabulates the DEGs detected in both tables, or only in one table
# with total, and directional partitions. It returns a four column
# data frame. 

vennIntersections<-function(
  df,                       # Data frame with a joined pair of results
  id_col="MGI.symbol",      # Unique Identifier for this gene
  Contrast_1 = "LE",        # Name of the first contrast in df
  Contrast_2 = "PCO",       # Name of the second contrast in df
  dg1="dg1",                # Prefix for contrast 1
  dg2="dg2",                # Prefix for contrast 2
  lfc="logFC",              # Column with log 2 fold change values
  pvl="p_value",            # Column with p value for pairwise test
  fdr="FDR",                # Column with FDR values
  g1 = "Group_1",           # Column with Group_1 label
  g2 = "Group_2",           # Column with Group_2 label
  a1 = "Avg1",              # Column with average values for Group_1
  a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
  lfcmin = 1,               # Fold Change threshold
  minDiff = 2,              # Minimum difference for biological significance
  minAvg  = 2               # Minimum Abundance for biologuical significance
){
  # Standardize column headers
  cols<-c(lfc=lfc, pvl=pvl, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
  df<-dfSubname(df, cols, pf="")
  
  prf<-c(dg1=dg1, dg2=dg2)
  df<-dfSubname(df, prf, po="")
  
  df <- df %>% 
    filter(abs(dg1.a1 - dg1.a2) > minDiff) %>%
    filter(dg1.a1 > minAvg | dg1.a2 > minAvg) %>%
    filter(abs(dg2.a1 - dg2.a2) > minDiff) %>%
    filter(dg2.a1 > minAvg | dg2.a2 > minAvg)
  
  Venns<-data.frame(
    Partition = c('Total', 'Increased', 'Decreased'),
    C1 = c(
      nrow(df %>% filter((abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))),
      nrow(df %>% filter((dg1.lfc > lfcmin & dg1.fdr < 0.05) & !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))),
      nrow(df %>% filter((dg1.lfc < 0-lfcmin & dg1.fdr < 0.05) & !(abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05)))),
    Both = c(
      nrow(df %>% filter((abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))), 
      nrow(df %>% filter((dg1.lfc > lfcmin & dg1.fdr < 0.05) & (dg2.lfc > lfcmin & dg2.fdr < 0.05))),
      nrow(df %>% filter((dg1.lfc < 0-lfcmin & dg1.fdr < 0.05) & (dg2.lfc < 0-lfcmin & dg2.fdr < 0.05)))),
    C2 = c(
      nrow(df %>% filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (abs(dg2.lfc) > lfcmin & dg2.fdr < 0.05))),
      nrow(df %>% filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (dg2.lfc > lfcmin & dg2.fdr < 0.05))),
      nrow(df %>% filter(!(abs(dg1.lfc) > lfcmin & dg1.fdr < 0.05) & (dg2.lfc < 0-lfcmin & dg2.fdr < 0.05))))
  )
  names(Venns)[grep('C1', names(Venns))]<-Contrast_1
  names(Venns)[grep('C2', names(Venns))]<-Contrast_2
  return(Venns)
}