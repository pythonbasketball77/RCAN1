library("DESeq2")
library('sva')
library("ggpubr")
library("tidyverse")
library(plotly)

pre_processing <- function(remove_failed ) {
   # loading data
   outdir2 <-  "Aug3_2021/Inputs/"
   single_counts <- read.csv(paste0(outdir2, "single.counts.csv"), header=TRUE, row.names = 1)
   paired_counts <- read.csv(paste0(outdir2, "paired.counts.csv"), header=TRUE, row.names = 1) 
   metadata <- read.csv(paste0(outdir2, "metadata.csv"), header=TRUE)
   
   # preparing counts data (single and paired end) 
   combined_counts <- cbind(single_counts, paired_counts)
   
   # preparing metadata
   metadata$seq <- factor(metadata$seq)
   metadata$type_tested <- factor(metadata$type_tested, levels  = c("Dp16", "Dp16dsRcan1", "WT",  "Rcan1het"))
   metadata$mousesample_id <- as.factor(metadata$mousesample_id)
   
   # counts table exclusing failed samples
   # making sure ncol(countData) == nrow(colData) as required for dds 
   metadata_no_failed <- metadata %>% filter(metadata$failed == "n")
   combined_counts_no_failed <- combined_counts[, metadata_no_failed$specificlabel]
   stopifnot(all(metadata_no_failed$specificlabel == colnames(combined_counts_no_failed)))
   
   if (remove_failed ) {
     dds <-  DESeqDataSetFromMatrix(countData = combined_counts_no_failed , colData = metadata_no_failed, design = ~ type_tested)
   } else {
     dds <- DESeqDataSetFromMatrix(countData = combined_counts , colData = metadata , design = ~ type_tested)
   }
   return (dds)
 }



 
 
 