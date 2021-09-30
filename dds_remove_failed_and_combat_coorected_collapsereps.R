source("Aug3_2021/Scripts/Utils.R")
dds <- pre_processing(remove_failed = TRUE)

# setting directory for products 
outdir <- "/Users/sebastianvelez/Desktop/R_Projects/CU_Boulder_RCAN1_Project/Aug3_2021/dds_remove_failed_and_combat_corrected_collapereps/"
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# PCA
vsd <- vst(dds, blind = TRUE) 
PCA_values <- plotPCA(vsd, intgroup =  "failed", returnData = TRUE)
PCA_creation <- lapply(colnames(colData(dds))[-c(1,2,3,5,11,14, 15, 17,19, 20,21,22,23,24,25,27, 28)],
                       function(metanow) {
                         plotPCA(vsd, intgroup = metanow,) + 
                           labs(title = metanow) +
                           coord_cartesian()
                       }) 
PCA_grid_display <- ggarrange(plotlist = PCA_creation, align = "hv")
ggsave(filename = "PCA1.jpg", plot = PCA_grid_display, width = 40 , height = 25 , limitsize = F)

# Batch corrected for single/paired end 
combat_counts_corrected <- ComBat_seq(as.matrix(assay(dds)), dds$singleorpaired)
combat_counts_corrected <- ComBat_seq(as.matrix(combat_counts_corrected), dds$gender)

dds <-  DESeqDataSetFromMatrix(countData = combat_counts_corrected , colData = metadata_no_failed, design = ~ type_tested)


# technical replicate 
dds2 <- collapseReplicates(dds,groupby = dds$RNAprep ,run = dds$RNAprep)

# PCA 
vsd <- vst(dds2, blind = TRUE) 
PCA_values <- plotPCA(vsd, intgroup =  "failed", returnData = TRUE)
PCA_creation <- lapply(colnames(colData(dds))[-c(1,2,3,5,11,14, 15, 17,19, 20,21,22,23,24,25,27, 28)],
                       function(metanow) {
                         plotPCA(vsd, intgroup = metanow,) + 
                           labs(title = metanow) +
                           coord_cartesian()
                       }) 
PCA_grid_display <- ggarrange(plotlist = PCA_creation, align = "hv")
ggsave(filename = "PCA2.jpg", plot = PCA_grid_display, width = 40 , height = 25 , limitsize = F)


# DESeq
DEdds <- DESeq(dds2)

#  MA_plot
res <- results(DEdds, contrast = c("type_tested","Dp16","Dp16dsRcan1"))
jpeg(paste0(outdir,"/",'MA_Plot_Removing_Failed_and_Combat_Corrected_for_SingleOrPairedEnd_Collapse_Replicates.jpg', sep=""))
MA_plot <- plotMA(res, ylim=c(-2,2), main = "Main Plot")
dev.off()

# res based ranked on padj
sink("Kobe.txt")
res2 <- res[ order( res$padj ), ]
sink()
# res summary 
sink("Res_Summary_Table.txt")
summary(res)
sink()

# plotDisp
jpeg(paste0(outdir,"/",'Plot_Disp_Removing_Failed_and_Combat_Corrected_for_SingleOrPairedEnd_Collapse_Replicates.jpg', sep=""))
plot_disp <- plotDispEsts(DEdds, main = "Dp16 vs Dp16dsRcan1")
dev.off()


# for go and enrichr and gsea 
alpha_val <- 0.1
res_sig <- subset(res, padj < alpha_val)
res_expressed <- as.data.frame(res)
res_expressed <- res[!is.na(res_expressed$padj), ]
write.csv(rownames(res_expressed), file = paste0(outdir, "backgroundgenes.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(rownames(res_sig), file = paste0(outdir, "siggenes"), row.names = FALSE, col.names = FALSE, quote = FALSE)

rnkdf <- tibble(gene = rownames(res2), rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>% arrange(desc(rnk)) %>% drop_na()




          
          
          
          
          

