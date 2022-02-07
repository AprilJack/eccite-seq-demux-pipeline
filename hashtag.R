# demultiplex scRNA based on hashtags
# follow the tutorial https://satijalab.org/seurat/articles/hashing_vignette.html
# NOTE: this pipeline assumes that all cite-seq data will be used as hashtags. If it is not ture, we need to manually adjust the HTO obj
# NOTE: the demultiplexed dataset contains raw input without any filtering.

library(Seurat) # may require Seurat v4
library(ggplot2)

########################################
# change the below part
# working directory
setwd("/gpfs/analyses/april/pathto/raw_data/")
# sequencing results 
dataDir<-c("sample1_FB_name/outs/filtered_feature_bc_matrix/",
           "sample2_FB_name/outs/filtered_feature_bc_matrix/")
# sample names
dataNames<-c("sample1_FB_name",
             "sample2_FB_name")
# how many hashtags included for each sample?
nTags<-c(3,
         3)
########################################


OBJ<-list()
for(i in seq(1,length(dataDir))) {
  tmpData<-Read10X(dataDir[i])
  tmpObj<-CreateSeuratObject(tmpData$`Gene Expression`,min.cells = 0, min.features = 0)
  # do not apply any filtering for now because we need to do demultiplex first
  tmpObj[["HTO"]]<-CreateAssayObject(counts = tmpData$`Antibody Capture`)
  if(nrow(tmpObj[["HTO"]])!=nTags[i]) {
    stop("contains more hashtags than specified. The demultiplex pipeline needs manual adjustment. stopped.")
  }
  tmpObj <- NormalizeData(tmpObj, assay = "HTO", normalization.method = "CLR") 
  # default margin=1 normalize across features
  tmpObj <- HTODemux(tmpObj, assay = "HTO", positive.quantile = 0.99)
  Idents(tmpObj) <- "HTO_maxID"
  RidgePlot(tmpObj, assay = "HTO", features = rownames(tmpObj[["HTO"]]), ncol = 2)
  ggsave(paste0("ridgePlot_",dataNames[i],".jpeg"),width=20,height=ceiling(nTags[i]/2)*10,units="cm",dpi=600)
  tb.global<-table(tmpObj$HTO_classification.global)
  jpeg(paste0("bar_global_",dataNames[i],".jpeg"),width=12,height=10,units="cm",res=600)
  par(mar=c(10,5,2,2))
  mybar<-barplot(tb.global,ylim=c(0,max(tb.global)*1.5),las=2)
  text(mybar,tb.global+max(tb.global)*0.2, tb.global,cex=1) 
  dev.off()
  
  tb.final<-table(tmpObj$hash.ID)
  jpeg(paste0("bar_final_",dataNames[i],".jpeg"),width=12,height=10,units="cm",res=600)
  par(mar=c(10,5,2,2))
  mybar<-barplot(tb.final,ylim=c(0,max(tb.final)*1.5),las=2)
  text(mybar,tb.final+max(tb.final)*0.2,tb.final ,cex=1) 
  dev.off()
  
  # we use hash.ID as the final ID to demultiplex the samples
  obj.list <- SplitObject(tmpObj, split.by = "hash.ID")
  for(id in setdiff(unique(tmpObj$hash.ID),c("Negative","Doublet"))) {
    sampleName<-paste0(dataNames[i],"_",id)
    print(sampleName)
    saveRDS(obj.list[[id]],paste0(sampleName,"_demultiplxed.rds"))
    OBJ[[sampleName]]<-obj.list[[id]]
  }
}

save.image("hashtag.RData")

