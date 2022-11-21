require(devtools)
devtools::install_github("immunogenomics/harmony", force=TRUE)
devtools::install_github("powellgenomicslab/scPred")
library(scPred)
library(magrittr)
library(Seurat)

data("hg38_cytoband")
getwd()

# expression data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110499

counts=read.table("pool_filtered_full_DEG/matrix_ms_filtered.csv",header=T,sep=",",as.is=T,row.names=1)
#counts_matrix <- as.matrix(counts)
treatment <- read.table("pool_filtered_full_DEG/meta_Sample.txt", header=T, sep="\t", as.is=T, row.names=1)
  
data <- CreateSeuratObject(counts = counts, project = "CD9", min.cells = 0, min.features = 0)
data <- FindVariableFeatures(data, do.plot = T, nfeatures =20000,selection_method= "disp" )
data <- SetAssayData(data, assay = "RNA", slot = "scale.data", new.data = as.matrix(counts))
#data <- ScaleData(data)
data <- AddMetaData(data, metadata=treatment)

ref_counts <- read.table("pool_filtered_full_DEG/matrix149_ms_filtered.csv", header=T, sep=",", as.is=T, row.names=1)
ref <- CreateSeuratObject(counts = ref_counts, project = "SUM149_ref", min.cells = 0, min.features = 0)
ref <- FindVariableFeatures(ref, do.plot = T, 
                            selection_method= "disp", nfeatures = 4000)
#ref <- ScaleData(ref)
ref <- SetAssayData(ref, assay = "RNA", slot = "scale.data", new.data = as.matrix(ref_counts))
ref_subtype <- read.table("CD9_combine_DEG/SUM149_v_clusters.csv", header=T, sep=",", as.is=T, row.names=1)
Cell_ID<-ref@meta.data
row.names(ref_subtype) <- rownames(Cell_ID)
ref <- AddMetaData(ref, metadata=ref_subtype)
ref <- RunPCA(ref)
DimPlot(ref, group.by="v_clusters", reduction="pca")
ref <- getFeatureSpace(ref, "v_clusters",)
ref <- trainModel(ref,model="mda", metric="ROC")
ref <- trainModel(ref, model="mda", reclassify= c("cluster9","resistant_DMSO"))
get_scpred(ref)
plot_probabilities(ref)

#DefaultAssay(object = data) <- "RNA"
data <- scPredict(data, ref,   threshold = 0.95)
#data <- RunPCA(data)
DimPlot(data, group.by="scpred_prediction", reduction = "scpred")
DimPlot(data, group.by="Sample", reduction = "scpred")
scpred_predicted_id <- subset.data.frame(data.frame( data@meta.data), select = c( scpred_prediction))
write.csv(scpred_predicted_id,"scpred_score_ms.csv",row.names=T)


