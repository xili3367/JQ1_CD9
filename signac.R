require(devtools)
devtools::install_github("immunogenomics/harmony", force=TRUE)
devtools::install_github("powellgenomicslab/scPred")
BiocManager::install("Rsamtools")
BiocManager::install("EnsDb.Hsapiens.v75")
install.packages("Signac")
BiocManager::install("biovizBase")
BiocManager::install("motifmatchr")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("chromVAR")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVAR)
set.seed(1234)

#counts <- read.table("GSM4960035_peak_matrix.csv",header=T,sep=",",as.is=T,row.names=1)
counts <- read.table("K562_peak_matrix.csv",header=T,sep=",",as.is=T,row.names=1)
counts.t <- t(counts)



metadata <- read.csv(
  file = "K562_CD9_group.csv",
  header = TRUE,
  row.names = 1
)
row.names(metadata) <- colnames(counts.t)


metadata_10x <- read.csv(
  file = "GSM4960035_Snubar-96plex.cell_metainfo.csv",
  header = TRUE,
  row.names = 1
)
metadata_10x <- metadata_10x[row.names(metadata),]
metadata_10x$CD9_Group <- metadata$CD9_Group

chrom_assay <- CreateChromatinAssay(
  counts = counts.t,
  sep = c("_","_"),
  genome = 'hg19',
  fragments = 'GSM4960035_Snubar-96plex.fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)

data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata_10x
)

granges(data)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(data) <- annotations

# compute nucleosome signal score per cell
data <- NucleosomeSignal(object = data)

# compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = FALSE)

data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters *100
data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

sdata <- subset(
  x = data,
  subset = peak_region_fragments > 3000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)


sdata$high.tss <- ifelse (sdata$TSS.enrichment > 2, "High","Low")
TSSPlot(sdata, group.by = 'high.tss') + NoLegend()
sdata$nucleosome_group <- ifelse (sdata$nucleosome_signal > 4, "NS>4","NS<4")
FragmentHistogram(object=sdata, group.by = "nucleosome_group")


sdata <- RunTFIDF(sdata)
sdata <- FindTopFeatures(sdata, min.cutoff = 'q0')
sdata <- RunSVD(sdata)

sdata <- RunUMAP(object = sdata, reduction = 'lsi', dims = 2:30)
sdata <- FindNeighbors(object = sdata, reduction = 'lsi', dims = 2:30)
sdata <- FindClusters(object = sdata, verbose = FALSE, algorithm = 3)
DimPlot(object = sdata, group.by = "ident") + NoLegend()

gene.activities <- GeneActivity(sdata)
sdata[['RNA']] <- CreateAssayObject(counts = gene.activities)
sdata <- NormalizeData(
  object = sdata,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sdata$nCount_RNA)
)

DefaultAssay(sdata) <- 'RNA'

FeaturePlot(
  object = sdata,
  features = c('CD9'),
  pt.size = 1,
  max.cutoff = 'q95',
  ncol = 3
)

DefaultAssay(sdata) <- 'peaks'

Idents(sdata) <- sdata@meta.data$CD9_Group
DimPlot(object = sdata, group.by = "ident", order=c("CD9_lo","CD9_hi")) 
da_peaks <- FindMarkers(
  object = sdata,
	logfc.threshold = 0.2,
  ident.1 = "CD9_hi",
  ident.2 = "CD9_lo",
  min.pct = 0.05,
  test.use = "LR",
  latent.vars = 'peak_region_fragments'
)

plot1 <- VlnPlot(
  object = sdata,
  features = rownames(da_peaks)[2],
  pt.size = 0.1,
  idents = c("CD9_hi","CD9_lo")
)
plot2 <- FeaturePlot(
  object = sdata,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
sdata <- AddMotifs(sdata, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = pwm)

top.da.peak <- rownames(da_peaks[da_peaks[da_peaks$avg_log2FC>0,]$p_val < 0.005, ])

enriched.motifs <- FindMotifs(
  object = sdata,
  features = top.da.peak 
)

sdata <- Footprint(
  object = sdata,
  motif.name = c("BATF::JUN"),
  genome = BSgenome.Hsapiens.UCSC.hg19
)
PlotFootprint(sdata, features = c("BATF::JUN"))

sdata <- RunChromVAR(
  object = sdata,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

DefaultAssay(sdata) <- 'chromvar'

# look at the activity of Mef2c
FeaturePlot(
  object = sdata,
  features = "CD9",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

cell="K562"

sdata <- readRDS( paste0(cell,"_signac.rds"))


VlnPlot(
  object = sdata,
  features = c("MA1558.1","MA1123.2","MA0511.2","MA0684.2","MA0103.3","MA0143.4"),
  pt.size = 0.,
  idents = c("CD9_lo","CD9_hi"),
	cols = c('CD9_hi' = 'salmon', 'CD9_lo' = '#96ceca'),

)
enriched.motifs[c("MA1558.1","MA1123.2","MA0511.2","MA0684.2","MA0103.3","MA0143.4"),]

write.csv(data.frame( sdata@meta.data)["seurat_clusters"],"K562_seurat_clusters.csv")


motif_activity <- data.frame(t(sdata@assays$chromvar@data))[c(
	"MA1558.1","MA1123.2","MA0511.2","MA0684.2","MA0103.3","MA0143.4")]
row.names(motif_activity) <- make.names(sdata@active.ident, unique=T)

write.csv(motif_activity,
	paste0(cell,"_motif_activity.csv"))


saveRDS(sdata, paste0(cell, "_signac.rds"))


 