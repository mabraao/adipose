rm(list = ls())

#BiocInstaller::biocLite("GreenleafLab/chromVAR")
#BiocManager::install("JASPAR2020")
#BiocManager::install("TFBSTools")

library(ggplot2)
library(dplyr)

library(Signac)
library(Seurat)


library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)

library("Matrix")
library("data.table")

set.seed(1234)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'

#file path for peak index
ygi = "data/scATAC/CARE_PORTAL/all.merged.ygi"
# file path for barcode index
xgi = "data/scATAC/CARE_PORTAL/all.index"
# mfile file path for matrix file in coo
mfile = "data/scATAC/CARE_PORTAL/all.coo.gz"

fragments_path = "data/scATAC/CARE_PORTAL/fragments.tsv.gz"

peak_ids = fread(
  ygi,
  sep="\t", header = FALSE)

cell_ids = fread(
  xgi,
  sep="\t", header = FALSE)

dat = fread(mfile, sep="\t")

mat = sparseMatrix(i = dat$V1 + 1, j =dat$V2 + 1, x=dat$V3)

mat <- mat %>%
  magrittr::set_rownames(cell_ids$V1) %>%
  magrittr::set_colnames(paste0(peak_ids$V1, ":", paste0(peak_ids$V2, "-", peak_ids$V3)))


chrom_assay <- CreateChromatinAssay(
  counts = t(mat),
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = fragments_path,
  validate.fragments=TRUE,
  annotation = annotations
#    min.cells = 10,
#    min.features = 200
)

sn <-  CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)	

meta <- read.table("data/scATAC/CARE_PORTAL/Merged_metadata.tsv", sep = "\t", header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)
cell_annot = read.delim("data/scATAC/CARE_PORTAL/all.cluster", header = F)

meta = merge(meta, cell_annot, by.x="row.names", by.y="V1")
rownames(meta) <- meta$Row.names
meta = meta[,-1]
colnames(meta)[9] <- "cell_type"

meta <- meta[colnames(sn), ]
sn <- AddMetaData(sn, metadata = meta)

head(sn@meta.data)
combined = sn

combined <- NucleosomeSignal(object = combined)
#GetTSSPositions(Annotation(combined))
# # compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined, fast = FALSE)
 
# # add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
## need to be added
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(combined, group.by = 'high.tss') + NoLegend()


combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
FragmentHistogram(object = combined, group.by = 'nucleosome_group')

 
VlnPlot(
  object = combined,
  features = c('TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

Annotation(combined)
genome(combined)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)

combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
DimPlot(object = combined, label = TRUE) + NoLegend()

DimPlot(combined, group.by = "cell_type", label = T)

saveRDS(combined, "data/scATAC/CARE_PORTAL/combined_cell_seurat.rds")


# 
# gene.activities <- GeneActivity(combined)
# 
# 
# pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
# pbmc <- NormalizeData(
#   object = combined,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(pbmc$nCount_RNA)
# )

# da_peaks <- FindMarkers(
#   object = combined,
#   ident.1 = "CD4 Naive",
#   ident.2 = "CD14+ Monocytes",
#   test.use = 'LR',
#   latent.vars = 'nCount_peaks'
# )
# GeneActivity(combined)
# 
# activity.matrix <- CreateGeneActivityMatrix(peak.matrix = combined[['peaks']],
#                                             annotation.file = "../data/Homo_sapiens.GRCh37.82.gtf",
#                                             seq.levels = c(1:22, "X", "Y"),
#                                             upstream = 2000,
#                                             verbose = TRUE)

# 
# # sn$dataset <- basename(mex_dir_path)
# # combined = sn
# combined[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

# 
# combined = readRDS("data/scATAC/combined.samples.rds")
# head(combined@meta.data)
# annotations = read.delim("data/scATAC/all.cluster.txt", header = F)
# head(annotations)


#.combined <- RenameCells(object = .combined, new.names = gsub("*_.", "", colnames(x = .combined[["peaks"]])))

# samples = unique(combined@meta.data$dataset)
# annotations$V1 = gsub("^.{0,14}", "", annotations$V1)
# annotations$V1 = gsub("\\+", "", annotations$V1)
# 
# 
# meta.index = match(annotations$V1, colnames(.combined))
# anno.index = which(!is.na(meta.index))
# 
# .combined@meta.data$cell.type = NA
# .combined@meta.data$cell.type[meta.index[!is.na(meta.index)]] <- annotations$V2[anno.index]
# 
# head(.combined@meta.data)
# table(.combined@meta.data$cell.type)
# 
# .combined@meta.data$cell.type[is.na(.combined@meta.data$cell.type)] <- "Unkown"
gene.activities <- GeneActivity(combined)


# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)


DefaultAssay(combined) <- 'RNA'

FeaturePlot(
  object = combined,
  features = c('MYL3','TBX5', 'CD3D', 'ACTA2', 'CD44', 'PDGFRB', 'LYZ', "PRPSAP1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)



#############################################################################################

endo = subset(combined, subset = cell_type == "Endothelial")
DimPlot(endo, group.by = "cell_type")

FeaturePlot(
  object = endo,
  features = c('SPHK1', 'CD3D', 'LEF1','TREM1', 'LYZ', 'EGFL7', 'MYH7'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# VlnPlot(endo, features = "SPHK1")
# 
# DefaultAssay(.combined) <- 'peaks'
# 
# DimPlot(object = .combined, label = TRUE) + NoLegend()

DefaultAssay(endo) <- 'peaks'
Idents(endo) <- endo@meta.data$cell_type

CoveragePlot(
  object = endo,
  region = "SPHK1",
  extend.upstream = 4000,
  extend.downstream = 2000
)


Idents(combined) <- combined@meta.data$cell_type
# 
CoveragePlot(
  object = combined,
  region = "SPHK1",
  extend.upstream = 4000,
  extend.downstream = 2000,
  annotation = TRUE,
  peaks = TRUE,
  tile = TRUE,
  links = TRUE
)



# 
# Idents(endo) <- endo@meta.data$cell.type
# 
# CoveragePlot(
#   object = endo,
#   region = "SPHK1",
#   extend.upstream = 40000,
#   extend.downstream = 20000
# )

# DimPlot(object = .combined, label = TRUE) + NoLegend()


#######################################################
library(BSgenome.Hsapiens.UCSC.hg19)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
# 
DefaultAssay(combined) <- "peaks"
combined
# 

# # 
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)
keep.peaks <- which(as.character(seqnames(granges(combined))) %in% main.chroms)
combined[["peaks"]] <- subset(combined[["peaks"]], features = rownames(combined[["peaks"]])[keep.peaks])
combined[["peaks"]] <- keepStandardChromosomes(combined[["peaks"]])
# 
# add motif information
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = pfm
)

######

library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library("BSgenome.Hsapiens.UCSC.hg19")

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(combined),
  pwm = pfm,
  genome = 'hg19'
)

Motifs(combined)

pfm_tmp <- GetMotifData(object = motif, slot = "pwm")


####
# # 
# Motifs(combined) <- motif
# # 
da_peaks <- FindMarkers(
  object = combined,
  ident.1 = 'Endothelial',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

####
DefaultAssay(combined) <- "peaks"
open.peaks <- AccessiblePeaks(combined, idents = c("Endothelial"), )
ch17.peak = open.peaks[grep("chr17", open.peaks)]

overlap_peaks <- findOverlaps(combined, StringToGRanges("chr17-74360000-74385000"))
peaksoi = as.character(granges(combined)[queryHits(overlap_peaks)])
peaksoi = gsub(":", "-", peaksoi)

DefaultAssay(endo) <- "peaks"
# test enrichment
enriched.motifs <- FindMotifs(
  object = endo,
  features = as.character(peaksoi)
)


p1 = DimPlot(combined, group.by = "cell_type", label = T)
p1
### motific activity
# 

combined

combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg19 
)

DefaultAssay(combined) <- 'chromvar'
enriched.motifs[enriched.motifs$motif.name %in% c("PRPSAP1", "MYH7", "EGFL7", "SPHK1"),]

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = e,
  features = "MA0066.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

# 
# 
# 
# # get top differentially accessible peaks
# top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
# 
# # test enrichment
# enriched.motifs <- FindMotifs(
#   object = combined,
#   features = top.da.peak
# )
# 
# MotifPlot(
#   object = combined,
#   motifs = head(rownames(enriched.motifs))
# )
# 
# 

p1 = DimPlot(combined, group.by = "cell_type", label = T)
p1
### motific activity
# 

combined

combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(combined) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = combined,
  features = "SPHK1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2


p2 <- FeaturePlot(
  object = combined,
  features = "PRPSAP1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2
# 
# 
# differential.activity <- FindMarkers(
#   object = combined,
#   ident.1 = 'Pvalb',
#   ident.2 = 'Sst',
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff"
# )
# 
# MotifPlot(
#   object = combined,
#   motifs = head(rownames(differential.activity)),
#   assay = 'peaks'
# )


endo = subset(combined, subset = cell_type == "Endothelial")
DimPlot(endo, group.by = "cell_type")

p1 = DimPlot(endo, group.by = "cell_type", label = T)
p1
### motific activity
# 
#############################################
endo
DefaultAssay(endo) <- "peaks"
endo[['peaks']]

endo <- RunChromVAR(
  object = endo,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

#enriched.motifs[enriched.motifs$motif.name %in% c("PRPSAP1", "MYH7", "EGFL7", "SPHK1"),]

de.motif <- rownames(endo)[grep("chr17", rownames(endo))]
bg.peaks <- tail(rownames(combined))

motifs = FindMotifs(
  object = endo,
  features = de.motif,
  background = bg.peaks
)

motifs[motifs$motif.name %in% c("PRPSAP1", "MYH7", "EGFL7", "SPHK1"),]

DefaultAssay(endo) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = combined,
  features = "MA0019.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2


p2 <- FeaturePlot(
  object = endo,
  features = "PRPSAP1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

p2 <- FeaturePlot(
  object = endo,
  features = "EGFL7",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

p2 <- FeaturePlot(
  object = endo,
  features = "MYH7",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2


DefaultAssay(combined) <- "peaks"

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(combined, idents = c("Endothelial"), )

# match the overall GC content in the peak set
meta.feature <- GetAssayData(combined, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = combined,
  features = top.da.peak
)

enriched.motifs[enriched.motifs$motif.name %in% c("PRPSAP1", "MYH7", "EGFL7", "SPHK1"),]

MotifPlot(
  object = combined,
  motifs = head(rownames(enriched.motifs))
)

FeaturePlot(
  object = combined,
  features = "KLF15",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)




