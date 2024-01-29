
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

devtools::install_github("omarwagih/ggseqlogo")

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("motifmatchr")
#BiocManager::install("chromVAR")

set.seed(1234)

combined = readRDS("data/scATAC/CARE_PORTAL/combined_cell_seurat.rds")
DimPlot(combined, group.by = "cell_type", label = T)

## I dont know why it complained about the fragment index
# added and removed
Fragments(combined@assays$peaks) <- NULL
fragments <- CreateFragmentObject(path = "data/scATAC/CARE_PORTAL/fragments.tsv.gz", cells = colnames(combined), validate.fragments = TRUE)
Fragments(combined@assays$peaks) <- fragments

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

DefaultAssay(endo) <- 'peaks'
Idents(endo) <- endo@meta.data$cell_type

CoveragePlot(
  object = endo,
  region = "SPHK1",
  extend.upstream = 4000,
  extend.downstream = 2000
)

DefaultAssay(combined) <- 'peaks'
Idents(combined) <- combined@meta.data$cell_type

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

####################################

#library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(chromVAR)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

#DefaultAssay(combined) <- "peaks"
#combined

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
table(as.character(seqnames(granges(combined))))
keep.peaks <- which(as.character(seqnames(granges(combined))) %in% main.chroms)
combined[["peaks"]] <- subset(combined[["peaks"]], features = rownames(combined[["peaks"]])[keep.peaks])
combined[["peaks"]] <- keepStandardChromosomes(combined[["peaks"]])
table(as.character(seqnames(granges(combined))))

# add motif information
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
combined
Motifs(combined)

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
top.da.peak

# test enrichment
enriched.motifs <- FindMotifs(
  object = combined,
  features = top.da.peak
)
enriched.motifs
table(enriched.motifs$motif.name %in% "PRPSAP1")

MotifPlot(
  object = combined,
  motifs = head(rownames(enriched.motifs))
)

#### Computing motif activities
p1 = DimPlot(combined, group.by = "cell_type", label = T)
p1

combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38 
)

DefaultAssay(combined) <- 'chromvar'

# look at the activity of Pparg::Rxra
p2 <- FeaturePlot(
  object = combined,
  features = "MA0065.2",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

#################

DefaultAssay(combined) <- "peaks"
open.peaks <- AccessiblePeaks(combined, idents = c("Endothelial"), )
ch17.peak = open.peaks[grep("chr17", open.peaks)]

overlap_peaks <- findOverlaps(combined, StringToGRanges("chr17-74360000-74385000"))
peaksoi = as.character(granges(combined)[queryHits(overlap_peaks)])
peaksoi = gsub(":", "-", peaksoi)


################################################ Endo specifc


endo
DefaultAssay(endo) <- "peaks"
genome(endo) <- "hg38"
endo[['peaks']]

# main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
# table(as.character(seqnames(granges(endo))))
# keep.peaks <- which(as.character(seqnames(granges(endo))) %in% main.chroms)
# endo[["peaks"]] <- subset(endo[["peaks"]], features = rownames(endo[["peaks"]])[keep.peaks])
# endo[["peaks"]] <- keepStandardChromosomes(endo[["peaks"]])
# table(as.character(seqnames(granges(endo))))

# add motif information
endo <- AddMotifs(
  object = endo,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

###
Motifs(endo)
# A Motif object containing 746 motifs in 287415 regions

endo <- RunChromVAR(
  object = endo,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

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

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

