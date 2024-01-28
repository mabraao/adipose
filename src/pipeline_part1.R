rm(list = ls())

#BiocInstaller::biocLite("GreenleafLab/chromVAR")
#BiocManager::install("JASPAR2020")
#BiocManager::install("TFBSTools")
#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

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
## need to be added
#combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
## need to be added
#combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

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
