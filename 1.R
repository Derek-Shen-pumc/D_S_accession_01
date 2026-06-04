library(Seurat)
library(tidyr)
library(ggplot2)
library(dplyr)
library(sctransform)
library(scDblFinder)

sce <- function(sc_data){
  sce <- as.SingleCellExperiment(sc_data) %>% scDblFinder(dbr=0.1)
  a<- CreateSeuratObject(
    counts = assay(sce, "counts"),  # 提取基因表达矩阵
    assay = "RNA",                  # 命名数据层为RNA
    meta.data = as.data.frame(colData(sce))  # 携带质检结果
  ) %>% 
    subset(scDblFinder.class == "singlet")
}

C1 <- CreateSeuratObject(counts = Read10X('data/C1/'), project = "C1", min.cells = 3, min.features = 200)
C1 <- sce(C1)
C2 <- CreateSeuratObject(counts = Read10X('data/C2/'), project = "C2", min.cells = 3, min.features = 200)
C2 <- sce(C2)
C3 <- CreateSeuratObject(counts = Read10X('data/C3/'), project = "C3", min.cells = 3, min.features = 200)
C3 <- sce(C3)
Control <- merge(x = C1, y = c(C2, C3))
Control$group <- Control$orig.ident
Control$stim <- 'Control'

P1 <- CreateSeuratObject(counts = Read10X('data/P1/'), project = "P1", min.cells = 3, min.features = 200)
P1 <- sce(P1)
P2 <- CreateSeuratObject(counts = Read10X('data/P2/'), project = "P2", min.cells = 3, min.features = 200)
P2 <- sce(P2)
P3 <- CreateSeuratObject(counts = Read10X('data/P3/'), project = "P3", min.cells = 3, min.features = 200)
P3 <- sce(P3)
Preg <- merge(x = P1, y = c(P2, P3))
Preg$group <- Preg$orig.ident
Preg$stim <- 'Pregnant'

psoriasis <- merge(x = Preg, y = Control)
psoriasis[["percent.mt"]] <- PercentageFeatureSet(psoriasis, pattern = "^mt-")
VlnPlot(psoriasis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(psoriasis, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = T)
plot2 <- FeatureScatter(psoriasis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = T)
plot1 + plot2
psoriasis <- subset(psoriasis, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10)

all_genes <- rownames(psoriasis)
mito_genes <- grep("^mt-", all_genes, value = TRUE)
ribo_genes <- grep("^Rp[sl]", all_genes, value = TRUE)
genes_to_keep <- setdiff(all_genes, c(ribo_genes, mito_genes))
psoriasis <- subset(psoriasis, features = genes_to_keep)
psoriasis <- NormalizeData(psoriasis) %>% FindVariableFeatures() %>% ScaleData()
psoriasis <- RunPCA(psoriasis)
ElbowPlot(psoriasis, ndims = 50)
psoriasis <- RunUMAP(psoriasis, dims = 1:30, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(psoriasis, reduction = 'umap.unintegrated', split.by = 'stim')
saveRDS(object = psoriasis, file = 'data1/psoriasis_normalized.rds', compress = F)
