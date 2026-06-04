library(Seurat)
library(tidyr)
library(ggplot2)
library(dplyr)
library(sctransform)
library(future)
library(clustree)
library(openxlsx)
library(RColorBrewer)

immune_cell <- readRDS('data1/immune_cell_round0.rds')
immune_cell[['RNA']] <- split(immune_cell[['RNA']], f = immune_cell$orig.ident)
immune_cell <- NormalizeData(immune_cell) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(immune_cell, ndims = 50)
immune_cell <- RunUMAP(immune_cell, dims = 1:30, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(immune_cell, reduction = 'umap.unintegrated', split.by = 'stim', group.by = 'orig.ident')
immune_cell <- IntegrateLayers(object = immune_cell, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                             verbose = T)
immune_cell[['RNA']] <- JoinLayers(immune_cell[['RNA']])
immune_cell <- RunUMAP(immune_cell, dims = 1:30, reduction = 'integrated.cca', reduction.name = 'umap.cca')
DimPlot(immune_cell, reduction = "umap.cca")

immune_cell <- FindNeighbors(immune_cell, reduction = 'integrated.cca', dims = 1:30)
a <- seq(from = 0.1, to = 1.5, by = 0.1)
immune_cell <- FindClusters(immune_cell, resolution = a)
clustree(immune_cell)
Idents(immune_cell) <- 'RNA_snn_res.0.8'
immune_cell$seurat_clusters <- immune_cell$RNA_snn_res.0.8
cluster_colors <- c(
  "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
  "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",
  "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000",
  "#AAFFC3"
)
DimPlot(immune_cell, reduction = 'umap.cca', label = T, label.size = 5, cols = cluster_colors)
immune.cell.markers <- FindAllMarkers(immune_cell, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
all.markers.list <- split(immune.cell.markers, f = immune.cell.markers$cluster)
openxlsx::write.xlsx(all.markers.list, file = "data1/all_marker_list_0.8_immune.xlsx")
saveRDS(immune.cell.markers, file = 'data1/0.8_immnue_markers.rds')
top30 <- immune.cell.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/0.8_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data1/top30_0.8_immune.xlsx")
saveRDS(immune_cell, file = 'data1/immune_cell_round1.rds', compress = F)

immune_cell <- readRDS('data1/immune_cell_round1.rds')
FeaturePlot(immune_cell, features = 'Klrk1', reduction = 'umap.cca')
celltype <- c('gdT_cell', 'M', 'gdT_cell', 'NEU', 'NEU', 'KC', 'DC', 'DC2', 'NK_cell', 'KC',
              'DC', 'M', 'Th_cell', 'DC', 'Mast_cell')
names(celltype) <- levels(immune_cell)
immune_cell <- RenameIdents(immune_cell, celltype)
immune_cell$celltype_1 <- immune_cell@active.ident
immune_cell <- subset(immune_cell, idents = 'KC', invert = T)
DimPlot(immune_cell, reduction = 'umap.cca', label = T, label.size = 5, cols = cluster_colors,
        group.by = 'celltype_1')

DC2 <- subset(immune_cell, idents = 'DC2')
DimPlot(DC2)
FeaturePlot(DC, features = 'Cd3e')
DC2[['RNA']] <- split(DC2[['RNA']], f = DC2$orig.ident)
DC2 <- NormalizeData(DC2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(DC2, ndims = 50)
DC2 <- RunUMAP(DC2, dims = 1:15, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(DC2, reduction = 'umap.unintegrated', split.by = 'stim')
DC2[['RNA']] <- JoinLayers(DC2[['RNA']])
DC2 <- FindNeighbors(DC2, reduction = 'pca', dims = 1:15)
DC2 <- FindClusters(DC2, resolution = 0.1)
DC2.markers <- FindAllMarkers(DC2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
saveRDS(DC2.markers, file = 'data1/DC2_markers.rds')
top30 <- DC2.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/DC2_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data1/top30_DC2_immune.xlsx")
celltype <- c('gdT_cell', 'DC')
names(celltype) <- levels(DC2)
DC2 <- RenameIdents(DC2, celltype)
DC2$DC2_cluster <- Idents(DC2)

Idents(immune_cell, cells = colnames(DC2)) <- Idents(DC2)

NEU <- subset(immune_cell, idents = 'NEU')
DimPlot(NEU)
NEU[['RNA']] <- split(NEU[['RNA']], f = NEU$orig.ident)
NEU <- NormalizeData(NEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(NEU, ndims = 20)
NEU <- RunUMAP(NEU, dims = 1:10, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(NEU, reduction = 'umap.unintegrated', split.by = 'stim')
NEU <- IntegrateLayers(object = NEU, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                               verbose = T)
NEU[['RNA']] <- JoinLayers(NEU[['RNA']])
NEU <- RunUMAP(NEU, dims = 1:10, reduction = 'integrated.cca', reduction.name = 'umap.cca')
DimPlot(NEU, split.by = 'stim')
NEU <- FindNeighbors(NEU, reduction = 'integrated.cca', dims = 1:10)
NEU <- FindClusters(NEU, resolution = 0.2)
DimPlot(NEU, reduction = 'umap.unintegrated', label = T, label.size = 5, cols = cluster_colors,
        group.by = 'RNA_snn_res.0.2',split.by = 'stim')
NEU.markers <- FindAllMarkers(NEU, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
saveRDS(NEU.markers, file = 'data1/NEU_markers.rds')
top30 <- NEU.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/NEU_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data1/top30_NEU_immune.xlsx")
celltype <- c('Ccrl2+NEU', 'Ly6g+NEU')
names(celltype) <- levels(NEU)
NEU <- RenameIdents(NEU, celltype)
NEU$NEU_cluster <- Idents(NEU)

M <- subset(psoriasis, idents = 'M')
M[['RNA']] <- split(M[['RNA']], f = M$orig.ident)
M <- NormalizeData(M) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(M, ndims = 40)
M <- RunUMAP(M, dims = 1:20, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(M, reduction = 'umap.unintegrated', split.by = 'stim')
M <- IntegrateLayers(object = M, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = T)

M <- RunUMAP(M, dims = 1:20, reduction = 'integrated.harmony', reduction.name = 'umap.harmony')
DimPlot(M,reduction = 'umap.cca')
M <- FindNeighbors(M, reduction = 'integrated.harmony', dims = 1:20)
a <- seq(from = 0.1, to = 1, by = 0.1)
M <- FindClusters(M, resolution = a)
clustree(M)
Idents(M) <- 'RNA_snn_res.0.2'
DimPlot(M, reduction = 'umap.harmony', label = T, label.size = 5, cols = cluster_colors,
        group.by = 'RNA_snn_res.0.2')
M[['RNA']] <- JoinLayers(M[['RNA']])
M.markers <- FindAllMarkers(M, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
saveRDS(M.markers, file = 'data1/M_markers.rds')
top30 <- M.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/NEU_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data1/top30_M_immune.xlsx")
celltype <- c('Ccrl2+NEU', 'Ly6g+NEU')
names(celltype) <- levels(NEU)
NEU <- RenameIdents(NEU, celltype)
NEU$NEU_cluster <- Idents(NEU)

VlnPlot(M, features = c('Mrc1', 'Cxcl13', 'Cd163','Il1b','Adgre1','Cxcl2'))

exper_matrix <- as(as.matrix(M@assays$RNA$counts),'sparseMatrix')
pdata <- M@meta.data
pdata$celltype <- M@active.ident
fdata <- data.frame(gene_short_name = rownames(M), row.names = rownames(M))
pd <- new('AnnotatedDataFrame', data = pdata)
fd <- new('AnnotatedDataFrame', data = fdata)
cds <- newCellDataSet(exper_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
express_gene <- subset(M.markers, p_val_adj < 0.05)$gene
cds <- setOrderingFilter(cds, express_gene)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "celltype") +  # 替换为你的基因名
  theme_bw()


DC <- subset(psoriasis, idents = 'DC')
DC[['RNA']] <- split(DC[['RNA']], f = DC$orig.ident)
DC <- NormalizeData(DC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(DC, ndims = 40)
DC <- RunUMAP(DC, dims = 1:20, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(DC, reduction = 'umap.unintegrated', split.by = 'stim')
DC <- IntegrateLayers(object = DC, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony",
                     verbose = T)

DC <- RunUMAP(DC, dims = 1:20, reduction = 'integrated.harmony', reduction.name = 'umap.harmony')
DimPlot(DC,reduction = 'umap.harmony')
DC <- FindNeighbors(DC, reduction = 'integrated.harmony', dims = 1:20)
a <- seq(from = 0.1, to = 1, by = 0.1)
DC <- FindClusters(DC, resolution = a)
clustree(DC)
Idents(DC) <- 'RNA_snn_res.0.1'
DimPlot(DC, reduction = 'umap.harmony', label = T, label.size = 5, cols = cluster_colors,
        group.by = 'RNA_snn_res.0.1')
DC[['RNA']] <- JoinLayers(DC[['RNA']])
DC.markers <- FindAllMarkers(DC, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
saveRDS(DC.markers, file = 'data1/DC_markers.rds')
top30 <- DC.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/NEU_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data1/top30_DC_immune.xlsx")
celltype <- c('Ccrl2+NEU', 'Ly6g+NEU')
names(celltype) <- levels(NEU)
NEU <- RenameIdents(NEU, celltype)
NEU$NEU_cluster <- Idents(NEU)

VlnPlot(DC, features = c('Cd207', 'Il4i1', 'Tnf','Il1b','Tgfb1','Itgae','Ifna7'))

exper_matrix <- as(as.matrix(DC@assays$RNA$counts),'sparseMatrix')
pdata <- DC@meta.data
pdata$celltype <- DC@active.ident
fdata <- data.frame(gene_short_name = rownames(DC), row.names = rownames(DC))
pd <- new('AnnotatedDataFrame', data = pdata)
fd <- new('AnnotatedDataFrame', data = fdata)
cds <- newCellDataSet(exper_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
express_gene <- subset(M.markers, p_val_adj < 0.05)$gene
cds <- setOrderingFilter(cds, express_gene)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "celltype") +  # 替换为你的基因名
  theme_bw()






Idents(immune_cell, cells = colnames(NEU)) <- Idents(NEU)
immune_cell$celltype_2 <- immune_cell@active.ident
DimPlot(immune_cell, reduction = 'umap.cca', label = T, label.size = 5, cols = cluster_colors,
        group.by = 'celltype_2')

immune_cell$cluster_names <- Idents(immune_cell)
saveRDS(immune_cell,file = 'data1/immune_cell_round2.rds', compress = F)


immune_cell <- subset(psoriasis, idents =c('gdT_cell','DC', 'M', 'NEU', 'NK_cell', 'Th_cell', 'Mast_cell') )

celltype_ratio <- DC@meta.data %>%
  group_by(stim, RNA_snn_res.0.1) %>%#分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
celltype_ratio$RNA_snn_res.0.1 <- factor(celltype_ratio$RNA_snn_res.0.1)
ggplot(celltype_ratio, aes(x=RNA_snn_res.0.1, y=relative_freq)) +
  geom_col(aes(fill=stim), color="black", position="dodge") +
  ylab("Relative Frequency") +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text = element_text(size=11),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.y=element_text(size=10, color='black'),
        axis.text.x=element_text(size=10, color='black', angle=45, hjust=1),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
openxlsx::write.xlsx(celltype_ratio, file = "data1/celltype_ratio_immune.xlsx")

Idents(psoriasis, cells = colnames(immune_cell)) <- Idents(immune_cell)
psoriasis$celltype_2 <- psoriasis@active.ident
DimPlot(psoriasis, reduction = 'umap.cca', group.by= 'celltype_2',raster = F, cols = cluster_colors)
saveRDS(psoriasis, 'data1/round1.rds', compress = F)
