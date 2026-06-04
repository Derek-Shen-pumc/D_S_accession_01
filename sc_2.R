library(Seurat)
library(tidyr)
library(ggplot2)
library(dplyr)
library(sctransform)
library(future)
library(clustree)
library(openxlsx)
library(monocle)


plan(multicore, workers = 4)

psoriasis <- readRDS('data1/psoriasis_normalized.rds')
options(future.globals.maxSize= 30*1024*1024^2)
psoriasis <- IntegrateLayers(object = psoriasis, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony",
                             verbose = T)
psoriasis <- IntegrateLayers(object = psoriasis, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                             verbose = T)
psoriasis[['RNA']] <- JoinLayers(psoriasis[['RNA']])
psoriasis <- RunUMAP(psoriasis, dims = 1:30, reduction = 'integrated.cca', reduction.name = 'umap.cca')
DimPlot(psoriasis, reduction = "umap.cca", split.by = 'stim')
psoriasis <- RunUMAP(psoriasis, dims = 1:40, reduction = 'integrated.harmony', reduction.name = 'umap.harmony')
DimPlot(psoriasis, reduction = "umap.cca", split.by = 'stim')
saveRDS(psoriasis,file = 'data1/psoriasis_integrated.rds',compress = F)

psoriasis <- readRDS('data1/psoriasis_integrated.rds')
psoriasis <- FindNeighbors(psoriasis, reduction = 'integrated.cca', dims = 1:30)
a <- seq(from = 0.1, to = 0.5, by = 0.1)
psoriasis <- FindClusters(psoriasis, resolution = a)
clustree(psoriasis)
Idents(psoriasis) <- 'RNA_snn_res.0.3'
DimPlot(psoriasis, reduction = 'umap.cca', group.by = 'RNA_snn_res.0.3', label = T, label.size = 5)
psoriasis.markers <- FindAllMarkers(psoriasis, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
all.markers.list <- split(psoriasis.markers, f = psoriasis.markers$cluster)
openxlsx::write.xlsx(all.markers.list, file = "data3/all_marker_list.xlsx")
saveRDS(psoriasis.markers, file = 'data1/psoriasis_markers.rds', compress = F)
top30 <- psoriasis.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/top30.rds', compress = F)
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data3/top30.xlsx")
saveRDS(psoriasis, file = 'data1/round0.rds', compress = F)
FeaturePlot(psoriasis,features = c('Acta2'),reduction = 'umap.cca')
FeaturePlot(psoriasis,features = 'Adipoq',reduction = 'umap.cca')
celltype_1 <- c('Stratum spinosum', 'Stratum spinosum','Stratum spinosum','Stratum basale','Stratum basale','Stratum basale','T cell/NK cell','Mono-Macro cell/DC','Fibroblast','Endothelial cell','Fibroblast','Neutrophil','Smooth muscle cell','dl','dl','Endothelial cell','dl')
names(celltype_1) <- levels(psoriasis)
psoriasis <- RenameIdents(psoriasis, celltype_1)
psoriasis$celltype_1 <- Idents(psoriasis)
DimPlot(psoriasis, reduction = 'umap.cca', group.by= 'celltype_1', label = T,
        label.size = 5, raster = T)
psoriasis <- subset(psoriasis, idents = 'dl', invert = T)
DimPlot(psoriasis, reduction = 'umap.cca', group.by= 'celltype_1', label = T,
        label.size = 5, raster = T)


  
  
  

cell_info <- psoriasis@meta.data %>%
  select(celltype_0, stim)
cell_prop <- cell_info %>%
  group_by(stim, celltype_0) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))
ggplot(cell_prop, aes(x = "", y = proportion, fill = celltype_0)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +  # 转换为饼图
  facet_wrap(~ stim) +  # 按组分开
  scale_fill_manual(values = cluster_colors) +
  labs(fill = "Cell Type") +  # 图例标题
  ggtitle("Cell Type Proportion by Group")  # 图表标题
saveRDS(psoriasis, file = 'data/round1.rds')


psoriasis <- readRDS('data1/round1.rds')
FB <- subset(psoriasis, idents = 'FB')
FB[['RNA']] <- split(FB[['RNA']], f = FB$orig.ident)
FB <- NormalizeData(FB) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(FB, ndims = 50)
FB <- RunUMAP(FB, dims = 1:20, reduction = 'pca', reduction.name = 'umap.unintegrated')
FB <- IntegrateLayers(object = FB, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony",
                      verbose = T)
FB[['RNA']] <- JoinLayers(FB[['RNA']])
FB <- RunUMAP(FB, dims = 1:20, reduction = 'integrated.harmony', reduction.name = 'umap.harmony')
DimPlot(FB, reduction = "umap.harmony")
FB <- FindNeighbors(FB, reduction = 'integrated.harmony', dims = 1:20)
a <- seq(from = 0.1, to = 1.0, by = 0.1)
FB <- FindClusters(FB, resolution = 0.1)
clustree(FB)
Idents(FB) <- 'RNA_snn_res.0.1'
DimPlot(FB, reduction = 'umap.harmony', label = T, label.size = 5, cols = cluster_colors, group = 'RNA_snn_res.0.1')
FB <- subset(FB, idents = c('1','3', '0'))
FB.markers <- FindAllMarkers(FB, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
all.markers.list <- split(FB.markers, f = FB.markers$cluster)
openxlsx::write.xlsx(all.markers.list, file = "data1/all_marker_list_0.1_FB.xlsx")
saveRDS(immune.cell.markers, file = 'data1/0.8_immnue_markers.rds')
top30 <- FB.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/0.1_FB_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data1/top30_0.1_FB.xlsx")
VlnPlot(FB, features = c('Fgf7', 'Fgf10', 'Cxcl12','Igf1','Tnfsf12','Tgfb2','Timp1','Fn1'))

exper_matrix <- as(as.matrix(FB@assays$RNA$counts),'sparseMatrix')
pdata <- FB@meta.data
pdata$celltype <- FB@active.ident
fdata <- data.frame(gene_short_name = rownames(FB), row.names = rownames(FB))
pd <- new('AnnotatedDataFrame', data = pdata)
fd <- new('AnnotatedDataFrame', data = fdata)
cds <- newCellDataSet(exper_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
express_gene <- subset(FB.markers, p_val_adj < 0.05)$gene
cds <- setOrderingFilter(cds, express_gene)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "celltype") +  # 替换为你的基因名
  theme_bw()


immune_cell <- subset(psoriasis, idents = c('NEU','M', 'DC', 'NK_cell/abT_cell', 'gdT_cell'))
saveRDS(immune_cell, file = 'data1/immune_cell_round0.rds', compress = F)
