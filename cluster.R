library(Seurat)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tibble)
library(sctransform)
library(future)
library(clustree)
library(openxlsx)
library(RColorBrewer)
library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpubr)
library(pathview)
library(enrichplot)
library(CellChat)


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
DimPlot(psoriasis, reduction = 'umap.cca', group.by= 'celltype_1', label = F,
        label.size = 5, raster = T)
saveRDS(psoriasis, file = 'data3/psoriasis0.rds', compress = F)

psoriasis <-readRDS('data3/psoriasis0.rds')
immune_cell <- subset(psoriasis, idents = c('T cell/NK cell','Mono-Macro cell/DC','Neutrophil'))

immune_cell[['RNA']] <- split(immune_cell[['RNA']], f = immune_cell$orig.ident)
immune_cell <- NormalizeData(immune_cell) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(immune_cell, ndims = 50)
immune_cell <- RunUMAP(immune_cell, dims = 1:30, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(immune_cell, reduction = 'umap.unintegrated', group.by = 'orig.ident')
immune_cell <- IntegrateLayers(object = immune_cell, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                               verbose = T)
immune_cell[['RNA']] <- JoinLayers(immune_cell[['RNA']])
immune_cell <- RunUMAP(immune_cell, dims = 1:30, reduction = 'integrated.cca', reduction.name = 'umap.cca')
DimPlot(immune_cell, reduction = "umap.cca")

immune_cell <- FindNeighbors(immune_cell, reduction = 'integrated.cca', dims = 1:30)
a <- seq(from = 0.1, to = 1, by = 0.1)
immune_cell <- FindClusters(immune_cell, resolution = a)
clustree(immune_cell)
DimPlot(immune_cell, reduction = "umap.cca",group.by = 'RNA_snn_res.0.6',label = T,
        label.size = 5)
FeaturePlot(immune_cell, features = 'Fcer1a',reduction = 'umap.cca')
Idents(immune_cell) <- 'RNA_snn_res.0.6'

immune.cell.markers <- FindAllMarkers(immune_cell, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
all.markers.list <- split(immune.cell.markers, f = immune.cell.markers$cluster)
openxlsx::write.xlsx(all.markers.list, file = "data3/all_marker_list_0.6_immune.xlsx")
saveRDS(immune.cell.markers, file = 'data1/0.8_immnue_markers.rds')
top30 <- immune.cell.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/0.8_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data3/top30_0.6_immune.xlsx")

celltype <- c('Neutrophil', 'Monocyte/Macrophage', 'γδT cell', 'γδT cell', 'dl', 'Dendritic cell', 'Dendritic cell', 'NK cell', 'Dendritic cell', 'Monocyte/Macrophage',
              'αβT cell', 'Neutrophil', 'Dendritic cell', 'dl')
names(celltype) <- levels(immune_cell)
immune_cell <- RenameIdents(immune_cell, celltype)
immune_cell$celltype_1 <- immune_cell@active.ident
immune_cell <- subset(immune_cell, idents = 'dl', invert = T)
DimPlot(immune_cell, reduction = 'umap.cca', label = T,
        label.size = 5, raster = T, group.by = 'celltype_1')

DC <- subset(immune_cell, idents = 'Dendritic cell')
DimPlot(DC,reduction = 'umap.cca')
DC[['RNA']] <- split(DC[['RNA']], f = DC$orig.ident)
DC <- NormalizeData(DC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(DC, ndims = 30)
DC <- RunUMAP(DC, dims = 1:15, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(DC, reduction = 'umap.unintegrated', split.by = 'stim')
DC <- IntegrateLayers(object = DC, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                      verbose = T, k.weight = 50)
DC[['RNA']] <- JoinLayers(DC[['RNA']])
DC <- RunUMAP(DC, dims = 1:15, reduction = 'integrated.cca', reduction.name = 'umap.cca')
DimPlot(DC, reduction = "umap.cca")

DC <- FindNeighbors(DC, reduction = 'integrated.cca', dims = 1:15)
DC <- FindClusters(DC, resolution = 0.1)
DimPlot(DC, reduction = "umap.cca",group.by = 'RNA_snn_res.0.1',label = T,
        label.size = 5)
FeaturePlot(DC, features = 'Il12b',reduction = 'umap.cca')
Idents(DC) <- 'RNA_snn_res.0.1'
DC.markers <- FindAllMarkers(DC, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
all.markers.list <- split(DC.markers, f = DC.markers$cluster)
openxlsx::write.xlsx(all.markers.list, file = "data3/all_marker_list_0.1_DC.xlsx")
saveRDS(immune.cell.markers, file = 'data1/0.8_immnue_markers.rds')
top30 <- DC.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/0.8_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data3/top30_0.1_DC.xlsx")

celltype <- c('Dendritic cell','Dendritic cell','dl','Dendritic cell','Dendritic cell')
names(celltype) <- levels(DC)
DC <- RenameIdents(DC, celltype)
DC$celltype_DC <- DC@active.ident
Idents(immune_cell, cells = colnames(DC)) <- Idents(DC)

M <- subset(immune_cell,idents = 'Monocyte/Macrophage')
M[['RNA']] <- split(M[['RNA']], f = M$orig.ident)
M <- NormalizeData(M) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(M, ndims = 40)
M <- RunUMAP(M, dims = 1:20, reduction = 'pca', reduction.name = 'umap.unintegrated')
DimPlot(M, reduction = 'umap.unintegrated')
M <- IntegrateLayers(object = M, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony",
                     verbose = T)
M[['RNA']] <- JoinLayers(M[['RNA']])
M <- RunUMAP(M, dims = 1:20, reduction = 'integrated.harmony', reduction.name = 'umap.harmony')
DimPlot(M, reduction = "umap.harmony")
M <- FindNeighbors(M, reduction = 'integrated.harmony', dims = 1:20)
a <- seq(from = 0.1, to = 1, by = 0.1)
M <- FindClusters(M, resolution = a)
DimPlot(M, reduction = 'umap.harmony', label = T, label.size = 5, 
        group.by = 'RNA_snn_res.0.1')
FeaturePlot(M, features = '', reduction = 'umap.harmony')
Idents(M) <- 'RNA_snn_res.0.1'

M.markers <- FindAllMarkers(M, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
saveRDS(M.markers, file = 'data1/M_markers.rds')
top30 <- M.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/NEU_immune_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data3/top30_M_immune.xlsx")
celltype <- c('Monocyte/Macrophage','Monocyte/Macrophage','Monocyte/Macrophage','Monocyte/Macrophage')
names(celltype) <- levels(M)
M <- RenameIdents(M, celltype)
M$celltype_M <- M@active.ident
Idents(immune_cell, cells = colnames(M)) <- Idents(M)

immune_cell$celltype_2 <- immune_cell@active.ident
immune_cell <- subset(immune_cell, idents = 'dl', invert = T)
saveRDS(immune_cell, file = 'data3/immune1.rds',compress = F)
DimPlot(immune_cell, reduction = 'umap.cca', label = T, label.size = 5, 
        group.by = 'celltype_2')
a <- subset(immune_cell, idents = 'dl', invert = T)
DimPlot(a, reduction = 'umap.cca', label = T, label.size = 5, 
        group.by = 'celltype_1')

celltype_ratio <- a@meta.data %>%
  group_by(stim, celltype_M) %>%#分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
celltype_ratio$celltype_M <- factor(celltype_ratio$celltype_M)
ggplot(celltype_ratio, aes(x=celltype_M, y=relative_freq)) +
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
openxlsx::write.xlsx(celltype_ratio, file = "data3/celltype_ratio_immune.xlsx")

sample_table <- as.data.frame(table(a@meta.data$ident,a@meta.data$celltype_1))
names(sample_table) <- c("Samples","celltype","CellNumber")
plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  #scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample


Idents(psoriasis, cells = colnames(immune_cell)) <- Idents(immune_cell)
psoriasis$celltype_2 <- psoriasis@active.ident
psoriasis <- subset(psoriasis, idents = 'dl', invert = T)
DimPlot(psoriasis, reduction = 'umap.cca', group.by= 'celltype_2',raster = T, cols = brewer.pal(11,'Paired'),label = F, label.size = 5)
saveRDS(psoriasis, 'data3/round1.rds', compress = F)

psoriasis <- readRDS('data3/round1.rds')
celltype_ratio <- psoriasis@meta.data %>%
  group_by(stim, celltype_2) %>%#分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))
celltype_ratio$celltype_2 <- factor(celltype_ratio$celltype_2)
ggplot(celltype_ratio, aes(x=celltype_2, y=relative_freq)) +
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
openxlsx::write.xlsx(celltype_ratio, file = "data3/celltype_ratio.xlsx")

sample_table <- as.data.frame(table(psoriasis@meta.data$ident,psoriasis@active.ident))
names(sample_table) <- c("Samples","celltype","CellNumber")
plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values= brewer.pal(12,'Paired')) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
plot_sample

FeaturePlot(psoriasis, reduction = 'umap.cca', features = 'Acta2',raster = T)
marker <- c('Adgre1','Itgam','H2-Ea','H2-Aa','Csf3r','S100a9','Trdc','Trgc1','Ncr1','Nkg7',
            'Cd28','Trac','Krt10','Krt1','Krt5','Krt14','Col1a1','Pdgfra',
            'Pecam1','Cdh5','Acta2','Tagln')
marker <- as.character(marker)
VlnPlot(psoriasis, features = marker, stack = TRUE, sort = F, flip = T) +
  theme(legend.position = "none") + ggtitle("Identity on y-axis")
DotPlot(psoriasis, features = marker)
