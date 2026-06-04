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
library(ComplexHeatmap)
library(monocle)

psoriasis <- readRDS('data3/round1.rds')

FB <- subset(psoriasis, idents = 'Fibroblast')
DimPlot(FB,reduction = 'umap.cca')
FB[['RNA']] <- split(FB[['RNA']], f = FB$orig.ident)
FB <- NormalizeData(FB) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
ElbowPlot(FB, ndims = 30)
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
DimPlot(FB, reduction = 'umap.harmony', label = F, label.size = 5, group = 'celltype_FB')
FB.markers <- FindAllMarkers(FB, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
all.markers.list <- split(FB.markers, f = FB.markers$cluster)
openxlsx::write.xlsx(all.markers.list, file = "data3/all_marker_list_0.1_FB.xlsx")
saveRDS(immune.cell.markers, file = 'data1/0.8_immnue_markers.rds')
top30 <- FB.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
saveRDS(top30, file = 'data1/0.1_FB_top30.rds')
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data3/top30_0.1_FB.xlsx")
FB <- subset(FB, idents = c('0','1','4'))
celltype <- c('1','2','3')
names(celltype) <- levels(FB)
FB <- RenameIdents(FB, celltype)
FB$celltype_FB <- FB@active.ident
saveRDS(FB, file = 'data3/FB.rds', compress = F)
FeaturePlot(FB, features = c('Esr1','Pgr', 'Prlr','Ghr'),reduction = 'umap.harmony')
VlnPlot(FB, features = c('Fgf7', 'Fgf10', 'Cxcl12','Igf1','Tnfsf12','Tgfb2','Timp1','Fn1','Gas6'))

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
plot_cell_trajectory(cds, color_by = "celltype_FB") +  # 替换为你的基因名
  theme_bw()+plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")

plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
saveRDS(cds, file = 'data3/FBcds.rds',compress = F)

pdf('figure/hormonetime.pdf', height = 6, width = 4.5)
plot_genes_branched_pseudotime(cds[c('Esr1','Pgr', 'Prlr','Ghr'),],
                               branch_point = 1,
                               color_by = 'celltype_FB',
                               ncol = 1)
dev.off()

FB$celltype.stim <- paste(FB$celltype_2, FB$stim, sep = "_")
Idents(FB) <- "celltype.stim"
allDE <- FindMarkers(FB, ident.1 = "Fibroblast_Pregnant", ident.2 = "Fibroblast_Control", verbose = FALSE)

DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data3/DE_FB.xlsx")
### 获得基因列表
gene <- a %>% dplyr::select(gene, avg_log2FC)
id <- bitr(gene$gene, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(gene, `gene` %in% id$SYMBOL == T)
gene <- gene$avg_log2FC
names(gene) <- id$ENTREZID
up <- names(gene[gene > 0])
down <- names(gene[gene <0])
dels <- list(up = up, down = down)
go <- compareCluster(dels,  # 分析的基因列表
                     fun = 'enrichGO',
                     OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
                     ont = "BP",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
                     pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
                     pvalueCutoff = 0.05,       # p值阈值
                     qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
                     readable = TRUE)          # 基因ID转换为基因名
barplot(go)
dotplot(go)
a <- as.data.frame(go)
a$GeneRatio=apply(a,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  GeneRatio
})

dat <- a
dat$Cluster <-  ifelse(dat$Cluster =='up',1,-1)
dat <- dat %>% group_by(Cluster) %>% top_n(n=20, wt = GeneRatio)
dat$GeneRatio <- -log10(dat$GeneRatio)
dat$GeneRatio <- dat$GeneRatio*dat$Cluster
dat=dat[order(dat$GeneRatio,decreasing = F),]

ggplot(dat, aes(x=reorder(Description,order(GeneRatio, decreasing = F)), 
                y=GeneRatio,fill=Cluster)) + 
  geom_bar(stat="identity",width = 0.78,position = position_dodge(0.7)) + 
  scale_fill_gradient(low="#3685af",high="#E64B35B2",guide = FALSE)+
  scale_x_discrete(name ="GO Terms") +
  #scale_y_discrete(labels =Description )
  scale_y_continuous(name ="-log10GeneRatio",limits = c(-2,2)) +
  coord_flip() +  # 翻转坐标
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  theme(axis.text.y = element_blank())+ #刻度标签空白
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())+ #边框为空白
  theme(axis.text=element_text(face = "bold",size = 15),
        axis.title = element_text(face = 'bold',size = 15), #坐标轴字体
        plot.title = element_text(size = 20,hjust = 0.3), 
        legend.position = "top",panel.grid = element_line(colour = 'white'))+
  geom_text(aes(label=Description),size=3.5, hjust = ifelse(dat$GeneRatio>0,-0.05,1.1))+ # 两侧加上标签文字
  theme(
    panel.grid.major = element_blank(), # 去除主网格线
    panel.grid.minor = element_blank(), # 去除次网格线
    panel.background = element_blank(), # 去除背景色
    plot.background = element_blank() # 去除绘图区域背景色
  )


M <- subset(psoriasis, idents = 'Monocyte/Macrophage')
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
        group.by = 'RNA_snn_res.0.2')
FeaturePlot(M, features = 'Fcgr1', reduction = 'umap.harmony')
Idents(M) <- 'RNA_snn_res.0.2'

M.markers <- FindAllMarkers(M, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
saveRDS(M.markers, file = 'data1/M_markers.rds')
top30 <- M.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
top30.list <- split(top30, f = top30$cluster)
openxlsx::write.xlsx(top30.list, file = "data3/top30_M_immune.xlsx")

celltype <- c('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5')
names(celltype) <- levels(M)
M <- RenameIdents(M, celltype)
M$celltype_M <- M@active.ident
saveRDS(M, file = 'data3/M.rds', compress = F)

DimPlot(M, reduction = 'umap.harmony', label = F, label.size = 5, 
        group.by = 'celltype_M')

a <- M
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
plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")


plot_genes_branched_pseudotime(cds[c('Ear2','Adgre1','Mrc1','Arg1','Cd163','Cxcl2','Cd80','Cd86','Clec4e','Il1b','Mif','Osm','Tgfb1','Tgfb2',
                                     'Grn','Mertk','Cd200r1','Vsir','Entpd1','Cd48'),],
                               branch_point = 1,
                               color_by = 'celltype',
                               ncol = 1)

