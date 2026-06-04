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
psoriasis$celltype.stim <- paste(psoriasis$celltype_2, psoriasis$stim, sep = "_")
Idents(psoriasis) <- "celltype.stim"
allDE <- FindMarkers(psoriasis, ident.1 = "Fibroblast_Pregnant", ident.2 = "Fibroblast_Control", verbose = FALSE)
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
openxlsx::write.xlsx(a, file = "data3/GO_fibroblast.xlsx")
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
openxlsx::write.xlsx(a, file = "data2/go_gdT_down.xlsx")

allDE <- FindMarkers(psoriasis, ident.1 = "Monocyte/Macrophage_Pregnant", ident.2 = "Monocyte/Macrophage_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data3/DE_M.xlsx")
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
openxlsx::write.xlsx(a, file = "data3/GO_monocyte_macrophage.xlsx")
a$GeneRatio=apply(a,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  GeneRatio
})

dat <- a
dat$Cluster <-  ifelse(dat$Cluster =='up',1,-1)
dat <- dat %>% group_by(Cluster) %>% top_n(n=10, wt = GeneRatio)
dat$GeneRatio <- -log10(dat$GeneRatio)
dat$GeneRatio <- dat$GeneRatio*dat$Cluster
dat=dat[order(dat$GeneRatio,decreasing = F),]
up<- dat[which(dat$Cluster > 0),]
down<- dat[which(dat$Cluster < 0),]

ggplot(dat,
       aes(y = reorder(Description,order(GeneRatio, decreasing = F)), x = GeneRatio,fill=Cluster)) + #数据映射
  geom_col()+ #绘制添加条形图
  theme_bw()+
  scale_fill_gradient(low="#3685af",high="#E64B35B2",guide = FALSE)+
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_line()
  )+
  geom_text(data = up,
            aes(x = -0.01, y = Description, label = Description),
            size= 3.5,
            hjust= 1)+ #标签右对齐
  geom_text(data = down,
            aes(x = 0.01, y = Description, label = Description),
            size= 3.5,
            hjust= 0)+ #标签左对齐
  labs(x = '-log10GeneRatio', y = 'GO Terms') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整

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
openxlsx::write.xlsx(a, file = "data2/go_gdT_down.xlsx")

allDE <- FindMarkers(psoriasis, ident.1 = "Dendritic cell_Pregnant", ident.2 = "Dendritic cell_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data3/DE_DC.xlsx")
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

allDE <- FindMarkers(psoriasis, ident.1 = "Stratum spinosum_Pregnant", ident.2 = "Stratum spinosum_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data3/DE_SS.xlsx")
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
openxlsx::write.xlsx(a, file = "data3/GO_stratum_spinosum_keratinocyte.xlsx")
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
up<- dat[which(dat$Cluster > 0),]
down<- dat[which(dat$Cluster < 0),]

ggplot(dat,
       aes(y = reorder(Description,order(GeneRatio, decreasing = F)), x = GeneRatio,fill=Cluster)) + #数据映射
  geom_col()+ #绘制添加条形图
  theme_bw()+
  scale_fill_gradient(low="#3685af",high="#E64B35B2",guide = FALSE)+
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_line()
  )+
  geom_text(data = up,
            aes(x = -0.01, y = Description, label = Description),
            size= 3.5,
            hjust= 1)+ #标签右对齐
  geom_text(data = down,
            aes(x = 0.01, y = Description, label = Description),
            size= 3.5,
            hjust= 0)+ #标签左对齐
  labs(x = '-log10GeneRatio', y = 'GO Terms') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整



allDE <- FindMarkers(psoriasis, ident.1 = "Stratum basale_Pregnant", ident.2 = "Stratum basale_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data3/DE_SB.xlsx")
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
openxlsx::write.xlsx(a, file = "data3/GO_stratum_basale_keratinocyte.xlsx")
a$GeneRatio=apply(a,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  GeneRatio
})

dat <- a
dat$Cluster <-  ifelse(dat$Cluster =='up',1,-1)
dat <- dat %>% group_by(Cluster) %>% top_n(n=10, wt = GeneRatio)
dat$GeneRatio <- -log10(dat$GeneRatio)
dat$GeneRatio <- dat$GeneRatio*dat$Cluster
dat=dat[order(dat$GeneRatio,decreasing = F),]
up<- dat[which(dat$Cluster > 0),]
down<- dat[which(dat$Cluster < 0),]

ggplot(dat,
       aes(y = reorder(Description,order(GeneRatio, decreasing = F)), x = GeneRatio,fill=Cluster)) + #数据映射
  geom_col()+ #绘制添加条形图
  theme_bw()+
  scale_fill_gradient(low="#3685af",high="#E64B35B2",guide = FALSE)+
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_line()
  )+
  geom_text(data = up,
            aes(x = -0.01, y = Description, label = Description),
            size= 3.5,
            hjust= 1)+ #标签右对齐
  geom_text(data = down,
            aes(x = 0.01, y = Description, label = Description),
            size= 3.5,
            hjust= 0)+ #标签左对齐
  labs(x = '-log10GeneRatio', y = 'GO Terms') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 14)) #主标题居中、字号调整








Idents(psoriasis) <- 'celltype_2'
KC <- subset(psoriasis, idents =c('Stratum spinosum','Stratum basale'))
KC$celltype.stim <- paste(KC$celltype_2, KC$stim, sep = " ")
Idents(KC) <- KC$celltype.stim
Idents(KC) <- factor(Idents(KC), levels = c('Stratum basale Pregnant','Stratum basale Control',
                                            'Stratum spinosum Pregnant','Stratum spinosum Control'))
KC <- AverageExpression(KC, group.by = 'celltype.stim',return.seurat = F)
marker <- c('Krt5','Krt14','Krt10','Krt1','Efna1','Efna4','Epha1',
            'Epha2','Ephb2','Sema3e','Il34',
            'Thbs1','Tgfbr1','Tgfbr2','Acvr1b','Acvr1','Acvr2a',
            'Grn','Gas6','Cd200','Mif','Prlr','Ghr')
a <- as.data.frame(KC$RNA)
b <- a[marker,]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(b, 1, cal_z_score))
data_subset_norm <- na.omit(data_subset_norm)
pheatmap(data_subset_norm, display_numbers = F,cluster_cols = F, cluster_rows=T
         ,   colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(512)
         ,number_format = "%.0f"
         , clustering_method = "complete",
         scale = "row")

pdf('figure/KCvln.pdf', height = 9, width = 5.5)
VlnPlot(KC, features = marker, stack = T, flip = T,raster = T, fill.by = 'ident')
dev.off()



FB <-readRDS('data3/FB.rds')
celltype <- c('Cluster 1','Cluster 2','Cluster 3')
names(celltype) <- levels(FB)
FB <- RenameIdents(FB, celltype)
FB$celltype_FB <- FB@active.ident

FB <- AverageExpression(FB, group.by = 'celltype_FB',return.seurat = F)
marker <- c('Esr1','Pgr','Prlr','Ar','Ghr','Sema3c','Thbs2','Tgfb2','Grn','Gas6','Mif','Rarres2',
            'Fgf7','Fgf10','Angptl2','Inhba','Cd47','Ccl11')
a <- as.data.frame(FB$RNA)
b <- a[marker,]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(b, 1, cal_z_score))
data_subset_norm <- na.omit(data_subset_norm)
pheatmap(data_subset_norm, display_numbers = F,cluster_cols = F, cluster_rows=T
         ,   colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(512)
         ,number_format = "%.0f"
         , main = "FB Average expression",border_color=T ,
         clustering_method = "complete",
         scale = "row")

pdf('figure/FBvln.pdf', height = 9, width = 4)
VlnPlot(FB, features = marker, stack = T, flip = T, fill.by = 'ident')
dev.off()

marker <- c('Ear2','Adgre1','Mrc1','Arg1','Cd163','Cxcl2','Cd86','Clec4e','Il1b','Mif','Osm',
  'Grn','Mertk','Cd200r1','Vsir','Entpd1','Cd48')
pdf('figure/Mvln.pdf', height = 9, width = 4)
VlnPlot(M, features = marker, stack = T, flip = T, fill.by = 'ident')
dev.off()
