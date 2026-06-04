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

psoriasis <- readRDS('data1/round1.rds')
gdT_cell <- subset(psoriasis, idents = c('gdT_cell'))
pseudo <- AggregateExpression(gdT_cell, assays = "RNA", return.seurat = T, group.by = c('stim','group'))
tail(Cells(pseudo))
pseudo$stim_group <- paste(pseudo$stim, pseudo$group, sep = "_")
Idents(pseudo) <- "stim_group"
allDE <- FindMarkers(object = pseudo, 
                     ident.2 = c("Control_C1","Control_C2" ,"Control_C3"), 
                     ident.1 = c("Pregnant_P1","Pregnant_P2","Pregnant_P3"),
                     test.use = "DESeq2")
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data1/DE_gdT.xlsx")
### 获得基因列表
gene <- a %>% dplyr::select(gene, avg_log2FC)
#基因名称转换，返回的是数据框
id <- bitr(gene$gene, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(gene, `gene` %in% id$SYMBOL == T)
gene <- gene$avg_log2FC
names(gene) <- id$ENTREZID
EGG <- enrichKEGG(gene = names(gene),
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
barplot(EGG)
dotplot(EGG)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data1/EGG_gdT.xlsx")

up <- gene[gene > 0]
down <- gene[gene <0]
go <- enrichGO(gene = names(down),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "ALL",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
barplot(go)
dotplot(go)
a <- as.data.frame(go)
openxlsx::write.xlsx(a, file = "data1/go_gdT_down.xlsx")


for (pathway.id in EGG@result$ID){
  print(pathway.id)
  pathview(gene,
           pathway.id = pathway.id,
           species    = "mmu"
  )
}




DC <- subset(psoriasis, idents = c('DC'))
pseudo <- AggregateExpression(DC, assays = "RNA", return.seurat = T, group.by = c('stim','group'))
tail(Cells(pseudo))
pseudo$stim_group <- paste(pseudo$stim, pseudo$group, sep = "_")
Idents(pseudo) <- "stim_group"
allDE <- FindMarkers(object = pseudo, 
                     ident.2 = c("Control_C1","Control_C2" ,"Control_C3"), 
                     ident.1 = c("Pregnant_P1","Pregnant_P2","Pregnant_P3"),
                     test.use = "DESeq2")
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data1/DE_DC.xlsx")
### 获得基因列表
gene <- a %>% dplyr::select(gene, avg_log2FC)
#基因名称转换，返回的是数据框
id <- bitr(gene$gene, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(gene, `gene` %in% id$SYMBOL == T)
gene <- gene$avg_log2FC
names(gene) <- id$ENTREZID
EGG <- enrichKEGG(gene = names(gene),
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
barplot(EGG)
dotplot(EGG)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data1/EGG_DC.xlsx")

up <- gene[gene > 0]
down <- gene[gene <0]
go <- enrichGO(gene = names(down),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "ALL",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
barplot(go)
dotplot(go)
a <- as.data.frame(go)
openxlsx::write.xlsx(a, file = "data1/go_DC_down.xlsx")


for (pathway.id in EGG@result$ID){
  print(pathway.id)
  pathview(gene,
           pathway.id = pathway.id,
           species    = "mmu"
  )
}



M <- subset(psoriasis, idents = c('M'))
pseudo <- AggregateExpression(M, assays = "RNA", return.seurat = T, group.by = c('stim','group'))
tail(Cells(pseudo))
pseudo$stim_group <- paste(pseudo$stim, pseudo$group, sep = "_")
Idents(pseudo) <- "stim_group"
allDE <- FindMarkers(object = pseudo, 
                     ident.2 = c("Control_C1","Control_C2" ,"Control_C3"), 
                     ident.1 = c("Pregnant_P1","Pregnant_P2","Pregnant_P3"),
                     test.use = "DESeq2")
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data1/DE_M.xlsx")
### 获得基因列表
gene <- a %>% dplyr::select(gene, avg_log2FC)
#基因名称转换，返回的是数据框
id <- bitr(gene$gene, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(gene, `gene` %in% id$SYMBOL == T)
gene <- gene$avg_log2FC
names(gene) <- id$ENTREZID
EGG <- enrichKEGG(gene = names(gene),
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
barplot(EGG)
dotplot(EGG)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data1/EGG_M.xlsx")

up <- gene[gene > 0]
down <- gene[gene <0]
go <- enrichGO(gene = names(down),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "ALL",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
barplot(go)
dotplot(go)
a <- as.data.frame(go)
openxlsx::write.xlsx(a, file = "data1/go_M_down.xlsx")


for (pathway.id in EGG@result$ID){
  print(pathway.id)
  pathview(gene,
           pathway.id = pathway.id,
           species    = "mmu"
  )
}


KC <- subset(psoriasis, idents = c('KC'))
pseudo <- AggregateExpression(KC, assays = "RNA", return.seurat = T, group.by = c('stim','group'))
tail(Cells(pseudo))
pseudo$stim_group <- paste(pseudo$stim, pseudo$group, sep = "_")
Idents(pseudo) <- "stim_group"
allDE <- FindMarkers(object = pseudo, 
                     ident.2 = c("Control_C1","Control_C2" ,"Control_C3"), 
                     ident.1 = c("Pregnant_P1","Pregnant_P2","Pregnant_P3"),
                     test.use = "DESeq2")
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data1/DE_KC.xlsx")
### 获得基因列表
gene <- a %>% dplyr::select(gene, avg_log2FC)
#基因名称转换，返回的是数据框
id <- bitr(gene$gene, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(gene, `gene` %in% id$SYMBOL == T)
gene <- gene$avg_log2FC
names(gene) <- id$ENTREZID
EGG <- enrichKEGG(gene = names(gene),
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
barplot(EGG)
dotplot(EGG)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data1/EGG_KC.xlsx")

up <- gene[gene > 0]
down <- gene[gene <0]
go <- enrichGO(gene = names(down),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "ALL",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
barplot(go)
dotplot(go)
a <- as.data.frame(go)
openxlsx::write.xlsx(a, file = "data1/go_KC_down.xlsx")




my_comparisons <- list(c("Control", "Pregnant"))

VlnPlot(psoriasis, features = 'Ccl18',add.noise = F, idents = 'M', group.by = 'stim')&
  theme_bw()&
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(color = 'black',face = "bold", size = 12), 
        axis.text.y = element_text(color = 'black', face = "bold"),       
        axis.title.y = element_text(color = 'black', face = "bold", size = 15),     
        panel.grid.major = element_blank(),     
        panel.grid.minor = element_blank(),      
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),     
        panel.spacing = unit(0.12, "cm"),      
        plot.title = element_text(hjust = 0.5, face = "bold.italic"),    
        legend.position = 'none')&  
  stat_compare_means(method="t.test",hide.ns = F,   
                     comparisons = my_comparisons,                   
                     label="p.signif",                   
                     bracket.size=0.8,                  
                     tip.length=0, 
                     size=6)&  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)
  ))
