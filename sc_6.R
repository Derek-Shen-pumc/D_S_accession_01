psoriasis$celltype.stim <- paste(psoriasis$celltype_2, psoriasis$stim, sep = "_")
Idents(psoriasis) <- "celltype.stim"
allDE <- FindMarkers(psoriasis, ident.1 = "gdT_cell_Pregnant", ident.2 = "gdT_cell_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data2/DE_gdT.xlsx")
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
openxlsx::write.xlsx(a, file = "data2/EGG_gdT.xlsx")

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
openxlsx::write.xlsx(a, file = "data2/go_gdT_down.xlsx")



allDE <- FindMarkers(psoriasis, ident.1 = "DC_Pregnant", ident.2 = "DC_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data2/DE_DC.xlsx")
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
openxlsx::write.xlsx(a, file = "data2/EGG_DC.xlsx")

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
openxlsx::write.xlsx(a, file = "data2/go_gdT_down.xlsx")



allDE <- FindMarkers(psoriasis, ident.1 = "M_Pregnant", ident.2 = "M_Control", min.pct=0.1, verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data2/DE_M.xlsx")
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
openxlsx::write.xlsx(a, file = "data2/EGG_M.xlsx")

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
openxlsx::write.xlsx(a, file = "data2/go_M_down.xlsx")


allDE <- FindMarkers(psoriasis, ident.1 = "KC_Pregnant", ident.2 = "KC_Control", verbose = FALSE)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data2/DE_KC.xlsx")
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
openxlsx::write.xlsx(a, file = "data2/EGG_KC.xlsx")

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
openxlsx::write.xlsx(a, file = "data2/go_KC_down.xlsx")


for (pathway.id in EGG@result$ID){
  print(pathway.id)
  pathview(gene,
           pathway.id = pathway.id,
           species    = "mmu"
  )
}
pathviewu(gene,
         pathway.id = 'mmu04723',
         species    = "mmu"
)


psoriasis$celltype.stim <- paste(psoriasis$celltype_2, psoriasis$stim, sep = "_")
Idents(psoriasis) <- "celltype.stim"
allDE <- FindMarkers(psoriasis, ident.1 = "FB_Pregnant", ident.2 = "FB_Control", verbose = FALSE, min.pct = 0.25)
DE <- allDE %>% 
  rownames_to_column(var = 'gene') %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(abs(avg_log2FC) >0.5) %>%
  arrange(desc(avg_log2FC)) %>% 
  distinct(gene,.keep_all = T) %>% 
  column_to_rownames(var = 'gene')
a <- rownames_to_column(DE, var = 'gene')
openxlsx::write.xlsx(a, file = "data2/DE_FB.xlsx")
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
                  pvalueCutoff = 0.05)
barplot(EGG)
dotplot(EGG)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data2/EGG_gdT.xlsx")

up <- gene[gene > 0]
down <- gene[gene <0]
go <- enrichGO(gene = names(up),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "ALL",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
barplot(go)
dotplot(go)
a <- as.data.frame(go)
openxlsx::write.xlsx(a, file = "data2/go_FB_up.xlsx")
