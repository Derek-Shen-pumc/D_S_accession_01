library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(pheatmap)
library(ggsankey)
library(cols4all)
library(cowplot)
library(circlize)
library(colorspace)
library(pheatmap)


gset <- read.xlsx('data/count_data.xlsx')
cts <- gset[,1:10] %>%   
  ## 重新排列
  select(id,M1, M2, M3, I1, I2, I3) %>%  
  ## rowMeans求出行的平均数(这边的.代表上面传入的数据)
  ## .[,-1]表示去掉出入数据的第一列，然后求行的平均值
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  ##排除所有样本均不表达的基因
  filter(rowMean >= 10 ) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，symbol留下第一个
  distinct(id,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  select(-rowMean) %>% 
  ## 列名转行名
  column_to_rownames("id")
coldata <- data.frame(
  group = factor(c(rep('IMQ',3),rep('Pregnant',3))),  # 分组，对应示例的 Treatment、Control
  row.names = colnames(cts)  # 样本名，和计数矩阵行名一致
)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ group)
dds <- DESeq(dds)
rld <- rlog(dds)



expr_matrix <- assay(rld)
pca_data <- prcomp(t(expr_matrix), scale. = TRUE)
pca_df <- data.frame(
  Sample = colnames(dds),  # 样本名
  PC1 = pca_data$x[, 1],   # 第一主成分
  PC2 = pca_data$x[, 2],   # 第二主成分
  Group = colData(dds)$group  # 假设 colData(dds) 里有分组信息列叫 group，替换成你实际的分组列名
)
percent_var <- round(100 * summary(pca_data)$importance[2, 1:2], 1)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) + 
  geom_point(size = 3) +  # 画散点
  xlab(paste0("PC1 (", percent_var[1], "%)")) +  # x 轴标题带解释度
  ylab(paste0("PC2 (", percent_var[2], "%)")) +  # y 轴标题带解释度
  theme_bw() +  # 简洁背景主题
  ggtitle("PCA of DESeq2 Normalized Counts") +  # 标题
  scale_color_manual(values = c("IMQ" = "orange", "Pregnant" = "purple")) +  # 自定义颜色，按需改
  scale_shape_manual(values = c(17, 19))  # 自定义点形状，按需改



normalized_counts <- counts(dds, normalized = TRUE)
plotData <- normalized_counts %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene") %>% 
  filter(gene%in%c('Esr1','Pgr','Prlr','Il17a','Il17f','Tnf','Il23a','Il1b','Il6') == T) %>% 
  tidyr::pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "expression"
  ) %>% 
  left_join(
    coldata %>% 
      tibble::rownames_to_column("sample"),
    by = "sample"
  )


ggplot(filter(plotData, gene == 'Il17a'),aes(x = group, y = expression))+
  geom_boxplot(width = 0.6, outlier.shape = NA,color= c('red','blue')) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = c(rep('red',3),rep('blue',3))) +
  stat_pvalue_manual(
    statTest <- compare_means(
      expression ~ group, 
      data = filter(plotData, gene == 'Il17a'),
      method = "t.test",  # 非参数检验
    ),
    y.position = max(plotData$expression)*0.02,
    label = "p.signif",  # 使用校正后的p值标记
    tip.length = 0.01,
  )+
  # 自定义配色（三组颜色）
  labs(
    x = NULL,
    y = "Normalized Expression",
    title = paste0("Il17a")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+ylim(0,150)

res <- results(dds, contrast = c('group', 'Pregnant', 'IMQ'))
plotMA(res, ylim=c(-5,5))
res <- res %>% as.data.frame() %>% rownames_to_column('gene_id')
DE <- res %>% filter(abs(log2FoldChange) > log2(1.5)) %>% filter(padj < 0.05)
openxlsx::write.xlsx(DE, file = "data/DEGs.xlsx")

res$diffexpressed <- "Not sig"
res$diffexpressed[res$log2FoldChange > log2(1.5) & res$padj < 0.05] <- "Up"
res$diffexpressed[res$log2FoldChange < -log2(1.5) & res$padj < 0.05] <- "Down"
res$diffexpressed <- factor(res$diffexpressed, levels = c("Not sig", "Down", "Up"))
colors <- c("Not sig" = "#CCCCCC", "Down" = "#377EB8", "Up" = "#E41A1C")

volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = diffexpressed), alpha = 0.6, size = 2) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "DESeq2 Differential Expression Results",
    subtitle = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted p-value)",
    color = "Expression"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "#EEEEEE"),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "gray60")

# 标记top基因（可选）
# 选择最显著的基因进行标记
top_genes <- res %>%
  filter(padj < 0.05) %>%
  arrange(abs(log2FoldChange)) %>%
  tail(20)  # 选择log2FC最高的10个基因

# 添加标签
if(nrow(top_genes) > 0) {
  volcano_plot <- volcano_plot +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene_id),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.5
    )
}
pdf('data/volcanoplot.pdf', height = 5.5, width = 4.5)
volcano_plot
dev.off()
a <- c('Esr1','Pgr','Prlr','Il17a','Il17f','Tnf','Il23a','Il1b','Il6','Krt10',
       'Krt16','Krt10','Cxcl10','Tnf','Traf1','Cxcl5','Mmp9','Ccl5',
       'Tnfrsf12a','Krt14','Krt5','Mmp2','Krt7','Tnfsf12','Oxtr',
       'Traf5','Ccl2','Ghr','Il17c','Cldn17','Il22','Tgfb1','Gas6')
b <- c()
for (i in a) {
  if (i %in% DE$gene_id == T) {
    b <- c(b,i)
  }
}
heatdata <- assay(rld)[rownames(assay(rld))%in%a==T,]
annotation_col <- data.frame(group=coldata$group)
rownames(annotation_col) <- rownames(coldata)
pheatmap(heatdata)
pdf('data/heatmap.pdf', height = 6.5, width = 4.5)
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = TRUE,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
)
dev.off()

gene <- DE %>% dplyr::select(gene_id,log2FoldChange)
up <- filter(gene, log2FoldChange > 0)
down <- filter(gene, log2FoldChange < 0)
id <- bitr(up$gene_id, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(up, `gene_id` %in% id$SYMBOL == T)
gene <- gene$log2FoldChange
names(gene) <- id$ENTREZID
EGG <- enrichKEGG(gene = names(gene),
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data/upreguated_KEGG_bulk.xlsx")
dotplot(EGG,title = 'KEGG pathways enriched by upregulated DEGs')

id <- bitr(down$gene_id, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(down, `gene_id` %in% id$SYMBOL == T)
gene <- gene$log2FoldChange
names(gene) <- id$ENTREZID
EGG <- enrichKEGG(gene = names(gene),
                  organism = 'mmu',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
a <- as.data.frame(EGG)
openxlsx::write.xlsx(a, file = "data/downreguated_KEGG_bulk.xlsx")
dotplot(EGG,title = 'KEGG pathways enriched by downregulated DEGs')


id <- bitr(up$gene_id, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(up, `gene_id` %in% id$SYMBOL == T)
gene <- gene$log2FoldChange
names(gene) <- id$ENTREZID

go <- enrichGO(gene = names(gene),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "BP",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
a <- as.data.frame(go)
openxlsx::write.xlsx(a, file = "data/upreguated_GO_bulk.xlsx")
dotplot(go, showCategory=10, title = 'GO terms enriched by upregulated DEGs')
barplot(go)

id <- bitr(down$gene_id, 
           fromType="SYMBOL", 
           toType="ENTREZID", 
           OrgDb="org.Mm.eg.db")
gene <- filter(down, `gene_id` %in% id$SYMBOL == T)
gene <- gene$log2FoldChange
names(gene) <- id$ENTREZID

go <- enrichGO(gene = names(gene),  # 分析的基因列表
               OrgDb = 'org.Mm.eg.db',   # 物种对应的注释数据库，如org.Hs.eg.db
               ont = "BP",              # ont可为BP、MF和CC，CC细胞组件，MF分子功能，BP生物学过程
               pAdjustMethod = "BH",     # p值校正方法，Benjamini-Hochberg
               pvalueCutoff = 0.05,       # p值阈值
               qvalueCutoff = 0.05,       # 调整后的p值（q值）阈值
               readable = TRUE)          # 基因ID转换为基因名
openxlsx::write.xlsx(a, file = "data/downreguated_GO_bulk.xlsx")
dotplot(go, showCategory=10, title = 'GO terms enriched by upregulated DEGs')
dotplot(go, showCategory=10, title = 'GO terms enriched by upregulated DEGs')
