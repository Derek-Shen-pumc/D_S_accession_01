library(openxlsx)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggfortify)
library(ggpubr)

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
  filter(gene == 'Prlr') %>% 
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
statTest <- compare_means(
  expression ~ group, 
  data = plotData,
  method = "t.test",  # 非参数检验
)

ggplot(plotData,aes(x = group, y = expression))+
  geom_boxplot(width = 0.6, outlier.shape = NA,color= c('red','blue','green')) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = c(rep('red',3),rep('blue',3),rep('green',3))) +
  stat_pvalue_manual(
    statTest,
    y.position = max(plotData$expression) * c(1.1, 1.2, 1.3),  # 错开标记位置
    label = "p.adj",  # 使用校正后的p值标记
    tip.length = 0.01
  )+
# 自定义配色（三组颜色）
labs(
  x = NULL,
  y = "Normalized Expression (CPM)",
  title = paste0(geneName, " Expression in Control vs IMQ vs Pregnant")
) +
theme_classic(base_size = 14) +
theme(
  legend.position = "none",
  plot.title = element_text(hjust = 0.5, face = "bold")
)
    


res_PM <- results(dds, contrast = c('group', 'Pregnant', 'IMQ'))
plotMA(res_PM, ylim=c(-5,5))



a <- read.xlsx('epdermathickness.xlsx')
b <- c('group','epidermis.thickness')
colnames(a) <- b
ggplot(a,aes(x=group,y=epidermis.thickness))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = NULL)+#图例位置
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("IMQ","IMQ+pregnancy"),c("IMQ","vehicle"),c("IMQ+pregnancy","vehicle")),#设置需要比较的组
              test = t.test, ##计算方法
              size=0.8,color="black",
              y_position = max(a$epidermis.thickness) * c(1.1, 1.2, 1.3),
              map_signif_level = T)

data <- read.xlsx('pasi.xlsx')
# 计算每个分组-时间点的均值和标准差（用于绘图）
stats_data <- data %>%
  group_by(Group, Day) %>%
  summarise(
    mean_score = mean(Total.score),  # 均值
    sd_score = sd(Total.score),      # 标准差
    n = n(),                         # 样本量
    .groups = "drop"
  )

# 定义时间点和分组
days <- unique(data$Day)
groups <- unique(data$Group)

# 存储t检验结果
t_test_results <- map_dfr(days, function(day) {
  # 提取当前时间点的两组数据
  group1 <- data %>% filter(Day == day, Group == groups[1]) %>% pull(Total.score)
  group2 <- data %>% filter(Day == day, Group == groups[2]) %>% pull(Total.score)
  
  # 独立样本t检验（不假设方差齐性）
  t_result <- t.test(group1, group2, var.equal = FALSE)
  
  # 判断显著性水平
  sig_level <- case_when(
    t_result$p.value < 0.001 ~ "***",
    t_result$p.value < 0.01 ~ "**",
    t_result$p.value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  # 返回结果
  tibble(
    Day = day,
    p_value = t_result$p.value,
    significance = sig_level
  )
})

# 绘制折线图
p <- ggplot(stats_data, aes(x = factor(Day), y = mean_score, color = Group, group = Group)) +
  # 绘制折线
  geom_line(size = 1.2, 
            aes(linetype = Group)) +  # 用线型区分分组
  # 绘制数据点
  geom_point(size = 3.5, 
             aes(shape = Group)) +   # 用形状区分分组
  # 添加误差棒（标准差）
  geom_errorbar(aes(ymin = mean_score - sd_score, 
                    ymax = mean_score + sd_score),
                width = 0.15, 
                size = 0.8) +
  # 添加显著性标记
  geom_text(data = t_test_results,  # 排除不显著结果
            aes(x = factor(Day), 
                y = max(stats_data$mean_score) + max(stats_data$sd_score) * 1.1,  # 标记位置
                label = significance),
            color = "black", 
            size = 6, 
            inherit.aes = FALSE) +
  # 设置坐标轴标签和标题
  labs(x = "Day", 
       y = "Total score", 
       color = "Group",
       linetype = "Group",
       shape = "Group") +
  # 设置主题
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.background =element_blank(), #背景
    axis.line=element_line(),
    panel.grid = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    )

# 显示图形
print(p)
