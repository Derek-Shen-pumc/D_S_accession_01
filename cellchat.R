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

psoriasis <- readRDS('data3/round1.rds') %>% subset(idents = 'Smooth muscle cell', invert = T)
pregnant <- subset(psoriasis, subset = stim == 'Pregnant')
control <- subset(psoriasis, subset = stim == 'Control')
rm(psoriasis)
gc()


cellchat <- createCellChat(object = pregnant, group.by = "ident")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB

table_genes <- CellChatDB.use$interaction

# Replace old gene names in that table
table_genes$ligand <- gsub("H2-BI", "H2-Bl", table_genes$ligand)
table_genes$receptor <- gsub("H2-BI", "H2-Bl", table_genes$receptor)
table_genes$ligand <- gsub("H2-Ea-ps", "H2-Ea", table_genes$ligand)
table_genes$receptor <- gsub("H2-Ea-ps", "H2-Ea", table_genes$receptor)

# Put back
CellChatDB.use$interaction <- table_genes

# Assign updated DB to cellchat
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
options(future.globals.maxSize= 30*1024*1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, type =  "truncatedMean",trim = 0.1,raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net.p <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.p <- netAnalysis_computeCentrality(cellchat, 
                                                slot.name = "netP")
saveRDS(cellchat.p, 'data3/cellchatp.rds', compress = F)
saveRDS(df.net.p, 'data3/dfnetp.rds', compress = F)
df.net.p <- df.net.p %>% group_by(source)
df.net.p <- split(df.net.p, f = df.net.p$source)
openxlsx::write.xlsx(df.net.p, file = "data3/dfnetp.xlsx")





cellchat <- createCellChat(object = control, group.by = "ident")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB

table_genes <- CellChatDB.use$interaction

# Replace old gene names in that table
table_genes$ligand <- gsub("H2-BI", "H2-Bl", table_genes$ligand)
table_genes$receptor <- gsub("H2-BI", "H2-Bl", table_genes$receptor)
table_genes$ligand <- gsub("H2-Ea-ps", "H2-Ea", table_genes$ligand)
table_genes$receptor <- gsub("H2-Ea-ps", "H2-Ea", table_genes$receptor)

# Put back
CellChatDB.use$interaction <- table_genes

# Assign updated DB to cellchat
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
options(future.globals.maxSize= 30*1024*1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, type =  "truncatedMean",trim = 0.1,raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net.c <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat.c <- netAnalysis_computeCentrality(cellchat, 
                                            slot.name = "netP")
saveRDS(cellchat.c, 'data3/cellchatc.rds', compress = F)
saveRDS(df.net.c, 'data3/dfnetc.rds', compress = F)
df.net.c <- df.net.c %>% group_by(source)
df.net.c <- split(df.net.c, f = df.net.c$source)
openxlsx::write.xlsx(df.net.c, file = "data3/dfnetc.xlsx")


object.list <- list(c = cellchat.c, p = cellchat.p)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Stratum basale",signaling.exclude = "MIF")

pdf('figure/flowchangeM.pdf', height = 6.5, width = 4.5)
rankNet(cellchat, mode = "comparison", sources.use = 'Monocyte/Macrophage',stacked = T, do.stat = TRUE)
dev.off()

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

pdf("data3/inpattern.pdf", width=4.5, height=6.5)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 30, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 30, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

netVisual_bubble(cellchat, sources.use = 'Fibroblast', targets.use = c('Monocyte/Macrophage','Dendritic cell'),comparison = c(1, 2), angle.x = 45)


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "p"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "p",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "c",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 'Fibroblast',  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 'Fibroblast',  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


pathways.show <- c("ITGB2") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("c", "p")) # set factor level
plotGeneExpression(cellchat, signaling = "GDF", split.by = "datasets", colors.ggplot = T)



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "p"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "p",ligand.logFC = 0.2, receptor.logFC = 0.2)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "c",ligand.logFC = -0.2, receptor.logFC = -0.2)
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 'Stratum basale',  slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 'Fibroblast',  slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
computeEnrichmentScore(net.up, species = 'mouse')

pathways.show <- c("GRN") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
