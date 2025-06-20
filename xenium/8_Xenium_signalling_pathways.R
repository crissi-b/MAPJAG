library(Seurat)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(dplyr)
library(magrittr)
library(grid)
library(gridExtra)
library(ggExtra)
library(stringr)
library(RColorBrewer)
library(splitstackshape)
library(pheatmap)
library(CellChat)

library(ggalluvial)
library(patchwork)
library(NMF)
library(ggpubr)

# Load single-cell tissue object
load(file="/rds/projects/c/croftap-mapjagdata/MAPJAGv2/2306/Global/sc-tissue.RData")

#################################################################

# Subset out the cell types in the scRNA-seq dataset that were identifiable in the Xenium dataset
iden1 <- as.character(unique(tissue$clusters2312))
iden <- grep("MerTK|Cycling F|LL F|SPP1|S1008A|Cycling M|LAMP3|DC2s|bright NK|IFNG|FOXP3|central mem|GZMK|KLRB1|Cycling T|Memory B|Plasmacytoid|Plasma|CXCL12|POSTN|SOX5|CD34|Venou|Arter|Capil|Peri|Lympha", iden1)
Idents(tissue) <- "clusters2312"
niches <- subset(tissue, idents=iden1[iden])
unique(niches$clusters2312)

# Create vectors of the groups of cells belonging to each niche
iden2 <- grep("MerTK|Cycling F|LL F|SPP1|S1008A|Cycling M|LAMP3|DC2s", iden1)
iden3 <- grep("bright NK|IFNG|FOXP3|central mem|GZMK|KLRB1|Cycling T|Memory B|Plasmacytoid", iden1)
iden4 <- grep("lasma cells", iden1)
iden5 <- grep("CXCL12|POSTN|SOX5|CD34", iden1)
iden6 <- grep("Venou|Arter|Capil|Peri|Lympha", iden1)

# Make a vector to relabel the cell types as niches instead
niche2 <- c(iden2, iden3, iden4, iden5, iden6)
names <- c(rep("Lining layer", length(iden2)), 
           rep("Myelo-lymphoid", length(iden3)), 
           rep("Plasma cell", length(iden4)),
           rep("Sublining stroma", length(iden5)),
           rep("Vascular", length(iden6)))

# Rename the cells of the scRNA-seq to give a new ident of niches
name_mapping <- setNames(names, iden1[niche2])
Idents(niches) <- "clusters2312"
niches <- RenameIdents(niches, name_mapping)
niches$groups <- niches@active.ident
table(niches$groups, niches$clusters2312)

#################################################################

# CellChat pipeline

data.input <- GetAssayData(niches, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(niches)
meta_data_seurat <- niches@meta.data
meta_data_seurat <- meta_data_seurat %>% select(sample, groups)
cellchat <- createCellChat(object = data.input, coordinates=NULL, meta = meta_data_seurat, group.by = "groups", datatype="RNA")

levels(cellchat@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 30)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#################################################################

# Extract aggregate scores for signalling pathways across the niches

object <- cellchat
slot.name <- "netP"

centr <- slot(object, slot.name)$centr
outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(object@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)

for (i in 1:length(centr)) {
  outgoing[,i] <- centr[[i]]$outdeg
  incoming[,i] <- centr[[i]]$indeg
}

heatcol <- rev(brewer.pal(8, "RdBu"))

# Select the top 25 pathways per row
get_top_5_indices <- function(row) {
  order(row, decreasing = TRUE)[1:25]
}

# Incoming signals
out <- incoming %>% as.data.frame()
outfb <- out[, colSums(out) !=0]
out5 <- apply(outfb, 1, get_top_5_indices)
out5 <- unique(as.vector(out5))
filtered_outfb <- outfb[c("Lining layer", "Myelo-lymphoid", "Plasma cell","Sublining stroma", "Vascular"), out5]
pheatmap(t(filtered_outfb), scale = "row", col = heatcol, cluster_cols = F)

# Outgoing signals
out <- outgoing %>% as.data.frame()
outfb <- out[, colSums(out) !=0]
out5 <- apply(outfb, 1, get_top_5_indices)
out5 <- unique(as.vector(out5))
filtered_outfb <- outfb[c("Lining layer", "Myelo-lymphoid", "Plasma cell","Sublining stroma", "Vascular"), out5]
pheatmap(t(filtered_outfb), scale = "row", col = heatcol, cluster_cols = F)


#################################################################
#################################################################

# Load Xenium object

load("/rds/projects/c/croftap-mapjagx1/MAPJAG-Xenium/xenium_obj.rds")

# Create a table of gene expression per cell
check <- GetAssayData(xen, assay="RNA", slot="data")

# Create vector of genes present in Xenium probes also present in the signalling pathways identified from scRNA-seq, as EGFR is more likely to intercept a signal than send an outgoing signal it is not included

genes <- c("CXCL10", "CXCL2", "CXCL6", "CXCL9", "PECAM1", "PTPRC",
           "CCL5", "CCL19", "CCL27", "CLEC10A", "CLEC14A", "CLEC1", "SELL", "THBS2", 
           "THY1", "PTN", "SEMA3C", "FGFBP1", "FGFBP2", "ANGPT2", "CD34", "COL5A2", "HLA-DQB2",
           "TNC")

# Create a dataframe with each cell barcode and it’s designated niche
niche_id <- pt_table[,c("cell", "niches_named")]

# Select out the relevant gene names and format so each cell barcode is labelled by it’s niche identity
# Remove adipose-rich niche as no scRNA-seq data for it
check2 <- check[rownames(check) %in% genes,]
check2 <- t(check2) %>% as.data.frame()
check2$niche <- rownames(check2)
check2 <- merge(check2, niche_id, by.x="niche", by.y="cell")
check2 <- check2[!is.na(check2$niches_named), ]
check2 <- check2[!(check2$niches_named %in% "Adipose-rich"), ]

# Identify average gene expression per niche per donor
donor_id <- xen[["patient"]]
donor_id$cell <- rownames(donor_id) 
check2 <- merge(check2, donor_id, by.x="niche", by.y="cell")
check3 <- check2[,-1]
check3 <- check3 %>% mutate(niches_named = recode(niches_named, "Vascular"="Vascular/Perivascular", "Perivascular"="Vascular/Perivascular")) # simplify Vascular and Perivascular niche to be one
result <- aggregate(. ~ patient + niches_named, data = check3, FUN = mean, na.rm = TRUE)

#################################################################

# Rename any niche except the niche in question as “Other” and plot average expression 
# Outgoing signals for vascular niche
colcol <- alpha(c("blue", "grey"), 0.4)

vresult <- result %>% mutate(niches_named = recode(niches_named, "Lining layer"="Other", "Myeloid / Lymphocyte"="Other", "Plasma cell-rich"="Other", "Sublining stroma"="Other"))
vresult$niches_named <- factor(vresult$niches_named, levels=c("Vascular/Perivascular", "Other"))
setv <- theme(axis.text.x = element_blank(), axis.title = element_blank())

p1 <- ggplot(vresult, aes(x=niches_named, y=PECAM1), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="PECAM1") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p2 <- ggplot(vresult, aes(x=niches_named, y=CD34),) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="CD34") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p4 <- ggplot(vresult, aes(x=niches_named, y=ANGPT2)) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="ANGPT2") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv

#################################################################

# Outgoing signals for stromal sublining niche

colcol <- alpha(c("brown", "grey"), 0.4)

sresult <- result %>% mutate(niches_named = recode(niches_named, "Lining layer"="Other", "Myeloid / Lymphocyte"="Other", "Plasma cell-rich"="Other", "Vascular/Perivascular"="Other"))
sresult$niches_named <- factor(sresult$niches_named, levels=c("Sublining stroma", "Other"))

p5 <- ggplot(sresult, aes(x=niches_named, y=COL5A2), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="COL5A2") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p6 <- ggplot(sresult, aes(x=niches_named, y=PTN),) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="PTN") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p7 <- ggplot(sresult, aes(x=niches_named, y=THY1)) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="THY1") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p8 <- ggplot(sresult, aes(x=niches_named, y=TNC)) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="TNC") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p9 <- ggplot(sresult, aes(x=niches_named, y=THBS2)) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="THBS2") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p10 <- ggplot(sresult, aes(x=niches_named, y=SEMA3C)) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="SEMA3C") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv

#################################################################

# Outgoing signals from the myeloid niche
colcol <- alpha(c("red", "grey"), 0.4)

mresult <- result %>% mutate(niches_named = recode(niches_named, "Lining layer"="Other", "Sublining stroma"="Other", "Plasma cell-rich"="Other", "Vascular/Perivascular"="Other"))
mresult$niches_named <- factor(mresult$niches_named, levels=c("Myeloid / Lymphocyte", "Other"))

p11 <- ggplot(mresult, aes(x=niches_named, y=PTPRC), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="PTPRC") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv
p12 <- ggplot(mresult, aes(x=niches_named, y=SELL),) + geom_boxplot() + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="SELL") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv

#################################################################

# Outgoing signals from the lining layer niche

colcol <- alpha(c("yellow", "grey"), 0.4)

lresult <- result %>% mutate(niches_named = recode(niches_named, "Myeloid / Lymphocyte"="Other", "Sublining stroma"="Other", "Plasma cell-rich"="Other", "Vascular/Perivascular"="Other"))
lresult$niches_named <- factor(lresult$niches_named, levels=c("Lining layer", "Other"))

p13 <- ggplot(lresult, aes(x=niches_named, y=`HLA-DQB2`), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="HLA-DQB2") + stat_compare_means(method = "wilcox.test", label = "p.signif") & setv


#################################################################

grid.arrange(p11,p12,p1,p2,p4,p5,p6,p7,p8,p9,p10,p13, nrow=1)

#################################################################
#################################################################

# The gene pathways that were not validated

colcol <- alpha(c("yellow", "red", "grey"), 0.4)

lresult <- result %>% mutate(niches_named = recode(niches_named, "Sublining stroma"="Other", "Plasma cell-rich"="Other", "Vascular/Perivascular"="Other"))
lresult$niches_named <- factor(lresult$niches_named, levels=c("Lining layer", "Myeloid / Lymphocyte", "Other"))

pa <- ggplot(lresult, aes(x=niches_named, y=CXCL10), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="CXCL10")  & setv
pb <- ggplot(lresult, aes(x=niches_named, y=CXCL9), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="CXCL9")  & setv
pc <- ggplot(lresult, aes(x=niches_named, y=CLEC10A), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="CLEC10A")  & setv
pd <- ggplot(lresult, aes(x=niches_named, y=CLECL1), fill=niches_named) + geom_boxplot(fill=colcol) + geom_jitter(size=0.5) + theme_classic() + labs(title="CLECL1")  & setv

grid.arrange(pa, pb, pc, pd, nrow=1)
