
library(Seurat)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(dplyr)
library(grid)
library(gridExtra)
library(ggExtra)
library(stringr)
library(RColorBrewer)
library(SeuratObject)
library(pheatmap)


# Load tissue object and integrated samples
load(file="/rds/projects/c/croftap-1003568/analysis/setup_repeat_2/mt10_20-nc-protein-10000.rdata")
load(file="/rds/projects/c/croftap-mapjagdata/MAPJAGv2/2306/Global/sc-tissue.RData")
heatcol <- rev(brewer.pal(8, "RdBu"))

# Genes associated with extension based on a fold change above 2 and significant P value
ext <- c('CHI3L1', 'TIMD4', 'CCL18', 'TMEM176B', 'C1QB', 'C2', 'IFI27', 'Y276A02', 
         'LOC646470', 'LOC100130429', 'CXCL9', 'DHX57', 'DOCK9', 'C1QC', 'FN1', 
         'SLC2A5', 'SETDB2', 'FERMT2', 'RIPK5', 'OSBPL1A', 'SLC8A1', 'PLTP', 
         'RNASE1', 'ERAP2', 'C1QA', 'MARCO', 'NUPR1', 'ZNF770', 'BRUNOL5', 
         'IGF2BP3', 'IL31RA', 'FABP3', 'ACTR2', 'NR4A2', 'SFRS12', 'SERPING1', 
         'MEGF9', 'LOC339260')

# Expression examined in SFMC to recapitulate the analysis of the original paper in single cell
# Remove clusters with < 300 cells

Idents(PBMC1) <- "TYPE"
test <- subset(PBMC1, idents="SF")
ptnum <- as.data.frame(table(test$clusters2312))
ptnum <- ptnum[ptnum$Freq < 300, ] %>% pull(Var1) %>% unique()
Idents(test) <- "clusters2312"
test <- subset(test, idents=ptnum, invert=T)

# Take average expression of genes, reformat names of rows and columns, make heatmap of gene expression
Idents(test) <- "clusters2312"
glob <- as.data.frame(AverageExpression(object = test, assay="RNA", features=ext)) %>% t() %>% as.matrix()
rownames(glob) <- gsub("RNA.", "", rownames(glob))
rownames(glob) <- gsub("\\.\\.", "+ ", rownames(glob))
rownames(glob) <- gsub("\\.", " ", rownames(glob))
pheatmap(glob,col= heatcol,scale="column")


####################################################################################
####################################################################################

# SOX5 analyses

# GWAS loci gene expression
# Create vector of GWAS genes
genes <- c('PTPN22', 'STAT4', 'PTPN2', 'ANKRD55', 'ATXN2', 'IL6R', 'IL2RA', 'AHI1', 
           'CCR3', 'TNFSF11', 'FOXP1', 'TYK2', 'ERAP2', 'UBE2L3', 'C5orf56', 'RUNX1', 
           'ZFP36L1', 'IL2RB', 'FAS', 'LTBR', 'IL6', 'TNFSF8', 'AFF3', 'ACTA2', 
           'COG6', '13q14', 'CASC3', 'PRR5L', 'RUNX3', 'PTPRD', 'RMI2', 'JAZF1', 
           'CACNG3', 'USP12', 'RGS14', 'CD247', 'PTPN2', 'TMEM108', 'GLB1', 'CEP70', 
           'PDE9A', 'IL12A-AS1', 'FRMD4A', 'SLC47A1', 'CCR9', 'DDX60L', 'ZC3H12C', 'PTPRM')

# Subset out tissue data to only include clusters of at least 300 cells

ptnum <- as.data.frame(table(tissue$clusters2312))
ptnum <- ptnum[ptnum$Freq < 300, ] %>% pull(Var1) %>% unique()
Idents(tissue) <- "clusters2312"
smaller <- subset(tissue, idents=ptnum, invert=T)

# Take average expression of genes, reformat names of rows and columns, make heatmap of gene expression

fbBRES <- as.data.frame(AverageExpression(object = smaller, assay="RNA", features=genes)) %>% t() %>% as.matrix()
rownames(fbBRES) <- gsub("RNA.", "", rownames(fbBRES))
rownames(fbBRES) <- gsub("\\.\\.", "+ ", rownames(fbBRES))
rownames(fbBRES) <- gsub("\\.", " ", rownames(fbBRES))
fbBRES <- fbBRES %>% t()
pheatmap(fbBRES, scale="row", col= heatcol)

####################################################################################
####################################################################################

# Modules of genes in stromal cells
# Load stromal object
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/stromal-object")
Idents(seurat_obj) <- "clusters2"
DefaultAssay(seurat_obj) <- "RNA"

Idents(seurat_obj) <- "TYPE"
seurat_obj <- subset(seurat_obj, idents="Tissue")
Idents(seurat_obj) <- "clusters2"

seurat_obj -> v3
v3$sampleTYPE <-paste(v3$clusters2, v3$patient, sep=",")  
Idents(v3) <- "sampleTYPE"

#AggregateExpression will sum counts when slot is set to "counts"
cts_v3 <- AggregateExpression(v3, group.by = c("sampleTYPE"), assays = "RNA", slot = "counts", return.seurat = F)
cts_v3<-cts_v3$RNA
meta_data=colnames(cts_v3)
cts_v3<-as.data.frame(cts_v3)
meta_data<-as.data.frame(meta_data)
meta_data$to_split<-meta_data$meta_data

meta_data<-cSplit(meta_data, splitCols = "to_split", sep=",")
colnames(meta_data)<-c("all", "TYPE", "sample")
all(colnames(cts_v3) == meta_data$all)
meta_data$sample <- factor(meta_data$sample)
meta_data$all <- factor(meta_data$all)
meta_data$TYPE <- factor(meta_data$TYPE)

dds <- DESeqDataSetFromMatrix(countData = cts_v3,
                              colData = meta_data,
                              design = ~1)

dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

mingenecount <- quantile(rowSums(counts(dds)), 0.5)
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))

# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

#############################################################################

design(dds) <- formula(~ TYPE)
print(design(dds))
dds <- DESeq(dds, test = "Wald")

print(resultsNames(dds))
targetvar <- "TYPE"
comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
head(comps1)

ress <- apply(comps1, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
  res
})

res1 <- Reduce(rbind, ress)

head(res1)      

res1 %>% 
  filter(padj < 0.01) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres

subres_up<-subres[subres$log2FoldChange > 1,]
subres_down<-subres[subres$log2FoldChange < -1,]
subres<-rbind(subres_up, subres_down)

head(subres)
dim(subres)
length(unique(subres$gene))

library(ComplexHeatmap)

deseq2VST <- vst(dds, blind=T)
feats <- unique(subres$gene)
print(length(feats))

# Sub-set matrix to relevant features
deseq2VST <- assay(deseq2VST)
deseq2VST<-as.matrix(deseq2VST)
sub_vst_mat <- deseq2VST[rownames(deseq2VST) %in% feats, ]
scale_sub_vst_mat <- t(scale(t(sub_vst_mat)))
head(scale_sub_vst_mat)
dim(scale_sub_vst_mat)

ss_sm <- meta_data[, c("sample", "TYPE")]

scol <- colorRampPalette(brewer.pal(8, "Greys"))(length(unique(ss_sm$sample)))
names(scol) <- unique(ss_sm$sample) 
col_ann <- HeatmapAnnotation(df = ss_sm, col= list(sample = scol))
nrow(scale_sub_vst_mat)
topedges <- 0.05
ggnet <- Rfast::cora(t(scale_sub_vst_mat))
ggnet_gather <- reshape2::melt(ggnet, id.vars = "V1")
print(dim(ggnet_gather))

ggnet_gather %>%
  filter(value > 0) %>%
  filter(Var1 != Var2) -> ggnet_gather

print(dim(ggnet_gather))

if(nrow(ggnet_gather) > 0) {
  
  edges <- arrange(ggnet_gather, -value)
  print(head(edges))
  print(tail(edges))
  edges <- edges[seq(1, nrow(edges), by = 2), ]
  print(head(edges))
  print(tail(edges))
  qtop <- quantile(edges$value, 1-topedges)
  print(qtop)
  
  edges %>%
    filter(value > qtop) -> top_edges
  
  print(dim(top_edges))
  
  top_edges %>%
    mutate('idx' = paste(Var1, Var2, sep = "_")) -> top_edges# %>%
  #dplyr::distinct(idx, .keep_all = TRUE) -> top_edges
  
  colnames(top_edges) <- c("Source", "Target", "Weight", "Id")
  #top_edges %>%
  #  select(c(Id, Source, Target, Weight)) -> top_edges
  print(head(top_edges))
  print(dim(top_edges))
  
  # - community detection
  # - igraph definition
  g <- igraph::graph_from_data_frame(top_edges[, c("Source", "Target")],
                                     directed = FALSE)
  g <- igraph::set_edge_attr(g, "weight", value = top_edges$Weight)
  g <- igraph::set_edge_attr(g, "name", value = top_edges$Id)  
  # leiden
  leiden_mod <- igraph::cluster_leiden(g, objective_function = "modularity")
  mods <- data.frame(cbind(igraph::V(g)$name, leiden_mod$membership))
  colnames(mods) <- c("Id", "leiden_gene_cluster")
  imods <- names(table(mods$leiden_gene_cluster)[table(mods$leiden_gene_cluster) > 10])
  print(imods)
  
  if(length(imods) > 0) {
    mods <- filter(mods, leiden_gene_cluster %in% imods)
    head(mods)
    
    row_cl <- data.frame('gene' = mods$Id,
                         'gene_cluster' = paste0("K", mods$leiden_gene_cluster))
    
    row_cl %>%
      arrange(gene_cluster) %>%
      data.frame -> row_cl
    
    rownames(row_cl) <- row_cl$gene
    
  }
}


if(exists("row_cl")) {
  cpal <- rev(brewer.pal(length(unique(row_cl$gene_cluster)), "Set1"))
  names(cpal) <- unique(row_cl$gene_cluster)
  row_ann <- rowAnnotation(df = row_cl[, -1], col=list(df = cpal))    
  
  p1 <- draw(
    Heatmap(scale_sub_vst_mat[rownames(row_cl),], 
            top_annotation = col_ann, 
            right_annotation = row_ann,
            col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
            row_names_gp = gpar(fontsize = 4), 
            cluster_columns = T,
            cluster_rows = FALSE,
            show_row_names = F,
            show_column_names = F))
}

# Find Gene Ontology annotations of modules
index <- match(res1$gene, row_cl$gene)
res1$cluster <- row_cl$gene_cluster[index]
res1_cleaned<-na.omit(res1)

index <- match(res1_cleaned$gene, annotation_gs$gene_name)
res1_cleaned$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- unique(res1$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)
res1_cleaned<-na.omit(res1_cleaned)

go.results <- runGO.all(results=res1_cleaned,
                        background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                        species = "hs")
go.results <- filterGenesets(go.results)
go.results2 <- go.results[go.results$ontology =="BP",]
go.results.top <- go.results2 %>% group_by(cluster) %>% top_n(n=40, -p.val)

# Select informative terms for visualisation
goterm <- unique(unlist(str_split("GO:0035329 GO:0038202 GO:0007223 GO:0002479 GO:0002040 GO:0030199 GO:0009152 GO:0032210 GO:0140014",pattern=" ")))
include <- go.results.top[go.results.top$geneset_id %in% goterm,]
include$geneset_id <- factor(include$geneset_id, levels=c(goterm))
include <- arrange(include, cluster)
sampleEnrichmentDotplot(include, selection_col = "description", selected_genesets = unique(include$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

# Visualise expression of gene modules across clusters
seurat_obj -> PBMC1
row_cl_K1<-row_cl[row_cl$gene_cluster == "K1",]
row_cl_K2<-row_cl[row_cl$gene_cluster == "K2",]
row_cl_K3<-row_cl[row_cl$gene_cluster == "K3",]
row_cl_K4<-row_cl[row_cl$gene_cluster == "K4",]
row_cl_K5<-row_cl[row_cl$gene_cluster == "K5",]
row_cl_K6<-row_cl[row_cl$gene_cluster == "K6",]
row_cl_K7<-row_cl[row_cl$gene_cluster == "K7",]

PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K1$gene), name = "K1_mod" )
PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K2$gene), name = "K2_mod" )
PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K3$gene), name = "K3_mod" )
PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K4$gene), name = "K4_mod" )
PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K5$gene), name = "K5_mod" )
PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K6$gene), name = "K6_mod" )
PBMC1<-AddModuleScore(PBMC1, features =list(row_cl_K7$gene), name = "K7_mod" )

library(viridis)
Idents(PBMC1) <- 'clusters2'

ord <- rev(c("K7_mod1","K6_mod1", "K5_mod1","K4_mod1", "K2_mod1", "K1_mod1"))
dotplot <- DotPlot(PBMC1, feature=c(ord), cluster.idents=T) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)  +
  scale_color_gradientn(colors=c("blue", "white", "red")) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) & RotatedAxis()

dotplot <- dotplot$data
dotplot <- dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp)
dotplot <- dotplot %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)
dotplot <- dotplot[dotplot$features.plot %in% ord,]
row.names(dotplot) <- dotplot$features.plot  
row.names(dotplot) <- gsub("_mod1", "_module", row.names(dotplot))
dotplot <- dotplot[,-1] %>% as.matrix()

Heatmap(dotplot,cluster_rows = F, cluster_columns=T,   
        show_column_names = T, column_names_rot = 90)

heatcol <- rev(brewer.pal(8, "RdBu"))
pheatmap(t(dotplot), scale="column", col= heatcol, cluster_cols = F)

####################################################################################
####################################################################################

# Examine fibroblast signalling pathways from ligand receptor analysis

# Use the ligand receptor analysis summarised in the ‘subcluster visualisation script’
# Load cellchat object of all clusters
heatcol <- rev(brewer.pal(8, "RdBu"))

include <- c("POSTN+ Fibroblasts", "CD34+ Fibroblasts", "SOX5+ Fibroblasts", "CXCL12+ Fibroblasts", "Pericytes", "Lymphatics", "Venous", "KDR+ Arterial", "FLI1 + Capillary")
mini <- subsetCellChat(cellchat, idents.use=include)

pattern <- "incoming"
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

# function to isolate top 10
get_top_5_indices <- function(row) {
  order(row, decreasing = TRUE)[1:10]
}

# Top 10 incoming signals per cluster
out <- incoming %>% as.data.frame()
outfb <- out[, colSums(out) !=0]
out5 <- apply(outfb, 1, get_top_5_indices)
out5 <- unique(as.vector(out5))
pheatmap(t(out5), scale = "row", col = heatcol, cluster_cols = F)


# Top 20 outgoing signals per cluster
get_top_5_indices <- function(row) {
  order(row, decreasing = TRUE)[1:20]
}

out <- outgoing %>% as.data.frame()
outfb <- out[, colSums(out) !=0]
out5 <- apply(outfb, 1, get_top_5_indices)
out5 <- unique(as.vector(out5))
pheatmap(t(out5), scale = "row", col = heatcol, cluster_cols = F)

####################################################################################
####################################################################################

# Co-clustering of adult and paediatric samples with re-annotation of merged dataset from scratch

# Load MAPJAG tissue and AMP2 samples of DMARD-naive patient samples undergoing knee-biopsy who had a disease duration < 1 year
load(file="amp2-knee-1yr")
load(file="jia-knee-1yr")

# Subset out stromal clusters
Idents(tissue) <- "global1"
kfibs <- subset(tissue, idents=c("Fibroblasts", "Pericytes", "Endothelial", "Lymphatic"))

Idents(amp2) <- "cell_type"
afibs <- subset(amp2, idents=c("Stromal cell", "Endothelial cell"))

afibs <- afibs %>% NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) 
  
kfibs <- kfibs %>% NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = c(afibs, kfibs), dims = 1:50, reduction="rpca")
fibs <- IntegrateData(anchorset = anchors, dims = 1:50)

DefaultAssay(fibs) <- "RNA"
fibs <- fibs %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

#identify principal componant in which subsequent PCs account for <0.01% of variance
pct <- fibs[["pca"]]@stdev / sum(fibs[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
npc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

fibs <- RunHarmony(fibs, group.by.vars = "orig.ident")
fibs <- RunUMAP(fibs, reduction = "harmony", dims=1:npc)
fibs <- FindNeighbors(fibs, reduction= "harmony", dims = 1:npc, verbose = T)
fibs <- FindClusters(fibs, resolution = c(0.2,0.4, 0.5, 0.6, 0.65), graph.name = "integrated_snn")

FeaturePlot(fibs, features=c("VWF", "COL1A1", "ACTA2"), order=T) & sparseplot

# Examine resolutions
plot.list <- list()
res <- c(0.3,0.4, 0.5, 0.6)
for (i in 1:length(res)) {
  resg <- paste0("integrated_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(fibs, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=2) 

FeaturePlot(fibs, features=c("VWF", "COL1A1", "ACTA2"), order=T) & sparseplot

Idents(fibs) <- "integrated_snn_res.0.4"
fibs <- RenameIdents(fibs, "0"="CXCL12+ Fibroblasts", "1"="POSTN+ Fibroblasts", "2"="Venous", "3"="LL Fibroblasts", "4"="FLI1+ Capillary",
  "5"="CD34+ Fibroblasts", "6"="Pericytes", "7"="KDR+ Arterial", "8"="SOX5+ Fibroblasts", "9"="Lymphatics", "10"="Pericytes")
fibs$fbclusters <- fibs@active.ident

FeaturePlot(fibs, features=c("CXCL12", "CD34", "SOX5", "CDH11","POSTN", "PRG4"), ncol=3) & sparseplot

DimPlot(fibs, group.by = "fbclusters", split.by="cohort") & scale_color_manual(values=colurs)

####################################################################################

# Visualise original cluster names per dataset
combinedclust <- afibs[["cluster_name"]]
colnames(combinedclust) <- "combined_orig"
combined2 <- kfibs[["clusters2312"]]
colnames(combined2) <-  "combined_orig"
combined <- rbind(combined2, combinedclust)
fibs <- AddMetaData(fibs, combined, col.name='combined_orig')

DimPlot(fibs, group.by="combined_orig", split.by="cohort") & scale_color_manual(values=colurs)

####################################################################################

# Visualise split of clusters by origin of data
ptnum <- as.data.frame(table(fibs$fbclusters, fibs$cohort))
ptnum <- ptnum %>% group_by(Var2) %>% mutate(percent =  Freq*100/sum(Freq))
ptnum$Var1 <- factor(ptnum$Var1, levels=c("LL Fibroblasts", "CXCL12+ Fibroblasts",  "CD34+ Fibroblasts", "POSTN+ Fibroblasts","SOX5+ Fibroblasts", "KDR+ Arterial", "Venous", 
   "FLI1+ Capillary","Pericytes","Lymphatics"))

ggplot(ptnum) + aes(x = Var1, y = percent, fill=Var2) + 
  geom_bar(position="fill", stat="identity") +
  labs(x="\n\n", y = "\nProportion of Specimen Type by Cluster\n") + 
  scale_fill_manual(values=c("pink", "blue")) +
  ggtitle("\n") + theme_bw() & RotatedAxis() 
  


