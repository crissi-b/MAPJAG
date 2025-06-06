

library(Seurat)
library(ggplot2)
library(knitr)
library(stringr)
library(tidyverse)
library(dplyr)

# QC

#################################################################

# Set graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
sparseplot <- list(NoAxes(), NoLegend())
colurs <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Paired"))


# Finer annotation xenium subclusters
xendat <- readRDS(file="2404-xenium-global.rds")

#################################################################

# Chose a global clustering level to work from

#identify principal componant in which subsequent PCs account for <0.01% of variance
pct <- xendat[["pca"]]@stdev / sum(xendat[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
npc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

xendat <- FindNeighbors(xendat, reduction= "harmony", dims = 1:npc, verbose = T)
xendat <- FindClusters(xendat, resolution = c(0.1), graph.name = "RNA_snn")

DimPlot(xendat, group.by="RNA_snn_res.0.1", raster=FALSE, cols = colurs, label=T, repel=T) & NoLegend()

xendat <- FindClusters(xendat, resolution = c(0.1, 0.12, 0.14, 0.16, 0.18, 0.2), graph.name = "RNA_snn")
plot.list <- list()
res <- c(0.1, 0.12, 0.14, 0.16, 0.18, 0.2)
for (i in 1:length(res)) {
  resg <- paste0("RNA_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(xendat, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=3) 

xendat <- FindClusters(xendat, resolution = c(0.1,0.13,0.15), graph.name = "RNA_snn")
DimPlot(xendat, group.by="RNA_snn_res.0.13", raster=FALSE, cols = colurs, label=T, repel=T) & NoLegend()

Idents(xendat) <- “RNA_snn_res.0.13”
check <- FindAllMarkers(xendat, only.pos = T, logfc.threshold = 3)
xendat <- RenameIdents(xendat, "0"="Fibroblasts", "1"="Myeloid", "2"="Endothelial cells", 
                       "3"="Plasma cells", "4"="Pericytes", "5"="T/NK cells", "6"="Lymphatics", 
                       "7"="B cells/ pDCs", "8"="Cycling cells", "9"="Adipose", "10"="Granulocytes / progenitors")
xendat$global13 <- xendat@active.ident 

Idents(xendat) <- "global13"
levels(xendat) <- rev(c("Fibroblasts", "Pericytes", "Endothelial cells", "Lymphatics", "Adipose", "Myeloid", "T/NK cells", "B cells/ pDCs", "Plasma cells", "Granulocytes / progenitors", "Cycling cells"))
gen2=c("PRG4", "PDGFRA", "ACTA2", "VWF","ACKR1", "PROX1","LYVE1", "ADIPOQ", "PLIN4","MARCO", "FCN1", "CD247","CD8A","GNLY", "MS4A1","TCL1A","MZB1","PRDM1", "MS4A2", "SLC18A2", "TOP2A", "MKI67", "PLD4")
DotPlot(xendat, features=gen2) + labs(title="Key markers in global cell clusters\n") & formplot

#################################################################

# Annotate T cells

tcells <- subset(xendat, idents="T/NK cells")     

DefaultAssay(tcells) <- "RNA"
tcells <- tcells %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

#identify principal componant in which subsequent PCs account for <0.01% of variance
pct <- tcells[["pca"]]@stdev / sum(tcells[["pca"]]@stdev) * 100
npc1 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
npc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

tcells <- RunHarmony(tcells, group.by.vars = "orig.ident")
tcells <- RunUMAP(tcells, reduction = "harmony", dims=1:npc1)
tcells <- FindNeighbors(tcells, reduction= "harmony", dims = 1:npc1, verbose = T)
tcells <- FindClusters(tcells, resolution = c(0.1, 0.2), graph.name = "RNA_nn")

DimPlot(tcells, group.by="RNA_snn_res.0.1", raster=FALSE, cols = colurs, label=T, repel=T) & NoLegend()
FeaturePlot(tcells, features=c("KLRC1","CD4", "CD8A", "PRG4", "FOXP3", "MS4A1", "KLRB1", "MZB1","MKI67"))

tcells <- FindClusters(tcells, resolution = c(0.14,0.16, 0.2, 0.24, 0.26, 0.28), graph.name = "RNA_snn")

plot.list <- list()
res <- c(0.22, 0.24, 0.26)
for (i in 1:length(res)) {
  resg <- paste0("RNA_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(tcells, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=3) 


Idents(tcells) <- "RNA_snn_res.0.26"
markers <- FindAllMarkers(tcells, only.pos=T, logfc.threshold = 1.5)
top4uniq <- unique(markers %>% 
                     group_by(cluster) %>% 
                     top_n(n=6, wt = avg_log2FC) %>% 
                     dplyr::pull(gene))
DotPlot(tcells, features = top4uniq, cluster.idents = T) + formplot

Idents(tcells) <- "RNA_snn_res.0.26"
tcells <- RenameIdents(tcells,"0"="T cells/B/Myeloid", "1"= "T cell/Fibroblast", "2"="CD4+ KLRB1+ T", "3"="GZMK+ CD8+ T", 
                       "4"="Vascular/T cell", "5"="NK cells/ILCs", "6"="FOXP3+ Tregs", "7"="T cell/Plasma cell",
                       "8"="GZMK+ CD8+ T", "9"="Cycling T", "10"="T cells/LL",
                       "11"="Vascular/T cell")
tcells$anno2 <- tcells@active.ident

levels(tcells) <- rev(c("CD4+ KLRB1+ T", "FOXP3+ Tregs","GZMK+ CD8+ T", "NK cells/ILCs","Cycling T","T cells/LL",
                    "T cell/Fibroblast","T cells/B/Myeloid","T cell/Plasma cell","Vascular/T cell"))
keym <- unlist(str_split("KLRB1,IL7R,CCR2,CD4,CTLA4,FOXP3,IL2RA,CD8A,GZMK,CCL5,PRF1,GNLY,KLRD1,KLRC1,PRG4,PDGFRA,FBLN1,CD1C,MS4A1,FCER1A,MZB1,PRDM1,VWF,MYH11", pattern=","))
DotPlot(tcells, features = keym) + formplot + labs(title="Key markers\n\n")



#################################################################
#################################################################

# Annotate B cells

Idents(xendat) <- "global13"
bcells <- subset(xendat, idents=c("B cells/ pDCs", "Plasma cells"))     

DefaultAssay(bcells) <- "RNA"
bcells <- bcells %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

#identify principal componant in which subsequent PCs account for <0.01% of variance
pct <- bcells[["pca"]]@stdev / sum(bcells[["pca"]]@stdev) * 100
npc1 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
npc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

bcells <- RunHarmony(bcells, group.by.vars = "orig.ident")
bcells <- RunUMAP(bcells, reduction = "harmony", dims=1:npc)

bcells <- FindNeighbors(bcells, reduction= "harmony", dims = 1:npc, verbose = T)
bcells <- FindClusters(bcells, resolution = c(0.05,0.1,0.2), graph.name = "RNA_snn")
FeaturePlot(bcells, features=c("BANK1", "CXCR4", "PLD4", "MS4A1”, "SERPINB9", "TCF4", "SELL", "TCL1A"))

Idents(bcells) <- "RNA_snn_res.0.2"
markers <- FindAllMarkers(bcells, only.pos=T)

top4uniq <- unique(markers %>% 
                     group_by(cluster) %>% 
                     top_n(n=6, wt = avg_log2FC) %>% 
                     dplyr::pull(gene))
DotPlot(bcells, features = top4uniq, cluster.idents = T) + formplot

bcells <- RenameIdents(bcells, "0" ="Plasma cells", "1"="Plasma cells/Fibroblasts", 
                       "2"="Plasma cells/Myeloid", "3"="B cells", "4"="Vascular plasma cells",
                       "5"="Plasma cell/Granulocyte","6"="Cycling plasma cells", "7"="pDCs")
bcells$anno2 <- bcells@active.ident

levels(bcells) <- rev(c("pDCs", "B cells", "Plasma cells", "Cycling plasma cells", "Plasma cells/Fibroblasts", "Plasma cells/Myeloid", "Plasma cell/Granulocyte", "Vascular plasma cells"))
keym <- unlist(str_split("PLD4,CD4,GZMB,MS4A1,BANK1,CD83,CD69,MZB1,PRDM1,SLAMF7,TOP2A,MKI67,PDGFRA,THBS2,SEMA3C,VSIG4,MARCO,FCN1,SLC18A2,MS4A2,GATA2,VWF,MYH11", pattern=","))
DotPlot(bcells, features = keym) + formplot + labs(title="Key markers\n\n")


#################################################################
#################################################################

# Annotating granulocytes- not proceeded with as lack of discriminatory markers within small cluster

Idents(xendat) <- "global13"
gran <- subset(xendat, idents=c("Granulocytes / progenitors"))     

DefaultAssay(gran) <- "RNA"
gran <- gran %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

#identify principal componant in which subsequent PCs account for <0.01% of variance
pct <- gran[["pca"]]@stdev / sum(gran[["pca"]]@stdev) * 100
npc1 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
npc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

gran <- RunHarmony(gran, group.by.vars = "orig.ident")
gran <- RunUMAP(gran, reduction = "harmony", dims=1:npc1)
gran <- FindNeighbors(gran, reduction= "harmony", dims = 1:npc1, verbose = T)
gran <- FindClusters(gran, resolution = c(0.05, 0.1, 0.2, 0.5, 0.8), graph.name = "RNA_snn")

DimPlot(gran, group.by="RNA_snn_res.0.08", raster=FALSE, cols = colurs, label=T, repel=T) & NoLegend()

Idents(gran) <- "RNA_snn_res.0.08"
markers <- FindAllMarkers(gran, only.pos=T)

Idents(gran) <- "RNA_snn_res.0.05”
DotPlot(gran, features=c("FCER1A","MS4A2", "SLC18A2","GATA2", "CPA3", "KIT", "CD69", "IL1RL1", "SOX18", "ELF5", "HES4", "ADIPOQ", "PTPRC", "PDGFRA", "MZB1", "GZMK")) & formplot & labs(title="Key markers\n\n")

#################################################################
#################################################################

# Finer clustering of cycling cells

Idents(xen) <- "named_sub_new2"
cycling <- subset(xen, idents="Cycling cells")
DefaultAssay(cycling) <- "RNA"
cycling <- cycling %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30, verbose = FALSE)
cycling <- cycling %>% FindNeighbors(dims = 1:30) 
cycling <- cycling %>% FindClusters(resolution = c(0.1), graph.name = "RNA_snn")

DimPlot(cycling, group.by="RNA_snn_res.0.1", label=T, repel=T) & scale_color_manual(values=colorRampPalette(colurs)(48)) & NoLegend()
DotPlot(cycling, features=c("CD8A", "THY1", "CD68", "PDPN", "KIT", "VWF", "MS4A2", "MS4A1", "PRDM1", "SLC18A2"))
Idents(cycling) <- "RNA_snn_res.0.1"
cycling <- RenameIdents(cycling, "0"="Cycling fibroblasts", "1"="Cycling T cells", "2"="Cycling myeloid", "3"="Cycling fibroblasts", "4"="Granulocytes/progenitors")
cycling$named_sub_new2 <- cycling@active.ident

Idents(xen) <- "named_sub_new2"
other <- subset(xen, idents="Cycling cells", invert=T)
clusters <- rbind(other[["named_sub_new2"]], cycling[["named_sub_new2"]])
xen <- AddMetaData(xen, clusters, col.name="temp2")

DimPlot(xen, group.by="temp2", label=T, repel=T) & scale_color_manual(values=colorRampPalette(colurs)(51)) & NoLegend()

#################################################################
#################################################################

# Cleaning up the annotations

# Group decision to label mixed populations by their main cell type for clarity
check <- as.data.frame(table(xen$temp2, xen$global13))
check <- check[check$Freq >0,]
result <- check %>%
  group_by(Var1) %>%
  filter(Freq == max(Freq)) %>%
  ungroup()

# Unchanged 'Adipose','pDCs','Pericytes', 'Lymphatics', 'Plasma cells', 'NK cells/ILCs', 'CD1C+ cDC2s','Endothelial cells'

Idents(xen) <- "temp2"
xen <- RenameIdents(xen, 'macs'="Fibroblasts", 'Plasma cells/Myeloid'="Plasma cells", 
                    'S1008A+ Monocytes'="S1008A+ monocytes", "Mixed fibroblasts/myeloid"="Myeloid cells", 
                    "Mixed myeloid/B cells"="Myeloid cells",
                    'Tcell_contam_myeloid'="Myeloid cells", 'Peri_contam_myeloid'="Myeloid cells", 'CD1C_DCs' ="pDCs",
                    'MertK+ macrophages'="MerTK+ macrophages", 'proliferating_myeloid'="Cycling myeloid", 
                    'B_cell_doublets'="Fibroblasts", 'LL'= "LL fibroblasts", 'Granulocytes / progenitors'= "Granulocytes",
                    'fib_contam_myeloid'="Myeloid cells", 'LAMP3 DCs'="LAMP3+ DCs", 'B_cell_contam_myeloid'="Myeloid cells", 
                    'Vascular plasma cells'="Plasma cells", "SPP1+ Macrophages"="SPP1+ macrophages",
                    'POSTN'= "POSTN+ fibroblasts", 'Fibroblasts_macs'="Fibroblasts", 'CD4+ KLRB1+ T'="CD4+ KLRB1+ T cells", 
                    'GZMK+ CD8+ T' = "GZMK+ CD8+ T cells", 'Tcell/Fibroblast'="Fibroblasts", 
                    'T cell/Fibroblast'="T cells",  'Plasma cell/Granulocyte'="Plasma cells", 
                    'Vascular/T cell'="T cells", 'CD34/MFAP5'= "CD34+ fibroblasts", 
                    'SOX5/CDH11'= "SOX5+ CDH11+ fibroblasts", 'mac/LL'="LL fibroblasts", 'SFRP/CXCL12'= "CXCL12+ fibroblasts",
                    'Tcell_contam'="Fibroblasts",  
                    'Plasma cells/Fibroblasts'="Plasma cells", 'Granulocytes/progenitors'= "Granulocytes", 
                    'T cells/LL'="Lining layer T cells", 
                    'T cell/Plasma cell'="T cells", 'Lamp3_DCs'="LAMP3+ DCs", 
                    "Cycling T cells"="Cycling T cells", "Cycling myeloid"="Cycling myeloid")
xen$named2407 <- xen@active.ident      

#################################################################
#################################################################

# Niche analysis


#################################################################
#################################################################

# Adding patient meta-data for Krenn scoring
# First add patient study ID
xen$patient <- stringr::str_extract(xen$orig.ident,"[^_]*_[^_]*")
xen$patient <- gsub("D251", "", xen$patient)
xen$patient <- gsub("_", "", xen$patient)
unique(xen$patient)

# Total Krenn score (LL hyperplasia and infiltrate)
Idents(xen) <- "patient"
xen <- RenameIdents(xen, "811"="2", "824"="3", "828"="3", "836"= "4", "840"="2", "904"=NA, "915"="1.5", "918"="3")
xen$krenn <- xen@active.ident

# Inflammatory Krenn score (infiltrate only)
Idents(xen) <- "patient"
xen <- RenameIdents(xen, "811"="1", "824"="2", "828"="1.5", "836"= "2.5", "840"="3", "904"=NA, "915"="1.5", "918"="1.5")
xen$infiltrate <- xen@active.ident

