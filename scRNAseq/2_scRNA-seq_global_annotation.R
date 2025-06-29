library(Seurat)
library(ggplot2)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(SeuratObject)
options(bitmapType='cairo')

# Load global object
load(file="global-object.rdata")

### set graph graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
thin <-   list(theme(plot.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y = element_blank()))

#################################################################
#################################################################

# Create UMAP split by specimen type

# First name idents- simplest labels
Idents(PBMC1) <- "integrated_snn_res.0.15"
PBMC1 <- RenameIdents(PBMC1, '0'='Myeloid cells', '9'='Endothelial cells', '13'='Pericytes', '1' = 'T cells', '3' = 'T cells', '5' = 'T cells',
                      '4'='NK cells', '8'='T cells', '2'='Fibroblasts', '12' = 'Cycling cells', '7'='B cells', '17'='Lymphatics', '6'='Innate lymphoid cells',
                      '11'='Plasma cells', '14'='Dendritic cells', '18'= 'Myeloid cells', '10'='Dendritic cells', '16'='Dendritic cells', '15'="Plasmacytoid DCs")
PBMC1$simple15 <- PBMC1@active.ident

# Naming idents- CD4 and CD8 T cell split, dendritic cells split
Idents(PBMC1) <- "integrated_snn_res.0.15"
PBMC1 <- RenameIdents(PBMC1, '0'='Myeloid cells', '9'='Endothelial cells', '13'='Pericytes', '1' = 'CD4+ T cells', '3' = 'CD8+ T cells', '5' = 'CD8+ T cells',
                      '4'='NK cells', '8'='yd T cells', '2'='Fibroblasts', '12' = 'Cycling cells', '7'='B cells', '17'='Lymphatics', '6'='Innate lymphoid cells',
                      '11'='Plasma cells', '16'='Dendritic cells (LAMP3+)', '10'='Dendritic cells (cDC2s)', '14'='Dendritic cells (cDC1s)', '15'="Plasmacytoid DCs")
PBMC1$label15s <- PBMC1@active.ident

# Set cluster order
Idents(PBMC1) <- "label15s"
mainOrder <- c('CD4+ T cells', 'CD8+ T cells', 'Innate lymphoid cells', 'NK cells', 'yd T cells', 
               'Cycling cells', 'Plasma cells', 'B cells',  'Plasmacytoid DCs', 
               'Myeloid cells', 'Dendritic cells (cDC2s)',  'Dendritic cells (LAMP3+)', 
               'Dendritic cells (cDC1s)', 'Fibroblasts', 'Endothelial cells', 'Pericytes',
               'Lymphatics')
levels(PBMC1) <- rev(mainOrder)
PBMC1$label15s <- factor(PBMC1$label15s, levels=mainOrder)


ordercol <- c( "blue","navy", "dodgerblue2","turquoise2", "aquamarine", "#A6CEE3", "deeppink", "plum1","blueviolet", "aquamarine4","yellow3","chartreuse", "darkgreen", 
               "#A65628","#E6AB02","#E31A1C", "black")

DimPlot(PBMC1, group.by="label15s", split.by="TYPE", raster=F) & scale_colour_manual(values=ordercol) & thin


#################################################################
#################################################################


# Create bar chart of total cell numbers per specimen types

ptnum <- as.data.frame(table(PBMC1$label15s, PBMC1$TYPE))
ptnum$Var2 <- factor(ptnum$Var2, levels=rev(c( "Tissue", "Blood", "SF")))
ptnum$Var1 <- factor(ptnum$Var1, levels=rev(mainOrder))
ggplot(ptnum) + aes(x = Var2, y = Freq, fill=Var1) + geom_bar(position="stack", stat="identity") +
  labs(x="\n\n", y = "\nAbsolute cell numbers\n") + 
  ggtitle("\n") +
  theme_bw() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=rev(ordercol)) + RotatedAxis() + FontSize("8") +
  coord_flip() + guides(fill=guide_legend(reverse=F))


#################################################################
#################################################################


# PCA of haematopoietic cells across samples
# Needs trouble-shooting

# Extract haematopoietic cells
Idents(PBMC1) <- "label15s"
haem <- subset(PBMC1, idents=c("Endothelial cells", "Lymphatics", "Fibroblasts", "Pericytes"), invert=T)


dat <- table(haem@meta.data$label15s,haem@meta.data$sample)
dat <- dat[,which(colSums(dat) > 0)]
dat <- dat[which(rowSums(dat) > 0),]
dat <- t(dat)
dat <- dat/rowSums(dat)
dat <- scale(dat)
pca <- prcomp(dat)
# Extract the variance explained by each PC
explained_variance <- (pca$sdev)^2
exp_vr <- explained_variance / sum(explained_variance)
# Variance labels
vPC1 <- paste0("PC1 (", round(exp_vr[1]*100, digits=0), "% variance)")
vPC2 <- paste0("PC2 (", round(exp_vr[2]*100, digits=0), "% variance)")
pca_coords <- as.data.frame(pca$x)[,1:2]
pca_coords$ID <- row.names(pca_coords)
pca_coords$ID <- gsub("[0-9]+", "", pca_coords$ID)
pca_coords$ID <- factor(pca_coords$ID, levels=c("T", "B", "SF"))
ggplot(pca_coords, aes(x=PC1, y=PC2, label=ID, color=ID)) + geom_point(size=4) + 
  scale_color_manual(values=c("burlywood3", "red",  "steelblue3")) +
  theme_bw() + theme(panel.grid = element_blank()) + labs(x=vPC1, y =vPC2)


#################################################################
#################################################################


# Create dotplot of key markers for Total-seq protein assay and scRNA-seq 

Idents(PBMC1) <- "label15s"
mainOrder <- c('CD4+ T cells','Innate lymphoid cells', 'CD8+ T cells', 'yd T cells', 'NK cells',
               'Plasmacytoid DCs',  'Plasma cells', 'B cells','Dendritic cells (cDC2s)', 'Dendritic cells (cDC1s)',
               'Dendritic cells (LAMP3+)','Myeloid cells',  'Cycling cells', 'Lymphatics',
               'Endothelial cells','Fibroblasts', 'Pericytes')
PBMC1$label15s <- factor(PBMC1$label15s, levels=rev(mainOrder))

Idents(PBMC1) <- "label15s"
DefaultAssay(PBMC1) <- "RNA"
keyy <-unlist(str_split("CD3G,IL7R,PRKCH,NKG7,GNLY,PLD4,JCHAIN,CD79A,MS4A1,CLEC10A,CLEC9A,LAMP3,CST3,LYZ,S100A8,MKI67,LYVE1,VWF,COL3A1,ACTA2",pattern=','))
DotPlot(PBMC1, features=keyy) & list(RotatedAxis(), FontSize("10"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())


# Protein assay of blood and SF

Idents(PBMC1) <- "TYPE"
sfb <- subset(PBMC1, idents=c("SF", "Blood"))
DefaultAssay(sfb) <- "Protein"

Idents(sfb) <- "label15s"
sfb <- subset(sfb, idents=c("Endothelial cells", "Fibroblasts"), invert=T)
prot <- unlist(str_split('CD4.1,IL7R.1,CD8A.1,KLRG1.1,huTCRVg2,IL2RB.1,NCAM1.1,CLEC4C.1,IL3RA.1,CD38.1,TFRC.1,CD19.1,MS4A1.1,CD1C.1,ITGAX.1,ANPEP.1,CD86.1,FCGR1A.1,ITGAM.1', pattern=','))
levels(sfb) <- rev(mainOrder)
DotPlot(sfb, features=prot) & formplot

#################################################################
