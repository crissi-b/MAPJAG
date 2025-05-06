
library(Seurat)
library(dplyr)
library(sctransform)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(harmony)
library(gridExtra)

##################################################################################

# Matrices extracted from QuPath
# Read in segmented matrices, normalise and merge all samples together

# Set graph graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
thin <- list(theme(axis.title.x=element_blank(), axis.title.y=element_blank()))
bright <- c("bisque3", "#FF00CC","lightpink", "forestgreen",  "darkorange","aquamarine3" , "steelblue3", "lightcyan",  "lemonchiffon",
            "gold2", "yellow4", "yellow", "cyan", "royalblue4" ,"brown",  "mediumpurple", "green", "red", "blue","purple4",
            "honeydew3","deeppink4","#826D9B", "seashell","black",  "lightskyblue1","yellowgreen", "plum1", "cyan3", "palevioletred",
            "lavenderblush2" ,"darkgrey","mediumorchid3","salmon2", "brown" ,"mediumspringgreen")
formplotpink <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2(low = "pink", mid = "white",high = "blue", midpoint = 0))

                                                                                                                   
##################################################################################

setwd("/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2/2401_segmentation")

results <- dir("./", pattern = "*txt", 
               full.names = TRUE)

# Create pt vector from file names
pt <- gsub("[^_]*_[^_]*_", "", results)
pt <- gsub("_", "", pt)
pt <- gsub(".txt", "", pt)

results <- gsub(".//", "", results)

data = list()
for (i in 1:length(results)) {
  data[[i]] <- read.delim(results[[i]])
}

# Format spatial coordinates into a dataframe
data_x_y = list()
for (i in 1:length(results)) {
  data_x_y[[i]] <- read.delim(results[[i]])
  rownames(data_x_y[[i]])<-seq_along(data_x_y[[i]][,1])
  data_x_y[[i]] <- data_x_y[[i]] %>% select(c('Centroid.X.µm', 'Centroid.Y.µm'))
}

names(data_x_y) <- pt
names(data) <- pt

# Remove markers where staining failed
genesRemove <- c("CD31", "CD138")

# Create Seurat object from expression matrices
for (i in 1:length(data)) {
  data[[i]]$coord <- paste0(data[[i]]$Centroid.X.µm,data[[i]]$Centroid.Y.µm) # Create a meta-data tag of combined X and Y coordinates
  Cell.coord <- as.data.frame(data[[i]][["coord"]])
  data[[i]] <- data[[i]] %>% select(contains('.mean')) # select mean expression of marker per cell
  data[[i]] <- data[[i]] %>% select(contains('cell'))
  colnames(data[[i]]) <- gsub("Cell..", "", colnames(data[[i]])) # reformat column name for clarity
  colnames(data[[i]]) <- gsub(".mean", "", colnames(data[[i]])) # reformat column name for clarity
  rownames(data[[i]])<-seq_along(data[[i]][,1])
  data[[i]]<-as.data.frame(t(data[[i]]))
  data[[i]]<- filter(data[[i]], !(rownames(data[[i]]) %in% genesRemove)) # get rid of markers that failed
  data[[i]]<-CreateSeuratObject(data[[i]])
  data[[i]]<- subset(data[[i]], subset=nCount_RNA>0)
  data[[i]]<-SCTransform(data[[i]])
  data[[i]] <- AddMetaData(data[[i]], metadata=Cell.coord, col.name="coordinates") # add the coordinates label made earlier 
}

head(data[[i]])

#rename
for (i in 1:length(results)) {
  data[[i]]$orig.ident<- names(data[i]) 
  }

##################################################################################

# To improve data quality only include cells that had a good DAPI signal at the start and end of staining, ie the tissue wasn’t dislodged and the nuclear staining wasn’t lost 
# Example of where the cells were detected in sample 1
i <- 1
test <- data_x_y[[i]]
ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm)) + geom_point(size=0.2) + scale_x_reverse()

i <- 1
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 9 & `DAPI-INIT` > 8.8)
data[[i]]<-SCTransform(data[[i]])

i <- 2
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 7 & `DAPI-INIT` > 6 )
data[[i]]<-SCTransform(data[[i]])

i <- 3
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 7 & `DAPI-INIT` > 6)
data[[i]]<-SCTransform(data[[i]])

i <- 4
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 7 & `DAPI-INIT` > 7.5)
data[[i]]<-SCTransform(data[[i]])

i <- 5
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0)
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 6 & `DAPI-INIT` > 6)
data[[i]]<-SCTransform(data[[i]])

i <- 6
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 6 & `DAPI-INIT` > 6)
data[[i]]<-SCTransform(data[[i]])

i <- 7
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 6 & `DAPI-INIT` > 6)
data[[i]]<-SCTransform(data[[i]])

i <- 8
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 7.5 & `DAPI-INIT` > 7.5)
data[[i]]<-SCTransform(data[[i]])

i <- 9
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 6 & `DAPI-INIT` > 5.8)
data[[i]]<-SCTransform(data[[i]])

#visualise the DAPI at start and end per sample
library(cowplot)

qcplot <- list()
for (i in 1:length(pt)) {
  qcplot[[i]] <- VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
}
plot_grid(plotlist = qcplot)

# Filter coordinates for those meeting the DAPI filtering steps
for (i in 1:length(pt)) {
  data_x_y[[i]]$coord <- paste0(data_x_y[[i]]$Centroid.X.µm, data_x_y[[i]]$Centroid.Y.µm)
  data_x_y[[i]] <- filter(data_x_y[[i]], data_x_y[[i]]$coord %in% data[[i]]$coordinates)  
}

##################################################################################

# Merge data into one object
all <- merge(x = data[[1]], y = data[-1], merge.data=TRUE)
all <- FindVariableFeatures(all, assay="RNA")
all <- ScaleData(all, assay="RNA") 
check <- all@assays$RNA@var.features
check <- check[grep("DAPI", check, invert=T)]
VariableFeatures(all) <- check
all <- RunPCA(all, assay = 'SCT')

# Find principal component where subsequent PC has < 0.1% and < 0.05% difference in variance
pct <- all[["pca"]]@stdev / sum(all[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

# With harmony batch correction
library(harmony)
all <- RunHarmony(all, group.by.vars = "orig.ident", assay.use="SCT", reduction="pca")
all <- RunUMAP(all, reduction="harmony", dims=1:21)
all<-FindNeighbors(all, dims = 1:21, reduction= "harmony", assay="SCT")
all<-FindClusters(all, resolution = c(0.2, 0.24, 0.28), graph.name="SCT_snn")

# View distribution of clusters across samples
ptnum2 <- as.data.frame(table(all$SCT_snn_res.0.28, all$orig.ident))
ggplot(ptnum2) + aes(x = Var1, y = Freq, fill=Var2) + geom_bar(position="fill", stat="identity") + labs(x="", y = "") + 
  theme(legend.title=element_blank()) + scale_fill_manual(values=bright) + RotatedAxis()

# check distribution of markers across clusters
gen2 <- c("CD15", "CD3", "MCT", "CD68", "CLU", "COL1", "COL3A1", "CD20",
          "COL6A1", "COMP", "CTHRC1", "DKK3", "FABP4", "LYVE1", "MMP3", "PDPN", "POSTN","SPARC",  "CD146", "COL4A1",  "SMA")
DefaultAssay(all) <- "SCT"
Idents(all) <- "SCT_snn_res.0.28"
DotPlot(all, features =gen2, cluster.idents = T) & formplot

##################################################################################

# Review how it looks with clustering

# Amend i to any number between 1 and 8 to review the samples
i <- 3
nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$SCT_snn_res.0.28

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

# Zoom in and review

small <- test[test$Centroid.X.µm < 2000  &  test$Centroid.Y.µm < 2000,]

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7) + theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

##################################################################################

# Now go through each cluster and recluster where necessary to find the best fit of the cell types, cross-check with the original image, and review the samples with the new clustering- after trying different resolutions this was the best representation of the fluorescence images we could annotate

# COL1 high cluster
Idents(all) <- "SCT_snn_res.0.28"
CC <- subset(all, idents=0)
CC <-FindNeighbors(CC, dims = 1:21, reduction= "harmony", assay="SCT")
CC <-FindClusters(CC, resolution = c(0.1), graph.name="SCT_snn")
table(CC$orig.ident, CC@active.ident)
Idents(CC) <- "SCT_snn_res.0.1"
CC <- RenameIdents(CC, "0"="SL Fibroblasts", "1"="COL1-hi Fibroblasts", "2"="COL1-hi Fibroblasts", "3"="COL1-hi Fibroblasts")
CC$take2 <- CC@active.ident

# Annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=0, invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], CC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Different sized vessels
Idents(all) <- "SCT_snn_res.0.28"
AAa <- subset(all, idents=1)

# Deal with CD15 high populations later as abnormal appearance to stain
AAf <- subset(AAa, CD15 > 4.8)
AAf$take2 <- "Fibrous"
AA <- subset(AAa, CD15 <= 4.8)

AA <-FindNeighbors(AA, dims = 1:21, reduction= "harmony", assay="SCT")
AA <-FindClusters(AA, resolution = c(0.12), graph.name="SCT_snn")
table(AA$orig.ident, AA@active.ident)
Idents(AA) <- "SCT_snn_res.0.12"
AA <- RenameIdents(AA,  "3"="Artefact", "2"="SMA-hi vessels")
AA$take2 <- AA@active.ident

# Split DKK3 fbs from vessels
Idents(AA) <- "SCT_snn_res.0.12"
AA5 <- subset(AA, idents=c("0"))
AA5a <- subset(AA5, PDPN > 7.5 & COL4A1 <= 8.3)
AA5a$take2 <- "DKK3+ Fibroblasts"
AA5b <- subset(AA5, PDPN > 7.5 & COL4A1 <= 8.3, invert=T)
AA5c <- subset(AA5b, PDPN > 6.8, invert=T)
AA5c$take2 <- "Endothelial cells-1"
AA5d <- subset(AA5b, PDPN > 6.8)
AA5d$take2 <- "Endothelial cells-2"

anno <- rbind(AA5a[["take2"]],AA5c[["take2"]], AA5d[["take2"]])
AA5 <- AddMetaData(AA5, anno, col.name="take2")

# Large vessel
Idents(AA) <- "SCT_snn_res.0.12"
AA1 <- subset(AA, idents=c("1"))
AA1 <-FindNeighbors(AA1, dims = 1:21, reduction= "harmony", assay="SCT")
AA1 <-FindClusters(AA1, resolution = c(0.2), graph.name="SCT_snn")
table(AA1$orig.ident, AA1@active.ident)
Idents(AA1) <- "SCT_snn_res.0.2"
AA1 <- RenameIdents(AA1, "3"="Arteriole")
AA1$take2 <- AA1@active.ident

Idents(AA1) <- "SCT_snn_res.0.2"
AAy <- subset(AA1, idents=c("3"))
AAz <- subset(AA1, idents=c("0", "1", "2"))
AAz1 <- subset(AAz, PDPN > 6.9)
AAz1$take2 <- "CD146-hi vessels"
AAz2 <- subset(AAz, PDPN > 6.9, invert=T)
AAz2$take2 <- "Vessels"

Idents(AA) <- "SCT_snn_res.0.12"
AA3 <- subset(AA, idents=c("3", "2"))
anno <- rbind(AAy[["take2"]], AAz2[["take2"]],AAz1[["take2"]], AA3[["take2"]], AAf[["take2"]], AA5[["take2"]])
AAa <- AddMetaData(AAa, anno, col.name="take2")

unique(AAa$take2)

# Annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("1"), invert=T)
#nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], AAa[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# CLU- high cluster
Idents(all) <- "SCT_snn_res.0.28"
FF <- subset(all, idents=2)
VlnPlot(FF, features=c("CD68", "SPARC", "COL4A1", "MCT", "PDPN", "CLU"), pt.size=0, ncol=3) & thin

FF <-FindNeighbors(FF, dims = 1:21, reduction= "harmony", assay="SCT")
FF <-FindClusters(FF, resolution = c(0.25), graph.name="SCT_snn")
table(FF$orig.ident, FF@active.ident)
Idents(FF) <- "SCT_snn_res.0.25"
DotPlot(FF, features=gen2) & formplot
FF <- RenameIdents(FF, "0"="Fibrous Macrophages", "1"="Fibrous", "2"="Low-staining Fibroblasts", "3"="Fibrous", "4"="Fibrous", "5"="Fibrous", "6"="CLU+ LL", "7"="Fibrous", "8"="Fibrous", "9"="Fibrous", "10"="Fibrous")
FF$take2 <- FF@active.ident

# Annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("2"), invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], FF[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# COMP high cluster
Idents(all) <- "SCT_snn_res.0.28"
TT1 <- subset(all, idents=3)

TTo <- subset(TT1, CD15 > 4.6 & `DAPI-INIT` > 9.2)
TTo$take2 <- "Neutrophils"
TT <- subset(TT1, CD15 > 4.6 & `DAPI-INIT` > 9.2, invert=T)

TT <-FindNeighbors(TT, dims = 1:21, reduction= "harmony", assay="SCT")
TT <-FindClusters(TT, resolution = c(0.05), graph.name="SCT_snn")
table(TT$orig.ident, TT@active.ident)
DotPlot(TT, features=gen2, cluster.idents=T) & formplot
Idents(TT) <- "SCT_snn_res.0.05"
TT <- RenameIdents(TT, "0"="COMP-hi Fibroblasts","1"="COMP-hi Fibroblasts",  "2"="Artefact")
TT$take2 <- TT@active.ident

anno <- rbind(TTo[["take2"]], TT[["take2"]])
TT1 <- AddMetaData(TT1, anno, col.name="take2")

Idents(TT1) <- "take2"
VlnPlot(TT1, features=c("CD15", "SPARC", "CLU", "CD68"), pt.size=0)

Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("3"), invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], TT1[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# T cells
Idents(all) <- "SCT_snn_res.0.28"
DKK1 <- subset(all, idents=4)

# Remove Fibrous Macrophages
DKK3 <- subset(DKK1, CLU > 8.8 & FABP4 > 4.5)
DKK3$take2 <- "Fibrous Macrophages"

DKK <- subset(DKK1, CLU > 8.8 & FABP4 > 4.5, invert=T)
DKK <-FindNeighbors(DKK, dims = 1:21, reduction= "harmony", assay="SCT")
DKK <-FindClusters(DKK, resolution = c(0.2), graph.name="SCT_snn")
table(DKK$orig.ident, DKK@active.ident)
Idents(DKK) <- "SCT_snn_res.0.2"
DKK <- RenameIdents(DKK, "0"="LL Macrophages", "1"="CD68+ Myeloid", "2"="T cells", "3"="T cells")
DKK$take2 <- DKK@active.ident

anno <- rbind(DKK[["take2"]], DKK3[["take2"]])
DKK1 <- AddMetaData(DKK1, anno, col.name="take2")

DKK2 <- subset(DKK, idents="T cells")
DKK2 <-FindNeighbors(DKK2, dims = 1:21, reduction= "harmony", assay="SCT")
DKK2 <-FindClusters(DKK2, resolution = c(0.15), graph.name="SCT_snn")
table(DKK2$orig.ident, DKK2@active.ident)
DKK2 <- RenameIdents(DKK2, "0"="CD3-low T cells", "1"="T cells", "2"="CD3-hi T cells")
DKK2$take2 <- DKK2@active.ident

Idents(DKK) <- "take2"
DKK4 <- subset(DKK, idents="T cells", invert=T)

anno <- rbind(DKK2[["take2"]], DKK3[["take2"]], DKK4[["take2"]])
DKK1 <- AddMetaData(DKK1, anno, col.name="take2")
table(DKK1$orig.ident, DKK1$take2)

# Annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("4"), invert=T)
#nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], DKK1[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Split out lymphatics from LL fibroblasts
Idents(all) <- "SCT_snn_res.0.28"
LLn <- subset(all, idents=5)
VlnPlot(LLn, features=c("COL6A1", "PDPN", "SPARC", "LYVE1", "CLU", "CD146"), pt.size=0)

LLn1 <- subset(LLn, LYVE1 > 6.2)
LLn1$take2 <- "Lymphatics"
LLn2 <- subset(LLn, LYVE1 <= 6.2 & CLU <= 9.5)
LLn2$take2 <- "LL Fibroblasts"
LLn3 <- subset(LLn, LYVE1 <= 6.2 & CLU > 9.5)
LLn3$take2 <- "Fibrous Fibroblasts"

anno <- rbind(LLn1[["take2"]], LLn2[["take2"]], LLn3[["take2"]])
LLn <- AddMetaData(LLn, anno, col.name="take2")
table(LLn$orig.ident, LLn$take2)

# Correct staining anomalies of lympahtics in sample 1003
Idents(LLn1) <- "orig.ident"
Lna <- subset(LLn1, ident= c("1003"), invert=T)
Ln3 <- subset(LLn1, ident= "1003")
VlnPlot(Ln3, features=c("DKK3", "PDPN", "SPARC", "LYVE1", "COL4A1", "MMP3"), pt.size=0)
Ln3a <- subset(Ln3, DKK3 <= 7.8)
Ln3a$take2 <- "LL Fibroblasts"
Ln3b <- subset(Ln3, DKK3 > 7.8)
 
# Correct staining anomalies of LL in sample 2002
Idents(LLn2) <- "orig.ident" 
Lnb <- subset(LLn2, ident= c("2002"), invert=T)
Ln2 <- subset(LLn2, ident= "2002")
VlnPlot(Lb4, features=c("COL6A1", "PDPN", "SPARC", "LYVE1", "CLU", "MMP3"), pt.size=0)
Lb3 <- subset(Ln2, SPARC <= 8.4 & PDPN > 7.4)
Lb3$take2 <- "Lymphatics"
Lb4 <- subset(Ln2,  SPARC <= 8.4 & PDPN > 7.4, invert=T)
Lb4$take2 <- "LL Fibroblasts"

anno <- rbind(Lna[["take2"]], Ln3a[["take2"]], Ln3b[["take2"]], Lnb[["take2"]], Lb3[["take2"]], Lb4[["take2"]],  LLn3[["take2"]])
LLn <- AddMetaData(LLn, anno, col.name="take2")

Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("5"), invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], LLn[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Explored and decided against this clustering

#Idents(all) <- "SCT_snn_res.0.28"
#CL <- subset(all, idents=6)
#VlnPlot(CL, features=c("PDPN", "MCT", "COL6A1", "CLU", "COL1", "SPARC"), pt.size=0)

#CL <-FindNeighbors(CL, dims = 1:21, reduction= "harmony", assay="SCT")
#CL <-FindClusters(CL, resolution = c(0.08), graph.name="SCT_snn")
#table(CL$orig.ident, CL@active.ident)
#DotPlot(CL, features=gen2, cluster.idents=T) & formplot
#Idents(CL) <- "SCT_snn_res.0.08"
#CL$take2 <- CL@active.ident
##CL <- RenameIdents(CL, "0"="Low staining Fibroblasts","1"="SL Fibroblasts")

#table(CL$take2, CL$orig.ident)
#Idents(all) <- "SCT_snn_res.0.28"
#nLL <- subset(all, idents=c("6"), invert=T)
#nLL$take2 <- nLL$SCT_snn_res.0.28
#anno <- rbind(nLL[["take2"]], CL[["take2"]])
#all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Splitting adipose from COL4 hi cluster
Idents(all) <- "SCT_snn_res.0.28"
AP <- subset(all, idents=7)

VlnPlot(AP, features=c("FABP4", "COL4A1", "SPARC", "CD15", "APOD"), pt.size=0)

AP <-FindNeighbors(AP, dims = 1:21, reduction= "harmony", assay="SCT")
AP <-FindClusters(AP, resolution = c(0.05), graph.name="SCT_snn")
table(AP$orig.ident, AP@active.ident)
DotPlot(AP, features=gen2, cluster.idents=T) & formplot
Idents(AP) <- "SCT_snn_res.0.05"
AP <- RenameIdents(AP, "0"="COL3-hi Fibroblasts","1"="Adipose", "2"="Adipose")
AP$take2 <- AP@active.ident

table(AP$take2, AP$orig.ident)
Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("7"), invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], AP[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Subcluster POSTN
Idents(all) <- "SCT_snn_res.0.28"
PF <- subset(all, idents=8)

PF <-FindNeighbors(PF, dims = 1:21, reduction= "harmony", assay="SCT")
PF <-FindClusters(PF, resolution = c(0.1), graph.name="SCT_snn")
table(PF$orig.ident, PF@active.ident)
DotPlot(PF, features=gen2, cluster.idents=T) & formplot
Idents(PF) <- "SCT_snn_res.0.08"
PF <- RenameIdents(PF, "0"="POSTN-hi Fibroblasts","1"="COL1-hi Fibroblasts", "2"="Fibrous")
PF$take2 <- PF@active.ident

Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("8"), invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], PF[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Split mast cells from fibrous material
Idents(all) <- "SCT_snn_res.0.28"
MC <- subset(all, idents=14)

VlnPlot(MC, features=c("COL6A1", "COL3A1", "COL1", "MCT", "CD15", "CLU"), pt.size=0)

MCn <- subset(MC, MCT > 7)
MCn$take2 <- "Mast cells"
MCm <- subset(MC, MCT <= 7)
MCm$take2 <- "Fibrous"

anno <- rbind(MCn[["take2"]], MCm[["take2"]])
MC <- AddMetaData(MC, anno, col.name="take2")
table(MC$orig.ident, MC$take2)

Idents(all) <- "SCT_snn_res.0.28"
nLL <- subset(all, idents=c("14"), invert=T)
nLL$take2 <- nLL$SCT_snn_res.0.28
anno <- rbind(nLL[["take2"]], MC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Annotate the rest of the clusters
Idents(all) <- "take2"
all <- RenameIdents(all, "6"="Low-staining cells", "9"="COL6-hi Fibroblasts", "10"="Blood/Fibrous",
                   "11"="Mast cells", "12"="LYVE1-hi Macrophages", "13"="B/T cell aggregates",
                   "15"="LL Macrophages", 
                   "16"="Fibrous", "17"="Artefact", "18"="Artefact", "19"="Artefact") 
all$take2 <- all@active.ident

##################################################################################

# Refine clustering further
Idents(all) <- "take2"
DKt <- subset(all, ident = "T cells")

DKt1 <- subset(DKt, CD3 > 7.8)
DKt1$take2 <- "CD3-low T cells"
DKt2 <- subset(DKt, CD3 <= 7.8)
DKt2$take2 <- "Low-staining Fibroblasts"

anno <- rbind(DKt1[["take2"]], DKt2[["take2"]])
DKt <- AddMetaData(DKt, anno, col.name="take2")

Idents(all) <- "take2"
nLL <- subset(all, ident = "T cells", invert=T)
anno <- rbind(nLL[["take2"]], DKt[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Refine clustering of fibrin / fibrous material

Idents(all) <- "take2"
fib <- subset(all, idents=c("Fibrous", "Fibrous Fibroblasts"))

fibm <- subset(fib, CD15 > 4.5)
fibm$take2 <- "Neutrophil-rich Fibrous"
fibo <- subset(fib, CD15 <= 4.5)
fibo$take2 <- "Fibrous"

anno <- rbind(fibm[["take2"]], fibo[["take2"]])
fib <- AddMetaData(fib, anno, col.name="take2")

Idents(fib) <- "take2"
VlnPlot(fib, features=c("CD15", "SPARC"), pt.size=0)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("Fibrous"), invert=T)
anno <- rbind(nLL[["take2"]], fib[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Refine fibrous clustering

Idents(all) <- "take2"
clu <- subset(all, idents=c("CLU+ LL"))

Idents(clu) <- "orig.ident"
#VlnPlot(clu, features=c("PDPN", "COL3A1", "CD15", "MCT", "COL1", "COL6A1"), pt.size=0)

clua <- subset(clu, PDPN > 6 & COL3A1 > 6.8 & CD15 <= 4.6)
clua$take2 <- "CLU+ LL"
club <- subset(clu, PDPN > 6 & COL3A1 > 6.8 & CD15 <= 4.6, invert=T)
club$take2 <- "Fibrin"
anno <- rbind(clua[["take2"]], club[["take2"]])
clu <- AddMetaData(clu, anno, col.name="take2")

table(clu$orig.ident, clu$take2)

# Annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("CLU+ LL"), invert=T)
anno <- rbind(nLL[["take2"]], clu[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Refine neutrophils in messy tissue
Idents(all) <- "take2"
fib3 <- subset(all, idents=c("Neutrophils", "Low-staining cells"))

fibA <- subset(fib3, CD3 > 8.6)
fibA$take2 <- "Neutrophil/T cell-rich"
fibB <- subset(fib3, CD3 <= 8.6 & CD15 > 4.5)
fibB$take2 <- "Neutrophil-rich"
fibC <- subset(fib3, CD3 <= 8.6 & CD15 <= 4.5)
fibC$take2 <- "Low-staining cells"

anno <- rbind(fibA[["take2"]], fibB[["take2"]], fibC[["take2"]])
fib3 <- AddMetaData(fib3, anno, col.name="take2")

Idents(all) <- "take2"
nLL <- subset(all, idents=c("Neutrophils", "Low-staining cells"), invert=T)
anno <- rbind(fib3[["take2"]], nLL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Refine macrophage populations
Idents(all) <- "take2"
LLM <- subset(all, idents=c("LL Macrophages"))

LLMa <- subset(LLM, SMA > 7)
LLMa$take2 <- "Fibrous Macrophages"
LLMb <- subset(LLM, SMA <= 7)

anno <- rbind(LLMa[["take2"]], LLMb[["take2"]])
LLM <- AddMetaData(LLM, anno, col.name="take2")

table(LLM$take2, LLM$orig.ident)

# Annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("LL Macrophages"), invert=T)
anno <- rbind(nLL[["take2"]], LLM[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# B / T cells
Idents(all) <- "take2"
BB <- subset(all, ident = "B/T cell aggregates")
VlnPlot(BB, features=c("CD20", "CD3"), pt.size=0)

BB1 <- subset(BB, CD20 <= 7)
BB1$take2 <- "T cells-1"
BB2 <- subset(BB, CD20 > 7)
BB2$take2 <- "B/T cell aggregates"

anno <- rbind(BB1[["take2"]], BB2[["take2"]])
BB <- AddMetaData(BB, anno, col.name="take2")

Idents(all) <- "take2"
nLL <- subset(all, ident = "B/T cell aggregates", invert=T)
anno <- rbind(nLL[["take2"]], BB[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Split LYVE1+ macrophages from lymphatics
Idents(all) <- "take2"
CF1 <- subset(all, idents="LYVE1-hi Macrophages")

CF11 <- subset(CF1, PDPN > 8.2)
CF11$take2 <- "Lymphatics"
CF12 <- subset(CF1, PDPN <= 8.2)
CF12$take2 <- "LYVE1-hi Macrophages"

anno <- rbind(CF11[["take2"]], CF12[["take2"]])
CF1 <- AddMetaData(CF1, anno, col.name="take2")

table(CF1$orig.ident, CF1$take2)

Idents(all) <- "take2"
nLL <- subset(all, idents="LYVE1-hi Macrophages", invert=T)
anno <- rbind(nLL[["take2"]], CF1[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Clarify naming further

Idents(all) <- "take2"
all <- RenameIdents(all, "CD3-hi T cells"="T cells", "CD3-low T cells"="T cells", "T cells-1"="T cells", "Fibrous"="Fibrin", 
                    "Blood/Fibrous"="Blood/Fibrin", "Fibrous Macrophages"="Fibrin Macrophages", 
                    "Neutrophil-rich Fibrous"="Neutrophil-rich Fibrin")
all$take2 <- all@active.ident

##################################################################################

Idents(all) <- "take2"
neut <- subset(all, idents=c("Neutrophil-rich Fibrin"))

neut1 <- subset(neut, CD15 > 5 & CD68 <= 7.5)
neut1$take2 <- "Neutrophils"
neut2 <- subset(neut, CD15 > 5 & CD68 <= 7.5, invert=T)
neut2$take2 <- "Fibrin-2"

Idents(all) <- "take2"
nLL <- subset(all, idents="Neutrophil-rich Fibrin", invert=T)
anno <- rbind(nLL[["take2"]], neut1[["take2"]], neut2[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

##################################################################################

# Check markers 
DefaultAssay(all) <- "SCT"
Idents(all) <- "take2"
iden <- unique(all$take2)
DotPlot(all, features =gen2) & formplot

##################################################################################

# Inspect each sample for artefacts or staining issues, select coordinates of cells that need amending
# Amendments corrected later in code

i <- 1
nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

# Artefactual coordinates will be dealt with later in code
small <- test[test$Centroid.Y.µm <1300 & test$named1 %in% "T cells",]
exclude1 <- small$coord 

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$named1 %in% "Mast cells",]
exclude1a <- small$coord

####################################################################################

i <- 2

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm > 3600  & test$Centroid.Y.µm >2800 ,]
small <- test[test$Centroid.X.µm > 4600  & test$Centroid.Y.µm <4100 ,]
small <- test[test$Centroid.X.µm > 4600  & test$Centroid.X.µm < 4750  & test$Centroid.Y.µm < 4100 & test$Centroid.Y.µm > 3700 & test$named1 %in% c("SMA-hi vessels", "Blood/Fibrin", "Endothelial cells-1"),]
exclude2e <- small$coord

small <- test[test$Centroid.X.µm > 2000  & test$Centroid.X.µm < 3200  &  test$Centroid.Y.µm < 3801 &  test$Centroid.Y.µm > 2500,]

exclude2a <- small[small$named1 %in% c("B/T cell aggregates", "T cells", "COMP-hi Fibroblasts", "POSTN-hi Fibroblasts", "LL Fibroblasts"),]
exclude2a <- exclude2a$coord

small <- test[test$Centroid.X.µm < 3850  & test$Centroid.X.µm > 2700  & test$Centroid.Y.µm > 3500 &  test$Centroid.Y.µm < 4800,]
exclude2b <- small[small$named1 %in% c("B/T cell aggregates", "T cells", "COMP-hi Fibroblasts", "POSTN-hi Fibroblasts","LL Fibroblasts"),]
exclude2b <- exclude2b$coord

#small <- test[test$Centroid.X.µm < 2950  & test$Centroid.X.µm > 2600  & test$Centroid.Y.µm > 3800 &  test$Centroid.Y.µm < 4100,]
#exclude2d <- small[small$named1 %in% c("Neutrophils"),]
#exclude2d <- exclude2d$coord

small <- test[test$Centroid.X.µm > 2000  &  test$Centroid.Y.µm < 4800 &  test$Centroid.Y.µm > 1800,]
small <- test[test$Centroid.X.µm > 2000  &  test$Centroid.Y.µm < 4800 &  test$Centroid.Y.µm > 1800 & test$named1 %in% c("Fibrin-2", "Neutrophils", "Neutrophil-rich", "Neutrophil/T cell-rich", "Fibrin"),]

small <- test[test$Centroid.X.µm > 3000  &  test$Centroid.X.µm < 4000  & test$Centroid.Y.µm > 3200 &  test$Centroid.Y.µm < 4800 & test$named1 %in% c("Neutrophils"),]
exclude2f <- small$coord

small <- test[test$Centroid.X.µm > 2400  &  test$Centroid.Y.µm < 4800 &  test$Centroid.Y.µm > 2000 & test$named1 %in% c("LL Fibroblasts","B/T cell aggregates", "T cells", "SL T cells"),]
fibrous2 <- small$coord

small <- test[test$Centroid.X.µm > 3900  &  test$Centroid.Y.µm < 1300 & test$named1 %in% c("Fibrin", "Neutrophils"),]
exclude2g <- small$coord

####################################################################################

i <- 3

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 1600  &  test$Centroid.Y.µm > 1350 & test$named1 %in% c("POSTN-hi Fibroblasts", "COL1-hi Fibroblasts"),]
exclude3a <- small$coord                

small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.X.µm < 3100  & 
                test$Centroid.Y.µm < 600 &  test$Centroid.Y.µm > 400 & test$named1 %in% "POSTN-hi Fibroblasts",]
exclude3b <- small$coord  

small <- test[test$Centroid.X.µm > 1200  &  test$Centroid.X.µm < 1600 &  test$Centroid.Y.µm > 1000 & 
               test$named1 %in% "Lymphatics",]
LLc2 <- small$coord  
              
small <- test[test$named1 %in% c("Fibrin Macrophages"),]
mac3 <- small$coord

####################################################################################

i <- 4

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm > 600  &  test$Centroid.X.µm < 900  & 
                test$Centroid.Y.µm < 3100 &  test$Centroid.Y.µm > 2700  & test$named1 %in% c("POSTN-hi Fibroblasts","COL1-hi Fibroblasts"),]
exclude4 <- small$coord  

small <- test[test$named1 %in% "Adipose",]
col4 <- small$coord

small <- test[test$named1 %in% c("Fibrin Macrophages"),]
mac4 <- small$coord

small <- test[test$named1 %in% c("Fibrin"),]
collagen4a <- small$coord

####################################################################################

i <- 5

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.Y.µm > 3000  &  test$Centroid.Y.µm < 5500 & test$Centroid.X.µm > 3300 & test$named1 =="Artefact",]
vessel5 <- small$coord

####################################################################################

i <- 6

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & coord_flip() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.Y.µm > 1900  &  test$Centroid.Y.µm < 3800 & test$Centroid.X.µm > 6000 ,]
small <- test[test$Centroid.Y.µm > 1900  &  test$Centroid.Y.µm < 3800 & test$Centroid.X.µm > 6000 & 
                test$named1 %in% c("Neutrophils"),]
fibrous6a <- small$coord

small <- test[test$Centroid.Y.µm > 1900  &  test$Centroid.Y.µm < 3800 & test$Centroid.X.µm > 6000 & 
                test$named1 %in% c("Neutrophils"),]
fibrous6a2 <- small$coord

small <- test[test$Centroid.Y.µm > 6800 & test$named1 %in% c("Neutrophils"),]
fibrous6b <- small$coord

small <- test[test$Centroid.X.µm > 1800  &  test$Centroid.X.µm < 3200 & test$Centroid.Y.µm > 3800 & test$Centroid.Y.µm < 5300 & test$named1 %in% "Lymphatics",]
LLc <- small$coord

small <- test[test$named1 %in% c("Fibrin Macrophages"),]
mac6a <- small$coord

####################################################################################

i <- 7 

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500 & 
                test$named1 %in% c("LYVE1-hi Macrophages", "LL Fibroblasts", "Lymphatics", "SMA-hi vessels"),]

small <- test[test$Centroid.Y.µm > 4000  &  test$Centroid.Y.µm < 5000  & 
                test$Centroid.X.µm < 3000,]
small <- test[test$Centroid.Y.µm > 4000  &  test$Centroid.Y.µm < 5000  & 
                test$Centroid.X.µm < 3000 & test$named1 %in% "Artefact",]
vessel7 <- small$coord

small <- test[test$Centroid.Y.µm > 6500  & 
                test$Centroid.X.µm < 7000 & test$Centroid.X.µm > 5000,]
small <- test[test$Centroid.Y.µm > 7500  & test$Centroid.Y.µm < 7800  &
                test$Centroid.X.µm < 5550 & test$Centroid.X.µm > 5300 & test$named1 %in% "Artefact",]
vessel7a <- small$coord

small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 4000 & test$named1 %in% c("Fibrin Macrophages"),]
LLmac7a <- small$coord

small <- test[test$Centroid.X.µm < 1000  &  test$Centroid.Y.µm > 7000 & test$named1 %in% c("Fibrin Macrophages"),]
LLmac7b <- small$coord

####################################################################################

small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500,]
small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500 & test$named1 %in% "LYVE1-hi Macrophages",]

#VlnPlot(check, features=c("CD15", "LYVE1", "PDPN", "COMP", "COL6A1", "CD3"), pt.size=0)
check <- subset(all, coordinates %in% small$coord)

ly1 <- subset(check,  LYVE1 <= 8.5 & COL6A1 > 8.3, invert=T)
ly1$take2 <- "Lymphatics"
ly2b <- subset(check, LYVE1 <= 8.5 & COL6A1 > 8.3)
ly2b$take2 <- "LYVE1-hi Macrophages"

anno <- rbind(ly2b[["take2"]], ly1[["take2"]])
check <- AddMetaData(check, anno, col.name="take2")

table(check$take2, check$orig.ident)

# Annotate with new splits from clustering
nLL <- subset(all, coordinates %in% small$coord, invert=T)
anno <- rbind(nLL[["take2"]], check[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

i <- 8

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.7)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.Y.µm > 1900  &  test$Centroid.Y.µm < 2500 & test$Centroid.X.µm > 1800 & test$named1 %in% c("LL Fibroblasts", "T cells"),]
exclude8 <- small$coord

####################################################################################

# Renaming annotation based on fine examination of samples

exclude.cells <- c(exclude1, exclude1a, exclude2a, exclude2b, exclude2e, exclude2f, exclude2g, exclude3a, exclude3b, exclude4, exclude8)
exclude <- subset(all, coordinates %in% exclude.cells)
exclude$take2 <- "Artefact"

fibrous <- subset(all, coordinates %in% c(fibrous2, fibrous6a, fibrous6a2, fibrous6b))
fibrous$take2 <- "Fibrin"

vessel <- subset(all, coordinates %in% c(vessel5, vessel7, vessel7a))
vessel$take2 <- "SMA-hi vessels"

col4 <- subset(all, coordinates %in% col4)
col4$take2 <- "COL3-hi Fibroblasts"

mac <- subset(all, coordinates %in% c(mac3, mac4, mac6a)) 
mac$take2 <- "CD68+ Myeloid"
 
LLc1 <- subset(all, coordinates %in% c(LLc, LLc2, collagen4a)) 
LLc1$take2 <- "LL Fibroblasts"

macLL <- subset(all, coordinates %in% c(LLmac7a, LLmac7b)) 
macLL$take2 <- "LL Macrophages"

total <- c(exclude$coordinates, vessel$coordinates, fibrous$coordinates, col4$coordinates, mac$coordinates, LLc1$coordinates, check$coordinates, macLL$coordinates)
include <- subset(all, coordinates %in% total, invert=T)

anno <- rbind(vessel[["take2"]], include[["take2"]], exclude[["take2"]], check[["take2"]],
              col4[["take2"]], fibrous[["take2"]], mac[["take2"]], LLc1[["take2"]], macLL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

individual.objects <- c(exclude1, exclude1a, exclude2a, exclude2b, exclude2e, exclude2f, exclude2g, exclude3a, 
                        exclude3b, exclude4, exclude8, fibrous2, fibrous6a, fibrous6a2, fibrous6b, 
                        vessel5, vessel7, vessel7a,col4,mac3, mac4, mac6a,LLc, LLc2, collagen4a,
                        LLmac7a, LLmac7b)

####################################################################################

# Clarify cell names
Idents(all) <- "take2"
all <- RenameIdents(all, "Fibrin-2"="Fibrin", "Arteriole"="SMA-hi vessels", 
                    "Vessels"="SMA-hi vessels","CLU+ LL"="CLU+ Fibroblasts","SL Fibroblasts"="COL1-hi Fibroblasts")
all$take2 <- all@active.ident

####################################################################################

# Further refinements
Idents(all) <- "take2"
CH <- subset(all, idents=c("COMP-hi Fibroblasts"))

Idents(CH) <- "orig.ident"

CH <-FindNeighbors(CH, dims = 1:21, reduction= "harmony", assay="SCT")
CH <-FindClusters(CH, resolution = c(0.08), graph.name="SCT_snn")
table(CH$orig.ident, CH@active.ident)

Idents(CH) <- "SCT_snn_res.0.08"
CH <- RenameIdents(CH, "0"="Sparse cells", "1"="Sparse cells", "2"="COMP-hi fibroblasts")
CH$take2 <- CH@active.ident

# Annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents="COMP-hi Fibroblasts", invert=T)
anno <- rbind(nLL[["take2"]], CH[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "take2"
RC <- subset(all, idents=c("Endothelial cells-1"))

Idents(RC) <- "orig.ident"
RCa <- subset(RC,  SMA < 7.5 & APOD < 6.8)
RCa$take2 <- "RBCs"
RCb <- subset(RC, SMA < 7.5 & APOD < 6.8, invert=T)
RCb$take2 <- "Endothelial cells-1"

anno <- rbind(RCb[["take2"]], RCa[["take2"]])
RC <- AddMetaData(RC, anno, col.name="take2")
table(RC$orig.ident,RC$take2)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("Endothelial cells-1"), invert=T)
anno <- rbind(nLL[["take2"]], RC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "take2"
RC <- subset(all, idents=c("Endothelial cells-2"))

Idents(RC) <- "orig.ident"

RCa <- subset(RC,  APOD < 6)
RCa$take2 <- "RBCs"
RCb <- subset(RC,  APOD < 6, invert=T)
RCb$take2 <- "Endothelial cells-2"

anno <- rbind(RCb[["take2"]], RCa[["take2"]])
RC <- AddMetaData(RC, anno, col.name="take2")
table(RC$orig.ident,RC$take2)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("Endothelial cells-2"), invert=T)
anno <- rbind(nLL[["take2"]], RC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "take2"
all <- RenameIdents(all, "DKK3+ Fibroblasts"="RBCs","Neutrophils"="Neutrophil-rich", "Blood/Fibrin"="RBCs”, "Endothelial cells-1"="Endothelial cells",  "Endothelial cells-2"="Endothelial cells", "Low-staining Fibroblasts"="Low-staining cells", "COMP-hi fibroblasts" = "COMP-hi Fibroblasts"
all$take2 <- all@active.ident


####################################################################################

# Create object containing missing cells
#make seperate objects of coordinates
data_x_y -> data_x_y_2

for (i in 1:length(pt)) {
  data_x_y[[i]]$coord <- paste0(data_x_y[[i]]$Centroid.X.µm, data_x_y[[i]]$Centroid.Y.µm)
  data_x_y[[i]] <- data_x_y[[i]] %>% filter(!(coord %in% data_x_y_2[[i]]$coord)) 
  data_x_y[[i]]$named1 <- "Artefact"
}

data_x_y -> missing.coord
data_x_y_2 -> data_x_y


# Review all samples after annotation

names <- names(data_x_y)
xy_list <- list()
for (i in 1:length(data_x_y)) {
  xy_list[[i]] <- data_x_y[[names[[i]]]]
}

all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()


for (i in 1:length(pt)) {
  s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
  xy_list[[i]]$named1 <- s_obj_meta[[i]]$take2;
  test <- rbind(missing.coord[[i]], xy_list[[i]])
  ggplot_ls[[i]] <- ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
    geom_point(size=0.4)+theme_classic()+scale_color_manual(values = bright) + scale_x_reverse() &
    guides(colour = guide_legend(override.aes = list(size=5)));
  print(ggplot_ls[[i]])
}

for (i in 1:length(names)) {
  s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
  xy_list[[i]]$named1 <- s_obj_meta[[i]]$take2;
  test <- rbind(missing.coord[[i]], xy_list[[i]])
  ggplot_ls[[i]] <- ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
    geom_point(size=0.2)+theme_classic()+scale_color_manual(values = group.color) + scale_x_reverse() &
    NoLegend() & thin
}




####################################################################################
####################################################################################

# Adult seropositive RA sample CellDIVE annotation

# Set graph graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
thin <- list(theme(axis.title.x=element_blank(), axis.title.y=element_blank()))

# Load segmentation files and create vector of sample names 
# This dataset was annotated as part of a larger set of slides of adult arthritis/synovitis and then subset for comparisons to only include 
treatment/ disease duration matched samples from those with a confirmed RA diagnosis

results <- dir("segmentation-files/", pattern = "*txt", 
               full.names = TRUE)

pt <- gsub("segmentation-files//", "", results)
sample <- substr(pt, 6,12)
pt <- substr(pt, 6,8)

data = list()
for (i in 1:length(results)) {
  data[[i]] <- read.delim(results[[i]])
}

# Format spatial, centroid x and y are your coordinates of each point
data_x_y = list()
for (i in 1:length(results)) {
  data_x_y[[i]] <- read.delim(results[[i]])
  rownames(data_x_y[[i]]) <- seq_along(data_x_y[[i]][,1])
  data_x_y[[i]] <- data_x_y[[i]] %>% select(c('Centroid.X.µm', 'Centroid.Y.µm'))
}

names(data_x_y) <- pt # assign names
names(data) <- pt

# Note markers here that failed that you don't want to include in the clustering
genesRemove <- c("CD31", "CD138")

# Take mean expression values for cell and create Seurat object
for (i in 1:length(data)) {
  data[[i]]$coord <- paste0(data[[i]]$Centroid.X.µm,data[[i]]$Centroid.Y.µm) # Create a unique tag for each cell based on it’s X and Y coordinates
  Cell.coord <- as.data.frame(data[[i]][["coord"]])
  data[[i]] <- data[[i]] %>% select(contains('.mean'))
  data[[i]] <- data[[i]] %>% select(contains('cell'))
  colnames(data[[i]]) <- gsub("Cell..", "", colnames(data[[i]])) # Format column names better
  colnames(data[[i]]) <- gsub(".mean", "", colnames(data[[i]]))
  rownames(data[[i]])<-seq_along(data[[i]][,1])
  data[[i]]<-as.data.frame(t(data[[i]]))
  data[[i]]<- filter(data[[i]], !(rownames(data[[i]]) %in% genesRemove))
  data[[i]]<-CreateSeuratObject(data[[i]])
  data[[i]]<- subset(data[[i]], subset=nCount_RNA>0)
  data[[i]]<-SCTransform(data[[i]])
  data[[i]] <- AddMetaData(data[[i]], metadata=Cell.coord, col.name="coordinates") # Re-assign the coordinate tag
}

for (i in 1:length(results)) {
  data[[i]]$orig.ident<- names(data[i]) 
}

####################################################################################

# Subset out the points that remain in the same place (continuing to have a DAPI stain at the start and end)

i <- 1
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 3.5 & `DAPI-FINAL` < 8.5 & `DAPI-INIT` > 5.8 & `DAPI-INIT` < 8.3)
data[[i]]<-SCTransform(data[[i]])

i <- 2
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 6 & `DAPI-INIT` > 6.4 & `DAPI-INIT` < 9.5)
data[[i]]<-SCTransform(data[[i]])

i <- 3
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL`> 6.3 & `DAPI-INIT` < 10)
data[[i]]<-SCTransform(data[[i]])

i <- 4
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL` < 10.3  & `DAPI-FINAL` > 6  & `DAPI-INIT` > 6 & `DAPI-INIT` < 10.3)
data[[i]]<-SCTransform(data[[i]])

i <- 5
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL` < 10  & `DAPI-FINAL` > 6.5  & `DAPI-INIT` > 6 & `DAPI-INIT` < 10)
data[[i]]<-SCTransform(data[[i]])

i <- 6
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-INIT` > 5.5 & `DAPI-INIT` < 10.3  & `DAPI-FINAL` > 5.5)
data[[i]]<-SCTransform(data[[i]])

i <- 7
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-INIT` > 6 & `DAPI-INIT` < 10  & `DAPI-FINAL` > 6)
data[[i]]<-SCTransform(data[[i]])

i <- 8
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]]<- subset(data[[i]], subset= `DAPI-FINAL` > 6  & `DAPI-INIT` > 6.5 & `DAPI-INIT` < 10.2)
data[[i]]<-SCTransform(data[[i]])

i <- 9
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0)
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 5 & `DAPI-INIT` > 5.8 & `DAPI-INIT` < 9.8)
data[[i]]<-SCTransform(data[[i]])

i <- 10
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) 
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 5 & `DAPI-INIT` > 5.5 & `DAPI-INIT` < 9)
data[[i]]<-SCTransform(data[[i]])

i <- 11
VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0)
data[[i]] <- subset(data[[i]], subset= `DAPI-FINAL`> 5.8 & `DAPI-FINAL` < 10)
data[[i]]<-SCTransform(data[[i]])

library(cowplot)

qcplot <- list()
for (i in 1:length(pt)) {
  qcplot[[i]] <- VlnPlot(data[[i]], features=c("DAPI-INIT", "DAPI-FINAL"), pt.size=0) & NoLegend()
}
plot_grid(plotlist = qcplot)

####################################################################################

# Filter the spatial coordinates by those that made the QC step above
for (i in 1:length(pt)) {
  data_x_y[[i]]$coord <- paste0(data_x_y[[i]]$Centroid.X.µm, data_x_y[[i]]$Centroid.Y.µm)
  data_x_y[[i]] <- filter(data_x_y[[i]], data_x_y[[i]]$coord %in% data[[i]]$coordinates)  
}

# Merge data into one object
all <- merge(x = data[[1]], y = data[-1], merge.data=TRUE)
VariableFeatures(all) <- check
all <- RunPCA(all, assay = 'SCT')

# Find principal componant with < 0.1% variance between it and the subsequent principal component to use for next steps
pct <- all[["pca"]]@stdev / sum(all[["pca"]]@stdev) * 100
npc <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

library(harmony)
all <- RunHarmony(all, group.by.vars = "orig.ident", assay.use="SCT", reduction="pca")
all <- RunUMAP(all, reduction="harmony", dims=1:npc)
all<-FindNeighbors(all, dims = 1:npc, reduction= "harmony", assay="SCT")
all<-FindClusters(all, resolution = c(0.1, 0.12, 0.14, 0.16, 0.18, 0.2), graph.name="SCT_snn")

# Included markers
gen2 <- c("SMA","COL4A1","CD146","COL6A1","DKK3", "COL3A1","COMP","FABP4","COL1","POSTN","SPARC", "MMP3", "PDPN","LYVE1", 
          "CD68",  "CD3", "MCT","CD20","CD15","CLU", "APOD", "DAPI-INIT", "DAPI-FINAL")

# Examine clusters for unique marker expression profiles
Idents(all) <- "SCT_snn_res.0.2"
DotPlot(all, features =gen2, cluster.idents = T) & formplot

# 0.16 resolution has no overlapping of marker profiles across clusters (identifies distinct cllusters without too much redundancy)
Idents(all) <- "SCT_snn_res.0.16"
DotPlot(all, features =gen2, cluster.idents = T) & formplot

####################################################################################

# Review how it looks with clustering

bright <- c( "plum1", "green", "yellow", "purple4", "darkorange",  "cyan","aquamarine3","deeppink4","#826D9B","bisque3","red",
             "lightpink", "forestgreen",  "steelblue3", "lightcyan",  "lemonchiffon", "blue",
            "gold2", "yellow4", "royalblue4" ,  "mediumpurple", "seashell","black",  "lightskyblue1","yellowgreen", "cyan3", "palevioletred",
            "lavenderblush2" ,"darkgrey","mediumorchid3","salmon2","mediumspringgreen","khaki" ,"brown")

i <- 11
nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$SCT_snn_res.0.16

ggplot(bright, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

# Zoom in on subsection
small <- test[test$Centroid.Y.µm < 3500  & test$Centroid.Y.µm > 2000 &  test$Centroid.X.µm > 3600 & test$Centroid.X.µm < 6100,]

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

####################################################################################

all$take2 <- all$SCT_snn_res.0.16

# Finer subclustering

Idents(all) <- "SCT_snn_res.0.16"
CO <- subset(all, idents=0)

#VlnPlot(CO, features=c("CD15", "CD3", "CD68", "COMP", "SPARC", "COL1"), pt.size=0)

COd <- subset(CO, CD68 <= 7.5)
COd$take2 <- "COMP-hi fibroblasts"
COc <- subset(CO, CD68 > 7.5)
COc$take2 <- "Macrophages"

anno <- rbind(COc[["take2"]], COd[["take2"]])
CO <- AddMetaData(CO, anno, col.name="take2")
table(CO$orig.ident, CO$take2)

Idents(CO) <- "take2"
VlnPlot(CO, features=c("CD15", "CD3", "CD68", "COMP", "SPARC", "MCT"), pt.size=0)
DotPlot(CO, features=gen2, cluster.idents=T) & formplot

# Annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=0, invert=T)
anno <- rbind(nLL[["take2"]], CO[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")


####################################################################################

# Tried subclustering more finely with FindClusters(). etc, but no major artefacts or subgroups, tried res 0.03, 0.4, 0.5, 0.8

Idents(all) <- "SCT_snn_res.0.16"
CC <- subset(all, idents=1)

CCa <- subset(CC, CD3 > 7.5)
CCa$take2 <- "T cells"
CCb <- subset(CC, CD3 <= 7.5 & POSTN > 8)
CCb$take2 <- "POSTN-hi fibroblasts"
CCc <- subset(CC, CD3 <= 7.5 & POSTN <= 8)
CCc$take2 <- "Low-staining cells"

anno <- rbind(CCa[["take2"]], CCb[["take2"]], CCc[["take2"]])
CC <- AddMetaData(CC, anno, col.name="take2")
table(CC$orig.ident, CC$take2)

#CC <-FindNeighbors(CC, dims = 1:21, reduction= "harmony", assay="SCT")
#CC <-FindClusters(CC, resolution = c(0.05), graph.name="SCT_snn")
#table(CC$orig.ident, CC@active.ident)

# Annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=1, invert=T)
anno <- rbind(nLL[["take2"]], CC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

# Subclustering unable to distinguish at res of 0.05, 0.08, 0.1, 0.12, 0.15
Idents(all) <- "SCT_snn_res.0.16"
AD <- subset(all, idents=2)

ADa <- subset(AD, CD68 > 7.4)
ADa$take2 <- "LL macrophages"
ADb <- subset(AD, CD68 <= 7.4)

ADc <- subset(ADb, PDPN > 7.5 & MMP3 > 7.8)
ADc$take2 <- "LL fibroblasts"
ADd <- subset(ADb, PDPN > 7.5 & MMP3 > 7.8, invert=T)

ADb1 <- subset(ADd, FABP4 <= 5.5 & COL1 > 5)
ADb1$take2 <- "COL3-hi fibroblasts"
ADb2 <- subset(ADd, FABP4 <= 5.5 & COL1 > 5, invert=T)
ADb2$take2 <- "Adipose"

anno <- rbind(ADa[["take2"]],ADc[["take2"]],ADb1[["take2"]],ADb2[["take2"]])
AD <- AddMetaData(AD, anno, col.name="take2")
table(AD$orig.ident, AD$take2)

Idents(AD) <- "take2"
VlnPlot(AD, features=c("FABP4", "COL3A1", "SPARC", "MMP3", "CD68", "PDPN"), pt.size=0)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=2, invert=T)
anno <- rbind(nLL[["take2"]], AD[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

# Subclustering not helpful for distinguishing tissue/cell types here

Idents(all) <- "SCT_snn_res.0.16"
CL <- subset(all, idents=3)

CLa <- subset(CL, CLU > 9.6 & SPARC > 7.5)
CLa$take2 <- "Fibrin"
CLb <- subset(CL, CLU > 9.6 & SPARC > 7.5, invert=T)

CLb1 <- subset(CLb, APOD > 7.2 & CD20 < 4 & CD68 < 6)
CLb1$take2 <- "Muscle"
CLc <- subset(CLb, APOD > 7.2 & CD20 < 4 & CD68 < 6, invert=T)

CLd1 <- subset(CLc, PDPN > 6.5 & CD68 <= 7.5)
CLd1$take2 <- "LL fibroblasts"
CLd2 <- subset(CLc, PDPN > 6.5 & CD68 > 7.5)
CLd2$take2 <- "LL macrophages"

CLe <- subset(CLc, PDPN <= 6.5 & FABP4 <= 5.8)
CLe$take2 <- "CLU-hi fibroblasts"
CLf <- subset(CLc, PDPN <= 6.5 & FABP4 > 5.6)
CLf$take2 <- "Adipose"

anno <- rbind(CLa[["take2"]], CLb1[["take2"]], CLe[["take2"]], CLd1[["take2"]], CLd2[["take2"]], CLf[["take2"]])
CL <- AddMetaData(CL, anno, col.name="take2")
table(CL$orig.ident, CL$take2)

Idents(CL) <- "take2"
VlnPlot(CL, features=c("SPARC", "PDPN", "CD3", "FABP4", "APOD", "CD68"), pt.size=0)

Idents(CL) <- "take2"
check <- subset(CL, idents="b")
Idents(check) <- "orig.ident"
DotPlot(check, features=gen2, cluster.idents = T) & formplot
VlnPlot(check, features=c("SPARC", "COL3A1", "CD146", "APOD", "CD20", "CD68"), pt.size=0)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=3, invert=T)
anno <- rbind(nLL[["take2"]], CL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "SCT_snn_res.0.16"
DK <- subset(all, idents=4)

DKa <- subset(DK, POSTN <= 5.5)
DKa$take2 <- "RBCs"
DKi <- subset(DK, POSTN > 5.5 & DKK3 > 8)
DKp <- subset(DK, POSTN > 5.5 & DKK3 < 8)
DKp$take2 <- "SMA-hi cells"

DKc <- subset(DKi, SPARC <= 7.6 & SMA > 7.5)
DKc$take2 <- "SMA-hi cells"
DKc1 <- subset(DKi, SPARC <= 7.6 & SMA <= 7.5)
DKc1$take2 <- "RBCs"
DKb <- subset(DKi, SPARC > 7.6 & CD146 > 4.8)
DKb$take2 <- "Small vessels"
DKb2 <- subset(DKi, SPARC > 7.6 & CD146 <= 4.8)
DKb2$take2 <- "COL4+ fibroblasts"

anno <- rbind(DKa[["take2"]], DKp[["take2"]], DKc[["take2"]], DKc1[["take2"]], DKb[["take2"]],DKb2[["take2"]])
DK <- AddMetaData(DK, anno, col.name="take2")
table(DK$orig.ident, DK$take2)

Idents(DK) <- "take2"
VlnPlot(DK, features=c("DKK3", "COL1", "COL6A1", "COL3A1"), pt.size=0)
table(DK$orig.ident, DK$take2)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=4, invert=T)
anno <- rbind(nLL[["take2"]], DK[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "SCT_snn_res.0.16"
SM <- subset(all, idents=5)

SMa <- subset(SM, CD146 <= 4.5 & SMA <= 7)
SMa$take2 <- "RBCs"
SMd <- subset(SM, CD146 <= 4.5 & SMA <= 7, invert=T)

SMc<- subset(SMd, DKK3 < 8 & CD146 < 4)
SMc$take2 <- "SMA-hi cells"
SMb <- subset(SMd, DKK3 < 8 & CD146 < 4, invert=T)

SMb <-FindNeighbors(SMb, dims = 1:21, reduction= "harmony", assay="SCT")
SMb <-FindClusters(SMb, resolution = c(0.05), graph.name="SCT_snn")
table(SMb$orig.ident, SMb@active.ident)

Idents(SMb) <- "SCT_snn_res.0.05"
SMb <- RenameIdents(SMb, "0"="SMA-hi vessels", "1"="Small vessels", "2"="SMA-hi vessels")
SMb$take2 <- SMb@active.ident

anno <- rbind(SMa[["take2"]], SMb[["take2"]], SMc[["take2"]])
SM <- AddMetaData(SM, anno, col.name="take2")
table(SM$orig.ident, SM$take2)

Idents(SM) <- "take2"
VlnPlot(SM, features=c("SMA", "CD146", "DAPI-INIT", "DKK3", "CD3", "COL4A1"), pt.size=0)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=5, invert=T)
anno <- rbind(nLL[["take2"]], SM[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

# Unable to find lymphatics despite clustering and using LYVE1 cut offs- manually examined slides, very few obvious vessels compared to paediatric slides
Idents(all) <- "SCT_snn_res.0.16"
LL <- subset(all, idents=6)

#VlnPlot(LL, features=c("SPARC", "PDPN", "MMP3", "CD68", "CD3", "CD15"), pt.size=0)

LLt <- subset(LL, CD68 <= 7.5)
LLt$take2 <- "LL fibroblasts"
LLm <- subset(LL, CD68 > 7.5)
LLo <- subset(LLm, PDPN <= 6 & COL3A1 <= 8)
LLo$take2 <- "Macrophages"
LLn <- subset(LLm, PDPN <= 6 & COL3A1 <= 8, invert=T)
LLn$take2 <- "LL macrophages"

anno <- rbind(LLt[["take2"]], LLn[["take2"]], LLo[["take2"]])
LL <- AddMetaData(LL, anno, col.name="take2")
table(LL$orig.ident, LL$take2)

Idents(LL) <- "take2"
VlnPlot(LL, features=c("SPARC", "PDPN", "COL3A1", "CD68", "MMP3", "CD15"), pt.size=0)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=6, invert=T)
anno <- rbind(nLL[["take2"]], LL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "SCT_snn_res.0.16"
C1 <- subset(all, idents=7)

C1 <-FindNeighbors(C1, dims = 1:21, reduction= "harmony", assay="SCT")
C1 <-FindClusters(C1, resolution = c(0.05), graph.name="SCT_snn")
table(C1$orig.ident, C1@active.ident)

Idents(C1) <- "SCT_snn_res.0.05"
C1 <- RenameIdents(C1, "0"="COL1-hi fibroblasts", "1"="POSTN+COMP+ fibroblasts")
C1$take2 <- C1@active.ident

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=7, invert=T)
anno <- rbind(nLL[["take2"]], C1[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "SCT_snn_res.0.16"
C6 <- subset(all, idents=9)

#VlnPlot(C6, features=c("CD3", "CD20", "MCT", "DKK3"), pt.size=0)

C6a <- subset(C6, DKK3 > 8)
C6a$take2 <- "RBCs"
C6c <- subset(C6, DKK3 <= 8 & MCT > 6)
C6c$take2 <- "Mast cells"
C6b <- subset(C6, DKK3 <= 8 & MCT <= 6)

C6b <-FindNeighbors(C6b, dims = 1:21, reduction= "harmony", assay="SCT")
C6b <-FindClusters(C6b, resolution = c(0.03), graph.name="SCT_snn")
table(C6b$orig.ident, C6b@active.ident)

Idents(C6b) <- "SCT_snn_res.0.03"
C6b <- RenameIdents(C6b, "0"="T cells", "1"="B/T aggregates", "2"="B/T aggregates")
C6b$take2 <- C6b@active.ident

anno <- rbind(C6a[["take2"]], C6b[["take2"]], C6c[["take2"]])
C6 <- AddMetaData(C6, anno, col.name="take2")
table(C6$orig.ident, C6$take2)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=9, invert=T)
anno <- rbind(nLL[["take2"]], C6[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")


####################################################################################

# Tried subclustering cluster 10, no major distintions, just DAPI mainly

Idents(all) <- "SCT_snn_res.0.16"
LY <- subset(all, idents=12)

LYb <- subset(LY, CD15 > 5)
LYb$take2 <- "Artefact"
LYc <- subset(LY, CD15 <= 5)
LYc$take2 <- "SPARC-hi fibroblasts"

anno <- rbind(LYc[["take2"]], LYb[["take2"]])
LY <- AddMetaData(LY, anno, col.name="take2")
table(LY$orig.ident, LY$take2)

#annotate with new splits from clustering
Idents(all) <- "SCT_snn_res.0.16"
nLL <- subset(all, idents=12, invert=T)
anno <- rbind(nLL[["take2"]], LY[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "SCT_snn_res.0.16"
all <- RenameIdents(all, "8"="COL6-hi fibroblasts", "10"= "LYVE1+ Macrophages", "11"="Mast cells")
all$take2 <- all@active.ident

####################################################################################

# Further annotation

Idents(all) <- "take2"
MC <- subset(all, idents="Mast cells")
#VlnPlot(MC, features=c("CD15", "DKK3", "CD68", "COMP", "SMA", "COL1"), pt.size=0)

MC <-FindNeighbors(MC, dims = 1:21, reduction= "harmony", assay="SCT")
MC <-FindClusters(MC, resolution = c(0.15), graph.name="SCT_snn")
table(MC$orig.ident, MC@active.ident)

Idents(MC) <- "SCT_snn_res.0.15"
MC <- RenameIdents(MC, "0"="RBCs", "1"="Mast cells", "2"="Mast cells", "3"="Mast cells", "4"="RBCs", "5"="Adipose", "6"="Artefact", "7"="Mast cells")
MC$take2 <- MC@active.ident

Idents(all) <- "take2"
nLL <- subset(all, idents="Mast cells", invert=T)
anno <- rbind(nLL[["take2"]], MC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "take2"
LC <- subset(all, idents=c("COMP-hi fibroblasts"))
VlnPlot(LC, features=c("FABP4", "LYVE1", "SMA", "PDPN", "CD68", "COMP"), pt.size=0)

LCb <- subset(LC, LYVE1 > 6  & PDPN > 6 & CD68 < 4.5 & FABP4 < 4.5)
LCb$take2 <- "Lymphatics"
LCc <- subset(LC, LYVE1 > 6  & PDPN > 6 & CD68 < 4.5 & FABP4 < 4.5, invert=T)
LCc$take2 <- "COMP-hi fibroblasts"

anno <- rbind(LCc[["take2"]], LCb[["take2"]])
LC <- AddMetaData(LC, anno, col.name="take2")
table(LC$orig.ident, LC$take2)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("COMP-hi fibroblasts"), invert=T)
anno <- rbind(nLL[["take2"]], LC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "take2"
EC <- subset(all, idents=c("RBCs"))
VlnPlot(EC, features=c("APOD", "CD146", "SPARC", "DAPI-INIT", "COL4A1", "DAPI-FINAL"), pt.size=0)

ECb <- subset(EC, SPARC > 7.5 & SMA > 6.5 & APOD < 8)
ECb$take2 <- "Small vessels"
ECc <- subset(EC, SPARC > 7.5 & SMA > 6.5 & APOD < 8, invert=T)
ECc$take2 <- "RBCs"

anno <- rbind(ECc[["take2"]], ECb[["take2"]])
EC <- AddMetaData(EC, anno, col.name="take2")
table(EC$orig.ident, EC$take2)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("RBCs"), invert=T)
anno <- rbind(nLL[["take2"]], EC[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################
####################################################################################

# Correct staining abnormalities with visual inspection of each fragment and slide
# Extract coordinates of regions requiring reassignment of identity per slide and correct at the end after going through all slides

i <- 1

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.Y.µm < 2200  & test$Centroid.Y.µm > 1600 &  test$Centroid.X.µm > 3600 & test$Centroid.X.µm < 4300,]
exclude1a <- small$coord

small <- test[test$Centroid.X.µm > 4200 & test$named1 %in% c("COL6-hi fibroblasts", "Muscle", "COL1-hi fibroblasts", "B/T aggregates", "T cells", "Small vessels"),]
exclude1b <- small$coord

small <- test[test$Centroid.Y.µm < 1500 & test$named1 %in% c("Muscle", "COL1-hi fibroblasts","B/T aggregates"),]
exclude1c <- small$coord

####################################################################################

i <- 2

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.Y.µm > 2200  & test$Centroid.Y.µm < 2900  & test$Centroid.X.µm > 6700 & test$Centroid.X.µm < 7140 & test$named1 %in% "RBCs",]
vessel2 <- small$coord

small <- test[test$Centroid.Y.µm > 2000  & test$Centroid.Y.µm < 2400  & test$Centroid.X.µm > 7100 & test$Centroid.X.µm < 7500 & test$named1 %in% c("SMA-hi vessels"),]
rbc2 <- small$coord

small <- test[test$Centroid.Y.µm > 3700  & test$Centroid.X.µm > 5500,]
exclude2a <- small$coord

small <- test[test$Centroid.Y.µm > 3150  & test$Centroid.Y.µm < 3700  & test$Centroid.X.µm > 5500 & test$Centroid.X.µm < 6480,]
exclude2b <- small$coord

small <- test[test$Centroid.Y.µm > 2000 & test$Centroid.Y.µm < 3600 & test$Centroid.X.µm < 3200,]

small <- test[test$Centroid.Y.µm < 3200 & test$Centroid.Y.µm > 3000 & test$Centroid.X.µm < 2550 & test$Centroid.X.µm > 2380 & test$named1 %in% "T cells",]
exclude2c <- small$coord

small <- test[test$Centroid.Y.µm > 4000 & test$Centroid.X.µm < 5000,]
small <- test[test$Centroid.Y.µm > 5500 & test$Centroid.X.µm < 3000 & test$named1 %in% "B/T aggregates",]
exclude2d <- small$coord

####################################################################################

i <- 3

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.Y.µm > 4900 & test$Centroid.Y.µm < 7000 & test$Centroid.X.µm < 2500,]
small <- test[test$Centroid.Y.µm > 5900 & test$Centroid.Y.µm < 7000 & test$Centroid.X.µm > 2000 & test$Centroid.X.µm < 2300 & test$named1 %in% "Mast cells",]
exclude3a <- small$coord

small <- test[test$Centroid.Y.µm > 6000 & test$Centroid.X.µm > 2000 & test$Centroid.X.µm < 4000,]
small <- test[test$Centroid.Y.µm > 6500 & test$Centroid.Y.µm < 6830 & test$Centroid.X.µm > 2300 & test$Centroid.X.µm < 3550 & test$named1 %in% "Small vessels",]
exclude3b <- small$coord
small <- test[test$Centroid.Y.µm > 6830 & test$Centroid.Y.µm < 7500 & test$Centroid.X.µm > 3300 & test$Centroid.X.µm < 3550 & test$named1 %in% "Small vessels",]
exclude3c <- small$coord

small <- test[test$Centroid.X.µm > 4000,]
small <- test[test$Centroid.X.µm > 4000 & test$Centroid.Y.µm < 6000 & test$named1 %in% "T cells",]
exclude3d <- small$coord

small <- test[test$Centroid.X.µm > 4000  & test$Centroid.Y.µm > 6000,]
small <- test[test$Centroid.Y.µm < 4500  & test$Centroid.X.µm < 4000,]

small <- test[test$Centroid.Y.µm > 1600  & test$Centroid.Y.µm < 1900  & test$Centroid.X.µm > 2100 & test$named1 %in% c("T cells", "B/T aggregates", "COL1-hi fibroblasts", "Macrophages", "LL macrophages", "LL fibroblasts"),]
exclude3e <- small$coord

small <- test[test$Centroid.Y.µm < 1550   & test$Centroid.X.µm < 2700 & test$named1 %in% c("T cells", "Small vessels", "B/T aggregates", "Mast cells"),]
fibrin3 <- small$coord

small <- test[test$Centroid.X.µm > 4500  & test$Centroid.Y.µm < 6000 & test$named1 %in% c("COL4-hi fibroblasts", "COMP-hi fibroblasts"),]
rbc3 <- small$coord
small <- test[test$Centroid.X.µm > 5800 & test$Centroid.Y.µm > 4250 & test$Centroid.Y.µm < 6000 & test$named1 %in% c("T cells"),]
fibrin3b <- small$coord
small <- test[test$Centroid.X.µm > 4500  & test$Centroid.X.µm < 5800 & test$Centroid.Y.µm > 4500 & test$Centroid.Y.µm < 6000 & test$named1 %in% c("T cells"),]
fibrin3b <- small$coord

small <- test[test$named1 %in% "Endothelial cells",]
rbcs3a <- small$coord

####################################################################################

i <- 4

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 4000  & test$Centroid.Y.µm < 5000,]
small <- test[test$Centroid.X.µm < 1600  & test$Centroid.Y.µm < 5000 & test$Centroid.Y.µm > 3000 & test$named1 %in% "SPARC-hi fibroblasts",]
LLfibro4 <- small$coord

small <- test[test$Centroid.X.µm > 4000  & test$Centroid.Y.µm > 2500 & test$Centroid.Y.µm < 6500,]
small <- test[test$Centroid.X.µm > 8000  & test$Centroid.Y.µm < 3750 & test$Centroid.Y.µm > 3000 & test$named1 %in% "Adipose",]
CLU4 <- small$coord

small <- test[test$Centroid.X.µm > 6000  & test$Centroid.Y.µm < 2000,]
small <- test[test$Centroid.X.µm > 6500 & test$Centroid.X.µm < 8000  & test$Centroid.Y.µm < 1100 & test$named1 %in% "T cells",]
exclude4a <- small$coord

small <- test[test$Centroid.X.µm > 6500 & test$Centroid.X.µm < 8600  & test$Centroid.Y.µm < 1000 & test$named1 %in% "COL3-hi fibroblasts",]
LLfibro4a <- small$coord

small <- test[test$Centroid.X.µm < 4500  & test$Centroid.Y.µm > 7400,]
small <- test[test$Centroid.X.µm < 4500  & test$Centroid.Y.µm > 7400 & test$named1 %in% c("COL3-hi fibroblasts"),]

########################## 

relabel <- subset(all, coordinates %in% small$coord)
relabel <-FindNeighbors(relabel, dims = 1:21, reduction= "harmony", assay="SCT")
relabel <-FindClusters(relabel, resolution = c(0.08), graph.name="SCT_snn")
table(relabel$orig.ident, relabel@active.ident)
relabel <- RenameIdents(relabel, "0"="LL fibroblasts", "1"="COL3-hi fibroblasts")
relabel$take2 <- relabel@active.ident

nLL <- subset(all, coordinates %in% small$coord, invert=T)
anno <- rbind(relabel[["take2"]], nLL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")


####################################################################################

i <- 5

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 1100  & test$Centroid.Y.µm < 1300,]
exclude5a <- small$coord

small <- test[test$Centroid.X.µm > 1500  & test$Centroid.X.µm < 2600  & test$Centroid.Y.µm < 1000,]
exclude5b <- small$coord

small <- test[test$Centroid.X.µm > 2800  & test$Centroid.Y.µm < 3500,]

small <- test[test$Centroid.X.µm > 2800  & test$Centroid.Y.µm < 3500 & test$Centroid.Y.µm > 2000,]
small <- test[test$Centroid.X.µm > 4100  & test$Centroid.X.µm < 4250  & test$Centroid.Y.µm < 2750 & test$Centroid.Y.µm > 2480 & test$named1 %in% "Mast cells",]
rbc5 <- small$coord
small <- test[test$Centroid.X.µm > 3000  & test$Centroid.X.µm < 3250  & test$Centroid.Y.µm < 2900 & test$Centroid.Y.µm > 2500 & test$named1 %in% "Mast cells",]
rbc5b <- small$coord
small <- test[test$Centroid.X.µm > 3200  & test$Centroid.X.µm < 4000  & test$Centroid.Y.µm > 1500 & test$Centroid.Y.µm < 2160 & test$named1 %in% "Mast cells",]
rbc5c <- small$coord
small <- test[test$Centroid.X.µm > 3100  & test$Centroid.X.µm < 3600  & test$Centroid.Y.µm > 2100 & test$Centroid.Y.µm < 2280 & test$named1 %in% c("Mast cells", "T cells", "Macrophages"),]
rbc5d <- small$coord
small <- test[test$Centroid.X.µm > 3000  & test$Centroid.X.µm < 3200  & test$Centroid.Y.µm > 2500 & test$Centroid.Y.µm < 2900 & !(test$named1 %in% "COL1-hi fibroblasts"),]
rbc5e <- small$coord

small <- test[test$Centroid.X.µm > 3950  & test$Centroid.X.µm < 4050  & test$Centroid.Y.µm > 2700 & test$Centroid.Y.µm < 2750 & test$named1 %in% "T cells",]
exclude5c <- small$coord

small <- test[test$Centroid.X.µm > 2800 & test$Centroid.Y.µm < 2000,]
small <- test[test$Centroid.X.µm > 3270 & test$Centroid.Y.µm < 1350 & test$Centroid.Y.µm > 1000 & test$named1 %in% "COL1-hi fibroblasts",]
exclude5d <- small$coord

small <- test[test$Centroid.X.µm > 1000 & test$Centroid.X.µm < 2800 & test$Centroid.Y.µm < 2800 & test$Centroid.Y.µm > 1000,]
small <- test[test$Centroid.X.µm > 1750 & test$Centroid.X.µm < 2100 & test$Centroid.Y.µm > 2350 & test$Centroid.Y.µm < 2700 & test$named1 %in% "Mast cells",]
exclude5e <- small$coord

small <- test[test$Centroid.X.µm < 1500 & test$Centroid.Y.µm > 1600,]
small <- test[test$Centroid.X.µm > 1500 & test$Centroid.Y.µm > 3200,]

small <- test[test$Centroid.X.µm > 2250 & test$Centroid.X.µm < 2550 & test$Centroid.Y.µm > 3200 & test$Centroid.Y.µm < 3700 & test$named1 %in% c("Mast cells", "SMA-hi cells", "SMA-hi vessels", "T cells"),]
rbc5f <- small$coord

small <- test[test$named1 %in% "Endothelial cells",]
rbcs5g <- small$coord

small <- test[test$named1 %in% "Macrophages",]
LLmac <- small$coord

####################################################################################

i <- 6

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 3000  & test$Centroid.Y.µm < 6000,]
small <- test[test$Centroid.X.µm < 3000 & test$Centroid.X.µm > 1500  & test$Centroid.Y.µm < 3000,]
small <- test[test$Centroid.X.µm > 8000  & test$Centroid.Y.µm < 3000,]

small <- test[test$Centroid.X.µm > 9500  & test$Centroid.X.µm < 9700  & test$Centroid.Y.µm < 1200 & test$Centroid.Y.µm > 1000 & test$named1 %in% "T cells",]
exclude6a <- small$coord

small <- test[test$Centroid.X.µm > 2000  & test$Centroid.Y.µm > 6000,]

####################################################################################

i <- 7

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

########################## 

small <- test[test$named1 %in% "T cells",]
relabel2 <- subset(all, coordinates %in% small$coord)
relabel2 <-FindNeighbors(relabel2, dims = 1:21, reduction= "harmony", assay="SCT")
relabel2 <-FindClusters(relabel2, resolution = c(0.15), graph.name="SCT_snn")
table(relabel2$orig.ident, relabel2@active.ident)

Idents(relabel2) <- "SCT_snn_res.0.15"
relabel2 <- RenameIdents(relabel2, "0"="T cells", "1"="POSTN-hi fibroblasts", "2"="T cells", "3"="B/T aggregates")
relabel2$take2 <- relabel2@active.ident

nLL <- subset(all, coordinates %in% small$coord, invert=T)
anno <- rbind(relabel2[["take2"]], nLL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

########################## 

small <- test[test$Centroid.Y.µm < 1700,]
exclude7a <- small$coord

small <- test[test$Centroid.X.µm > 5000  & test$Centroid.X.µm < 7500  & test$Centroid.Y.µm <3000 & test$Centroid.Y.µm > 1700,]
small <- test[test$Centroid.X.µm > 5500 & test$Centroid.Y.µm > 3200 & test$Centroid.Y.µm < 8000,]

small <- test[test$Centroid.X.µm > 8250 & test$Centroid.Y.µm > 5250 & test$Centroid.Y.µm < 7500,]
exclude7b <- small$coord

small <- test[test$Centroid.X.µm > 8000 & test$Centroid.Y.µm > 7500 & test$Centroid.Y.µm < 7560,]
exclude7c <- small$coord

small <- test[test$Centroid.X.µm > 5500 & test$Centroid.Y.µm > 7500,]
small <- test[test$Centroid.X.µm > 8350 & test$Centroid.Y.µm > 7500 & test$Centroid.Y.µm < 7800 & test$named1 %in% "POSTN-hi fibroblasts",]
exclude7d <- small$coord

small <- test[test$Centroid.X.µm < 5500 & test$Centroid.Y.µm > 2500 & test$Centroid.Y.µm < 7300,]

small <- test[test$Centroid.X.µm < 2600 & test$Centroid.Y.µm > 4000 & test$Centroid.Y.µm < 5000,]

small <- test[test$Centroid.X.µm < 2600 & test$Centroid.X.µm > 1900 & test$Centroid.Y.µm > 4370 & test$Centroid.Y.µm < 4550 & test$named1 %in% c("POSTN-hi fibroblasts", "COL1-hi fibroblasts", "POSTN+COMP+ fibroblasts", "Mast cells"),]
exclude7e <- small$coord

small <- test[test$Centroid.X.µm < 1500 & test$Centroid.X.µm > 700 & test$Centroid.Y.µm > 4150 & test$Centroid.Y.µm < 4500 & test$named1 %in% c("T cells", "B/T aggregates", "SMA-hi cells"),]
exclude7f <- small$coord

small <- test[test$Centroid.X.µm < 6000 & test$Centroid.X.µm > 5200 & test$Centroid.Y.µm > 4900 & test$Centroid.Y.µm < 6050 & test$named1 %in% c("T cells", "B/T aggregates", "SMA-hi cells", "COMP-hi fibroblasts"),]
exclude7g <- small$coord

small <- test[test$Centroid.X.µm > 2900 & test$Centroid.X.µm < 3500 & test$Centroid.Y.µm > 7700 & test$Centroid.Y.µm < 8000 & test$named1 %in% c("COL1-hi fibroblasts"),]
exclude7h <- small$coord

small <- test[test$Centroid.X.µm < 5000 & test$Centroid.Y.µm > 7000,]
small <- test[test$Centroid.X.µm > 2900 & test$Centroid.X.µm < 3500 & test$Centroid.Y.µm > 7800 & test$named1 %in% c("POSTN-hi fibroblasts"),]
exclude7i <- small$coord

small <- test[test$named1 %in% "Lymphatics",]
exclude7j <- small$coord

small <- test[test$named1 %in% "Macrophages",]
LLmac7a <- "LL macrophages"

####################################################################################

i <- 8

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 2300 & test$Centroid.Y.µm > 2000 & test$Centroid.Y.µm < 4200,]
small <- test[test$Centroid.X.µm < 2300 & test$Centroid.Y.µm > 4200,]

small <- test[test$Centroid.X.µm < 2300 & test$Centroid.Y.µm > 4200 & test$named1 %in% "COMP-hi fibroblasts",]
muscle8a <- small$coord

small <- test[test$Centroid.X.µm < 2300 & test$Centroid.Y.µm > 4200 & test$named1 %in% "LL fibroblasts",]
lymphatics8a <- small$coord

small <- test[test$Centroid.X.µm < 3700 & test$Centroid.Y.µm < 2000,]

small <- test[test$Centroid.X.µm > 3700 & test$Centroid.Y.µm < 2000 & test$Centroid.X.µm < 5500,]
small <- test[test$Centroid.X.µm > 3700 & test$Centroid.Y.µm < 2000 & test$Centroid.X.µm < 5500 & test$named1 %in% "Adipose",]
exclude8a <- small$coord

small <- test[test$Centroid.X.µm > 3500 & test$Centroid.Y.µm > 1700 & test$Centroid.Y.µm < 3300,]

small <- test[test$Centroid.X.µm > 5000 & test$Centroid.Y.µm > 3300,]

####################################################################################

i <- 9

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 2000 & test$Centroid.Y.µm < 2500,]
small <- test[test$Centroid.X.µm > 800 & test$Centroid.X.µm < 1550 & test$Centroid.Y.µm < 2000 & test$named1 %in% "B/T aggregates",]
exclude9a <- small$coord

small <- test[test$Centroid.X.µm < 2500 & test$Centroid.Y.µm < 3800 & test$Centroid.Y.µm > 2200,]

small <- test[test$Centroid.X.µm < 2500 & test$Centroid.Y.µm < 3800 & test$Centroid.Y.µm > 2200 & test$named1 %in% c("Mast cells","B/T aggregates"),]
rbc9a <- small$coord
small <- test[test$Centroid.X.µm < 2350 & test$Centroid.X.µm > 2000 & test$Centroid.Y.µm < 3500 & test$Centroid.Y.µm > 3000 & test$named1 %in% "LL fibroblasts",]
rbc9b <- small$coord

small <- test[test$Centroid.X.µm < 3000 & test$Centroid.Y.µm > 3800,]
small <- test[test$Centroid.X.µm < 3000 & test$Centroid.Y.µm > 3800 & test$named1 %in% "Muscle",]
exclude9a <- small$coord

small <- test[test$Centroid.X.µm > 2000 & test$Centroid.Y.µm <2000,]
small <- test[test$Centroid.X.µm > 2000 & test$Centroid.X.µm < 3200 & test$Centroid.Y.µm < 800 & test$named1 %in% c("Mast cells, COL1-hi fibroblasts", "Muscle", "B/T aggregates"),]
exclude9b <- small$coord

small <- test[test$Centroid.X.µm < 3300 & test$Centroid.X.µm > 2700 & test$Centroid.Y.µm < 750 & test$named1 %in% c("T cells"),]
exclude9c <- small$coord

small <- test[test$Centroid.X.µm > 2000 & test$Centroid.Y.µm >2000 & test$Centroid.Y.µm < 3800,]
small <- test[test$Centroid.X.µm > 2000 & test$Centroid.Y.µm >2000 & test$Centroid.Y.µm < 3800 & test$named1 %in% "Muscle", ]
exclude9d <- small$coord

small <- test[test$Centroid.X.µm > 3500 & test$Centroid.Y.µm >2500,]
small <- test[test$Centroid.X.µm > 3500 & test$Centroid.Y.µm >2500 & test$named1 %in% "Muscle",]
exclude9e <- small$coord

small <-test[test$named1 %in% c("Endothelial cells"),]
rbc9c <- small$coord

small <- test[test$Centroid.X.µm > 1950 & test$Centroid.X.µm < 2400 & test$Centroid.Y.µm > 3050 & test$Centroid.Y.µm < 3700 & !(test$named1 %in% "RBCs"),]
exclude9f <- small$coord

small <- test[test$Centroid.X.µm > 700 & test$Centroid.X.µm < 1500 & test$Centroid.Y.µm > 3600 & test$Centroid.Y.µm < 4350,]
exclude9g <- small$coord

####################################################################################

i <- 10

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 1000 & test$Centroid.Y.µm < 4200  & test$Centroid.Y.µm > 2200,]

small <- test[test$Centroid.X.µm > 650 & test$Centroid.X.µm < 1300 & test$Centroid.Y.µm < 3950  & test$Centroid.Y.µm > 3200 & !(test$named1 %in% c("RBCs")),]
exclude10a <- small$coord
small <- test[test$Centroid.X.µm > 300 & test$Centroid.X.µm < 1000 & test$Centroid.Y.µm < 4000  & test$Centroid.Y.µm > 3750 & !(test$named1 %in% c("RBCs")),]
exclude10b <- small$coord

small <- test[test$Centroid.X.µm < 2200 & test$Centroid.Y.µm > 4300,]
small <- test[test$Centroid.X.µm < 2600 & test$Centroid.Y.µm < 1500,]
small <- test[test$Centroid.X.µm > 2600 & test$Centroid.Y.µm < 1500,]

small <- test[test$Centroid.X.µm > 2000 & test$Centroid.X.µm < 3500 & test$Centroid.Y.µm < 3000 & test$Centroid.Y.µm > 800,]
small <- test[test$Centroid.X.µm > 2200 & test$Centroid.X.µm < 3500 & test$Centroid.Y.µm < 2600 & test$Centroid.Y.µm > 1300 & test$named1 %in% c("Muscle", "B/T aggregates"),]
exclude10c <- small$coord

small <- test[test$Centroid.X.µm > 1000 & test$Centroid.X.µm < 2300 & test$Centroid.Y.µm < 3600 & test$Centroid.Y.µm > 1800,]
small <- test[test$Centroid.X.µm > 1250 & test$Centroid.X.µm < 1600 & test$Centroid.Y.µm < 2900 & test$Centroid.Y.µm > 2450 & test$named1 %in% c("T cells", "B/T aggregates", "Mast cells"),]
exclude10d <- small$coord

small <- test[test$Centroid.X.µm > 1000 & test$Centroid.X.µm < 3700 & test$Centroid.Y.µm < 4700 & test$Centroid.Y.µm > 2700,]
small <- test[test$Centroid.X.µm > 2400 & test$Centroid.X.µm < 2700 & test$Centroid.Y.µm < 4300 & test$Centroid.Y.µm > 3700 & test$named1 %in% c("T cells", "B/T aggregates"),]
exclude10e <- small$coord

small <- test[test$Centroid.X.µm >3000 & test$Centroid.Y.µm > 1900,]

small <- test[test$Centroid.X.µm < 2500 & test$Centroid.X.µm > 2200 &  test$Centroid.Y.µm > 1500  & test$Centroid.Y.µm < 2200 & test$named1 %in% c("Macrophages"),]
LLmacs10 <- small$coord

small <- test[test$Centroid.X.µm < 2300 & test$Centroid.Y.µm < 1500,]
small <- test[test$Centroid.X.µm < 1350 & test$Centroid.X.µm > 1200 & test$Centroid.Y.µm > 600 & test$Centroid.Y.µm < 1000 & test$named1 %in% "Macrophages",]
LLmacs10b <- small$coord

small <- test[test$Centroid.X.µm < 2300 & test$Centroid.X.µm > 2050 & test$Centroid.Y.µm > 1000 & test$Centroid.Y.µm < 1500 & test$named1 %in% "Macrophages",]
LLmacs10c <- small$coord

####################################################################################

i <- 11

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1)+theme_classic() & scale_y_reverse() & scale_color_manual(values=bright) &  
  guides(colour = guide_legend(override.aes = list(size=5)))

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5)))

small <- test[test$Centroid.X.µm < 3500 & test$Centroid.Y.µm < 3000,]
small <- test[test$Centroid.X.µm < 2800 & test$Centroid.X.µm > 2000 & test$Centroid.Y.µm < 2200 & test$Centroid.Y.µm > 1200,]
small <- test[test$Centroid.X.µm < 2000 & test$Centroid.X.µm > 1650 & test$Centroid.Y.µm > 800 & test$Centroid.Y.µm < 1450 & test$named1 %in% c("COL1-hi fibroblasts", "POSTN+COMP+ fibroblasts"),]
exclude11a <- small$coord

small <- test[test$Centroid.X.µm < 1200 & test$Centroid.Y.µm > 2200 & test$Centroid.Y.µm < 4200,]
small <- test[test$Centroid.X.µm < 3000 & test$Centroid.Y.µm > 3200,]
small <- test[test$Centroid.X.µm > 3000 & test$Centroid.X.µm < 4800 & test$Centroid.Y.µm < 4000,]
small <- test[test$Centroid.X.µm > 4000 & test$Centroid.X.µm < 6000 & test$Centroid.Y.µm > 2200 & test$Centroid.Y.µm < 3500,]

small <- test[test$Centroid.X.µm > 5500 & test$Centroid.X.µm < 5700 & test$Centroid.Y.µm > 3100 & test$Centroid.Y.µm < 3250 & test$named1 %in% "Adipose",]
sma11a <- small$coord

small <- test[test$Centroid.X.µm > 6000 & test$Centroid.Y.µm < 3800 & test$Centroid.Y.µm > 1900,]

small <- test[test$Centroid.X.µm > 4000 & test$Centroid.X.µm < 6700 & test$Centroid.Y.µm > 3500,]

small <- test[test$Centroid.X.µm < 2050 & test$Centroid.Y.µm < 1550 & test$named1 %in% "Lymphatics",]
exclude11b <- small$coord

small <- test[test$Centroid.Y.µm < 2700 & test$Centroid.Y.µm > 1000 & test$Centroid.X.µm < 2000 & test$Centroid.X.µm > 1130 & test$named1 %in% c("POSTN+COMP+ fibroblasts", "COMP-hi fibroblasts"),]
exclude11c <- small$coord

####################################################################################


exclude.cells <- c(exclude1a, exclude1b, exclude1c, exclude2a, exclude2b, exclude2c, exclude2d, exclude3a, exclude3b, 
                   exclude3c, exclude3d, exclude3e, exclude4a, exclude5a, exclude5b, exclude5c, exclude5d, exclude5e, 
                   exclude6a,exclude7a, exclude7b, exclude7c, exclude7d,exclude7e, exclude7f, exclude7g, exclude7h, 
                   exclude7i, exclude7j, exclude8a, exclude9a, exclude9b, exclude9c,exclude9d, 
                   exclude9e, exclude9f,  exclude9g,exclude10a, exclude10b, exclude10c, exclude10d, exclude10e, 
                   exclude11a, exclude11b, exclude11c)
exclude <- subset(all, coordinates %in% exclude.cells)
exclude$take2 <- "Artefact"

fibrous <- subset(all, coordinates %in% c(fibrin3, fibrin3b))
fibrous$take2 <- "Fibrin"

clu.cells <- subset(all, coordinates %in% c(CLU4))
clu.cells$take2 <- "CLU-hi fibroblasts"

lining <- subset(all, coordinates %in% c(LLfibro4, LLfibro4a))
lining$take2 <- "LL fibroblasts"

vessel <- subset(all, coordinates %in% c(vessel2, sma11a))
vessel$take2 <- "SMA-hi vessels"

rbcs <- subset(all, coordinates %in% c(rbc2,rbc3,rbcs3a,rbc5,rbc5b,rbc5c,rbc5d,rbc5e,rbc5f,rbcs5g,rbc9a,rbc9b, rbc9c))
rbcs$take2 <- "RBCs"

lymph <- subset(all, coordinates %in% c(lymphatics8a))
lymph$take2 <- "Lymphatics"

mus <- subset(all, coordinates %in% c(muscle8a))
mus$take2 <- "Muscle"

LLmacs <- subset(all, coordinates %in% c(LLmac, LLmac7a, LLmacs10, LLmacs10b, LLmacs10c))
LLmacs$take2 <- "LL macrophages"

total <- c(exclude$coordinates, fibrous$coordinates, clu.cells$coordinates,  vessel$coordinates, 
           lining$coordinates, rbcs$coordinates, lymph$coordinates, mus$coordinates, LLmacs$coordinates)
include <- subset(all, coordinates %in% total, invert=T)

anno <- rbind(include[["take2"]], exclude[["take2"]], fibrous[["take2"]], clu.cells[["take2"]], lining[["take2"]],
              vessel[["take2"]], rbcs[["take2"]], lymph[["take2"]], mus[["take2"]], LLmacs[["take2"]])  
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

# Further refinement of subclustering

# Looking for lymphatics
Idents(all) <- "take2"
LP <- subset(all, idents=c("LL fibroblasts"))

LPb <- subset(LP, LYVE1 > 7 & CD68 < 5)
LPb$take2 <- "Lymphatics"
LPc <- subset(LP,  LYVE1 > 7 & CD68 < 5, invert=T)
LPc$take2 <- "LL fibroblasts"

anno <- rbind(LPc[["take2"]], LPb[["take2"]])
LP <- AddMetaData(LP, anno, col.name="take2")
table(LP$orig.ident, LP$take2)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("LL fibroblasts"), invert=T)
anno <- rbind(nLL[["take2"]], LP[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

Idents(all) <- "take2"
VL <- subset(all, idents=c("SMA-hi vessels"))

VlnPlot(VL, features=c("FABP4", "LYVE1", "SMA", "PDPN", "CD68", "COMP"), pt.size=0)

VLb <- subset(VL, LYVE1 > 6.5 & COL4A1 > 7.8 & POSTN > 8.5 & PDPN <= 5.5)
VLb$take2 <- "Lymphatics"
VLc <- subset(VL, LYVE1 > 6.5 & COL4A1 > 7.8 & POSTN > 8.5 & PDPN <= 5.5, invert=T)
VLc$take2 <- "SMA-hi vessels"

anno <- rbind(VLc[["take2"]], VLb[["take2"]])
VL <- AddMetaData(VL, anno, col.name="take2")
table(VL$orig.ident, VL$take2)

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents=c("SMA-hi vessels"), invert=T)
anno <- rbind(nLL[["take2"]], VL[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

# Remove COMP artefact

Idents(all) <- "take2"
CH <- subset(all, idents=c("COMP-hi fibroblasts"))

Idents(CH) <- "orig.ident"
VlnPlot(CH, features=c("CLU", "LYVE1", "SMA", "PDPN", "CD68", "COMP"), pt.size=0)

CH <-FindNeighbors(CH, dims = 1:21, reduction= "harmony", assay="SCT")
CH <-FindClusters(CH, resolution = c(0.1), graph.name="SCT_snn")
table(CH$orig.ident, CH@active.ident)

Idents(CH) <- "SCT_snn_res.0.1"
CH <- RenameIdents(CH, "0"="Sparse cells", "1"="Adipose", "2"="Artefact", "3"="Sparse cells", "4"="Collagen-dense region", "5"="Artefact", "6"="Artefact", "7"="Artefact")
CH$take2 <- CH@active.ident

#annotate with new splits from clustering
Idents(all) <- "take2"
nLL <- subset(all, idents="COMP-hi fibroblasts", invert=T)
anno <- rbind(nLL[["take2"]], CH[["take2"]])
all <- AddMetaData(all, anno, col.name="take2")

####################################################################################

# Create object containing missing cells / those lost to artefactual staining

data_x_y -> data_x_y_2
load(file="2407_all-coordinates.RData")

for (i in 1:length(pt)) {
  data_x_y[[i]]$coord <- paste0(data_x_y[[i]]$Centroid.X.µm, data_x_y[[i]]$Centroid.Y.µm)
  data_x_y[[i]] <- data_x_y[[i]] %>% filter(!(coord %in% data_x_y_2[[i]]$coord)) 
  data_x_y[[i]]$named1 <- "Artefact"
}

data_x_y -> missing.coord
data_x_y_2 -> data_x_y

# Visualise the cell types in each sample

group.color <- c("SMA-hi vessels"="cyan2", "SMA-low vessels"="blue", "Endothelial cells"="red","RBCs"="brown",
                 "CLU-hi fibroblasts"="aquamarine3", "COL1-hi fibroblasts"="lemonchiffon",  "COL6-hi Fibroblasts"="snow3",
                 "POSTN-hi fibroblasts"="salmon","COL6-hi fibroblasts"="steelblue","Sparse cells"= "darkseagreen1",
                 "LL fibroblasts"= "magenta", "Collagen-dense region"="lightsalmon",
                 "SPARC-hi fibroblasts"="yellow4", 
                 "LL macrophages"="turquoise4", "Macrophages"="olivedrab1", "LYVE1+ macrophages"="orchid4",
                 "T cells" ="mediumorchid1", "SMA-hi cells"="lightskyblue","B/T aggregates"="yellow",
                 "Adipose"="green2",  "Mast cells"="darkorange","Muscle"="grey97", 
                 "Lymphatics"="palevioletred",
                 "Artefact"="azure","Fibrin"="skyblue2")

names <- names(data_x_y)
xy_list <- list()
for (i in 1:length(data_x_y)) {
  xy_list[[i]] <- data_x_y[[names[[i]]]]
}

all_meta <- all@meta.data
s_obj_meta <- list()
ggplot_ls <- list()

for (i in 1:length(pt)) {
  s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
  xy_list[[i]]$named1 <- s_obj_meta[[i]]$take2;
  test <- rbind(missing.coord[[i]], xy_list[[i]])
  ggplot_ls[[i]] <- ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
    geom_point(size=0.4)+theme_classic()+scale_color_manual(values = group.color) &  scale_y_reverse() &
    guides(colour = guide_legend(override.aes = list(size=5)));
  print(ggplot_ls[[i]])
}

for (i in 1:length(names)) {
  s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
  xy_list[[i]]$named1 <- s_obj_meta[[i]]$take2;
  test <- rbind(missing.coord[[i]], xy_list[[i]])
  ggplot_ls[[i]] <- ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
    geom_point(size=0.2)+theme_classic()+scale_color_manual(values = bright) &
    NoLegend() & thin
}

grid.arrange(grobs = ggplot_ls, ncol=4)  


####################################################################################

# Subset to only include relevant samples
Idents(all) <- "orig.ident"
all <- subset(all, idents=c(“115”, “127”, “031”, “195”, “020”, “240”))

# Order clusters for visualisation
iden <- unique(all$take2)
iden2 <- c("B/T aggregates", "T cells")
iden3 <- as.character(iden[grep("LL", iden)])
iden4  <- c("SL macrophages", "LYVE1+ macrophages") 
iden5  <- as.character(iden[grep("POSTN-hi", iden)]) 
iden5c  <- as.character(iden[grep("SPARC", iden)]) 
iden5d  <- as.character(iden[grep("COL6", iden)]) 
iden5e  <- as.character(iden[grep("COL1", iden)]) 
iden5f  <- as.character(iden[grep("CLU", iden)]) 
iden6 <- as.character(iden[grep("vess|Lymph", iden)]) 
iden7 <- as.character(iden[grep("Mast", iden)]) 
iden8 <- c("Adipose", "Muscle", "Fibrin")
iden9 <- as.character(iden[grep("SMA-hi|RBCs|Collagen|Sparse", iden)]) 

orderclust <- unique(c(iden2, iden3,iden4,iden6,iden5,iden5c,iden5d, iden5e, iden5f,iden7,iden8, iden9))

Idents(all) <- "take2"
dotp <- subset(all, idents="Artefact", invert=T)
levels(dotp) <- rev(orderclust)

gen2 <- c("CD20", "CD3", "PDPN", "MMP3", "COL3A1", "CD68", "LYVE1", "COL4A1", "SMA", "DKK3", "CD146", 
          "POSTN",  "SPARC","COL6A1", "COL1","COMP", "CLU",  "MCT", "CD15", "FABP4", "APOD")

# Unclear what SMA-hi cells are (not vessels) so remove for less ambiguity, red blood cells (RBCs) are part of artefact
incl <- iden[grep("RBCs|SMA-hi cells|Sparse|Collagen-den”, iden, invert=T)]
DotPlot(dotp, features =gen2, idents=incl) & formplot


####################################################################################

####################################################################################

####################################################################################


