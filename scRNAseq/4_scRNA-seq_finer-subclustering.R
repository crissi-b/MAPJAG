
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
library(SeuratObject)
library(harmony)

#################################################################

# set graph graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
sparseplot <- list(NoAxes(), NoLegend())
thin <- list(theme(axis.title.x=element_blank(), axis.title.y=element_blank()), NoAxes())
colurs <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), brewer.pal(8, "Paired"))

#################################################################

# Stromal clustering
Idents(PBMC1)<-'cm_global'
stroma <-subset(x = PBMC1, idents = c("fibs", "pericytes", "Endo_cells", "lymphatic")   )
stroma <- PBMC1[,grepl("fibs|pericytes|Endo_cells|lymphatic", PBMC1$cm_global, ignore.case=TRUE)]


stroma <- stroma %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
DimPlot(stroma)

stroma<-FindNeighbors(stroma, dims=1:30)
stroma<-FindClusters(stroma, resolution = c(0.01, 0.05, 0.1, 0.2))
stroma<-FindClusters(stroma, resolution = c(0.15))
stroma<-FindClusters(stroma, resolution = c(0.12))

cols <- ArchR::paletteDiscrete(stroma@meta.data[, "RNA_snn_res.0.12"])


DimPlot(stroma, group.by = "RNA_snn_res.0.12", label = T, cols=cols)
DimPlot(stroma, group.by = "orig.ident")
FeaturePlot(stroma, features=c("APOD", "APOE", "APOC1", "LYZ"))
FeaturePlot(stroma, features=c("PDGFRA", "COL1A1"))


Idents(stroma)<-'RNA_snn_res.0.12'
straom_markers<-FindAllMarkers(stroma, only.pos = T, logfc.threshold = 1)

FeaturePlot(stroma, features = c("PDGFRA", "COL1A1", "THY1", "CD248"))
DimPlot(stroma, group.by = "orig.ident")
table(stroma$orig.ident)

Idents(stroma)<-'RNA_snn_res.0.12'
levels(stroma)

stroma_clean <-subset(x = stroma, idents = c("0" ,"1", "2", "3", "5", "6", "8")   )


stroma_clean <- stroma_clean %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
DimPlot(stroma_clean)

FeaturePlot(stroma_clean, features=c("RUNX1", "SOX5", "CDH11"))

FeaturePlot(stroma_clean, features = c("PDGFRA", "COL1A1", "THY1", "CD248"))

stroma_clean<-FindNeighbors(stroma_clean, dims=1:30)
stroma_clean<-FindClusters(stroma_clean, resolution = c(0.01, 0.05, 0.1, 0.2))

Idents(stroma_clean)<-'orig.ident'
DimPlot(stroma_clean, split.by = "orig.ident", label = F)

stroma_clean<-CellSelector(p1, stroma_clean, ident = "remove")
levels(stroma_clean)


stroma_clean <-subset(x = stroma_clean, idents = c("0" ,"1", "2", "3", "5", "6")   )


stroma_clean <- stroma_clean %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
DimPlot(stroma_clean, group.by = "orig.ident")

FeaturePlot(stroma_clean, features=c("SOX5", "CDH11"))

FeaturePlot(stroma_clean, features = c("PDGFRA", "COL1A1", "THY1", "CD248"))

stroma_clean_h<-RunHarmony.Seurat_CM(stroma_clean, group.by.vars = "orig.ident")
stroma_clean_h<-RunUMAP(stroma_clean_h, reduction="harmony", dims=1:30)
DimPlot(stroma_clean_h, group.by = "orig.ident")
FeaturePlot(stroma_clean_h, features=c("RUNX1", "SOX5", "CDH11"))
FeaturePlot(stroma_clean_h, features = c("PDGFRA", "COL1A1", "THY1", "MKI67"))

pt <- table(stroma_clean$orig.ident , stroma_clean$RNA_snn_res.0.2)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) +RotatedAxis()

stroma_clean_h<-FindNeighbors(stroma_clean_h, dims=1:30, reduction="harmony")
stroma_clean_h<-FindClusters(stroma_clean_h, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3))
stroma_clean_h<-FindClusters(stroma_clean_h, resolution = c(0.4))
stroma_clean_h<-FindClusters(stroma_clean_h, resolution = c(0.5, 0.6))
DimPlot(stroma_clean_h, group.by = "RNA_snn_res.0.5", label = T)
FeaturePlot(stroma_clean_h, features = c("PRG4", "TSPAN15"))
FeaturePlot(stroma_clean_h, features = c("CD34", "PI16"))
FeaturePlot(stroma_clean_h, features = c("POSTN", "COMP"))
FeaturePlot(stroma_clean_h, features = c("SFRP1", "CXCL12"))

Idents(stroma_clean_h)<-'RNA_snn_res.0.5'
markers_stroma_clean_h<-FindAllMarkers(stroma_clean_h, only.pos = T, logfc.threshold = 1)

cols <- ArchR::paletteDiscrete(stroma_clean_h@meta.data[, "cm_global"])
DimPlot(stroma_clean_h, group.by = "cm_global", label = F, cols=cols)


ggplot(stroma_clean_h@meta.data, aes(x=orig.ident, fill=orig.ident)) + geom_bar()+ theme_ArchR() +RotatedAxis()
       
ggplot(stroma_clean_h@meta.data, aes(x=cm_global, fill=cm_global)) + geom_bar()+ theme_ArchR()

ggplot(stroma_clean_h@meta.data, aes(x=cm_global, fill=orig.ident)) + geom_bar(position = "fill")+ theme_ArchR()


cols <- ArchR::paletteDiscrete(stroma_clean_h@meta.data[, "RNA_snn_res.0.5"])
DimPlot(stroma_clean_h, group.by = "RNA_snn_res.0.5", label = T, cols=cols)

stroma_clean_h$clusters<-stroma_clean_h$RNA_snn_res.0.5

current.sample.ids <- c( "0" , "1" , "2" , "3",  "4" , "5" , "6",  "7" , "8" , "9" , "10", "11", "12")
new.sample.ids <- c( "Venous" , "CD34/MFAP5" , "POSTN" , "SFRP/CXCL12",  "Pericytes" , "LL" , "CXCL14",  "SOX5_CDH11" , "KDR_Arterial" , "FLi1_Callilary" , "MMPs", "Lymphatic", "Prolif")

stroma_clean_h@meta.data[["clusters"]] <- plyr::mapvalues(x = stroma_clean_h@meta.data[["clusters"]], from = current.sample.ids, to = new.sample.ids)
cols <- ArchR::paletteDiscrete(stroma_clean_h@meta.data[, "clusters"])
DimPlot(stroma_clean_h, group.by = "clusters", cols=cols, label = T, repel = T, label.box =T) +NoAxes()+NoLegend()

library(gsfisher)
#repeat but with numbered clusters, re run findall markers agter reordering idents so that they are in order!
Idents(stroma_clean_h)<-'clusters'
levels(stroma_clean_h)

bg_genes<-unique(markers_stroma_clean_h$gene)

annotation_gs <- fetchAnnotation(species="hs", ensembl_version=NULL, ensembl_host=NULL)

index <- match(markers_stroma_clean_h$gene, annotation_gs$gene_name)
markers_stroma_clean_h$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- bg_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]

seurat_obj.res <- markers_stroma_clean_h
seurat_obj <- stroma_clean_h
seurat_obj.res <- seurat_obj.res[!is.na(seurat_obj.res$ensembl),]
ensemblUni <- na.omit(ensemblUni)

go.results <- runGO.all(results=seurat_obj.res,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=2, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T,fill_colors = c("#d9f0a3", "#41ab5d", "#49006a"))




#################################################################
#################################################################


# Myeloid clustering
Idents(PBMC1)<-'global1'
DimPlot(PBMC1)
myel <-subset(x = PBMC1, idents = c("Myeloid")   )


myel <- myel %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

myel<-RunUMAP(myel, dims=1:30,  verbose = FALSE)
  
Idents(myel)<-'TYPE'
DimPlot(myel, group.by="TYPE")
DimPlot(myel, group.by="orig.ident")
DimPlot(myel, split.by="TYPE")


Idents(myel)<-'patient'

DimPlot(myel, group.by="patient")

myel<-FindNeighbors(myel, dims=1:30)
myel<-FindClusters(myel, resolution = c(0.01, 0.05, 0.1, 0.2))
myel<-FindClusters(myel, resolution = c(0.15))
myel<-FindClusters(myel, resolution = c(0.12))

DimPlot(myel, group.by = "RNA_snn_res.0.12", label = T)
DimPlot(myel, group.by = "RNA_snn_res.0.2", label = T)

Idents(myel)<-"T_clonotype_id"
DimPlot(myel)+NoLegend()

Idents(myel)<-"B_clonotype_id"
DimPlot(myel)+NoLegend()


DefaultAssay(myel)<-'RNA'
Idents(myel)<-'RNA_snn_res.0.2'
Allmarkers_res0.2<-FindAllMarkers(myel, only.pos = T, logfc.threshold = 0.5)

Idents(myel)<-'RNA_snn_res.0.2'
myel <-subset(x = myel, idents = c("0" , "1" , "2",  "3" , "4" , "5",  "6" , "8" , "9" , "10" ,"11", "12")   )


myel <- myel %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

myel_h<-RunHarmony.Seurat_CM(myel, group.by.vars = "orig.ident")
myel_h<-RunUMAP(myel_h, reduction="harmony", dims=1:30)

DimPlot(myel_h, group.by = "orig.ident")+NoAxes()
DimPlot(myel_h, group.by = "TYPE")+NoAxes()

myel_h<-FindNeighbors(myel_h, reduction = "harmony", dims=1:30)
myel_h<-FindClusters(myel_h, resolution = c(0.05, 0.1, 0.2))
DimPlot(myel_h, group.by = "RNA_snn_res.0.05")+NoAxes()
DimPlot(myel_h, group.by = "RNA_snn_res.0.1")+NoAxes()
DimPlot(myel_h, group.by = "RNA_snn_res.0.2")+NoAxes()
myel_h<-FindClusters(myel_h, resolution = c(0.3, 0.5, 0.4))
DimPlot(myel_h, group.by = "RNA_snn_res.0.3")+NoAxes()
DimPlot(myel_h, group.by = "RNA_snn_res.0.4")+NoAxes()
DimPlot(myel_h, group.by = "RNA_snn_res.0.5")+NoAxes()

myel_h_type<-RunHarmony.Seurat_CM(myel, group.by.vars = "TYPE")
myel_h_type<-RunUMAP(myel_h_type, reduction="harmony", dims=1:30)

DimPlot(myel_h_type, group.by = "orig.ident")+NoAxes()
DimPlot(myel_h_type, group.by = "TYPE")+NoAxes()

FeaturePlot(myel_h_type, features = c("MERTK", "LYVE1"))

myel_h_type<-FindNeighbors(myel_h_type, reduction = "harmony", dims=1:30)
myel_h_type<-FindClusters(myel_h_type, resolution = c(0.05, 0.1, 0.2))
DimPlot(myel_h_type, group.by = "RNA_snn_res.0.05")+NoAxes()
DimPlot(myel_h_type, group.by = "RNA_snn_res.0.1")+NoAxes()
DimPlot(myel_h_type, group.by = "RNA_snn_res.0.2")+NoAxes()
myel_h_type<-FindClusters(myel_h_type, resolution = c(0.3,  0.4))
DimPlot(myel_h_type, group.by = "RNA_snn_res.0.3")+NoAxes()
DimPlot(myel_h_type, group.by = "RNA_snn_res.0.4")+NoAxes()

Idents(myel_h_type)<-'TYPE'
DimPlot(myel_h_type, split.by = "TYPE")+NoAxes()



pt <- table(myel_h_type$orig.ident, myel_h_type$RNA_snn_res.0.3)
pt <- as.data.frame(pt)
pt<-pt[!pt$Var2 == "stromal_contam", ]
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis()

pt <- table(myel_hs$orig.ident, myel_hs$named)
pt <- as.data.frame(pt)
pt<-pt[!pt$Var2 == "stromal_contam", ]
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis()

myel_h_type_orig<-RunHarmony.Seurat_CM(myel, group.by.vars = c("TYPE", "orig.ident"))
myel_h_type_orig<-RunUMAP(myel_h_type_orig, reduction="harmony", dims=1:30)

DimPlot(myel_h_type_orig, group.by = "orig.ident")+NoAxes()
DimPlot(myel_h_type_orig, group.by = "TYPE")+NoAxes()

Idents(myel_h_type_orig)<-'TYPE'
DimPlot(myel_h_type_orig, split.by = "TYPE")+NoAxes()


FeaturePlot(myel_h_type_orig, features = c("MERTK", "LYVE1"))

myel_h_type_orig<-FindNeighbors(myel_h_type_orig, reduction = "harmony", dims=1:30)
myel_h_type_orig<-FindClusters(myel_h_type_orig, resolution = c(0.05, 0.1, 0.2))
DimPlot(myel_h_type_orig, group.by = "RNA_snn_res.0.05")+NoAxes()
DimPlot(myel_h_type_orig, group.by = "RNA_snn_res.0.1")+NoAxes()
DimPlot(myel_h_type_orig, group.by = "RNA_snn_res.0.2")+NoAxes()
myel_h_type_orig<-FindClusters(myel_h_type_orig, resolution = c(0.3,  0.4))
DimPlot(myel_h_type_orig, group.by = "RNA_snn_res.0.3", label=T)+NoAxes()

myel_h_type_orig<-FindClusters(myel_h_type_orig, resolution = c(0.6,  0.7))

DimPlot(myel_h_type_orig, group.by = "RNA_snn_res.0.7", label=T)+NoAxes()

pt <- table(myel_h_type_orig$TYPE, myel_h_type_orig$RNA_snn_res.0.3)
pt <- as.data.frame(pt)
pt<-pt[!pt$Var2 == "stromal_contam", ]
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank())+RotatedAxis()

Idents(myel_h)<-'RNA_snn_res.0.3'
markers_res0.3<-FindAllMarkers(myel_h, only.pos = T, logfc.threshold = 1)
write.csv(markers_res0.3, "/rds/projects/c/croftap-mapjagdata/MAPJAGv1/step3_subclusters/Chris/afterSoupX_stroma/markers_res0.3.csv")

Idents(myel_h)<-'TYPE'
DE_across_all_tissues<-FindAllMarkers(myel_h, only.pos = T, logfc.threshold = 1)
write.csv(DE_across_all_tissues, "/rds/projects/c/croftap-mapjagdata/MAPJAGv1/step3_subclusters/Chris/afterSoupX_stroma/Myel_DE_across_all_tissues.csv")

myel_h$named <- myel_h@meta.data[["RNA_snn_res.0.3"]]
Idents(myel_h) <- 'named'
levels(myel_h)
current.sample.ids <- c("0" , "1" , "10", "2" , "3" , "4" , "5", "6",  "7" , "8" , "9" )
new.sample.ids <- c("SPP1_FABP5",  "S1008A_FCN1" , "S1008A_FCN1", "MERTK_SELEOP",  "Il1B" , "DC_CD1C" , "NEAT1_SLC8A1",  "DC_CLEC9A",  "DC_LAMP3",  "DC_TCF7L2" , "stromal_contam" )

myel_h@meta.data[["named"]] <- plyr::mapvalues(x = myel_h@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)
DimPlot(myel_h, group.by = "named", label = T)+NoAxes()
Idents(myel_h) <- 'named'
levels(myel_h)

myel_hs <-subset(x = myel_h, idents = c("SPP1_FABP5",  "S1008A_FCN1" , "S1008A_FCN1", "MERTK_SELEOP",  "Il1B" , "DC_CD1C" , "NEAT1_SLC8A1",  "DC_CLEC9A",  "DC_LAMP3",  "DC_TCF7L2"  )   )


myel_hs <- myel_hs %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE) 


myel_hs<-RunHarmony.Seurat_CM(myel_hs, group.by.vars = "orig.ident")
myel_hs<-RunUMAP(myel_hs, reduction="harmony", dims=1:30)
DimPlot(myel_hs, group.by = "named", label = T)+NoAxes()

Idents(myel_hs)<-'TYPE'
DimPlot(myel_hs,group.by = "orig.ident", label = F)+NoAxes()

markers_res0_3_f=markers_res0_3[!markers_res0_3$cluster == "9",]
library(gsfisher)

bg_genes<-unique(markers_res0_3_f$gene)

#annotation_gs <- fetchAnnotation(species="hs", ensembl_version=NULL, ensembl_host=NULL)

index <- match(markers_res0_3_f$gene, annotation_gs$gene_name)
markers_res0_3_f$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- bg_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
markers_res0_3_f<-na.omit(markers_res0_3_f)
ensemblUni<-na.omit(ensemblUni)


go.results <- runGO.all(results=markers_res0_3_f,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=2, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T,fill_colors = c("#d9f0a3", "#41ab5d", "#49006a"))



#################################################################
#################################################################


# select out T cells / NK cells / ILC- the cycling cluster strongly expressed CD8
load(file="global-object.rdata")

DimPlot(PBMC1, group.by = "integrated_snn_res.0.15", label = T, repel=T) & NoLegend()
FeaturePlot(PBMC1, features=c("CD3G", "CD247", "NKG7", "CD8", "TRDC"))

Idents(Tcell) <- "integrated_snn_res.0.15"
Tcell <- subset(PBMC1, idents= c(1, 3, 5, 6, 12, 8, 4))

DefaultAssay(Tcell) <- "integrated"
Tcell <- Tcell %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

# Find contaminating clusters with finer subclustering
Tcell <- FindNeighbors(Tcell, reduction= "pca", dims = 1:30)
Tcell <- FindClusters(Tcell, resolution = 3, graph.name = "integrated_snn")
DimPlot(Tcell, group.by = "integrated_snn_res.3", label = T, repel=T) & NoLegend()

#################################################################

# Extract cycling myeloid cells from the cycling cluster
load(file="t-cell-object")

# Check expression of other cell types
DefaultAssay(Tcell) <- "RNA"
mang = unique(unlist(str_split('CD4, CHST11, scrublet_score, CD8A, CD8B, CD68, MAFB, CD1C, LYZ, SPARC, APOE, CST3, PLA2G2A, IGHG1, IGKC, MZB1, IGHG3, JCHAIN, COL1A1, AIF1, FCGRT, CHST11, TOX, PRKCH, CD3G, CXCL13, CD8A, CD8B, MKI67, TOP2A, PPBP, CD34, KLRF1, XCL1, MS4A1, CD27,', pattern = (', '))))
DotPlot(Tcell, features = mang, cluster.idents=T) & formplot & NoLegend()
mang2 = unlist(str_split('MS4A1, IGHG1, MZB1, CD27, GZMB, CD4, CHST11, PRKCH, CD3G, CD8A, CD8B, TRDC, MKI67, TOP2A, KLRF1, XCL1, LYZ, CD68, CD14, CD1C, CLEC9A, LAMP3,  PPBP, LYVE1, VWF, ACTA2, TAGLN, COL1A1, COL3A1', pattern = (', ')))
DotPlot(Tcell, features= mang2, cluster.idents = T) & formplot & NoLegend()
check <- as.data.frame(AverageExpression(Tcell, features =c("MS4A1", "LYZ", "CD68", "COL1A1", "COL3A1"), assay="RNA")) %>% t() %>% as.data.frame()
mang <- c("CD68", "LYZ", "CD14", "CD3G", "MKI67", "TOP2A")


# Subset cycling myeloid cells and extract the cellular barcodes
Idents(Tcell) <- "integrated_snn_res.3"
cyclM <- subset(Tcell, idents= c(55,48,42))
cyclM <- RenameIdents(cyclM, "55"="Likely Doublet", "48"="Likely Doublet", "42"="Cycling Myeloid")
cyclM$named <- cyclM@active.ident
clusters <- cyclM[["named"]]
write.csv(clusters, file="CyMy_clusters")

# Subset out doublets / contaminants that are not T cells
Tcell <- subset(Tcell, idents= c(55,48,42), invert = T)

DotPlot(Tcell, features= c("CD4", "CHST11", "PRKCH", "CD8A", "CD8B", "TRDC", "TRDV1", "MKI67"), cluster.idents = T) & formplot
FeaturePlot(Tcell, features=c("CD3G", "AOAH", "CD4", "CD8A", "TRDC", "TRAV1-2"), ncol=3, order=T) & thin
FeaturePlot(Tcell, features=c("RORA", "CHST11", "TOX", "FOXP3"), order=T)
FeaturePlot(Tcell, features=c("KLRB1", "CXCL13", "IL7R", "SELL"), ncol=2) & thin
VlnPlot(Tcell, features=c("nCount_RNA", "percent.mt", "nFeature_RNA"), pt.size=0)

# Data not included in this analysis
#Tcell$abTCR <- !(is.na(Tcell$T_clonotype_id))
#DimPlot(Tcell, group.by = "abTCR")
DimPlot(Tcell, group.by = "orig.ident")

#################################################################

# Data well-integrated without use of batch correction 

DefaultAssay(Tcell) <- "integrated"
Tcell <- Tcell %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA()

# Identify principal component with variance < 0.1% difference
pct <- Tcell[["pca"]]@stdev / sum(Tcell[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()
DimHeatmap(Tcell, assay= "integrated", dims = 18:23, cells = 500, balanced = TRUE)

Tcell <- RunUMAP(Tcell, dims = 1:20)
Tcell <- FindNeighbors(Tcell, reduction= "pca", dims = 1:20)
Tcell <- FindClusters(Tcell, resolution = c(0.1, 0.4, 0.6, 0.8,1,1.2), graph.name = "integrated_snn")
DimPlot(Tcell, group.by = "integrated_snn_res.0.1", label = T, repel=T) & NoLegend()

plot.list <- list()
res <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4)
for (i in 1:length(res)) {
  resg <- paste0("integrated_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(Tcell, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=3)  

DimPlot(Tcell, group.by = "integrated_snn_res.0.8", label = T, repel=T) & NoLegend()

#################################################################

# Decide on naive T cell cluster merging based on the number of differentially-expressed genes between them
Idents(Tcell) <- "integrated_snn_res.0.8"
check <- FindMarkers(Tcell, ident.1 =15, ident.2=18, logfc.threshold = 0.58)
check <- FindMarkers(Tcell, ident.1 =15, ident.2=6, logfc.threshold = 0.58)
check <- FindMarkers(Tcell, ident.1 =15, ident.2=1, logfc.threshold = 0.58)

# Combine similar clusters with few differentially expressed genes
Idents(Tcell) <- "integrated_snn_res.0.8"
Tcell <- RenameIdents(Tcell, "15"="1", "18"="1", "6"="1")
Tcell$integrated_snn_res.0.8 <- Tcell@active.ident
DimPlot(Tcell, group.by = "integrated_snn_res.0.8", label = T, repel=T, split.by="TYPE") & NoLegend()

#################################################################

# Check differentially-expressed genes between ILC subclusters
Idents(Tcell) <- "integrated_snn_res.0.8"
check <- FindMarkers(Tcell, ident.1 =21, ident.2=20, logfc.threshold = 0.58)
check <- FindMarkers(Tcell, ident.1 =10, ident.2=13, logfc.threshold = 0.58)
check <- FindMarkers(Tcell, ident.1 =20, ident.2=13, logfc.threshold = 0.58)

#################################################################

Idents(Tcell) <- "integrated_snn_res.0.8"
Tcell <- RenameIdents(Tcell, "0"="CD4+ KLRB1+ memory T", "1"="CD4+ naive/central memory T", "2"="CD8+ naive T", "3"="CXCL13+ TpH", 
                      "4"="CD8+ GZMB+/GZMK+ memory T", "5"="Activated NK-like T", "7"="CD8+ GZMK+ memory T", 
                      "8"="Vdelta2 yd", "9"="IFNG+ NK", "10"="THEMIS+ IL7R+ ILC", "11"="CD56-dim CD16+ NK", "12"="CD4+ FOXP3+ Tregs",
                      "14"="CD56+bright NK", "16"="MAIT cells", "17"="Cycling T", "19"="NK-like ILC1", "20"="THEMIS+ IL7R+ ILC", "21"="THEMIS+ IL7R+ ILC")
Tcell$simple8 <- Tcell@active.ident

Idents(Tcell) <- "simple8"
idents8 <- c("CD4+ naive/central memory T", "CD4+ KLRB1+ memory T", "CXCL13+ TpH", "CD4+ FOXP3+ Tregs", "2"="CD8+ naive T", "Activated NK-like T",
                                                 "CD8+ GZMK+ memory T", "CD8+ GZMB+/GZMK+ memory T", "Cycling T", "MAIT cells", "Vdelta2 yd", "CD56+bright NK", "IFNG+ NK","CD56-dim CD16+ NK",
                                                 "NK-like ILC1","THEMIS+ IL7R+ ILC")

levels(Tcell) <- rev(idents8)
Tcell$simple8 <- factor(Tcell$simple8, levels=rev(idents8))
 clusters <- Tcell[["simple8"]]
write.csv(clusters, file="T_clusters")


#################################################################
#################################################################

# Add scoring for expression of tissue resident memory markers 

# Define TRM gene module
moduleTRM <- list(c("ITGAE", "ITGA1", "PDCD1", "CXCR6"))

# Calculate module score
Tcell <- AddModuleScore(Tcell, features = moduleTRM, name = "TRM")

FeaturePlot(Tcell, features = "TRM1", split.by = "TYPE", 
            min.cutoff = "q10", max.cutoff = "q90", raster = FALSE) & 
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank())


#################################################################
#################################################################



# Subclustering B cells

#select out B and plasma from the global object
#load(file="global-object.rdata")
Idents(PBMC1) <- "integrated_snn_res.0.15"
Bcell <- subset(PBMC1, idents= c(7, 11))

DimPlot(Bcell, group.by = "label15", label = T, repel=T) & NoLegend()

# Remove contaminating clusters
DefaultAssay(Bcell) <- "integrated"
Bcell <- Bcell %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

# Find contaminating clusters using high resolution clustering
Bcell <- FindNeighbors(Bcell, reduction= "pca", dims = 1:30)
Bcell <- FindClusters(Bcell, resolution = 3, graph.name = "integrated_snn")
DimPlot(Bcell, group.by = "integrated_snn_res.3", label = T, repel=T) & NoLegend()

# Create gene list to plot
DefaultAssay(Bcell) <- "RNA"
mang2 = unlist(str_split('MS4A1, IGHG1, MZB1, CD27, GZMB, CD4, CHST11, PRKCH, CD3G, CD8A, CD8B, TRDC, MKI67, TOP2A, KLRF1, XCL1, LYZ, CD68, CD14, CD1C, CLEC9A, LAMP3, VWF, ACTA2, TAGLN, COL1A1, COL3A1', pattern = (', ')))
#plot it
DotPlot(Bcell, features= mang2, cluster.idents = T) & formplot & NoLegend()

mang = unique(unlist(str_split('CD4, CHST11, scrublet_score, CD8A, CD8B, CD68, MAFB, CD1C, LYZ, SPARC, APOE, CST3, PLA2G2A, IGHG1, IGKC, MZB1, IGHG3, JCHAIN, COL1A1, AIF1, FCGRT, CHST11, TOX, PRKCH, CD3G, CXCL13, CD8A, CD8B, MKI67, TOP2A, PPBP, CD34, KLRF1, XCL1, MS4A1, CD27', pattern = (', '))))
DotPlot(Bcell, features = mang) & formplot & NoLegend()

# See which clusters highly express myeloid/T cell genes 
check <- as.data.frame(AverageExpression(Bcell, features =c("CD8A", "CD3G", "CD68"), assay="RNA")) %>% t() %>% as.data.frame()

# Subset out doublet / non-B cell contaminants using fine clustering resolutions
Idents(Bcell) <- "integrated_snn_res.3"
Bcell <- subset(Bcell, idents= c(25, 27, 30), invert = T)

# Need to repeat the process now you have removed cells
DefaultAssay(Bcell) <- "integrated"
Bcell <- Bcell %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE) 

pct <- Bcell[["pca"]]@stdev / sum(Bcell[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

Bcell <- RunUMAP(Bcell, dims = 1:20, verbose = FALSE)
Bcell <- FindNeighbors(Bcell, reduction= "pca", dims = 1:20)
Bcell <- FindClusters(Bcell, resolution = 0.3, graph.name = "integrated_snn")
DimPlot(Bcell, group.by = "integrated_snn_res.0.3", label = T, repel=T) & NoLegend()
p1 <- DimPlot(Bcell, group.by = "TYPE")
p2 <- DimPlot(Bcell, group.by = "patient", cols = colurs)
grid.arrange(grobs = list(p1, p2), ncol=2) 

Bcell <- FindClusters(Bcell, resolution = c(0.05, 0.1, 0.5, 0.8, 1), graph.name = "integrated_snn")
DimPlot(Bcell, group.by = "integrated_snn_res.1", label = T, repel=T) & NoLegend()

# Plot multiple resolutions for inspection
plot.list <- list()
res <- c(0.05, 0.1, 0.3, 0.5, 0.8, 1)
for (i in 1:length(res)) {
  resg <- paste0("integrated_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(Bcell, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=3) 

# B cell markers
DefaultAssay(Bcell) <- "RNA"
FeaturePlot(Bcell, features=c("XBP1","TCL1A", "CD27", "ITGAX", "NR4A1", "BACH2"), ncol=3) & sparseplot

# Compare clusters to see if there are some that don't have many DE genes between them
Idents(Bcell) <- "integrated_snn_res.0.5"
check <- FindMarkers(Bcell, ident.1= 6, ident.2=3, logfc.threshold = 0.58)
check2 <- FindMarkers(Bcell, ident.1= 6, ident.2=10, logfc.threshold = 0.58)

# Check well-integrated
DimPlot(Bcell, group.by = "sample")
DimPlot(Bcell, group.by = "TYPE", split.by="TYPE")

# Batch correction by sample as some clusters not well-integrated across sample types- some versions of harmony throw an error code that needs to be corrected, advice available on online forums

DefaultAssay(Bcell) <- "RNA"
Bcell <- Bcell %>% 
  FindVariableFeatures(nfeatures=2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

Bcell <- RunHarmony.Seurat(Bcell, group.by.vars = c("sample"), assay.use = "RNA", plot_convergence = T, max.iter.harmony = 30)
pct <- Bcell[["pca"]]@stdev / sum(Bcell[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

Bcell <- RunUMAP(Bcell, reduction="harmony", dims = 1:20)
Bcell <- FindNeighbors(Bcell, reduction= "harmony", dims = 1:20, verbose = T)

Bcell <- FindClusters(Bcell, resolution = c(0.05, 0.1, 0.3, 0.5, 0.7), graph.name = "RNA_snn")
p3 <- DimPlot(Bcell, group.by = "RNA_snn_res.0.3", label = T, repel=T) & NoLegend()
tplot <- DimPlot(Bcell, group.by = "TYPE")
p2 <- DimPlot(Bcell, group.by = "patient", cols =colurs)

grid.arrange(grobs = list(p3, p1, p2), ncol=3) 

# Re-visualise clustering resolutions
plot.list <- list()
res <- c(0.05, 0.1, 0.3, 0.5, 0.7)
for (i in 1:length(res)) {
  resg <- paste0("RNA_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(Bcell, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=3) 

Idents(Bcell) <- "RNA_snn_res.0.3"
mang = unique(unlist(str_split('CD4, CHST11, scrublet_score, CD8A, CD8B, CD68, MAFB, CD1C, LYZ, SPARC, APOE, CST3, PLA2G2A, IGHG1, IGKC, MZB1, IGHG3, JCHAIN, COL1A1, AIF1, FCGRT, CHST11, TOX, PRKCH, CD3G, CXCL13, CD8A, CD8B, MKI67, TOP2A, PPBP, CD34, KLRF1, XCL1, MS4A1, CD27', pattern = (', '))))
DotPlot(Bcell, features = mang) & formplot & NoLegend()

Bcell <- subset(Bcell, idents=c(9,11), invert=T)

# Second round of removing doublet / non-B cell contaminants

DefaultAssay(Bcell) <- "RNA"
Bcell <- Bcell %>% 
  FindVariableFeatures(nfeatures=2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

Bcell <- RunHarmony.Seurat(Bcell, group.by.vars = c("sample"), assay.use = "RNA", plot_convergence = T, max.iter.harmony = 30)
pct <- Bcell[["pca"]]@stdev / sum(Bcell[["pca"]]@stdev) * 100
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 %>% return()
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1 %>% return()

Bcell <- RunUMAP(Bcell, reduction="harmony", dims = 1:18)
Bcell <- FindNeighbors(Bcell, reduction= "harmony", dims = 1:18, verbose = T)
Bcell <- FindClusters(Bcell, resolution = 0.5, graph.name = "RNA_snn")
DimPlot(Bcell, group.by = "RNA_snn_res.0.5", label = T, repel=T) & NoLegend()

p3 <- DimPlot(Bcell, group.by = "TYPE")
p2 <- DimPlot(Bcell, group.by = "patient", cols =colurs)
Bcell <- FindClusters(Bcell, resolution = c(0.05,0.1,0.6, 0.8,1), graph.name = "RNA_snn")

grid.arrange(grobs = list(p3, p1, p2), ncol=3) 

plot.list <- list()
res <- c(0.05, 0.1, 0.3, 0.4, 0.5, 0.6)
for (i in 1:length(res)) {
  resg <- paste0("RNA_snn_res.", res[[i]])
  plot.list[[i]] <- DimPlot(Bcell, group.by = resg, label = T) & sparseplot
}
grid.arrange(grobs = plot.list, ncol=3) 

# Inspect variable genes between clusters to identify those that are very similar
Idents(Bcell) <- "RNA_snn_res.0.5"
check <- FindMarkers(Bcell, ident.1=12, ident.2=4, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=11, ident.2=0, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=10, ident.2=0, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=10, ident.2=5, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=3, ident.2=0, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=8, ident.2=2, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=9, ident.2=3, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=9, ident.2=0, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=11, ident.2=4, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=7, ident.2=4, logfc.threshold = 0.58)
check <- FindMarkers(Bcell, ident.1=11, ident.2=12, logfc.threshold = 0.58)


# Merge cluster 10 with cluster 0
Bcell$label5 <- Bcell$RNA_snn_res.0.5
Idents(Bcell) <- "label5"
Bcell <- RenameIdents(Bcell, "10"="0")
Bcell$label5 <- Bcell@active.ident

# Find top markers for RNA
Markers <- FindAllMarkers(Bcell, logfc.threshold = 0.58, only.pos = T)
top4uniq <- unique(Markers %>% 
                     group_by(cluster) %>% 
                     top_n(n=6, wt = avg_log2FC) %>% 
                     dplyr::pull(gene))
DotPlot(Bcell, features = top4uniq, cluster.idents = T) + formplot

# Find top protein markers
DefaultAssay(Bcell) <- "Protein"
prot <- as.data.frame(unique(row.names(Bcell)))
Markers2 <- FindAllMarkers(Bcell, logfc.threshold = 0.58, only.pos = T)
Markers2 <- read.csv(file="SNN8-B-allprot.csv")
top4uniq <- unique(Markers2 %>% 
                     group_by(cluster) %>% 
                     top_n(n=6, wt = avg_log2FC) %>% 
                     dplyr::pull(gene))
DotPlot(Bcell, features = top4uniq) + formplot

# Plot key protein markers in SFMC and PBMC for inspection
Idents(Bcell) <- "TYPE"
citesub <- subset(Bcell, idents=c("Blood", "SF"))
Idents(citesub) <- "RNA_snn_res.0.5"
FeaturePlot(citesub, features=c("IGHD.1", "IGHM.1", "IGKC.1", "ITGAM.1", "ITGAX.1", "CD1C.1"), max.cutoff="q90", min.cutoff="q10", ncol=3)
DefaultAssay(Bcell) <- "RNA"

Idents(PBMC1) <- "RNA_snn_res.0.5"
PBMC1 <- RenameIdents(PBMC1, "0"="IgD+IgM+ Naive B", "1"="Age-associated B cells", "2"="Memory B cells", "3"="Transitional B", "4"="Plasma cells", 
                      "5"="IgD+ GC-like B", "6"="GC-like Memory B", "7"="Activated Plasma cells", 
                      "8"="PAX5+ Pro-B-like", "9"="Plasma cells")
PBMC1$label5 <- PBMC1@active.ident


b.labels <- rev(c("Activated Plasma cells", "Plasma cells","GC-like Memory B",
                  "IgD+ GC-like B", "Age-associated B cells", "Memory B cells",
                  "PAX5+ Pro-B-like", "Transitional B", "IgD+IgM+ Naive B"))
levels(PBMC1) <- b.labels
