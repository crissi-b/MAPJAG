
#xenium_chrissy <- `2404-xenium-global`
rm(`2404-xenium-global`)
DimPlot(xenium_chrissy, group.by = "global13", raster = F)

Idents(xenium_chrissy) <- 'global13'
Myeloid_xenium <- subset(xenium_chrissy, idents=c("Myeloid"))

Myeloid_xenium <- Myeloid_xenium %>% 
    ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:40)

library(harmony)
Myeloid_xenium <- RunHarmony(Myeloid_xenium, group.by.vars = c("orig.ident"))
Myeloid_xenium <- RunUMAP(Myeloid_xenium, reduction="harmony",  dims=1:20)
Myeloid_xenium <- FindNeighbors(Myeloid_xenium, reduction="harmony", dims=1:20)

DimPlot(Myeloid_xenium, group.by="orig.ident", raster=FALSE)+NoLegend()

Myeloid_xenium <- FindClusters(Myeloid_xenium, resolution = c(0.01, 0.05, 0.1, 0.2))

Myeloid_xenium <- FindClusters(Myeloid_xenium, resolution = c(0.3, 0.4))
Myeloid_xenium <- FindClusters(Myeloid_xenium, resolution = c(0.5, 0.6))


DimPlot(Myeloid_xenium, group.by = "RNA_snn_res.0.6", raster=FALSE)

Idents(Myeloid_xenium) <- 'RNA_snn_res.0.5'
myeloid_0.5_markers <- FindAllMarkers(Myeloid_xenium, only.pos = T)

DimPlot(Myeloid_xenium, group.by = "RNA_snn_res.0.5", raster=FALSE, label = T)
#rename and identify contams
Myeloid_xenium$named_sub <- Myeloid_xenium@meta.data[["RNA_snn_res.0.5"]]
Idents(Myeloid_xenium) <- 'named_sub'
levels(Myeloid_xenium)
current.sample.ids <- c( "0" , "1" , "10", "11", "12", "13", "14", "2" , "3",  "4",  "5",  "6" , "7",  "8",  "9" )
new.sample.ids <- c("C0" , "C1" , "C10" ,"C11", "proliferating_myeloid", "Tcell_contam_myeloid", "C14", "fib_contam_myeloid" , "C3" , "C4",  "C5",  "B_cell_contam_myeloid",  "Peri_contam_myeloid",  "Tcell_contam_myeloid" , "fib_contam_myeloid")

Myeloid_xenium@meta.data[["named_sub"]] <- plyr::mapvalues(x = Myeloid_xenium@meta.data[["named_sub"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(Myeloid_xenium, group.by = "named_sub", raster=FALSE, label=T)+NoLegend()
Idents(Myeloid_xenium) <- 'named_sub'
myeloid_xenium_markers <- FindAllMarkers(Myeloid_xenium, only.pos = T)
myeloid_xenium_markers_f <- myeloid_xenium_markers %>% filter(p_val_adj < 0.05)
myeloid_markers <- myeloid_markers %>% filter(p_val_adj < 0.05)


myeloid_MAPJAG_marker_genes_f_xenium_genes <- myeloid_markers[myeloid_markers$gene %in% rownames(Myeloid_xenium),]


xenium_markers <- list()
xenium_clusters_overlap <- list()
df_list <- list()

for (i in 1:length(unique(myeloid_xenium_markers_f$cluster))){
xenium_markers[[i]] <- myeloid_xenium_markers_f %>% filter(cluster == unique(myeloid_xenium_markers_f$cluster)[[i]])
xenium_clusters_overlap[[i]] <- myeloid_MAPJAG_marker_genes_f_xenium_genes[myeloid_MAPJAG_marker_genes_f_xenium_genes$gene %in% xenium_markers[[i]]$gene,]

df <- table(xenium_clusters_overlap[[i]]$cluster) %>% as.data.frame()
df2 <- table(myeloid_MAPJAG_marker_genes_f_xenium_genes$cluster) %>% as.data.frame()
df <- merge(x = df, y = df2, by = "Var1", all = TRUE)
df$pct <- df$Freq.x/df$Freq.y*100
df_list[[i]] <- df %>% replace(is.na(.), 0)
df_list[[i]]$cluster <- unique(myeloid_xenium_markers_f$cluster)[[i]]
print(ggplot(df_list[[i]], aes(x=Freq.x, y=pct)) + 
  geom_point()+theme_ArchR()+xlim(0,17)+ geom_text_repel(aes(label = df$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_list[[i]]$Freq.y)/2, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = max(df_list[[i]]$Freq.x)/2, linetype="dotted", 
                color = "red", size=0.5)+ggtitle(unique(myeloid_xenium_markers_f$cluster)[[i]]))



}

names(xenium_clusters_overlap) <- unique(myeloid_xenium_markers_f$cluster)


#rename and identify contams
Myeloid_xenium$named_sub_new <- Myeloid_xenium@meta.data[["named_sub"]]
Idents(Myeloid_xenium) <- 'named_sub_new'
levels(Myeloid_xenium)
current.sample.ids <- c( "C0" , "C1" , "C10" ,"C11", "proliferating_myeloid", "Tcell_contam_myeloid", "C14", "fib_contam_myeloid" , "C3" , "C4",  "C5",  "B_cell_contam_myeloid",  "Peri_contam_myeloid",  "Tcell_contam_myeloid" , "fib_contam_myeloid")
new.sample.ids <- c("MertK+ macrophages" , "MertK+ macrophages" , "CD1C+ cDC2s" ,"LAMP3 DCs", "proliferating_myeloid", "Tcell_contam_myeloid", "CD1C+ cDC2s", "fib_contam_myeloid" , "SPP1+ Macrophages" , "S1008A+ Monocytes",  "SPP1+ Macrophages",  "B_cell_contam_myeloid",  "Peri_contam_myeloid",  "Tcell_contam_myeloid" , "fib_contam_myeloid")

Myeloid_xenium@meta.data[["named_sub_new"]] <- plyr::mapvalues(x = Myeloid_xenium@meta.data[["named_sub_new"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(Myeloid_xenium, group.by = "named_sub_new", raster=FALSE, label=T)+NoLegend()


DimPlot(tcells, group.by = "anno2")

Idents(tcells) <- 'anno2'
t_cell_markers <- FindAllMarkers(tcells, only.pos = T)

t_cell_markers_xenium_markers_f <- t_cell_markers %>% filter(p_val_adj < 0.05)
Tcell_markers <- Tcell_markers %>% filter(p_val_adj < 0.05)


Tcell_MAPJAG_marker_genes_f_xenium_genes <- Tcell_markers[Tcell_markers$gene %in% rownames(T_NK_cells_xenium),]


xenium_markers <- list()
xenium_clusters_overlap <- list()
df_list <- list()

for (i in 1:length(unique(t_cell_markers_xenium_markers_f$cluster))){
xenium_markers[[i]] <- t_cell_markers_xenium_markers_f %>% filter(cluster == unique(t_cell_markers_xenium_markers_f$cluster)[[i]])
xenium_clusters_overlap[[i]] <- Tcell_MAPJAG_marker_genes_f_xenium_genes[Tcell_MAPJAG_marker_genes_f_xenium_genes$gene %in% xenium_markers[[i]]$gene,]

df <- table(xenium_clusters_overlap[[i]]$cluster) %>% as.data.frame()
df2 <- table(Tcell_MAPJAG_marker_genes_f_xenium_genes$cluster) %>% as.data.frame()
df <- merge(x = df, y = df2, by = "Var1", all = TRUE)
df$pct <- df$Freq.x/df$Freq.y*100
df_list[[i]] <- df %>% replace(is.na(.), 0)
df_list[[i]]$cluster <- unique(t_cell_markers_xenium_markers_f$cluster)[[i]]
print(ggplot(df_list[[i]], aes(x=Freq.x, y=pct)) + 
  geom_point()+theme_ArchR()+xlim(0,17)+ geom_text_repel(aes(label = df$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_list[[i]]$Freq.y)/2, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = max(df_list[[i]]$Freq.x)/2, linetype="dotted", 
                color = "red", size=0.5)+ggtitle(unique(t_cell_markers_xenium_markers_f$cluster)[[i]]))



}

names(xenium_clusters_overlap) <- unique(t_cell_markers_xenium_markers_f$cluster)

#rename and identify contams
tcells$named_sub_new <- tcells@meta.data[["anno2"]]
Idents(tcells) <- 'named_sub_new'
levels(tcells)
current.sample.ids <- c( "T cells/B/Myeloid" , "T cell/Fibroblast",  "CD4+ KLRB1+ T"  ,    "GZMK+ CD8+ T"  ,     "Vascular/T cell" ,   "NK cells/ILCs"   ,  
  "FOXP3+ Tregs"   ,    "T cell/Plasma cell" ,"Cycling T"    ,      "T cells/LL"     )
new.sample.ids <- c("CD4+ naive/central memory T" , "T cell/Fibroblast",  "CD4+ KLRB1+ T"  ,    "GZMK+ CD8+ T"  ,     "Vascular/T cell" ,   "NK cells/ILCs"   ,    "FOXP3+ Tregs"   ,    "T cell/Plasma cell" ,"Cycling T", "T cells/LL" )

tcells@meta.data[["named_sub_new"]] <- plyr::mapvalues(x = tcells@meta.data[["named_sub_new"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(tcells, group.by = "named_sub_new", raster=FALSE, label=T)+NoLegend()

DimPlot(bcells, group.by = "anno2")

Idents(bcells) <- 'anno2'
b_cell_markers <- FindAllMarkers(bcells, only.pos = T)


b_cell_markers_xenium_f <- b_cell_markers %>% filter(p_val_adj < 0.05)
B_plasma_markers <- B_plasma_markers %>% filter(p_val_adj < 0.05)


Bcell_MAPJAG_marker_genes_f_xenium_genes <- B_plasma_markers[B_plasma_markers$gene %in% rownames(bcells),]


xenium_markers <- list()
xenium_clusters_overlap <- list()
df_list <- list()

for (i in 1:length(unique(b_cell_markers_xenium_f$cluster))){
xenium_markers[[i]] <- b_cell_markers_xenium_f %>% filter(cluster == unique(b_cell_markers_xenium_f$cluster)[[i]])
xenium_clusters_overlap[[i]] <- Bcell_MAPJAG_marker_genes_f_xenium_genes[Bcell_MAPJAG_marker_genes_f_xenium_genes$gene %in% xenium_markers[[i]]$gene,]

df <- table(xenium_clusters_overlap[[i]]$cluster) %>% as.data.frame()
df2 <- table(Bcell_MAPJAG_marker_genes_f_xenium_genes$cluster) %>% as.data.frame()
df <- merge(x = df, y = df2, by = "Var1", all = TRUE)
df$pct <- df$Freq.x/df$Freq.y*100
df_list[[i]] <- df %>% replace(is.na(.), 0)
df_list[[i]]$cluster <- unique(b_cell_markers_xenium_f$cluster)[[i]]
print(ggplot(df_list[[i]], aes(x=Freq.x, y=pct)) + 
  geom_point()+theme_ArchR()+xlim(0,17)+ geom_text_repel(aes(label = df$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_list[[i]]$Freq.y)/2, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = max(df_list[[i]]$Freq.x)/2, linetype="dotted", 
                color = "red", size=0.5)+ggtitle(unique(b_cell_markers_xenium_f$cluster)[[i]]))



}



#rename and identify contams
bcells$named_sub_new <- bcells@meta.data[["anno2"]]
Idents(bcells) <- 'named_sub_new'
levels(bcells)
current.sample.ids <- c("Plasma cells"    ,         "Plasma cells/Fibroblasts", "Plasma cells/Myeloid"   ,  "B cells"          ,        "Vascular plasma, cells"   , "Plasma cell/Granulocyte",  "Cycling plasma cells"  ,   "pDCs"   )
new.sample.ids <- c("Plasma cells"    ,         "Plasma cells/Fibroblasts", "Plasma cells/Myeloid"   ,  "Memory B cells"          ,        "Vascular plasma, cells"   , "Plasma cell/Granulocyte",  "Cycling plasma cells"  ,   "pDCs"  )

bcells@meta.data[["named_sub_new"]] <- plyr::mapvalues(x = bcells@meta.data[["named_sub_new"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(bcells, group.by = "named_sub_new", raster=FALSE, label=T)+NoLegend()





bcells_meta <- bcells@meta.data
tcells_meta <- tcells@meta.data
Myeloid_xenium_meta <- Myeloid_xenium@meta.data

write.csv(bcells_meta, "/rds/projects/c/croftap-mapjagb6/xenium2/bcells_meta.csv")
write.csv(tcells_meta, "/rds/projects/c/croftap-mapjagb6/xenium2/tcells_meta.csv")
write.csv(Myeloid_xenium_meta, "/rds/projects/c/croftap-mapjagb6/xenium2/myeloid_meta.csv")




```{r}

DimPlot(xenium_chrissy, group.by="global13", raster=FALSE)



```
```{r}

Idents(xenium_chrissy) <- 'global13'
fibs_X <- subset(xenium_chrissy, idents=c("Fibroblasts"))

fibs_X <- fibs_X %>% 
    ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:40)

library(harmony)
fibs_X <- RunHarmony(fibs_X, group.by.vars = c("orig.ident"))
fibs_X <- RunUMAP(fibs_X, reduction="harmony",  dims=1:20)
fibs_X <- FindNeighbors(fibs_X, reduction="harmony", dims=1:20)

DimPlot(fibs_X, group.by="orig.ident", raster=FALSE)+NoLegend()

fibs_X <- FindClusters(fibs_X, resolution = c(0.01, 0.05, 0.1, 0.2))

fibs_X <- FindClusters(fibs_X, resolution = c(0.3, 0.4))
fibs_X <- FindClusters(fibs_X, resolution = c(0.5, 0.6))


DimPlot(fibs_X, group.by = "RNA_snn_res.0.6", raster=FALSE)

Idents(fibs_X) <- 'RNA_snn_res.0.6'
fibs_0.6_markers <- FindAllMarkers(fibs_X, only.pos = T)

DimPlot(fibs_X, group.by = "RNA_snn_res.0.6", raster=FALSE, label = T)


```

```{r}
FeaturePlot(fibs_X, features="CD14", raster=F)
FeaturePlot(fibs_X, features="MARCO", raster=F)

FeaturePlot(fibs_X, features="PDGFRA", raster=F)

FeaturePlot(fibs_X, features="THY1", raster=F)

FeaturePlot(fibs_X, features="PECAM1", raster=F)
FeaturePlot(fibs_X, features="ACTA2", raster=F)

DimPlot(fibs_X, group.by = "RNA_snn_res.0.6", raster=FALSE, label = T)

```



```{r}
fibs_X$named_sub <- fibs_X@meta.data[["RNA_snn_res.0.6"]]
Idents(fibs_X) <- 'named_sub'
levels(fibs_X)
current.sample.ids <- c( "0" , "1" , "10" ,"11", "12", "2" , "3" , "4",  "5",  "6",  "7",  "8" , "9")
new.sample.ids <- c("IFD1_COL5A2" , "LL",  "EC_contam", "TNC1_FBN1" , "macs" , "LL" , "macs", "pericyte" , "endo?" , "MFAP5_CD34",  "mac/LL",  "endo_contam", "THY1_MFAP5")

fibs_X@meta.data[["named_sub"]] <- plyr::mapvalues(x = fibs_X@meta.data[["named_sub"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(fibs_X, group.by = "named_sub", raster=FALSE, label=T)+NoLegend()





load("/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/analysis.RData")
rm(list=ls()[! ls() %in% c("stroma_clean_h")])

Idents(stroma_clean_h)<-"clusters"
levels(stroma_clean_h)
fibs<-subset(stroma_clean_h, idents=c("CD34/MFAP5"   ,  "POSTN"     ,     "SFRP/CXCL12" ,        "LL" ,          
 "CXCL14"  ,       "SOX5_CDH11" ,    "MMPs"   ))

fibs <- fibs %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()

rm(stroma_clean_h)


#filter markers genes for avg_log2FC >0.5 adjP val < 0.05


markers_named <- read_csv("/rds/projects/c/croftap-mapjagdata/MAPJAGv2/Chris/markers_named.csv")

fibs_MAPJAG_marker_genes_f <- markers_named %>% filter(cluster %in% c("CD34/MFAP5"   ,  "POSTN"     ,     "SFRP/CXCL12" ,        "LL" ,          
 "CXCL14"  ,       "SOX5_CDH11" ,    "MMPs"   ) & avg_log2FC >0.5 & p_val_adj < 0.05)

hist(fibs_MAPJAG_marker_genes_f$avg_log2FC)
min(fibs_MAPJAG_marker_genes_f$avg_log2FC)
table(fibs_MAPJAG_marker_genes_f$cluster)
max(fibs_MAPJAG_marker_genes_f$p_val_adj)


#subset xenium to just fibs (remove contaminating clusters)

Idents(fibs_X) <- 'named_sub'
fibs_clean_X <- subset(fibs_X, idents=c("IFD1_COL5A2" , "PDGFRA",  "TNC1_FBN1" ,  "LL" ,  "MFAP5_CD34",  "THY1_MFAP5"))


Idents(fibs_clean_X) <- 'named_sub'
fibs_clean_X_markers <- FindAllMarkers(fibs_clean_X, only.pos = T)
fibs_clean_X_markers <- fibs_clean_X_markers %>% filter(p_val_adj < 0.05)
fibs_MAPJAG_marker_genes_f_xenium_genes <- fibs_MAPJAG_marker_genes_f[fibs_MAPJAG_marker_genes_f$gene %in% rownames(fibs_clean_X),]


fibs_clean_X_markers$cluster %>% unique()


xenium_markers <- list()
xenium_clusters_overlap <- list()
df_list <- list()

for (i in 1:length(unique(fibs_clean_X_markers$cluster))){
xenium_markers[[i]] <- fibs_clean_X_markers %>% filter(cluster == unique(fibs_clean_X_markers$cluster)[[i]])
xenium_clusters_overlap[[i]] <- fibs_MAPJAG_marker_genes_f_xenium_genes[fibs_MAPJAG_marker_genes_f_xenium_genes$gene %in% xenium_markers[[i]]$gene,]

df <- table(xenium_clusters_overlap[[i]]$cluster) %>% as.data.frame()
df2 <- table(fibs_MAPJAG_marker_genes_f_xenium_genes$cluster) %>% as.data.frame()
df <- merge(x = df, y = df2, by = "Var1", all = TRUE)
df$pct <- df$Freq.x/df$Freq.y*100
df_list[[i]] <- df %>% replace(is.na(.), 0)
df_list[[i]]$cluster <- unique(fibs_clean_X_markers$cluster)[[i]]
print(ggplot(df_list[[i]], aes(x=Freq.x, y=pct)) + 
  geom_point()+theme_ArchR()+xlim(0,17)+ geom_text_repel(aes(label = df$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_list[[i]]$Freq.y)/2, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = max(df_list[[i]]$Freq.x)/2, linetype="dotted", 
                color = "red", size=0.5)+ggtitle(unique(fibs_clean_X_markers$cluster)[[i]]))



}

names(xenium_clusters_overlap) <- unique(fibs_clean_X_markers$cluster)

library(data.table)
total_elements <- rbindlist(xenium_clusters_overlap)

fibs_MAPJAG_marker_genes_f_xenium_genes <- fibs_MAPJAG_marker_genes_f_xenium_genes %>% as.data.frame()
rownames(fibs_MAPJAG_marker_genes_f_xenium_genes) <- fibs_MAPJAG_marker_genes_f_xenium_genes$...1

sc_marker_genes_list <- list()
for (i in 1:length(unique(fibs_MAPJAG_marker_genes_f_xenium_genes$cluster))){
  sc_marker_genes_list[[i]] <- fibs_MAPJAG_marker_genes_f_xenium_genes %>% filter(cluster == unique(fibs_MAPJAG_marker_genes_f_xenium_genes$cluster)[[i]]) %>% rownames()
  }
names(sc_marker_genes_list) <- unique(fibs_MAPJAG_marker_genes_f_xenium_genes$cluster)

par(mar=c(0.2,0.2,0.2,0.2))
overlap_list <- list()
for (i in 1:length(unique(fibs_clean_X_markers$cluster))){
  overlap_list[["list"]] <- xenium_clusters_overlap[[6]]$gene
    names(overlap_list) <- unique(xenium_clusters_overlap$cluster)[[6]]
  overlap_list_comb <- c(overlap_list, sc_marker_genes_list)
  gom.obj <- newGOM(overlap_list_comb, genome.size=nrow(total_elements))
drawHeatmap(gom.obj, adj.p=TRUE, cutoff=1, # show all.
 grid.col="Blues", note.col="black")+ theme(axis.text = element_text(size = 2))
   }


library(data.table)
all_dfs <- rbindlist(df_list)

all_dfs <- all_dfs[,-c(2,3)] %>% pivot_wider(names_from = cluster, values_from = pct) %>% as.data.frame()
rownames(all_dfs) <- all_dfs$Var1
all_dfs$Var1 <- NULL
all_dfs %>% as.matrix() %>% Heatmap( col=colorRamp2(c(0, 60), c("white", "red")))


ggplot(all_dfs, aes(Var1, cluster, fill = pct)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


DimPlot(fibs_clean_X, raster=F, label=T)

fibs_clean_X$named_sub_new <- fibs_clean_X@meta.data[["named_sub"]]
Idents(fibs_clean_X) <- 'named_sub_new'
levels(fibs_clean_X)
current.sample.ids <- c("IFD1_COL5A2", "PDGFRA"  ,    "TNC1_FBN1" ,  "LL"    ,      "MFAP5_CD34" , "THY1_MFAP5" )
new.sample.ids <- c("CD34/MFAP5", "POSTN"  ,    "SOX5/CDH11" ,  "LL"    ,      "SFRP/CXCL12" , "CD34/MFAP5")

fibs_clean_X@meta.data[["named_sub"]] <- plyr::mapvalues(x = fibs_clean_X@meta.data[["named_sub"]], from = current.sample.ids, to = new.sample.ids)

DimPlot(fibs_clean_X, group.by = "named_sub", raster=FALSE, label=T)+NoLegend()
DimPlot(fibs_clean_X, group.by = "named_sub", raster=FALSE, label=T)+NoLegend()
DimPlot(fibs_X, group.by = "named_sub", raster=FALSE, label=T)+NoLegend()

DimPlot(fibs_X, raster=F, label=T, group.by = "named_sub")

fibs_X$named_sub_new <- fibs_X@meta.data[["named_sub"]]
Idents(fibs_X) <- 'named_sub_new'
levels(fibs_X)
current.sample.ids <- c("IFD1_COL5A2",  "PDGFRA"   ,    "Tcell_contam", "TNC1_FBN1"   , "macs"   ,      "LL"     ,      "pericyte" ,    "endo?"    ,    "MFAP5_CD34"  , "mac/LL"  ,     "endo_contam",  "THY1_MFAP5" )


new.sample.ids <- c("CD34/MFAP5",  "POSTN"   ,    "Tcell_contam", "SOX5/CDH11"   , "macs"   ,      "LL"     ,      "pericyte" ,    "EC"    ,    "SFRP/CXCL12"  , "mac/LL"  ,     "EC",  "CD34/MFAP5")

fibs_X@meta.data[["named_sub_new"]] <- plyr::mapvalues(x = fibs_X@meta.data[["named_sub_new"]], from = current.sample.ids, to = new.sample.ids)


DimPlot(fibs_X, group.by = "named_sub_new", raster=FALSE, label=T)+NoLegend()


fibs_X_meta <- fibs_X@meta.data

rm(list=ls()[! ls() %in% c("fibs_X_meta")])













