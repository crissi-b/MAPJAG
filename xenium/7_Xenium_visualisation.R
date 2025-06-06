remotes::install_github("korsunskylab/spatula")

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(tidyverse)
library(sf)
library(ArchR)
library(cowplot)
library(tidyverse)
library(pheatmap)
library(grid)
library(data.table)
library(devtools)
library(spatula)
library(furrr)
library(knitr)
library(magrittr)
library(stringr)
library(ggpubr)
library(rstatix)
library(DescTools)
library(textshape)
library(ComplexHeatmap)
library(splitstackshape)


load("/rds/projects/c/croftap-mapjagx1/MAPJAG-Xenium/xenium_obj.rds")

#################################################################
#################################################################

# Set graph graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())
thin <-   list(theme(plot.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y = element_blank()))


# Set order for dotplot

iden1 <- as.character(unique(xen$named2407))
i1 <- grep("FOXP3", iden1)
i2 <- grep("central", iden1)
i3 <- c("CD4+ KLRB1+ T cells", "GZMK+ CD8+ T cells", "NK cells/ILCs", "Cycling T cells", "T cells", "Lining layer T cells")
i4 <- grep("Memory B|lasma cells", iden1)
i5 <- grep("pDC", iden1)
i6 <- grep("CD1C|LAMP3", iden1)
i7 <- grep("Mer|SPP1", iden1)
i7a <- grep("S100|Myel|Cycling my", iden1)
i8 <- "Granulocytes"
i9 <- grep("LL f|SOX5", iden1)
i9a <- grep("POSTN|CD34|CXCL12", iden1)
i9b <- grep("ibroblast", iden1)
i10 <- grep("Peri|Endothel|Lymphat", iden1)
i11 <- "Adipose"

xen_order <- unique(c(iden1[i1],iden1[i2],i3,iden1[i4],iden1[i5],iden1[i6],iden1[i7],
               iden1[i7a],i8,iden1[i9],iden1[i9a],iden1[i9b],iden1[i10],i11))
#xen_order

levels(xen) <- rev(xen_order)
genes <- c("CD2", "CTLA4", "CD4",  "IL7R", "KLRB1", "CD69", 
            "CD8A","GZMK", "PRF1", "CCL5",  "NKG7","GNLY", "KLRD1", "MKI67", 
            "PRDM1","CD79A", "MS4A1",  "PLD4", "PTGDS", "IRF8","CLEC10A", "FCER1A", "CD1C", "CXCL9",
           "MARCO", "VSIG4","LYVE1", "MRC1",  "FCN1", "S100A12",
           "CPA3", "CTSG", "SLC18A2", 
           "PRG4","CFB", "THBS2", "PTN",
           "SFRP2", "FBLN1", "THY1", "MFAP5", "CD34", "FGFBP2",
           "ACTA2","EGFL7", "ANGPT2", "PLIN4", "ADIPOQ")
DotPlot(xen, features=genes) & formplot

#################################################################

# Visualise clusters

colcol <- c("#009E73", "#9A6324", "#FFE119",  "#42D4F4",  "#F58231", "#911EB4",
            "#BFEF45","#F032E6",  "#DBA206", "#3CB44B",  "#FABEBE","red", "#E6BEFF",
            "#8E0152", "#9ADE00", "#FFFAC8", "#A9A9A9", "black", "#FFD8B1", "#000075",
            "#808000", "#DCBEFF",  "#469990", "#0072B2", "#D55E00","#CC79A7","#4363D8",
             "#56B4E9","#E6194B", "#F0E442", "#D3D3D3", "#2B8CBE")

DimPlot(xen, group.by="named2407", label=T, repel=T) & scale_color_manual(values=colcol) & thin


#################################################################
#################################################################

# Visualise niches

load("/rds/projects/c/croftap-mapjagx1/analysis/niches_final/analysis_chrissy/analysis.RData")

cols <- c("deeppink", "navy", "yellow", "red", "lavender", "skyblue")
geoms_f_f[["D2518_36_Z3"]] %>% as.data.frame %>% 
  ggplot()+
  geom_sf(aes(geometry = geometry, fill = niches_new), alpha = 0.7,
          color = "black")+ scale_fill_manual(values = cols) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#################################################################
#################################################################


# Visualise main cell types within niches as heatmaps

heatcol <- rev(brewer.pal(8, "RdBu"))

# Make dataframe of cell types
named_clusters <- xen[["named2407"]]

# Add to table of niche identity per cell
rownames(named_clusters) -> named_clusters$cell
pt_table2 <- pt_table[,-1] 
pt_table2 <- pt_table2 %>% left_join(named_clusters, by="cell")
pt_table <- cbind(pt_table, pt_table2$named2407)  
colnames(pt_table)[10] <- "named2407"

# Tidy names up
pt_table <- pt_table %>% mutate(niches_named = recode(niches_new, "LL"="Lining layer", "Myelo_Tcell"="Myeloid / Lymphocyte", "Plasma_cell"="Plasma cell-rich", "SL"="Sublining stroma", "Adipose"="Adipose-rich", "Myelo_vascular"="Perivascular"))
# Simplify clusters
pt_table <- pt_table %>% mutate(global_mid = recode(global_mid, "LAMP3+ DCs"="Dendritic cells", "CD1C+ cDC2s"="Dendritic cells"))

predictions <- table(pt_table$niches_named, pt_table$global_mid)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions2<-predictions %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  as.data.frame() 
predictions2 <-predictions2 %>% remove_rownames %>% column_to_rownames("Var1")
predictions2 <- predictions2[c("Lining layer","Myeloid / Lymphocyte","Plasma cell-rich","Sublining stroma", "Perivascular", "Vascular", "Adipose-rich"),c("Fibroblasts", "Myeloid cells", "Dendritic cells", "CD8+ T cells", "CD4+ T cells", "NK cells/ILCs",  "pDCs","B cells",  "Plasma cells", "Granulocytes", "Pericytes","Lymphatics", "Endothelial cells",  "Adipose")]

predictions2 %>% scale() %>% 
  Heatmap(cluster_columns= F, cluster_rows = F, colorRamp2(c(-max(predictions2), 0, max(predictions2)), c("lightblue", "white", "darkred")), border=T)

#################################################################

# Visualise finer cell types within niches as heatmaps

predictions <- table(pt_table$niches_named, pt_table$named2407)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
predictions2<-predictions %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  as.data.frame() 
predictions2 <-predictions2 %>% remove_rownames %>% column_to_rownames("Var1")
predictions2 <- predictions2[c("Lining layer","Myeloid / Lymphocyte","Plasma cell-rich","Sublining stroma","Perivascular", "Vascular",  "Adipose-rich"),]

# Please note dendrogram is rotated about it’s branch points in final visualisation for clearer interpretation
predictions2 %>% scale() %>% 
  Heatmap( cluster_rows = F, colorRamp2(c(-max(predictions2), 0, max(predictions2)), c("lightblue", "white", "darkred")), cluster_columns= T, border=T)

#################################################################
#################################################################

# Visualise cells in tissue of niches

# Lining layer niche

region <- "D2518_36_Z3"
xmin <- 1220
xmax <- 1500
ymin <- 3250
ymax <- 3506
geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

# Lining layer cell colours
cellcols <- c("MerTK+ macrophages"= "navy","Cycling fibroblasts"= "orange", "LL fibroblasts" ="dodgerblue",
              "SPP1+ macrophages"=  "deeppink", "S1008A+ monocytes"= "#A6CEE3","Cycling myeloid"= "yellow",
              "LAMP3+ DCs"= "aquamarine", "Lining layer T cells" ="mediumpurple1", "CD1C+ cDC2s"="chartreuse",
             "Other"="grey87", "Pericytes"= "seashell", "Other"="#D0D1CF", "Endothelial cells"="seashell", "Lymphatics"="seashell")

filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") +
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(
    axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 

#################################################################

# Myeloid-lymphocyte niche

 Myeloid T cells
region <- "D2519_18_R16"
xmin <- 9080
xmax <- 9500
ymin <- 6850
ymax <- 7200
geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

cellcols <- c("Lining layer T cells"="lavender","CD1C+ cDC2s"="black", "NK cells/ILCs"="deeppink", 
              "FOXP3+ Tregs"="forestgreen", 
              "CD4+ naive/central memory T"="blue","GZMK+ CD8+ T cells"=  "green2","CD4+ KLRB1+ T cells"= "cyan",
              "Cycling T cells"= "orange", "T cells" = "aquamarine2", "Memory B cells" ="lightpink",
              "pDCs"="indianred3", "Myeloid cells"= "yellow",
              "Other"="#C1C1C1", "Pericytes"= "grey97",  "Endothelial cells"="grey97", "Lymphatics"="grey97",
              "Lining layer T cells"="darkgrey", "SPP1+ macrophages"="darkgrey",
              "LL fibroblasts"="darkgrey", "Cycling fibroblasts"= "darkgrey","LAMP3+ DCs"= "darkgrey",
              "Cycling myeloid"= "darkgrey")

filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") +
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(
    axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 

#################################################################

# Vascular niche

region <- "D2518_36_Z3"
xmin <- 1220
xmax <- 1650
ymin <- 3130
ymax <- 3506
geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

cellcols <- c("Lymphatics"= "brown","Endothelial cells" ="lightgreen", "Pericytes" ="blue", "Other"="#D0D1CF",              "Lining layer T cells"="#505050", "SPP1+ macrophages"="#505050",
              "LL fibroblasts"="#505050", "Cycling fibroblasts"= "#505050","LAMP3+ DCs"= "#505050","Cycling myeloid"= "#505050")

filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") +
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 

#################################################################

# Stromal sublining niche- NB image in figure is rotated to match lining layer orientation of other niches

region <- "D2518_24_J11"
xmin <- 2450
xmax <- 3000
ymin <- 6550
ymax <- 7100
geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

cellcols <- c("CXCL12+ fibroblasts" = "blue2", "CD34+ fibroblasts"= "forestgreen","SOX5+ CDH11+ fibroblasts" ="chartreuse",
              "POSTN+ fibroblasts" ="yellow", "Granulocytes"= "cyan", "Adipose"="orange",
              "Fibroblasts"= "mediumorchid","Other"="#D0D1CF", "Pericytes"= "lightblue",  "Endothelial cells"="lavender", 
              "Lining layer T cells"="khaki", "SPP1+ macrophages"="khaki",
              "LL fibroblasts"="khaki", "Cycling fibroblasts"= "khaki","LAMP3+ DCs"= "khaki",
              "Cycling myeloid"= "khaki")
              
filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") + scale_x_reverse() + 
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax))

#################################################################

# Plasma cell niche

region <- "D2519_18_R16"
xmin <- 9080
xmax <- 9500
ymin <- 6850
ymax <- 7200
geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

cellcols <- c("Myeloid cells"= "yellow", "Plasma cells"="deeppink", "Cycling plasma cells"="blue",
              "Granulocytes"="orange", 
              "Other"="grey97", "Pericytes"= "grey80",  "Endothelial cells"="grey80", "Lymphatics"="grey80",
              "Lining layer T cells"="darkgrey", "SPP1+ macrophages"="darkgrey",
              "LL fibroblasts"="darkgrey", "Cycling fibroblasts"= "darkgrey","LAMP3+ DCs"= "darkgrey",
              "Cycling myeloid"= "darkgrey")             

filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") +
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(
    axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 

#################################################################

# Adipose-rich niche

region <- "D2519_15_Z6"
xmin <- 2400
xmax <- 3100
ymin <- 11630
ymax <- 12200
geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

cellcols <- c("Adipose"="green", "Granulocytes"="orange", "CXCL12+ fibroblasts"="blue", "POSTN+ fibroblasts"="yellow",
              "Other"="#C1C1C1", "Pericytes"= "grey80",  "Endothelial cells"="grey80", "Lymphatics"="grey80",
              "Lining layer T cells"="darkgrey", "SPP1+ macrophages"="darkgrey",
              "LL fibroblasts"="darkgrey", "Cycling fibroblasts"= "darkgrey","LAMP3+ DCs"= "darkgrey",
              "Cycling myeloid"= "darkgrey")

filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") +
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(
    axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 


#################################################################
#################################################################

# Check correlations between Krenn scores and infiltrate scores of each patient
# Krenn scores added in 6_xenium_annotation script

# Ensure xenium seurat object has niches annotated per cell type
niche_tb <- pt_table[,c("niches_named", "cell")] %>% as.data.frame()
rownames(niche_tb) <- as.character(niche_tb$cell)
niche_tb$cell <- NULL
head(niche_tb)
xen <- AddMetaData(xen, niche_tb, col.name="niches_named")

# Examine correlations
Idents(xen) <- "niches_named"
ptnum2 <- as.data.frame(table(xen$patient, xen@active.ident))
ptnum2 <- group_by(ptnum2, Var1) %>% mutate(percent = Freq*100/sum(Freq))

# Replace xen$infiltrate with xen$krenn for the total Krenn score 
age <- as.data.frame(table(xen$patient, xen$infiltrate))  %>% distinct()
age <- age[age$Freq > 0,1:2]
ptnum2 <- merge(ptnum2, age, by="Var1")

summary <- list()
image.list <- list()
summary2 <- data.frame()

iden1 <- unique(ptnum2$Var2.x)

# Check correlation between Krenn infiltrate scores and proportion of niche
for (i in 1:length(iden1)) {
  iden <- iden1[[i]]
  selec <- ptnum2[ptnum2$Var2.x %in% iden, ]
  selec$Var2.y <- as.numeric(as.character(selec$Var2.y))
  coro <- cor.test(selec$Var2.y, selec$percent, method="pearson")
  labo <- paste0("\nr=",round(coro$estimate, digits=2), "\n", "P =", round(coro$p.value, digits=3))
  image.list[[i]] <- ggplot(selec, aes(x=Var2.y, y=percent)) + geom_point() + labs(x="Infiltrate score", y="Tissue cells (%)", title=paste0(iden, labo)) + geom_smooth(method="lm", se=FALSE) + theme(panel.border = element_rect(colour = "black", fill = NA),panel.grid = element_blank(),panel.background = element_blank())
  summary[i] <- paste0(iden1[i],"  ", labo)
  #summary[i] <- paste0(sig[i],"  ", labo)
  values <- as.data.frame(cbind(paste0(iden[1]), round(coro$estimate, digits=2),  round(coro$p.value, digits=3)))
  summary2 <- rbind(values, summary2)
}

grid.arrange(grobs = image.list[c(1,4,5,3,2,7,6)], ncol=4)  

#################################################################

# Repeat for individual cell types

Idents(xen) <- "named2407"

ptnum2 <- as.data.frame(table(xen$patient, xen@active.ident))
ptnum2 <- group_by(ptnum2, Var1) %>% mutate(percent = Freq*100/sum(Freq))

age <- as.data.frame(table(xen$patient, xen$infiltrate))  %>% distinct()
age <- age[age$Freq > 0,1:2]
ptnum2 <- merge(ptnum2, age, by="Var1")

summary <- list()
image.list <- list()
summary2 <- data.frame()

iden1 <- unique(ptnum2$Var2.x)

for (i in 1:length(iden1)) {
  iden <- iden1[[i]]
  selec <- ptnum2[ptnum2$Var2.x %in% iden, ]
  selec$Var2.y <- as.numeric(as.character(selec$Var2.y))
  coro <- cor.test(selec$Var2.y, selec$percent, method="pearson")
  labo <- paste0("\nr=",round(coro$estimate, digits=2), "\n", "P =", round(coro$p.value, digits=3))
  image.list[[i]] <- ggplot(selec, aes(x=Var2.y, y=percent)) + geom_point() + labs(x="Infiltrate score", y="Tissue cells (%)", title=paste0(iden, labo)) + geom_smooth(method="lm", se=FALSE) + theme(panel.border = element_rect(colour = "black", fill = NA),panel.grid = element_blank(),panel.background = element_blank())
  summary[i] <- paste0(iden1[i],"  ", labo)
  #summary[i] <- paste0(sig[i],"  ", labo)
  values <- as.data.frame(cbind(paste0(iden[1]), round(coro$estimate, digits=2),  round(coro$p.value, digits=3)))
  summary2 <- rbind(values, summary2)
}

# Plot cell types with a correlation > 0.5 or < -0.5
grid.arrange(grobs = image.list[c(26,23,17,31,18,29,14)], ncol=4)  
grid.arrange(grobs = image.list[c(5,12,7,6,24,28,1,2)], ncol=4)  


#################################################################
#################################################################


# Proximity analysis

geoms_f <- fread("/rds/projects/c/croftap-mapjagx1/MAPJAG-Xenium/x_y_cells.csv", select=c("x", "y", "cell", "donor"))

colnames(xen) -> xen$cell_id
nom <- unique(geoms_f$donor)
xen$donor <- stringr::str_extract(xen$orig.ident, "[^_]*_[^_]*")
seurat_orig_meta <- xen[[c("named2407", "donor")]]
rownames(seurat_orig_meta) -> seurat_orig_meta$cell_id
seurat_orig_meta <- split(seurat_orig_meta, seurat_orig_meta$donor)
ggsplit <- split(geoms_f, geoms_f$donor)

for (i in 1:length(ggsplit)){
  index <- match(ggsplit[[i]]$cell, seurat_orig_meta[[i]]$cell_id)
  ggsplit[[i]]$named <- seurat_orig_meta[[i]]$named2407[index]
}

ggsplit <- rbindlist(ggsplit)
ggsplit <- ggsplit[!is.na(ggsplit$named),]

#################################################################

# NN analysis functions
# It compares observed frequencies to a distribution of frequencies obtained through permutations,
# providing statistical measures of co-localization significance. The code is designed for a 
# specific analysis and likely fits into a broader context of understanding relationships within 
# a network represented by the adjacency matrix.

coloc_one_type = function(index_type, adj, y, nperm = 100, max_dist=30, compartments=NULL, verbose=TRUE) {
  if (verbose) message(index_type)
  #########    CREATE VARIABLES
  types = unique(y)
  # create new variable & assign values when values in the vector y are equal to the 
  # value of the index_type
  i_index = which(y == index_type)
  #creates set of indices not part of the i_index set (difference between the order of index from 1:length(y) and the i_index)
  i_shuffle = setdiff(seq_len(length(y)), i_index)
  ##########.   MATRIX MULTIPLICATION
  #selects row from the sparse matrix corresponding to the indices in i_index
  #. %*% performs matrix multiplication then 'sparse.model.matrix(~0+y)' converts y into 
  # a sparse matrix where categorical variables are assigned as binary outcomes for each
  # essentially, a matrix of value freuqnecies is generated
  X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+y) %>% as.matrix()
  colnames(X) = gsub('^y', '', colnames(X))
  #########.    FREQUENCY CALCULATION     --- frequency of each unique type in the modified matrix
  freq = (colSums(X) / nrow(X))[types]
  ###########.   PERMUTATION LOOP
  freq_perm = map(seq_len(nperm), function(i) {
    set.seed(i)
    yperm = y
    if (is.null(compartments)) {
      yperm[i_shuffle] = sample(y[i_shuffle])
    } else {
      ## shuffle inside compartments, to preserve total composition within compartment
      .x = split(i_shuffle, compartments[i_shuffle]) %>%
        map(function(.i) {
          ## CAUTION: if .i is a single number, sample will interpret it as 1:.i
          if (length(.i) == 1) {
            message('No shuffling is taking place, check code')
            res = .i
          } else {
            res = sample(.i) ## shuffle non-index cells inside hub
          }
          names(res) = .i
          return(res)
        }) %>%
        reduce(c)
      yperm[as.integer(names(.x))] <- y[.x]
    }
    X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+yperm) %>% as.matrix() #%>% prop.table(1)
    colnames(X) = gsub('^yperm', '', colnames(X))
    (colSums(X) / nrow(X))[types]
  }) %>%
    purrr::reduce(rbind2)
  ###########.   STORE STATS
  stats = tibble(
    type = types,
    freq,
    zscore = (freq - apply(freq_perm, 2, mean)) / apply(freq_perm, 2, sd),
    pval = exp(log(2) + (pnorm(-abs(zscore), log.p = TRUE, lower.tail = TRUE))), ## one-tailed
    fdr = p.adjust(pval)
  ) %>%
    cbind(dplyr::rename(data.frame(t(apply(freq_perm, 2, quantile, c(.025, .975)))), q025 = `X2.5.`, q975 = `X97.5.`)) %>% ## 95% CI
    subset(type != index_type) %>%
    dplyr::mutate(index_type = index_type) %>%
    dplyr::select(index_type, type, everything()) %>%
    arrange(fdr)
  return(stats)
}


# Loop through all cell types
coloc_all_types = function(index_types, coords, y, nperm = 100, nsteps=1, max_dist=30, compartments=NULL, parallel=TRUE, verbose=TRUE) {
  if (parallel & length(index_types) > 1) {
    plan(multicore)
  } else {
    plan(sequential)
  }
  ## Define neighbors
  ## NOTE: max_dist only refers to directly adjacent neighbors
  adj = spatula::getSpatialNeighbors(coords, return_weights = TRUE)
  adj@x[adj@x > max_dist] = 0
  adj = Matrix::drop0(adj)
  adj@x = rep(1, length(adj@x))
  ## If nsteps>1, consider not only your adjacent neighbors
  ##   but also your neighbor's neighbors etc.
  if (nsteps > 1) {
    adj = adj + Matrix::Diagonal(n = nrow(adj)) ## add self
    for (iter in seq_len(nsteps - 1)) {
      adj = adj %*% adj
    }
    ## Ignore weights. Only care if cell is a neighbor or not
    adj@x = rep(1, length(adj@x))
    ## Remove self as neighbor
    adj = adj - Matrix::Diagonal(n = nrow(adj))
    adj = Matrix::drop0(adj)
  }
  index_types %>%
    future_map(coloc_one_type, adj, y, nperm, max_dist, compartments, verbose, .options = furrr::furrr_options(seed = 1)) %>%
    rbindlist()  %>%
    identity
}


#################################################################

# Now run functions on code
coloc_res_coarse_global<- coloc_all_types(
  index_type = unique(ggsplit$named),
  coords = ggsplit[, c("x", "y")],
  y = ggsplit$named,
  compartments = NULL,
  max_dist = 40,
  nperm = 1000,
  parallel = TRUE
)

# Create proximity matrix
plt_df <- coloc_res_coarse_global %>%
  subset(pval < 0.05) %>%
  dplyr::select(index_type, type, zscore) %>%
  spread(type, zscore, fill = 0) %>%
  column_to_rownames('index_type') %>%
  as.matrix

# Plot proximity
plt_df %>%
  Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)

# Make proximity matrix symmetrical for easier visualisation
plt_df_use <- t(plt_df)+plt_df
plt_df_use %>%
  Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, 
          cluster_columns = T, cluster_rows = T) 


#################################################################
#################################################################


# S100A8 proximity analysis

look <- plt_df[, "S1008A+ monocytes"]
look <- look %>% as.data.frame()
colnames(look) <- "prox"
rownames(look) -> look$cells
top_5 <- look[order(-look$prox), ][1:5, ]

ggplot(top_5, aes(y = reorder(cells, prox), x = prox), fill=cells) +
  geom_bar(stat = "identity", fill=c("red", "blue", "green", "yellow", "purple"), alpha=0.6) + theme_minimal() 

# Visualise where S100A8+ monocytes are 

# Find slides with highest S100A8 cell numbers
count_data <- list()
# Loop through each table in the list
for (i in seq_along(geoms_f_f)) {
  count_data[[i]] <- table(geoms_f_f[[i]][["named2407"]])
  names(count_data)[i] <- names(geoms_f_f)[i]
}

result_df <- map_df(count_data, ~ as.data.frame(.), .id = "sample")
result_df <- dcast(result_df, sample ~ Var1)

# View results and order by S100A8 to find sample slides with highest number
check <- result_df[,grep("S100|LAMP3|MerTK|sample", colnames(result_df))]

# Visualise the top results, re-run the code below for each region value
region <- "D2519_18_S16"
region <- "D2519_18_J4"
region <- "D2519_18_R17"

cellcols <- c("MerTK+ macrophages"= "navy","Cycling fibroblasts"= "orange", "LL fibroblasts" ="dodgerblue",
              "SPP1+ macrophages"=  "deeppink", "S1008A+ monocytes"= "#A6CEE3","Cycling myeloid"= "yellow",
              "LAMP3+ DCs"= "aquamarine", "Lining layer T cells" ="mediumpurple1", "CD1C+ cDC2s"="chartreuse",
              "Other"="grey87", "Pericytes"= "red", "Endothelial cells"="red",
              "Lymphatics"="brown")

geoms_f_f[[region]] %>% as.data.frame %>% 
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7,
                     color = "black")+ scale_fill_manual(values = cellcols) +
  theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#################################################################
#################################################################


# SOX5 proximity analysis

# Extract highest colocalising values
# Convert the matrix into a data frame
df <- as.data.frame(plt_df)  # Replace 'matrix_values' with the name of your matrix

df1 <- df[order(df$`SOX5+ CDH11+ fibroblasts`, decreasing=T),] %>% head(6) 
df1 <- df1["SOX5+ CDH11+ fibroblasts"]
df1$Cell_Interactions <- rownames(df1)
colnames(df1)[1] <- "Proximity"
ggplot(df1, aes(x = Proximity, y = reorder(Cell_Interactions, Proximity), fill=Cell_Interactions)) +
  geom_bar(stat = "identity") +
  theme_minimal() + labs(y="", x="\nProximity score (Xenium)") & NoLegend()

#################################################################

# Visualisation of SOX5 high region
region <- "D2518_11_O7"

xmin <- 1480
xmax <- 1850
ymin <- 1600
ymax <- 1950

geoms_f_f_sf <- geoms_f_f[[region]] %>% as.data.frame() %>% st_as_sf()
bbox <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
filtered_geoms_f_f_sf <- geoms_f_f_sf %>%
  filter(st_intersects(geometry, bbox, sparse = FALSE))

cellcols <- c("CXCL12+ fibroblasts" = "cyan", "CD34+ fibroblasts"= "forestgreen","SOX5+ CDH11+ fibroblasts" ="blue2",
              "POSTN+ fibroblasts" ="yellow", "Granulocytes"= "orange", "LL fibroblasts"="deeppink", 
              "Cycling fibroblasts"="lavender",
              "Other"="#D0D1CF", "Pericytes"= "lightblue",  "Endothelial cells"="lightblue")

filtered_geoms_f_f_sf <- filtered_geoms_f_f_sf %>%
  mutate(fill_color = ifelse(named2407 %in% names(cellcols), named2407, "Other"))

# Plotting the filtered cells
filtered_geoms_f_f_sf %>%
  ggplot() + geom_sf(aes(geometry = geometry, fill = named2407), alpha = 0.7, color = "black") +
  scale_fill_manual(values = cellcols) + theme(panel.background = element_rect(fill = 'white')) +
  theme(
    axis.line = element_line(color = 'black'), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) 
