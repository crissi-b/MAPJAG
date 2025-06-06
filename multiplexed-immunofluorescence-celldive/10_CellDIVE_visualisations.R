remotes::install_github("korsunskylab/spatula")

options(bitmapType='cairo')
library(Seurat)
library(dplyr)
library(sctransform)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

library(devtools)
library(spatula)
library(furrr)

library(data.table)
library(splitstackshape)
library(ggpubr)
library(rstatix)
library(DescTools)
library(textshape)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)

# load CellDive objects

# Set graph graphics
formplot <- list(RotatedAxis(), FontSize("8"), theme(axis.title.x=element_blank(), axis.title.y=element_blank()), scale_color_gradient2())

####################################################################################

# Dotplot of markers per cell type

gen3 <- c("FABP4","COL4A1","SMA","CD146","DKK3","COMP","POSTN","COL3A1","COL6A1","COL1", "CLU", "SPARC", "MMP3","PDPN","LYVE1", 
          "CD68",  "CD3","CD20", "MCT","CD15", "DAPI-FINAL")

DefaultAssay(all) <- "SCT"
Idents(all) <- "take2"

iden <- unique(as.character(all$take2))
c1 <- iden[grep("Adip", iden)]
c2 <- iden[grep("SMA", iden)]
c4 <- iden[grep("146", iden)]
c5 <- c("Endothelial cells")
c6 <- iden[grep("ibrobl", iden)]
c7 <- iden[grep("CLU|LL F", iden)]
c8 <- iden[grep("Lymph", iden)]
c9 <- iden[grep("LYVE", iden)]
c10 <- iden[grep("LL M|CD68", iden)]
c11 <- "T cells"
c12 <- iden[grep("agg", iden)]
c13 <- iden[grep("Mast", iden)]
c14 <- iden[grep("Low-staining cells", iden)]
c15 <- iden[grep("Fibrin", iden)]
c16 <- c("Neutrophil-rich", "Neutrophil/T cell-rich")
c17 <- iden[grep("Artef|Fibrin", iden)]
c18 <- iden[grep("Sparse|RBCs", iden)]
idenCD <- unique(c(c18,c17,c14,c16,c15,c13,c12,c11,c10,c9,c8,c7,c6,c5,c4,c2,c1))

# Check all the cell types have been included
sort(idenCD) == sort(as.character(unique(all$take2)))

# Set the order for the dotplot
all$take2 <- factor(all$take2, levels=idenCD)
Idents(all) <- "take2"
levels(all) <- idenCD

Idents(all) <- "take2"
dotp <- subset(all, idents= c("Artefact","Low-staining cells"), invert=T)

# Create the order excluding artefact and low staining cells
c17 <- iden[grep("Fibrin|Sparse|RBC", iden)]
idenCD <- unique(c(c18,c17,c16,c15,c13,c12,c11,c10,c9,c8,c7,c6,c5,c4,c2,c1))

levels(dotp) <- idenCD
DotPlot(dotp, features =gen3) & formplot

####################################################################################

# Visualise niches on the from CellDive groups of cells

# Set colours of cells
group.color <- c("SMA-hi vessels"="steelblue1", "CD146-hi vessels"="green", 
                 "Endothelial cells"="brown","RBCs"="red",
                 "CLU+ Fibroblasts"="aquamarine3", "COL1-hi Fibroblasts"="lemonchiffon",  
                 "POSTN-hi Fibroblasts"="salmon","SL Fibroblasts"="palegoldenrod",
                 "LL Fibroblasts"= "magenta", "Sparse cells"= "darkseagreen1",
                 "COL6-hi Fibroblasts"="lavender", "COL3-hi Fibroblasts"="orchid4",
                 "LL Macrophages"="yellow",  "CD68+ Myeloid"="yellow3", 
                 "T cells" ="cyan", 
                 "B/T cell aggregates"=  "lavenderblush2",
                 "Adipose"="mediumspringgreen",  "Mast cells"="darkorange", 
                 "Lymphatics"="palevioletred", "LYVE1-hi Macrophages"="mediumpurple1", "Low-staining cells"="lavender",
                 "Artefact"="azure", "Fibrin Macrophages"="yellow4", "Fibrin"="skyblue2", "Blood/Fibrin"="burlywood",
                 "Neutrophil-rich"="plum2", "Neutrophil/T cell-rich"="steelblue", "Fibrin-2"="beige")

setv <- list(theme(panel.background = element_rect(fill = 'black'), 
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text = element_blank(),
                   legend.background = element_rect(fill = "black", color = NA),
                   legend.key = element_rect(color = "gray", fill = "black"),
                   legend.text = element_text(color = "white"))) 

# Visualise example
i <- 7 
nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2
iden <- unique(all$take2)

# Select small region
small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500,]

# All cells 
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

# All cells within the niches specified
small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500,]
iden1 <- iden[grep("LL F|LL M|CD68|T cells|agg|CLU|CD146|Lymph|CD146|SMA|Endo|COL3|COL6|POSTN|COL1|Adi|Sparse|Mast|LYVE", iden)]
small <- small[small$named1 %in% iden1,]
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

####################################################################################

# Visualise lining layer niche
small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500,]
iden1 <- iden[grep("LL F|LL M", iden)]
small <- small[small$named1 %in% iden1,]
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

# Visualise myeloid-lymphocyte niche
iden1 <- iden[grep("CD68|T cells|agg|CLU|CD146|Lymph", iden)]
small <- small[small$named1 %in% iden1,]
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

# Visualise vascular niche
iden1 <- iden[grep("CD146|SMA|Endo", iden)]
small <- small[small$named1 %in% iden1,]
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

# Visualise stromal niche
iden1 <- iden[grep("COL3|COL6|POSTN|COL1", iden)]
small <- small[small$named1 %in% iden1,]
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

# Visualise adipose niche
i <- 8
nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2
iden <- unique(all$take2)

small <- test[test$Centroid.X.µm < 1500  &  test$Centroid.Y.µm > 2000,]
iden1 <- iden[grep("Adi|Sparse|Mast|LYVE|SMA", iden)]
small <- small[small$named1 %in% iden1,]
ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) & scale_x_reverse() &
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & scale_y_reverse() &
  guides(colour = guide_legend(override.aes = list(size=5))) & setv


####################################################################################
####################################################################################

# Visualise myeloid fibrin example

i <- 6

nom <- pt[[i]]
Idents(all) <- "orig.ident"
M915 <- subset(all, idents=nom)
M915_meta <- M915@meta.data
test <- data_x_y[[nom]]
test$named1 <- M915$take2
iden <- unique(all$take2)

small <- test[test$Centroid.Y.µm > 3900 & test$Centroid.X.µm < 3000,]
iden1 <- iden[grep("Mast|CD68|Fibrin|LL F|Neut|ves", iden)]
small <- small[small$named1 %in% iden1,]

group.color <- c("SMA-hi vessels"="red", "CD146-hi vessels"="red", 
                 "LL Fibroblasts"= "steelblue3",  "CD68+ Myeloid"="yellow", 
                 "Mast cells"="darkorange",  "Fibrin"="lightpink", 
                 "Neutrophil-rich"="deeppink", "Neutrophil/T cell-rich"="deeppink")
                 
                 ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
                 geom_point(size=1) & scale_color_manual(values=group.color) &
                 guides(colour = guide_legend(override.aes = list(size=5))) & setv
                 
                 ####################################################################################
                 ####################################################################################
                 
                 # Visualise neutrophil fibrin example
                 
                 i <- 2
                 
                 nom <- pt[[i]]
                 Idents(all) <- "orig.ident"
                 M915 <- subset(all, idents=nom)
                 M915_meta <- M915@meta.data
                 test <- data_x_y[[nom]]
                 test$named1 <- M915$take2
                 iden <- unique(all$take2)
                 
                 small <- test[test$Centroid.X.µm < 2500  &  test$Centroid.Y.µm < 2000,]
                 iden1 <- iden[grep("CD68|Fibrin|Neut", iden)]
                 
                 group.color <- c("Fibrin "Macrophages= "yellow",  "CD68+ Myeloid"="yellow",  "Fibrin"="lightpink", 
                 "Neutrophil-rich"="deeppink", "Neutrophil/T cell-rich"="deeppink")

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=1) & coord_flip() & scale_color_manual(values=group.color) & scale_x_reverse() &
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

####################################################################################
####################################################################################

# Proximity analysis
# Read in Seurat object, filtered coordinates, vector of patients and artefactual coordinates

setwd("/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2/")

load(file="2408_kids-panel2.RData")
load(file="filtered_coordinates-2.RData")
load(file="2401-Inc-patients-noPRG4.RData")
load(file="2401_Artefact_coordinates.RData")

####################################################################################

#NN analysis- functions to use
###  It compares observed frequencies to a distribution of frequencies obtained through permutations,
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


#loop through all cell types
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


####################################################################################

#need to add a column with ident name

names <- names(data_x_y)
xy_list <- list()
for (i in 1:length(data_x_y)) {
  xy_list[[i]] <- data_x_y[[names[[i]]]]
}

all_meta <- all@meta.data
s_obj_meta <- list()
all_x_y <- list()

#need to ensure the labels and the artefactual cells are included
for (i in 1:length(pt)) {
  s_obj_meta[[i]] <- all_meta[all_meta$orig.ident == names[[i]],];
  xy_list[[i]]$named1 <- s_obj_meta[[i]]$take2;
  all_x_y[[i]] <- rbind(missing.coord[[i]], xy_list[[i]])
}

all_x_y <- rbindlist(all_x_y)

# Define neighbors
# NOTE: max_dist only refers to directly adjacent neighbors
# NOTE: Delauney triangulation draws triangles between points that are closest to each other
# expand list of triplets into symmetric list of pairs
# prune very large distances
# convert list of edges to sparse adjacency matrix and account for duplicate pairs
#take values greater than 40 and make them equal to 0
#drop 0s and covert any value above 0 as 1
#creating a diagonal matrix is a matrix where all the entries outside the main diagonal are 0
# Ignore weights. Only care if cell is a neighbor or not


coloc_res_coarse <- coloc_all_types(
  index_type = unique(all_x_y$named1),
  coords = all_x_y[, c("Centroid.X.µm", "Centroid.Y.µm")],
  y = all_x_y$named1,
  compartments = NULL,
  max_dist = 40,
  nperm = 1000,
  parallel = TRUE
)

####################################################################################

# Visualise proximity
plt_df <- coloc_res_coarse %>%
  subset(pval < 0.05) %>%
  dplyr::select(index_type, type, zscore) %>%
  spread(type, zscore, fill = 0) %>%
  column_to_rownames('index_type') %>%
  as.matrix

plt_df %>%
  Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)

# Remove artefact and low staining cells and make symmetrical
plt_df_use <- plt_df[-c(2,16),-c(2,16)]
plt_df_use2 <- plt_df_use+t(plt_df_use)
plt_df_use2 %>%
  Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)


