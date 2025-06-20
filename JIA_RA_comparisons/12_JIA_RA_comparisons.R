devtools::install_github("rpolicastro/scProportionTest")

library(Seurat)
library(ggplot2)
library(stringr)
library(tidyverse)
library(dplyr)
library(magrittr)
library(ggbiplot)
library(splitstackshape)
library(pheatmap)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library("scProportionTest")

# Include patients who underwent knee biopsies (not arthroscopies) and were DMARD-naive, with RA not OA (n = 12 AMP2 samples- Zhang et al)

# Load amp2 samples 
load(file="ampKNEE-RA-object")

ampknee <- c("BRI-413","BRI-415","BRI-421","BRI-462","BRI-475","BRI-503","BRI-542","BRI-562","BRI-605","BRI-623","BRI-625","BRI-566","BRI-401")
allAmp2 <- subset(amp2, sample %in% ampknee)
# Subset down to samples with the disease for < 1 year (n=12)
amp2 <- subset(allAmp2, idents="BRI-401", invert=T)

# Load MAPJAG tissue data & create an object of all the MAPJAG patients irrespective of disease duration
load(file="/rds/projects/c/croftap-mapjagdata/MAPJAGv2/2306/Global/sc-tissue.RData")
tissue -> allJIA

# Subset JIA tissue to be < 1 year duration (n=7)
Idents(tissue) <- "patient"
tissue <- subset(tissue, idents=c("98", "817", "84"), invert=T)

#################################################################

# Visualise main cluster proportions as a pie chart

# Rename main cell types so they share a common naming approach between both datasets
Idents(tissue) <- "integrated_snn_res.0.15"
tissue <- RenameIdents(tissue, '0'='Myeloid cells', '9'='Endothelial cells', '13'='Stromal cells', '1' = 'T cells', '3' = 'T cells', '5' = 'T cells',
                       '4'='NK cells / ILCs', '8'='T cells', '2'='Stromal cells', '12' = 'T cells', '7'='B cells', '17'='Endothelial cells', '6'='NK cells / ILCs',
                       '11'='Plasma cells', '14'='Myeloid cells', '10'='Myeloid cells', '16'='Myeloid cells', '15'="Myeloid cells")
tissue$cell_type <- tissue@active.ident

# Label plasma cells in AMP2 data and rename to give a common naming approach
Idents(amp2) <- "cluster_name"
rest <- subset(amp2, idents=c("B-2: IgG1+IgG3+ plasma", "B-6: IgM+ plasma"), invert=T)
plasma <- subset(amp2, idents=c("B-2: IgG1+IgG3+ plasma", "B-6: IgM+ plasma"))
plasma$cell_type <- "Plasma cells"
new <- rbind(plasma[["cell_type"]], rest[["cell_type"]])
amp2 <- AddMetaData(amp2, new, col.name="cell_type2")
Idents(amp2) <- "cell_type2"
amp2 <- RenameIdents(amp2, "T cell"="T cells", "NK"="NK cells / ILCs", "B cell/plasma cell"="B cells", "Myeloid cell"= "Myeloid cells", 
                     "Stromal cell"="Stromal cells", "Endothelial cell"="Endothelial cells")
amp2$cell_type2 <- amp2@active.ident

# Create pie charts 
ptnum <- as.data.frame(table(amp2$cell_type2))
ptnum$Var1 <- factor(ptnum$Var1, levels=c("T cell", "NK",  "Plasma cells", "B cell/plasma cell","Myeloid cell", "Stromal cell", "Endothelial cell"))
p1 <- ggplot(ptnum, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + scale_fill_manual(values=c("navy", "#A6CEE3","deeppink","grey", "forestgreen",   "#A65628","#E6AB02")) +
  theme_cowplot() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text = element_blank())

ptnum <- as.data.frame(table(tissue$cell_type))
ptnum$Var1 <- factor(ptnum$Var1, levels=c("T cells", "NK cells / ILCs", "Plasma cells", "B cells", "Myeloid cells", "Stromal cells", "Endothelial cells"))
p2 <- ggplot(ptnum, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + scale_fill_manual(values=c("navy", "#A6CEE3","deeppink","grey", "forestgreen",   "#A65628","#E6AB02"))  +
    theme_cowplot() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text = element_blank())

grid.arrange(p1, p2, nrow=1)

#################################################################

# Create bar charts of each patient and order by age

ptnum <- as.data.frame(table(allAmp2$orig.ident, allAmp2$cell_type2))
ptnum <- ptnum[ptnum$Freq >0,]
# Include patient with disease duration > 1 year but order them to be at the end of the line
ptnum$Var1 <- factor(ptnum$Var1, levels=c("BRI-625", "BRI-415", "BRI-542", "BRI-475","BRI-503", "BRI-421", "BRI-566", "BRI-462", "BRI-413", "BRI-605", "BRI-623", "BRI-562","BRI-401))
ptnum$Var2 <- factor(ptnum$Var2, levels=c("T cell", "NK",  "Plasma cells","B cell/plasma cell", "Myeloid cell", "Stromal cell", "Endothelial cell"))

p1 <- ggplot(ptnum) + aes(x = Var1, y = Freq, fill=Var2) + geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values=c("navy", "#A6CEE3", "deeppink","grey", "forestgreen",   "#A65628","#E6AB02")) + 
  RotatedAxis() + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_blank()) +
  labs(x="\n\n", y ="\n\n", title= "\n\n", fill="Cell Type")

ptnum2 <- as.data.frame(table(allJIA$patient, allJIA$cell_type))
ptnum2 <- ptnum2[ptnum2$Freq >0,]
# Including those with disease > 1 year but ordering so they are at the end of the line
ptnum2$Var1 <- factor(ptnum2$Var1, levels=c("99", "910","95","96","811", "814","88", "84", "97", "91","817","98"))
ptnum2$Var2 <- factor(ptnum2$Var2, levels=c("T cells", "NK cells / ILCs", "Plasma cells","B cells",  "Myeloid cells", "Stromal cells", "Endothelial cells"))

p2 <- ggplot(ptnum2) + aes(x = Var1, y = Freq, fill=Var2) + geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values=c("navy", "#A6CEE3","deeppink","grey", "forestgreen",   "#A65628","#E6AB02")) + 
  RotatedAxis() + theme(plot.title=element_blank(), axis.title=element_blank(), axis.text=element_blank()) +
  labs(x="\n\n", y ="\n\n", title= "\n")

grid.arrange(p1, p2, nrow=1)

#################################################################

# PCA of cell proportions- treatment and disease duration matched

ptnum <- as.data.frame(table(amp2$orig.ident,amp2$cell_type2))
ptnum <- ptnum[ptnum$Freq >0,]

ptnum2 <- as.data.frame(table(tissue$patient, tissue$cell_type))
ptnum2 <- ptnum2[ptnum2$Freq >0,]

combined <- rbind(ptnum, ptnum2)
dat <- dcast(combined, Var2 ~ Var1, value.var="Freq")
row.names(dat) <- dat$Var2
dat <- dat[,-1]
dat <- dat[,which(colSums(dat) > 0)]
dat <- t(dat)
dat <- dat/rowSums(dat)
dat <- scale(dat) %>% as.data.frame()

pca <- prcomp(dat)
dat$disease <- ifelse(grepl("^BRI", rownames(dat)), "RA", "JIA")
summary(pca)
pca$rotation

#stats
dat$disease <- as.factor(dat$disease)
lapply(dat[dat$disease == "JIA", -which(names(dat) == "disease")], shapiro.test)
boxM(as.matrix(dat[,-8]), dat$disease)

ggbiplot(pca, obs.scale = 1, var.scale = 1, 
         groups = dat$disease, ellipse = FALSE, 
         circle = FALSE) + geom_point(aes(color = dat$disease), size = 4) + 
  scale_color_manual(values = c("RA" = "purple", "JIA" = "turquoise")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        plot.margin = margin(t = 30, r = 10, b = 20, l = 10))

#################################################################

# Label transfer for comparing cluster proportions of the two datasets

tissue <- tissue %>% ScaleData(features = rownames(tissue)) %>%
  FindVariableFeatures() %>%
  RunPCA()
  
amp2 <- amp2 %>% ScaleData(features = rownames(amp2)) %>%
  FindVariableFeatures() %>%
  RunPCA()

# Transfer labels from JIA to RA- to compare using adult labels switch the 'reference' and 'query' objects in FindTransferAnchors() and the 'refdata' in TransferData()
transfer.anchors <- FindTransferAnchors(reference = tissue, query = amp2, dims=1:30, reduction = "rpca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = tissue$clusters2312, dims =1:30)
amp2 <- AddMetaData(amp2, metadata = celltype.predictions)

# Create one merged dataset of both JIA and RA data
tissue$predicted.id <- tissue$clusters2312
tissue$age <- "Paediatric"
amp2$age <- "Adult"
test2 <- merge(amp2, tissue)

# Subset to only visualise clusters that make up > 1% of total cells
ptnum <- as.data.frame(table(test2$predicted.id))
ptnum <- ptnum[ptnum$Freq < 680, ] %>% pull(Var1) %>% unique()
Idents(test2) <- "predicted.id"
test3 <- subset(test2, idents=ptnum, invert=T)
ptnum2 <- as.data.frame(table(test3$predicted.id, test3$age))

# Compare proportions
test <- sc_utils(test3)
prop.test <- permutation_test(test, cluster_identity = "predicted.id", 
                              sample_1="Paediatric", sample_2="Adult", sample_identity="age", n_permutations=10000)

permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T) + ylab("log2 fold difference in proportions") + xlab("")


#################################################################

# Label transfer for comparing cluster identities
# For concordance of adult and paediatric cell phenotypes we used all AMP2 RA samples and all MAPJAG tissue samples 

# Load all AMP2 samples & subset out the RA samples
load("/rds/projects/c/croftap-celldive01/amp2/amp2_full.RData")
diseaseID <- c('BRI-401', 'BRI-403', 'BRI-405', 'BRI-407', 'BRI-409', 'BRI-411', 'BRI-413', 'BRI-415', 'BRI-417', 'BRI-419',
               'BRI-421', 'BRI-425', 'BRI-427', 'BRI-429', 'BRI-431', 'BRI-436', 'BRI-440', 'BRI-458', 'BRI-460', 'BRI-462',
               'BRI-475', 'BRI-479', 'BRI-481', 'BRI-483', 'BRI-485', 'BRI-503', 'BRI-505', 'BRI-507', 'BRI-509', 'BRI-511',
               'BRI-515', 'BRI-525', 'BRI-527', 'BRI-534', 'BRI-536', 'BRI-538', 'BRI-540', 'BRI-542', 'BRI-544', 'BRI-546',
               'BRI-548', 'BRI-550', 'BRI-552', 'BRI-554', 'BRI-556', 'BRI-558', 'BRI-560', 'BRI-562', 'BRI-564', 'BRI-566',
               'BRI-570', 'BRI-581', 'BRI-583', 'BRI-589', 'BRI-601', 'BRI-603', 'BRI-605', 'BRI-607', 'BRI-611', 'BRI-613',
               'BRI-615', 'BRI-617', 'BRI-619', 'BRI-621', 'BRI-623', 'BRI-625', 'BRI-629', 'BRI-631', 'BRI-635')
Idents(amp2) <- "sample"
amp2 <- subset(amp2, sample %in% diseaseID)

# Load all MAPJAG tissue data
load(file="/rds/projects/c/croftap-mapjagdata/MAPJAGv2/2306/Global/sc-tissue.RData")

tissue <- tissue %>% ScaleData(features = rownames(tissue)) %>%
  FindVariableFeatures() %>% RunPCA()
  
amp2 <- amp2 %>% ScaleData(features = rownames(amp2)) %>%
  FindVariableFeatures() %>% RunPCA()

transfer.anchors <- FindTransferAnchors(reference = amp2, query = tissue, 
                                      features = VariableFeatures(object = amp2),
                                      reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = amp2$cluster_name,
                                     weight.reduction = tissue[["pca"]], dims = 1:30)

tissue <- AddMetaData(tissue, metadata = celltype.predictions)


#################################################################
#################################################################

# CellDive quantifications

# Load in adult RA and paediatric JIA CellDive coordinates, seurat objects of slides, vector of sample names
# Ensure only adult RA slides included in analysis
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2-adult/2407-panel2-object-2")
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2-adult/2407_filtered_coordinates-2.RData")
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2-adult/2407-Patients.RData")
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2-adult/2407_Artefact_coordinates.RData")

#Exclude non matched samples
include <- c("195", "127", "240", "115", "031", "020")
Idents(all) <- 'orig.ident'
adults <- subset(all, idents=include)
adult.missing <- missing.coord[names(missing.coord) %in% include]
adult.data <- data_x_y[names(data_x_y) %in% include]
adult.pt <- pt[pt %in% include]
save(adults, adult.missing, adult.data, adult.pt, file="/rds/projects/c/croftap-mapjagdata/CellDive/Comparisons/2407-adult-panel2-objects")

#################################################################

load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2/2408_kids-panel2.RData")
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2/filtered_coordinates-2.RData")
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2/2401-Inc-patients-noPRG4.RData")
load(file="/rds/projects/c/croftap-mapjagdata/CellDive/Panel-2/2401_Artefact_coordinates.RData")

kids.missing <- missing.coord
kids <- all
kid.data <- data_x_y
kpt <- pt

save(kids, kid.data, kids.missing, kpt, file="/rds/projects/c/croftap-mapjagdata/CellDive/Comparisons/2407-kids-panel2-objects")

#################################################################

# Exclude poor quality sample from analysis
Idents(kids) <- "orig.ident"
kids <- subset(kids, idents="1001", invert=T)
  
# Remove Fibrin regions from kids slides as these are not synovium and may have a large impact on cell proportions
Idents(kids) <- "take2"
kids <- subset(kids, idents=c("Fibrin", "Fibrin Macrophages",  "RBCs", "Neutrophil-rich", "Neutrophil/T cell-rich",
                               "Artefact", "Blood/Fibrin"), invert=T)

# Ensure common naming strategy
kids <- RenameIdents(kids, "Endothelial cells-1"="Vascular", "Endothelial cells-2"="Vascular", 
                     "CD146-hi vessels"="Vascular", "SMA-hi vessels"="Vascular", "CD68+ Myeloid"="Myeloid cells", 
                     "LL Macrophages"="Myeloid cells", "Macrophages-2"="Myeloid cells", 
                     "LYVE1-hi Macrophages"="Myeloid cells")
kids$take3 <- kids@active.ident

# Remove fibrin regions from adult slides
Idents(adults) <- "take2"
adults <- subset(adults, idents=c("Fibrin", "Artefact", "RBCs", "Muscle"), invert=T)

adults <- RenameIdents(adults, "SMA-hi vessels"="Vascular", "SMA-low vessels"="Vascular", 
                       "LYVE1+ macrophages"="Myeloid cells", "LL macrophages"="Myeloid cells", 
                       "SL macrophages"="Myeloid cells", "B/T aggregates"="B/T cell aggregates", 
                       "LL fibroblasts"="LL Fibroblasts")
adults$take3 <- adults@active.ident

# Proportion of cells per donor per slide
ptnum2 <- as.data.frame(table(kids$take3, kids$orig.ident))
ptnum <- as.data.frame(table(adults$take3, adults$orig.ident))

calculate_percentage <- function(data) {
  data %>% group_by(Var2) %>% mutate(Percentage = Freq / sum(Freq) * 100) }

ptnum2 <- calculate_percentage(ptnum2)
ptnum <- calculate_percentage(ptnum)


# Function to perform appropriate statistical test
perform_test <- function(data, group_col, value_col) {
  # Perform Shapiro-Wilk normality test for each group
  shapiro_kids <- shapiro.test(data[[value_col]][data[[group_col]] == "Kids"])
  shapiro_adults <- shapiro.test(data[[value_col]][data[[group_col]] == "Adults"])
  if (shapiro_kids$p.value > 0.05 && shapiro_adults$p.value > 0.05) {
    # Both groups are normally distributed, use t-test
    test_result <- t.test(as.formula(paste(value_col, "~", group_col)), data = data)
  } else {
    # At least one group is not normally distributed, use Wilcoxon test
    test_result <- wilcox.test(as.formula(paste(value_col, "~", group_col)), data = data)
  }
  return(test_result$p.value)
}

plots <- list()
cell.types <- c("Vascular", "Myeloid cells", "Lymphatics")

for (i in seq_along(cell.types)) { 
  vascular_kids <- ptnum2 %>% filter(Var1 == cell.types[i])
  vascular_adults <- ptnum %>% filter(Var1 == cell.types[i])
  combined_data <- rbind(data.frame(Group = "Kids", vascular_kids), data.frame(Group = "Adults", vascular_adults))
  p_value <- perform_test(combined_data, "Group", "Percentage")
  plots[[i]] <- ggplot(combined_data, aes(x = Group, y = Percentage, fill = Group)) +
    geom_boxplot(alpha=0.2) + geom_jitter(show.legend=FALSE) + theme_classic() + labs(x = "\nGroup",y = "Percentage\n", fill = "Group") +
    scale_fill_manual(values = c("Kids" = "green3", "Adults" = "purple")) + labs(title=paste0(unique(combined_data$Var1), "\n")) +
    geom_signif(comparisons = list(c("Kids", "Adults")), 
                annotations = paste("p =", format(p_value, digits = 3)), 
                textsize = 4, y_position = max(combined_data$Percentage) + 1) #+ ylim(min(combined_data$Percentage), (max(combined_data$Percentage) + max(combined_data$Percentage)/3))
}

grid.arrange(grobs=plots, ncol=3)

#################################################################

# Visualise JIA CellDive example with a reasonable lining layer shown

setv <- list(theme(panel.background = element_rect(fill = 'black'), 
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text = element_blank(),
                   legend.background = element_rect(fill = "black", color = NA),
                   legend.key = element_rect(color = "gray", fill = "black"),
                   legend.text = element_text(color = "white"))) 


i <- 7
nom <- kpt[[i]]
Idents(kids) <- "orig.ident"
M915 <- subset(kids, idents=nom)
M915_meta <- M915@meta.data
test <- kid.data[[nom]]
test$named1 <- M915$take2
iden <- unique(kids$take2)

ggplot(test, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.6) & scale_y_reverse() & scale_color_manual(values=bright) & 
  guides(colour = guide_legend(override.aes = list(size=5))) 

iden1 <- iden[grep("LL|vess|Lymph|LYVE|CD146|Endo|RBC|Mye", iden)]
small <- test[test$Centroid.X.µm > 2500  &  test$Centroid.Y.µm < 3500,]
small <- small[small$named1 %in% iden1,]

group.color2 <- c("SMA-hi vessels"="cyan", "CD146-hi vessels"="blue", "CD68+ Myeloid"="lemonchiffon",
                 "LL Fibroblasts"= "magenta", "LL Macrophages"="yellow", "LYVE1-hi Macrophages"="orange",
                 "RBCs"="red", "Endothelial cells-1"="palevioletred", "Endothelial cells-2"="palevioletred",
                 "Lymphatics"="burlywood1")

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.4) & scale_y_reverse() & scale_color_manual(values=group.color2) &
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

###################################################################### 

# Adult RA example fragment with reasonable lining layer

i <- 4
nom <- adult.pt[[i]]
Idents(adults) <- "orig.ident"
M915 <- subset(adults, idents=nom)
M915_meta <- M915@meta.data
test <- adult.data[[nom]]
#test$named1 <- M915$SCT_snn_res.0.16
test$named1 <- M915$take2
iden <- unique(adults$take2)

iden1 <- iden[grep("LL|vess|Lymph|LYVE|RBCs", iden)]
small <- test[test$Centroid.Y.µm > 7400  & test$Centroid.X.µm > 6000,]
small <- small[small$named1 %in% iden1,]

group.color <- c("SMA-hi vessels"="cyan", "SMA-low vessels"="blue", "SL macrophages"="olivedrab1",
                 "LL fibroblasts"= "magenta", "LL macrophages"="yellow", "LYVE1+ macrophages"="orange",
                 "Lymphatics"="palevioletred", "RBCs"="red")

ggplot(small, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named1)) +
  geom_point(size=0.4) & scale_y_reverse() & scale_color_manual(values=group.color) & 
  guides(colour = guide_legend(override.aes = list(size=5))) & setv

