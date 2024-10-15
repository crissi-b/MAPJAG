
load("/rds/projects/c/croftap-celldive01/xenium/analysis/analysis_step1.RData")
rm(list=ls()[! ls() %in% c("all_tx_f")])


all_csv <- dir(path="/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs", pattern=".csv")
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/", tx_csv, sep="")
samples_area<-paste("/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/", cell_stat_csvs, sep="")



data = list()
for (i in 1:length(samples)) {
    data[[i]] <- read.csv(samples[[i]]);
    data[[i]] <- data[[i]] %>% filter(is_noise == "false")
}


library(data.table)
data_plot<-rbindlist(data)

data_plot_f <- data_plot %>% filter(is_noise == "false")


tx_to_counts <- function(genes, cells, remove_bg = TRUE) {
    if (remove_bg) {
        idx <- which(cells != 0)
        cells <- cells[idx]
        genes <- genes[idx]
    }
    genes <- factor(genes)
    cells <- factor(cells)
    counts <- Matrix::sparseMatrix(
        i = as.integer(genes), 
        j = as.integer(cells), 
        x = rep(1, length(genes)),
        dims = c(length(levels(genes)), length(levels(cells)))
    )
    rownames(counts) <- levels(genes)
    colnames(counts) <- levels(cells)
    return(counts)
}




count_matrices=list()
for (i in 1:1) {
  count_matrices[[i]]<-tx_to_counts(genes =data[[i]]$gene, data[[i]]$cell)
  }


smaples_fov<-gsub(pattern=".csv", x=tx_csv, replacement="")
seurat_obj<-list()
for (i in 1:length(samples)) {
  seurat_obj[[i]]<-CreateSeuratObject(counts =count_matrices[[i]], project = smaples_fov[[i]])
  
}

all_merged<-merge(x=seurat_obj[[1]], y=seurat_obj[c(2:length(seurat_obj))])

median(all_merged$nCount_RNA)  #19
median(all_merged$nFeature_RNA)  #17
hist(all_merged$nCount_RNA)  
hist(all_merged$nFeature_RNA)

library(splitstackshape)
unique(all_merged$orig.ident)

df_plot <- all_merged$nCount_RNA %>% as.data.frame() %>% cbind(all_merged$orig.ident) %>% cSplit(splitCols ="all_merged$orig.ident", sep="_" )
colnames(df_plot) <- c("c1",                       "c2", "c3", "c4")
df_plot$donor <- paste(df_plot$c2, df_plot$c3, sep="_")
df_plot %>% 
ggplot(aes(x=c1,fill=donor)) + geom_histogram()+xlim(0,400)
  
  
df_plot <- all_merged$nFeature_RNA %>% as.data.frame() %>% cbind(all_merged$orig.ident) %>% cSplit(splitCols ="all_merged$orig.ident", sep="_" )
colnames(df_plot) <- c("c1",                       "c2", "c3", "c4")
df_plot$donor <- paste(df_plot$c2, df_plot$c3, sep="_")
df_plot %>% 
ggplot(aes(x=c1,fill=donor)) + geom_histogram()+xlim(0,100)
  
all_merged_f <- subset(all_merged, subset = nFeature_RNA > 15 & nCount_RNA > 20)

all_merged_f <- all_merged_f %>% 
  NormalizeData(scale.factor = median(all_merged_f$nCount_RNA)) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:40)

DimPlot(all_merged_f, raster=FALSE)+NoLegend()

data_plot[data_plot$fov_name %in% c("N6"),] %>% ggplot(aes(x=x,y=y)) + 
  geom_point(alpha=0.6, shape=".") +
  geom_point(data=data_plot[data_plot$fov_name %in% c("N6") & data_plot$gene %in% c("VWF", "ACTA2", "TNC", "THY1"),],
             aes(x=x,y=y, color=gene), 
                   size=0.6)+ scale_color_brewer(palette="Dark2")+xlim(1100,1250)+ylim(1300,1400)


data_plot[data_plot$fov_name %in% c("P6"),] %>% ggplot(aes(x=x,y=y)) + 
  geom_point(alpha=0.6, shape=".") +
  geom_point(data=data_plot[data_plot$fov_name %in% c("P6") & data_plot$gene %in% c("PRG4"),],
             aes(x=x,y=y, color=gene), 
                   size=0.6)+ scale_color_brewer(palette="Dark2")
}
library(harmony)
all_merged_f <- RunHarmony(all_merged_f, group.by.vars = c("orig.ident"))
all_merged_f <- RunUMAP(all_merged_f, reduction="harmony",  dims=1:50)
all_merged_f <- FindNeighbors(all_merged_f, reduction="harmony", dims=1:50)
all_merged_f <- FindClusters(all_merged_f, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4))
all_merged_f <- FindClusters(all_merged_f, resolution = c(0.6, 0.7, 0.8))
DimPlot(all_merged_f, group.by="orig.ident", raster=FALSE)+NoLegend()
DimPlot(all_merged_f, group.by="RNA_snn_res.0.8", raster=FALSE)+NoLegend()
DimPlot(all_merged_f, group.by="RNA_snn_res.0.4", label = T, raster=FALSE)+NoLegend()

DimPlot(all_merged_f, group.by="RNA_snn_res.0.1", raster=FALSE)+NoLegend()

FeaturePlot(all_merged_f, features = "PRG4", raster=FALSE)

FeaturePlot(all_merged_f, features = "PECAM1", raster=FALSE)
FeaturePlot(all_merged_f, features = "CD79A", raster=FALSE)



rownames(all_merged_f)

Idents(all_merged_f) <- 'RNA_snn_res.0.1'
markers_0.1 <- FindAllMarkers(all_merged_f, only.pos = T)


Idents(all_merged_f) <- 'RNA_snn_res.0.4'
markers_0.4 <- FindAllMarkers(all_merged_f, only.pos = T)

rownames(all_merged_f)

all_merged_f$named <- all_merged_f@meta.data[["RNA_snn_res.0.4"]]
Idents(all_merged_f) <- 'named'
levels(all_merged_f)
current.sample.ids <- c( "0" , "1",  "10", "11", "12", "13", "14" ,"15", "16", "17", "18" ,"19", "2",  "20", "21", "22" ,"23", "24", "25","26","27", "3" , "4" , "5" , "6" , "7",  "8",  "9" )
new.sample.ids <- c(  "macs" , "fibs",  "fibs", "blood_prog", "Bcells", "Tcells", "macs" ,"Tcells", "Bcells", "adipose", "fibs" ,"Bcells", "EC",  "CD1C_DCs", "Lamp3_DCs", "fibs" ,"fibs", "fibs", "fibs","26","27", "pericytes" , "Bcells" , "Tcells" , "LL" , "prolif",  "fibs",  "lymphatics" )

all_merged_f@meta.data[["named"]] <- plyr::mapvalues(x = all_merged_f@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)

cols <- ArchR::paletteDiscrete(all_merged_f@meta.data[, "named"])
DimPlot(all_merged_f, group.by="named", cols=cols, label=T, raster=FALSE)+NoAxes()

all_merged_f$named %>% 

DimPlot(all_merged_f, group.by="RNA_snn_res.0.4", raster=FALSE,label=T)+NoLegend()

all_merged_f_chrissy_check <- `2404-xenium-global`
cols <- ArchR::paletteDiscrete(all_merged_f_chrissy_check@meta.data[, "global13"])
DimPlot(all_merged_f_chrissy_check, group.by="global13", raster=FALSE,label=F, cols=cols)


library(furrr)
library(sfdct)
cellgeoms_baysor<-function(segfile){
    #system.time({ 
    #    plan(multicore, workers = availableCores(constraints = 'multicore') - 1)
        # delete transcripts that are noise
        #segfile<-segfile %>% filter(is_noise == FALSE)
        transcriptspercell<-furrr::future_map_dfr(.x = unique(segfile$cell), 
                                                  .f = ~ data.frame(
                                                      cell = .x, 
                                                      num_transcripts = sum(segfile$cell == .x)
                                                  ), 
                                                  .options = furrr_options(seed = TRUE)
        )
        cellidx <- transcriptspercell$cell[transcriptspercell$num_transcripts > 5]
        segfile.new <- furrr::future_map_dfr(.x = cellidx, function(.x) { 
            res <- st_as_sf(segfile[segfile$cell == .x, c('x', 'y')], coords = c('x', 'y')) %>%
                st_union() %>% #dont remove the union. It is needed here.
                ct_triangulate()
            resdf <- data.frame(cell = .x, geometry = res)
            return(resdf)
        }, .options = furrr_options(seed = TRUE))
        
        cellgeoms_final<-segfile.new$geometry %>% 
            furrr::future_map(purrr::reduce, st_union, .options = furrr_options(seed = TRUE)) %>%
            st_sfc() %>%
            as.data.frame()
        
        cellgeoms_final<-cellgeoms_final %>%
            cbind(transcriptspercell[transcriptspercell$cell %in% cellidx, ])
        
        return(cellgeoms_final)
        
    
    
}


library(sf)

all_tx_f_P6 <- all_tx_f %>% filter(fov_name == "P6")
all_tx_f_P6$cell <- all_tx_f_P6$cell_id_new
all_tx_f_P6$x <- all_tx_f_P6$x_location
all_tx_f_P6$y <- all_tx_f_P6$y_location

cellgeoms_baysorfov_P6 <- cellgeoms_baysor(all_tx_f_P6)
cellgeoms_baysorfov_P6$cell <- gsub(pattern = "_44", cellgeoms_baysorfov_P6$cell, replacement="")

all_merged_f_meta <- all_merged_f@meta.data %>% filter(orig.ident=="D2518_11_P6")

cellgeoms_baysorfov_P6$cell <- paste(cellgeoms_baysorfov_P6$cell, "_44", sep="")
index <- match(cellgeoms_baysorfov_P6$cell, rownames(all_merged_f_meta))
cellgeoms_baysorfov_P6$named_new <- all_merged_f_meta$named_coarse[index]


cellgeoms_baysorfov_P6 %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal()+ scale_fill_manual(values = cols$colors)#+ylim(2600,2900)+
   #geom_point(data=all_tx_f_P6,aes(x, y), shape = '.', size=0.0) +theme_minimal()+ylim(2600,2900)



cols <- as.data.frame(ArchR::paletteDiscrete(all_merged_f@meta.data[, "named_coarse"]))
colnames(cols)<-"colors"


cellgeoms_baysorfov_P6 %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal() +
  geom_point(data=data_plot[data_plot$fov_name %in% c("P6") & data_plot$gene %in% c("PRG4", "CD68"),],
             aes(x=x,y=y, color=gene), 
                   size=1)+ scale_fill_manual(values = cols$colors)#+xlim(1000,1200)+ylim(2700,2900)


cellgeoms_baysorfov_P6 %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.5,
              color = NA)+theme_minimal() +
  geom_point(data=data_plot[data_plot$fov_name %in% c("P6") & data_plot$gene %in% c("CD3E", "CCL5", "CD8A"),],
             aes(x=x,y=y, color=gene), 
                   size=1.6)+ scale_fill_manual(values = cols$colors)+xlim(1100,1200)+ylim(2750,2880)
library(readr)
tx_08 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-08__20231107__215839/transcripts.csv")
tx_11 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-11__20231107__215839/transcripts.csv")
tx_04 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2519-04__20231107__215839/transcripts.csv")

tx_18_24 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data_batch2/Xenium_data_for_Chris_Mahoney/output-XETG00150__0018573__2518-24__20240227__000410/transcripts.csv.gz")
tx_19_15 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data_batch2/Xenium_data_for_Chris_Mahoney/output-XETG00150__0018573__2519-15__20240227__000410/transcripts.csv.gz")

tx_18_28 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014095__2518-28__20240331__210931/transcripts.csv.gz")

tx_18_36 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014095__2518-36__20240331__210931/transcripts.csv.gz")
  
tx_18_40 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014097__2518-40__20240331__210931/transcripts.csv.gz")
  
tx_19_18 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014097__2519-18__20240331__210931/transcripts.csv.gz")

tx_08$slide <- '1'
tx_11$slide <- '2'
tx_04$slide <- '3'
tx_18_24$slide <- '4'
tx_19_15$slide <- '5'

tx_18_28$slide <- '6'

tx_18_36$slide <- '7'

tx_18_40$slide <- '8'

tx_19_18$slide <- '9'

#####


all_no_f_slide <- rbind(tx_08, tx_11,tx_04, tx_18_24, tx_19_15, tx_18_28, tx_18_36, tx_18_40,tx_19_18)

all_no_f_slide[all_no_f_slide$transcript_id %in% sample(all_no_f_slide$transcript_id, size=nrow(all_no_f_slide)*0.01),] %>%  ggplot(aes(x=x_location, y=y_location, color=fov_name))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()+facet_wrap("slide")+RotatedAxis()

all_no_f_slide[all_no_f_slide$transcript_id %in% sample(all_no_f_slide$transcript_id, size=nrow(all_no_f_slide)*0.01),] %>%  ggplot(aes(x=x_location, y=y_location))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()+facet_wrap("slide")+RotatedAxis()


library(readr)

tx_08 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-08__20231107__215839/transcripts.csv")
tx_11 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-11__20231107__215839/transcripts.csv")
tx_04 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2519-04__20231107__215839/transcripts.csv")

files <- list(tx_04, tx_08, tx_11)
donor <- c("D2518_08", "D2518_11", "D2519_04")

for (i in 1:length(donor)){
files[[i]]$donor <- donor[[i]]
files[[i]]$donor_FOV <- paste(files[[i]]$donor, files[[i]]$fov_name, sep="_")
  }
names(files) <- donor

library(data.table)
all_tx <- rbindlist(files)
head(all_tx)

all_tx_f <- all_tx %>% filter(feature_name != "BLANK" & feature_name != "NegControl")

all_tx_f <- all_tx %>% filter(!grepl('BLANK|NegControl', feature_name) & cell_id != "UNASSIGNED")

all_tx_f$cell_id_new <- as.numeric(as.factor(all_tx_f$cell_id))


tx_18_24 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data_batch2/Xenium_data_for_Chris_Mahoney/output-XETG00150__0018573__2518-24__20240227__000410/transcripts.csv.gz")
tx_19_15 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data_batch2/Xenium_data_for_Chris_Mahoney/output-XETG00150__0018573__2519-15__20240227__000410/transcripts.csv.gz")


files2 <- list(tx_18_24, tx_19_15)
donor <- c("D2518_24", "D2519_15")

for (i in 1:length(donor)){
files2[[i]]$donor <- donor[[i]]
files2[[i]]$donor_FOV <- paste(files2[[i]]$donor, files2[[i]]$fov_name, sep="_")
  }
names(files2) <- donor

library(data.table)
all_tx_b2 <- rbindlist(files2)
head(all_tx_b2)

all_tx_b2_f <- all_tx_b2 %>% filter(!grepl('BLANK|NegControl|UnassignedCodeword_', feature_name) & cell_id != "UNASSIGNED")

head(all_tx_b2_f)

all_tx_b2_f$cell_id_new <- as.numeric(as.factor(all_tx_b2_f$cell_id))


tx_18_28 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014095__2518-28__20240331__210931/transcripts.csv.gz")

tx_18_36 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014095__2518-36__20240331__210931/transcripts.csv.gz")
  
tx_18_40 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014097__2518-40__20240331__210931/transcripts.csv.gz")
  
tx_19_18 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014097__2519-18__20240331__210931/transcripts.csv.gz")



files3 <- list(tx_18_28, tx_18_36, tx_18_40, tx_19_18)
donor <- c("D2518_28", "D2518_36", "D2518_40", "D2519_18")

for (i in 1:length(donor)){
files3[[i]]$donor <- donor[[i]]
files3[[i]]$donor_FOV <- paste(files3[[i]]$donor, files3[[i]]$fov_name, sep="_")
  }
names(files3) <- donor

library(data.table)
all_tx_b3 <- rbindlist(files3)
head(all_tx_b3)

all_tx_b3_f <- all_tx_b3 %>% filter(!grepl('BLANK|NegControl|UnassignedCodeword_', feature_name) & cell_id != "UNASSIGNED")

all_tx_b3_f$cell_id_new <- as.numeric(as.factor(all_tx_b3_f$cell_id))

#basor cell assingment

data_plot


rm(list=ls()[! ls() %in% c("data_plot_f")])
data_plot_f


library(furrr)
#library(spatula)
library(sp)
library(sf)
library(sfdct)
FOVs <- list()
geoms_all <- list()
for (i in 1:length(unique(data_plot_f$donor_FOV))){
 FOVs[[i]] <-  data_plot_f %>% filter(donor_FOV == unique(data_plot_f$donor_FOV)[[i]])
    geoms_all[[i]] <- cellgeoms_baysor(FOVs[[i]])
}

names(geoms_all) <- unique(data_plot_f$donor_FOV)

length(geoms_all)
length(unique(all_merged_f@meta.data$orig.ident))



library(furrr)
#library(spatula)
library(sp)
library(sf)
library(sfdct)
FOVs <- list()
geoms_all <- list()
for (i in 1:length(unique(data_plot_f$donor))){
 FOVs[[i]] <-  data_plot_f %>% filter(donor == unique(data_plot_f$donor)[[i]])
    geoms_all[[i]] <- cellgeoms_baysor(FOVs[[i]])
}

names(geoms_all) <- unique(data_plot_f$donor)





rm(list=ls()[! ls() %in% c("data_plot_f", "geoms_all", "FOVs")])
gc()

fovs_df <- unique(all_merged_f@meta.data$orig.ident) %>% as.data.frame()
colnames(fovs_df) <- "all"
fov_order <- fovs_df$all
geoms_f <- geoms_all[c(fov_order)]

meta_list <- list()
cols <- as.data.frame(ArchR::paletteDiscrete(all_merged_f@meta.data[, "named"]))
colnames(cols)<-"colors"
ggplots <- list()
numbers <- 1:length(geoms_f)

for (i in 1:length(unique(all_merged_f@meta.data$orig.ident))){
   meta_list[[i]] <-   all_merged_f@meta.data %>% filter(orig.ident==unique(all_merged_f@meta.data$orig.ident)[[i]])
  # geoms_f[[i]]$cell <- gsub(pattern = "_", geoms_f[[i]]$cell, replacement="")
 #geoms_f[[i]]$cell <- paste(geoms_f[[i]]$cell, "_",numbers[[i]], sep="")
index <- match(geoms_f[[i]]$cell, rownames(meta_list[[i]]))
geoms_f[[i]]$named_new <- meta_list[[i]]$named[index]
ggplots[[i]] <- geoms_f[[i]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal()+ scale_fill_manual(values = cols$colors)+ggtitle(unique(all_merged_f@meta.data$orig.ident)[[i]])
print(ggplots[[i]])
}


geoms_f[["D2518_36_Y16"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = "black")+ scale_fill_manual(values = cols$colors) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

cols$cols2 <- c("grey", "#272E6A", "grey", "grey", "grey", "grey", "darkgreen", "grey" ,"grey", "red", "grey", "grey", "grey")

geoms_f[["D2518_36_Y16"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = "black")+ scale_fill_manual(values = cols$cols2) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())


geoms_f[["D2518_11_O7"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal() + scale_fill_manual(values = cols$colors)+xlim(1350,1800)+ylim(1700,2000)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


geoms_f[["D2518_11_O7"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.1,
              color = NA)+theme_minimal() +
  geom_point(data=data_plot[data_plot$donor_FOV %in% c("D2518_11_O7") & data_plot$gene %in% c("PRG4", "CD68"),],
             aes(x=x,y=y, color=gene), 
                   size=0.4)+ scale_fill_manual(values = cols$colors)+xlim(1350,1800)+ylim(1700,2000)


geoms_f[["D2518_11_O7"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.1,
              color = NA)+theme_minimal() +
  geom_point(color="#89C75F", data=data_plot[data_plot$donor_FOV %in% c("D2518_11_O7") & data_plot$gene %in% c("PRG4"),],
             aes(x=x,y=y), 
                   size=0.4)+ scale_fill_manual(values = cols$colors)+xlim(1350,1800)+ylim(1700,2000)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#rm(list=ls()[! ls() %in% c("fibs_X_meta")])
 #xenium_chrissy <- readRDS("/rds/projects/c/croftap-celldive01/xenium/analysis/chrissy_annotation/2404-xenium-global.rds")

bcells_meta <- read.csv("/rds/projects/c/croftap-mapjagb6/xenium2/bcells_meta.csv", row.names = 1)
tcells_meta <- read.csv("/rds/projects/c/croftap-mapjagb6/xenium2/tcells_meta.csv", row.names = 1)
myeloid_meta <- read.csv("/rds/projects/c/croftap-mapjagb6/xenium2/myeloid_meta.csv", row.names = 1)


xenium_chrissy_meta <- xenium_chrissy@meta.data

xenium_chrissy_meta <- xenium_chrissy_meta[!rownames(xenium_chrissy_meta) %in% rownames(fibs_X_meta),]

myeloid_meta_f <- myeloid_meta %>% dplyr::select(named_sub_new)
tcells_meta_f <- tcells_meta %>% dplyr::select(named_sub_new)
bcells_meta_f <- bcells_meta %>% dplyr::select(named_sub_new)

xenium_chrissy_meta <- xenium_chrissy_meta[!rownames(xenium_chrissy_meta) %in% rownames(myeloid_meta_f),]
xenium_chrissy_meta <- xenium_chrissy_meta[!rownames(xenium_chrissy_meta) %in% rownames(tcells_meta_f),]
xenium_chrissy_meta <- xenium_chrissy_meta[!rownames(xenium_chrissy_meta) %in% rownames(bcells_meta_f),]


xenium_chrissy_meta <- xenium_chrissy_meta %>% dplyr::select(global13)
colnames(xenium_chrissy_meta) <- 'named_sub_new'
fibs_X_meta_f <- fibs_X_meta %>% dplyr::select(named_sub_new)



meta_new <- rbind(xenium_chrissy_meta, fibs_X_meta_f, myeloid_meta_f, tcells_meta_f, bcells_meta_f)

xenium_chrissy <- Seurat::AddMetaData(xenium_chrissy, metadata=meta_new)

meta_list <- list()
cols <- as.data.frame(ArchR::paletteDiscrete(xenium_chrissy@meta.data[, "named_sub_new"]))
colnames(cols)<-"colors"

cols$colors <- c( "grey", "grey", "#208045", "grey", "#DE6C3E", "grey", "grey", "grey", "grey" ,"grey", "grey", "#B8A46C", "#CE7E96", "grey", "#87221E", "grey", "#5A5297", "grey", "grey", "#D0A164", "grey")

table(xenium_chrissy$named_sub_new, xenium_chrissy$orig.ident) %>% as.data.frame() %>% filter(Var1=="SOX5/CDH11")

ggplots <- list()
numbers <- 1:length(geoms_f)

for (i in 1:length(unique(xenium_chrissy@meta.data$orig.ident))){
   meta_list[[i]] <-   xenium_chrissy@meta.data %>% filter(orig.ident==unique(xenium_chrissy@meta.data$orig.ident)[[i]])
  # geoms_f[[i]]$cell <- gsub(pattern = "_", geoms_f[[i]]$cell, replacement="")
 #geoms_f[[i]]$cell <- paste(geoms_f[[i]]$cell, "_",numbers[[i]], sep="")
index <- match(geoms_f[[i]]$cell, rownames(meta_list[[i]]))
geoms_f[[i]]$named_new <- meta_list[[i]]$named_sub_new[index]
ggplots[[i]] <- geoms_f[[i]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal()+ scale_fill_manual(values = cols$colors)+ggtitle(unique(all_merged_f@meta.data$orig.ident)[[i]])
print(ggplots[[i]])
}


#work out labels from chrissy

xenium_chrissy$global13
xenium_chrissy$named_sub_new
xenium_chrissy$named

DimPlot(xenium_chrissy, group.by = "global13", raster=FALSE)

DimPlot(xenium_chrissy, group.by = "named_sub_new", raster=FALSE)


xenium_chrissy$named_sub_new_named <- paste(xenium_chrissy$named_sub_new, xenium_chrissy$named, sep="_")

Idents(xenium_chrissy) <- 'named_sub_new_named'

keep <- unique(xenium_chrissy$named_sub_new_named)[grep("Fibroblasts_", unique(xenium_chrissy$named_sub_new_named))]

subset_xenium <- subset(xenium_chrissy, idents=keep[c(1,2,3,5,6,7,10,12,14,15,17)])


DimPlot(subset_xenium, group.by = "named_sub_new_named", raster=FALSE, label = T)+NoLegend()

xenium_chrissy_meta <- xenium_chrissy@meta.data

xenium_chrissy_meta <- xenium_chrissy_meta[!rownames(xenium_chrissy_meta) %in% colnames(subset_xenium),]
xenium_chrissy_meta <- xenium_chrissy_meta %>% dplyr::select(named_sub_new)


Idents(subset_xenium) <- 'named_sub_new_named'
levels(subset_xenium) %>% unique()

levels(xenium_chrissy) %>% unique()


current.sample.ids <- c("Fibroblasts_Bcells" ,    "Fibroblasts_Tcells" ,    "Fibroblasts_macs",       "Fibroblasts_prolif"   , "Fibroblasts_blood_prog" ,"Fibroblasts_pericytes" , "Fibroblasts_CD1C_DCs" ,"Fibroblasts_Lamp3_DCs" , "Fibroblasts_lymphatics" ,"Fibroblasts_adipose",    "Fibroblasts_EC")
new.sample.ids <- c("B_cell_doublets" ,    "Tcell/Fibroblast" ,    "Fibroblasts_macs",       "Cycling cells"   , "Granulocytes/progenitors" ,"Pericytes" , "CD1C_DCs" ,"Lamp3_DCs" , "Lymphatics" ,"Adipose",    "Endothelial cells")

subset_xenium@meta.data[["named_sub_new_named"]] <- plyr::mapvalues(x = subset_xenium@meta.data[["named_sub_new_named"]], from = current.sample.ids, to = new.sample.ids)

subset_xenium_meta <- subset_xenium@meta.data
subset_xenium_meta <- subset_xenium_meta %>% dplyr::select(named_sub_new_named)
colnames(subset_xenium_meta) <- 'named_sub_new'

meta_new <- rbind(xenium_chrissy_meta, subset_xenium_meta)

xenium_chrissy <- Seurat::AddMetaData(xenium_chrissy, metadata=meta_new)


DimPlot(xenium_chrissy, group.by = "named_sub_new", raster=FALSE, label = T)+NoLegend()


Idents(xenium_chrissy) <- 'named_sub_new'
levels(xenium_chrissy)
current.sample.ids <- c("Endothelial cells"  ,         "Plasma cells"     ,           "Pericytes"    ,              
  "Lymphatics"       ,           "Cycling cells"    ,           "Adipose"  ,                  
  "Granulocytes / progenitors",  "CD34/MFAP5"       ,           "POSTN"   ,                   
 "Tcell_contam"         ,       "SOX5/CDH11"        ,          "macs"     ,                  
 "LL"                   ,       "pericyte"          ,          "EC"        ,                 
 "SFRP/CXCL12"          ,       "mac/LL"            ,          "SPP1+ Macrophages"  ,        
 "Tcell_contam_myeloid" ,       "S1008A+ Monocytes"  ,         "Peri_contam_myeloid" ,       
 "MertK+ macrophages"   ,       "proliferating_myeloid" ,      "CD1C+ cDC2s"          ,      
 "fib_contam_myeloid"   ,       "LAMP3 DCs"             ,      "B_cell_contam_myeloid" ,     
 "CD4+ naive/central memory T", "Vascular/T cell"       ,      "T cell/Fibroblast"      ,    
 "CD4+ KLRB1+ T"              , "FOXP3+ Tregs"          ,      "GZMK+ CD8+ T"  ,             
 "T cells/LL"                 , "NK cells/ILCs"          ,     "Cycling T"          )

new.sample.ids <- c("Endothelial cells"  ,         "Plasma cells"     ,           "Pericytes"    ,              
  "Lymphatics"       ,           "Cycling cells"    ,           "Adipose"  ,                  
  "Granulocytes / progenitors",  "CD34/MFAP5"       ,           "POSTN"   ,                   
 "Tcell_contam"         ,       "SOX5/CDH11"        ,          "macs"     ,                  
 "LL"                   ,       "Pericytes"          ,          "Endothelial cells"        ,                 
 "SFRP/CXCL12"          ,       "mac/LL"            ,          "SPP1+ Macrophages"  ,        
 "Tcell_contam_myeloid" ,       "S1008A+ Monocytes"  ,         "Peri_contam_myeloid" ,       
 "MertK+ macrophages"   ,       "proliferating_myeloid" ,      "CD1C+ cDC2s"          ,      
 "fib_contam_myeloid"   ,       "LAMP3 DCs"             ,      "B_cell_contam_myeloid" ,     
 "CD4+ naive/central memory T", "Vascular/T cell"       ,      "T cell/Fibroblast"      ,    
 "CD4+ KLRB1+ T"              , "FOXP3+ Tregs"          ,      "GZMK+ CD8+ T"  ,             
 "T cells/LL"                 , "NK cells/ILCs"          ,     "Cycling cells" )

xenium_chrissy@meta.data[["named_sub_new"]] <- plyr::mapvalues(x = xenium_chrissy@meta.data[["named_sub_new"]], from = current.sample.ids, to = new.sample.ids)



Idents(xenium_chrissy) <- 'named_sub_new'
levels(xenium_chrissy)

xenium_chrissy$named_sub_new2 <- xenium_chrissy$named_sub_new
current.sample.ids <- c("Endothelial cells"   ,        "Plasma cells"      ,          "Pericytes",                  
  "Lymphatics"              ,    "Cycling cells" ,              "Adipose"       ,             
  "Granulocytes / progenitors",  "CD34/MFAP5"     ,             "POSTN"          ,            
 "Tcell_contam"       ,         "SOX5/CDH11"       ,           "macs"             ,          
 "LL"                  ,        "pericyte"          ,          "SFRP/CXCL12"       ,         
 "mac/LL"               ,       "SPP1+ Macrophages"  ,         "Tcell_contam_myeloid" ,      
 "S1008A+ Monocytes"     ,      "Peri_contam_myeloid" ,        "MertK+ macrophages"    ,     
 "proliferating_myeloid"  ,     "CD1C+ cDC2s"          ,       "fib_contam_myeloid"     ,    
 "LAMP3 DCs"              ,     "B_cell_contam_myeloid" ,      "CD4+ naive/central memory T",
 "Vascular/T cell"        ,     "T cell/Fibroblast"     ,      "CD4+ KLRB1+ T"      ,        
 "FOXP3+ Tregs"           ,     "GZMK+ CD8+ T"          ,      "T cells/LL"          ,       
 "NK cells/ILCs"          ,     "T cell/Plasma cell"     ,     "pDCs"                 ,      
 "Plasma cells/Myeloid"   ,     "Vascular plasma cells"   ,    "Memory B cells"        ,     
 "Plasma cells/Fibroblasts",    "Cycling plasma cells"   ,     "Plasma cell/Granulocyte",    
 "B_cell_doublets"          ,   "Tcell/Fibroblast"        ,    "Fibroblasts_macs"     ,      
 "Granulocytes/progenitors"  ,  "CD1C_DCs"                ,    "Lamp3_DCs"               )

new.sample.ids <- c("Endothelial cells"   ,        "Plasma cells"      ,          "Pericytes",                  
  "Lymphatics"              ,    "Cycling cells" ,              "Adipose"       ,             
  "Granulocytes / progenitors",  "CD34/MFAP5"     ,             "POSTN"          ,            
 "Tcell_contam"       ,         "SOX5/CDH11"       ,           "macs"             ,          
 "LL"                  ,        "Pericytes"          ,          "SFRP/CXCL12"       ,         
 "mac/LL"               ,       "SPP1+ Macrophages"  ,         "Tcell_contam_myeloid" ,      
 "S1008A+ Monocytes"     ,      "Peri_contam_myeloid" ,        "MertK+ macrophages"    ,     
 "proliferating_myeloid"  ,     "CD1C+ cDC2s"          ,       "fib_contam_myeloid"     ,    
 "LAMP3 DCs"              ,     "B_cell_contam_myeloid" ,      "CD4+ naive/central memory T",
 "Vascular/T cell"        ,     "T cell/Fibroblast"     ,      "CD4+ KLRB1+ T"      ,        
 "FOXP3+ Tregs"           ,     "GZMK+ CD8+ T"          ,      "T cells/LL"          ,       
 "NK cells/ILCs"          ,     "T cell/Plasma cell"     ,     "pDCs"                 ,      
 "Plasma cells/Myeloid"   ,     "Vascular plasma cells"   ,    "Memory B cells"        ,     
 "Plasma cells/Fibroblasts",    "Cycling plasma cells"   ,     "Plasma cell/Granulocyte",    
 "B_cell_doublets"          ,   "Tcell/Fibroblast"        ,    "Fibroblasts_macs"     ,      
 "Granulocytes/progenitors"  ,  "CD1C_DCs"                ,    "Lamp3_DCs"               )

xenium_chrissy@meta.data[["named_sub_new2"]] <- plyr::mapvalues(x = xenium_chrissy@meta.data[["named_sub_new2"]], from = current.sample.ids, to = new.sample.ids)

xenium_chrissy$named_sub_new <- xenium_chrissy$named_sub_new2

meta_list <- list()
cols <- as.data.frame(ArchR::paletteDiscrete(xenium_chrissy@meta.data[, "named_sub_new"]))
colnames(cols)<-"colors"

ggplots <- list()
numbers <- 1:length(geoms_f)

for (i in 1:length(unique(xenium_chrissy@meta.data$orig.ident))){
   meta_list[[i]] <-   xenium_chrissy@meta.data %>% filter(orig.ident==unique(xenium_chrissy@meta.data$orig.ident)[[i]])
  # geoms_f[[i]]$cell <- gsub(pattern = "_", geoms_f[[i]]$cell, replacement="")
 #geoms_f[[i]]$cell <- paste(geoms_f[[i]]$cell, "_",numbers[[i]], sep="")
index <- match(geoms_f[[i]]$cell, rownames(meta_list[[i]]))
geoms_f[[i]]$named_new <- meta_list[[i]]$named_sub_new[index]
ggplots[[i]] <- geoms_f[[i]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal()+ scale_fill_manual(values = cols$colors)+ggtitle(unique(all_merged_f@meta.data$orig.ident)[[i]])
print(ggplots[[i]])
}



cols$colors <- c( rep("grey", 47))

rownames(cols)


cols$colors[8] <- 'blue'
cols$colors[11] <- 'yellow'
cols$colors[20] <- 'red'



table(xenium_chrissy$named_sub_new, xenium_chrissy$orig.ident) %>% as.data.frame() %>% filter(Var1=="SOX5/CDH11") %>% sort()

D2519_15_O5

D2519_15_P5

geoms_f[["D2519_15_P5"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = "black")+ scale_fill_manual(values = cols$colors) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())+NoLegend()+xlim(2100,2300)+ylim(4370,4550)




library(devtools)
library(spatula)
library(furrr)
 
xenium_chrissy_meta <- xenium_chrissy@meta.data
xenium_chrissy_meta <- xenium_chrissy_meta %>% dplyr::select("named_sub_new", "named")
xenium_chrissy_meta$cell_id <- rownames(xenium_chrissy_meta)
library(splitstackshape)

xenium_chrissy_meta <- cSplit(xenium_chrissy_meta, splitCols = "cell_id", sep="-")

xenium_chrissy_meta$cell_id_1 <- paste(xenium_chrissy_meta$cell_id_1, "-1", sep="")

index <- match(data_plot_f$cell, xenium_chrissy_meta$cell_id_1)
data_plot_f$named_sub_new <- xenium_chrissy_meta$named_sub_new[index]

table(data_plot_f$named_sub_new, useNA="ifANY")

index <- match(data_plot_f$cell, xenium_chrissy_meta$cell_id_1)
data_plot_f$named_sub_new <- xenium_chrissy_meta$named_sub_new[index]


#take list of all FOVs from baysor
data #433


table(xenium_chrissy$named_sub_new)
#split seurat obj into FOVs and extract meta.data for each one
seurat_orig <- SplitObject(xenium_chrissy, split.by = "orig.ident")

seurat_orig_meta <- list()

for(i in 1:length(seurat_orig)){
  seurat_orig_meta[[i]] <- seurat_orig[[i]]@meta.data
    seurat_orig_meta[[i]]$cell_id <- rownames(seurat_orig_meta[[i]])

  }
  
  
  
#  seurat_orig_meta[[i]] <- cSplit(seurat_orig_meta[[i]], splitCols = "cell_id", sep="-")
#seurat_orig_meta[[i]]$cell_id_1 <- paste(seurat_orig_meta[[i]]$cell_id_1, "-1", sep="")


#remove duplicates in baysor fov list

data_unique <- data

for (i in 1:length(data_unique)){
data_unique[[i]] <- data_unique[[i]][!duplicated(data_unique[[i]]$cell), ]
}

names(data_unique) <- gsub(".csv", "", tx_csv)

data_unique <- data_unique[c(fov_order)]


#add meta data from each seurat (matching)

for (i in 1:length(data_unique)){
  index <- match(data_unique[[i]]$cell, seurat_orig_meta[[i]]$cell_id)
data_unique[[i]]$named_sub_new <- seurat_orig_meta[[i]]$named_sub_new[index]
}

for (i in 1:length(data_unique)){
  index <- match(data_unique[[i]]$cell, seurat_orig_meta[[i]]$cell_id)
data_unique[[i]]$named <- seurat_orig_meta[[i]]$named[index]
}


data_unique[[100]]$named_sub_new %>% table()
seurat_orig_meta[[1]]$named_sub_new %>% table()


data_unique[[100]]$named %>% table()
seurat_orig_meta[[100]]$named %>% table()


#run proximity analysis
library(data.table)
data_unique_prox <- rbindlist(data_unique)



#need to add a column to all_x_y with the ident for each cell. in this case mine is called 'named'
 
library(devtools)
library(spatula)
library(furrr)

nrow(data_unique_prox)

data_unique_prox_noNA <- data_unique_prox[!is.na(data_unique_prox$named_sub_new),]

coloc_res_coarse_global<- coloc_all_types(
        index_type = unique(data_unique_prox_noNA$named),
        coords = data_unique_prox_noNA[, c("x", "y")],
        y = data_unique_prox_noNA$named,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )



coloc_res_coarse_sub<- coloc_all_types(
        index_type = unique(data_unique_prox_noNA$named_sub_new),
        coords = data_unique_prox_noNA[, c("x", "y")],
        y = data_unique_prox_noNA$named_sub_new,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )


library(data.table)
library(splitstackshape)
library(ggpubr)
library(rstatix)
library(DescTools)
library(tidyr)
library(tidyverse)
 
plt_df<-coloc_res_coarse_global %>%
    subset(pval < 0.05) %>%
    dplyr::select(index_type, type, zscore) %>%
    spread(type, zscore, fill = 0) %>%
    column_to_rownames('index_type') %>%
    as.matrix
 
mt2 <- plt_df +t(plt_df)

mt2%>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)




plt_df<-coloc_res_coarse_sub %>%
    subset(pval < 0.05) %>%
    dplyr::select(index_type, type, zscore) %>%
    spread(type, zscore, fill = 0) %>%
    column_to_rownames('index_type') %>%
    as.matrix



plt_df %>% as.data.frame() %>% write.csv("/rds/projects/c/croftap-mapjagb6/xenium2/proximity_subtypes.csv")


plt_df["MertK+ macrophages",] %>% sort() %>% as.data.frame()

plt_df["SOX5/CDH11",] %>% sort() %>% as.data.frame() %>% filter(. >0) %>%  tibble::rownames_to_column("VALUE") %>% filter(VALUE %in% c("MertK+ macrophages", "Granulocytes/progenitors", "SFRP/CXCL12", "LL", "POSTN", "CD34/MFAP5")) %>% ggplot(aes(x=reorder(VALUE, -.),y=.)) +
  geom_bar(stat="identity")+RotatedAxis()+coord_flip()+theme_ArchR()

 
mt2 <- plt_df +t(plt_df)

mt2%>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)


cols <- as.data.frame(ArchR::paletteDiscrete(xenium_chrissy@meta.data[, "named_sub_new"]))
colnames(cols)<-"colors"


cols$colors2 <- c("darkgreen", "grey", "red", "grey", "grey", "grey", "grey", "grey", "grey" ,"grey", "yellow", "grey", "grey", "grey", "grey" ,"grey", "grey", "grey", "grey" ,"grey", "grey", "grey", "grey" ,"grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey" ,"grey", "grey", "grey", "grey" ,"grey", "grey", "grey", "grey" ,"grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")



geoms_f[["D2518_11_O7"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal() + scale_fill_manual(values = cols$colors2)+xlim(1350,1800)+ylim(1700,2000)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+NoLegend()



geoms_f[["D2519_15_AD5"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = "black")+theme_minimal() + scale_fill_manual(values = cols$colors2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlim(2300,2400)+ylim(14850, 15000)



