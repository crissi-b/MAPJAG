


load("/rds/projects/c/croftap-celldive01/xenium/analysis/analysis_step1.RData")
rm(list=ls()[! ls() %in% c("all_tx_f")])



unique(all_tx_f$donor_FOV)[grep("D2518_08_AC6", unique(all_tx_f$donor_FOV))]


all_csv <- dir(path="/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs", pattern=".csv")
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/", tx_csv, sep="")
samples_area<-paste("/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/", cell_stat_csvs, sep="")



data = list()
for (i in 1:length(samples)) {
    data[[i]] <- read.csv(samples[[i]]);
}


#data_area = list()
#for (i in 1:length(samples_area)) {
#    data_area[[i]] <- read.csv(samples_area[[i]]);
#    data_area[[i]]$radius<-sqrt(data_area[[i]]$area/3.14)
#}

library(data.table)
data_plot<-rbindlist(data)


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
for (i in 1:length(samples)) {
  count_matrices[[i]]<-tx_to_counts(genes =data[[i]]$gene, data[[i]]$cell_id_new)
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

all_merged_f <- subset(all_merged, subset = nFeature_RNA > 10 & nCount_RNA > 10)

all_merged_f <- all_merged_f %>% 
  NormalizeData(scale.factor = median(all_merged_f$nCount_RNA)) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:40)

DimPlot(all_merged_f, raster=FALSE)+NoLegend()

table(data_plot$fov_name)

data_plot$gene["LGI4"]



#data_area_all <- rbindlist(data_area)

#hist(data_area_all$radius)


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

library(harmony)
all_merged_f <- RunHarmony(all_merged_f, group.by.vars = c("orig.ident"))
all_merged_f <- RunUMAP(all_merged_f, reduction="harmony",  dims=1:50)
all_merged_f <- FindNeighbors(all_merged_f, reduction="harmony", dims=1:50)
all_merged_f <- FindClusters(all_merged_f, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4))
all_merged_f <- FindClusters(all_merged_f, resolution = c(0.6, 0.7, 0.8))
DimPlot(all_merged_f, group.by="orig.ident", raster=FALSE)+NoLegend()
DimPlot(all_merged_f, group.by="RNA_snn_res.0.8", raster=FALSE)+NoLegend()
DimPlot(all_merged_f, group.by="RNA_snn_res.0.4", label = T, raster=FALSE)+NoLegend()


FeaturePlot(all_merged_f, features = "PECAM1", raster=FALSE)



Idents(all_merged_f) <- 'RNA_snn_res.0.1'
markers_0.1 <- FindAllMarkers(all_merged_f, only.pos = T)


Idents(all_merged_f) <- 'RNA_snn_res.0.4'
markers_0.4 <- FindAllMarkers(all_merged_f, only.pos = T)

all_merged_f$named <- all_merged_f@meta.data[["RNA_snn_res.0.4"]]
Idents(all_merged_f) <- 'named'
levels(all_merged_f)
current.sample.ids <- c( "0" , "1",  "10", "11", "12", "13", "14" ,"15", "16", "17", "18" ,"19", "2",  "20", "21", "22" ,"23", "24", "25", "3" , "4" , "5" , "6" , "7",  "8",  "9" )
new.sample.ids <- c( "macs" , "SL_fibs",  "Bcell", "SL_fibs", "Tcells", "Adipocytes", "SL_fibs" ,"monocytes", "CD19_Bcells", "SL_fibs", "SL_fibs" ,"SL_fibs", "NK_vasc?",  "macs", "SL_fibs", "SL_fibs" ,"macs", "??", "??", "Pericytes" , "Endothelial" , "Tcell" , "LL" , "Lymphatics",  "Prolif",  "Myelid_prog")

all_merged_f@meta.data[["named"]] <- plyr::mapvalues(x = all_merged_f@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)

cols <- ArchR::paletteDiscrete(all_merged_f@meta.data[, "named"])
DimPlot(all_merged_f, group.by="named", cols=cols, label=T, raster=FALSE)+NoAxes()


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

tx_08$slide <- '1'
tx_11$slide <- '2'
tx_04$slide <- '3'
tx_18_24$slide <- '4'
tx_19_15$slide <- '5'

unique(tx_08$fov_name)

all_no_f_slide <- rbind(tx_08, tx_11,tx_04, tx_18_24, tx_19_15)

all_no_f_slide <- all_no_f_slide[all_no_f_slide$transcript_id %in% all_tx_f$transcript_id,]
all_no_f_slide <- all_no_f_slide[all_no_f_slide$cell_id %in% all_tx_f$cell_id,]
all_no_f_slide <- all_no_f_slide[all_no_f_slide$feature_name %in% all_tx_f$feature_name,]

index <- match(all_no_f_slide$fov_name, all_tx_f$fov_name)
all_tx_f$slide <- all_no_f_slide$slide[index]

all_tx_f %>% filter(slide==5)  %>% ggplot(aes(x=x_location, y=y_location, color=fov_name))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()+facet_wrap(facets="slide")

cols <- as.data.frame(ArchR::paletteDiscrete(unique(all_tx_f$fov_name)))
colnames(cols)<-"colors"
all_tx_f %>% filter(fov_name== unique(tx_08$fov_name)) %>% ggplot(aes(x=x_location, y=y_location, color=fov_name))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+ scale_fill_manual(values = cols$colors)+theme_ArchR()+NoLegend()

all_tx_f %>% filter(fov_name== unique(tx_11$fov_name)) %>% ggplot(aes(x=x_location, y=y_location, color=fov_name))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+ scale_fill_manual(values = cols$colors)+theme_ArchR()+NoLegend()

all_tx_f %>% filter(fov_name == unique(tx_04$fov_name)) %>% ggplot(aes(x=x_location, y=y_location, color=fov_name))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+ scale_fill_manual(values = cols$colors)+theme_ArchR()+NoLegend()


P1 <- all_tx_f %>% filter(fov_name== unique(all_tx_f$fov_name)) %>% ggplot(aes(x=x_location, y=y_location))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()+facet_wrap("donor")

p2 <- all_tx_f %>% filter(fov_name== unique(tx_11$fov_name)) %>% ggplot(aes(x=x_location, y=y_location))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()

p3 <- all_tx_f %>% filter(fov_name == unique(tx_04$fov_name)) %>% ggplot(aes(x=x_location, y=y_location))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()


plot_grid(P1, p2, p3, ncol=3)



all_tx_f %>% filter(all_tx_f$fov_name == "O7") %>% ggplot(aes(x=x_location, y=y_location))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()+xlim(1600,1800)+ylim(1600,1800)


geoms_f[[1]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry), alpha = 0.9,
              color = "black")+theme_minimal()+
  geom_point(data=data_plot[data_plot$fov_name %in% c("AA10"),],
             aes(x=x,y=y), 
                   size=0.1)+ scale_color_brewer(palette="Dark2")+ylim(2800, 2900)


all_tx_f %>% filter(fov_name== unique(all_tx_f$fov_name)) %>% ggplot(aes(x=x_location, y=y_location, color=fov_name))+
  geom_point(alpha=0.6, shape=".")+
  coord_equal()+theme_ArchR()+NoLegend()+facet_wrap("donor")


unique(all_tx_f$donor)

all_tx_f$cell_id_new <- as.numeric(as.factor(all_tx_f$cell_id))

all_tx_f_2 <- all_tx_f %>% filter(donor != 'D2518_08')

unique(all_tx_f_2$donor_FOV) %>% length()


library(furrr)
#library(spatula)
library(sp)
library(sf)
library(sfdct)
FOVs <- list()
geoms_all <- list()
for (i in 1:length(unique(all_tx_f_2$donor_FOV))){
 FOVs[[i]] <-  all_tx_f_2 %>% filter(donor_FOV == unique(all_tx_f_2$donor_FOV)[[i]])
  FOVs[[i]]$cell <- FOVs[[i]]$cell_id_new
FOVs[[i]]$x <- FOVs[[i]]$x_location
FOVs[[i]]$y <- FOVs[[i]]$y_location
  geoms_all[[i]] <- cellgeoms_baysor(FOVs[[i]])
}

names(geoms_all) <- unique(all_tx_f_2$donor_FOV)

length(geoms_all)
length(unique(all_merged_f@meta.data$orig.ident))



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
 geoms_f[[i]]$cell <- paste(geoms_f[[i]]$cell, "_",numbers[[i]], sep="")
index <- match(geoms_f[[i]]$cell, rownames(meta_list[[i]]))
geoms_f[[i]]$named_new <- meta_list[[i]]$named[index]
ggplots[[i]] <- geoms_f[[i]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal()+ scale_fill_manual(values = cols$colors)
print(ggplots[[i]])
}








geoms_f[[1]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = named_new), alpha = 0.7,
              color = NA)+theme_minimal()+ scale_fill_manual(values = cols$colors)


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



