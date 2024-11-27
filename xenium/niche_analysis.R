

```{r}
#Step1- Voronoi tossolations

#packages and functions

library(spatula)
library(sfdct)
library(sf)
source("/rds/projects/c/croftap-mapjagx1/analysis/niches_int/repeat_clean/functions.txt")


```


```{r}

#Need geoms objs (cell boundaries for each FOV) and combine this in a list

geoms_f_f <- geoms_f
geoms_f_f <- do.call(rbind, geoms_f_f)


#Annotated Seurat Xenoum obj to obtain meta.data list
meta_list <- list()
for (i in 1:length(unique(xenium_chrissy@meta.data$orig.ident))){
   meta_list[[i]] <-   xenium_chrissy@meta.data %>% filter(orig.ident==unique(xenium_chrissy@meta.data$orig.ident)[[i]])
   geoms_f_f2 <- geoms_f_f[geoms_f_f$cell %in% rownames(meta_list[[i]]),]
  meta_list[[i]]$cell_centroid <- st_centroid(geoms_f_f2$geometry)
   meta_list[[i]]$cellID <- rownames(meta_list[[i]])
   }

names(meta_list) <- unique(xenium_chrissy@meta.data$orig.ident)


```


```{r}

# Next need tx counts fro each fov and cell (this uses baysor output)

#read in required data
all_csv <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern=".csv")
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs/", tx_csv, sep="")
rm(cell_stat_csvs)
rm(all_csv)


samplefov_list <- gsub(".csv", "",tx_csv)
tfiles_list <- samples

library(Matrix)

sample_list <- c("D2518_11", "D2518_24","D2518_28", "D2518_36",  "D2518_40","D2519_04", "D2519_15", "D2519_18")



#length(tfiles_list) prev
txfile <- list()
txfile <- lapply(seq_along(tfiles_list), function(i) {
  fread(tfiles_list[[i]])[, 
    .(x, y, z_location, gene, cell, SampleID = donor, FOV = fov_name)
  ][, 
    SampleFOV := paste(SampleID, FOV, sep = "_")
  ][, 
    `:=`(x = x * 1, y = y * 1) # Multiplication by 1 is redundant unless coercing
  ] %>% 
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE) %>% 
    setDT()
})


names <- tfiles_list
names <- gsub("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", "", names)
names <- gsub(".csv", "", names)
names <- gsub("/", "", names)
names(txfile) <- names



```


```{r}

#make sure cells that appear in all files are the same (in this case a subset of data has been taken)
txfile <- txfile[names(txfile) %in% names(meta_list)]

for (i in 1:length(txfile)){
txfile[[i]] <- txfile[[i]][txfile[[i]]$cell %in% meta_list[[i]]$cellID,]
geoms_f[[i]] <- geoms_f[[i]][geoms_f[[i]]$cell %in% txfile[[i]]$cell,]

}

#prepare column names for downstream functions

for (i in 1:length(txfile)){
  txfile[[i]]$z <- txfile[[i]]$z_location
        meta_list[[i]]$cell_center_x <- meta_list[[i]]$Pos_X
    meta_list[[i]]$cell_center_y <- meta_list[[i]]$Pos_Y
        meta_list[[i]]$bbox_tx <- meta_list[[i]]$cell_centroid
                meta_list[[i]]$SampleFOV <- meta_list[[i]]$orig.ident
                    meta_list[[i]]$SampleID <- meta_list[[i]]$donor
                              meta_list[[i]]$geometry    <-   geoms_f[[i]]$geometry
}



```


```{r}

###not used see below


voronoi_obj_all <- list() 

for (i in 1:length(txfile)){
    bbox_tx<-st_rectangle(min(txfile[[i]]$x), max(txfile[[i]]$x), min(txfile[[i]]$y), max(txfile[[i]]$y))
    message(" tessallate")
    cells_voronoi<-pts_to_voronoi(meta_list[[i]]$cell_center_x, meta_list[[i]]$cell_center_y, bbox = bbox_tx)
    # let us add cell ids to the Voronoi polygons so any changing of the rows of this dataframe is resistant to cellIDs
    cells_voronoi<-cells_voronoi %>% 
        as.data.frame %>% 
        cbind(meta_list[[i]] %>% dplyr::select(cellID, cell_centroid, SampleFOV, SampleID)) %>% 
        dplyr::rename(polygon = geometry) %>% 
        dplyr::mutate(tileID = cellID)
    # alpha shape cells
    cells_voronoi_alpha_shaped<-alpha_shape_cells(cells_voronoi)
    # generate gcmat
tfile_voronoi<-txfile[[i]][cell != 0] %>% 
        st_sf %>% 
        st_join(
            cells_voronoi_alpha_shaped %>% st_sf %>% 
                st_set_geometry("polygon") %>% 
                dplyr::select(cellID, tileID, polygon),
            .predicate = st_intersects
        ) 


#tfile_voronoi$cellID <- tfile_voronoi$cell
#tfile_voronoi$tileID <- tfile_voronoi$cell


# remove tx not belonging to any cells
    # message(" # tx not belonging to any voronoi polygon: ", data.table(tfile_voronoi_tissue)[is.na(tileID), .N])
    tfile_voronoi<-data.table(tfile_voronoi)[!is.na(tileID)]    
    gcmat_voronoi<-spatula::tx_to_counts(
        tfile_voronoi$gene, 
        tfile_voronoi$cellID, 
        remove_bg = TRUE
    )
    
    voronoi_obj <- list()
    voronoi_obj$metadata<-cells_voronoi_alpha_shaped  %>% as.data.frame
    voronoi_obj$counts<-gcmat_voronoi[, as.character(voronoi_obj$metadata$cellID)]
     #voronoi_obj$counts<-gcmat_voronoi[,colnames(gcmat_voronoi) %in% voronoi_obj$metadata$cellID]

    
    
    voronoi_obj$metadata<-voronoi_obj$metadata %>% 
        dplyr::mutate(region = "tissue")
    
    voronoi_obj_all[[i]] <- voronoi_obj

}



```



```{r}

#this worked and fast!

library(data.table)
library(sf)
library(dplyr)

process_single <- function(i) {
  bbox_tx <- st_rectangle(
    min(txfile[[i]]$x), 
    max(txfile[[i]]$x), 
    min(txfile[[i]]$y), 
    max(txfile[[i]]$y)
  )
  
  # Voronoi tessellation
  cells_voronoi <- pts_to_voronoi(
    meta_list[[i]]$cell_center_x, 
    meta_list[[i]]$cell_center_y, 
    bbox = bbox_tx
  ) %>% 
    as.data.frame() %>% 
    cbind(
      meta_list[[i]] %>% select(cellID, cell_centroid, SampleFOV, SampleID)
    ) %>% 
    rename(polygon = geometry) %>% 
    mutate(tileID = cellID)
  
  # Alpha shape for cells
  cells_voronoi_alpha_shaped <- alpha_shape_cells(cells_voronoi)
  
  # Generate gcmat
  tfile_voronoi <- st_as_sf(txfile[[i]][cell != 0]) %>% 
    st_join(
      st_as_sf(cells_voronoi_alpha_shaped) %>% 
        st_set_geometry("polygon") %>% 
        select(cellID, tileID, polygon), 
      join = st_intersects
    )
  
  # Remove transcripts not belonging to any cells
  tfile_voronoi <- as.data.table(tfile_voronoi)[!is.na(tileID)]
  
  gcmat_voronoi <- spatula::tx_to_counts(
    tfile_voronoi$gene, 
    tfile_voronoi$cellID, 
    remove_bg = TRUE
  )
  
  # Compile results into a voronoi object
  voronoi_obj <- list()
  voronoi_obj$metadata <- as.data.frame(cells_voronoi_alpha_shaped) %>% 
    mutate(region = "tissue")
  voronoi_obj$counts <- gcmat_voronoi[, as.character(voronoi_obj$metadata$cellID), with = FALSE]
  
  return(voronoi_obj)
}


library(parallel)

# Parallelize processing
num_cores <- detectCores() - 1 # Leave one core free
voronoi_obj_all <- mclapply(
  seq_along(txfile), 
  process_single, 
  mc.cores = num_cores
)


```



```{r}
library(rapportools)
library(purrr)


glass_gridded_ls<- list()


#working but proble with sample #7
  

for (i in 1:length(voronoi_obj_all)){
    #i=7
    #glass_obj_prelim<-check_glass(obj = voronoi_obj_all[[i]], polygon_col_name =  "polygon",tfile =  txfile[[i]], buffer_size = 30)
    
    
  obj=voronoi_obj_all[[i]]
    
    tfile=txfile[[i]]
    
    tissue_region<-obj$metadata %>% 
        st_sf %>% 
        st_set_geometry("polygon") %>% 
        st_union %>% 
        st_cast("MULTIPOLYGON") %>% 
        st_buffer(30)  
   
      fov_region<-st_bbox(c(xmin = min(tfile$x), xmax = max(tfile$x), ymin = min(tfile$y), ymax = max(tfile$y))) %>% 
            st_as_sfc 
    
    glass_region<-st_difference(fov_region, tissue_region) %>% 
        as.data.frame %>% 
        dplyr::mutate(region = "glass") 

    library(spatstat.geom)
    
    res<-!(spatstat.geom::is.empty(glass_region$geometry))
    obj$tissue_region<-tissue_region
    obj$glass_region<-glass_region
    obj$glass_exists<-res
    
    if(obj$glass_exists)  {
        message("Gridding glass:")
        ncounts_glass<-IQR(colSums(voronoi_obj_all[[i]]$counts)) * 2.5
        glass_gridded_df<-grid_glass(glass_region_input = obj$glass_region, 
            txfile_input = tfile , 
            cellgeoms_input = voronoi_obj_all[[i]]$metadata, 
            by_counts = TRUE, 
            ncounts = ncounts_glass
        ) 
       } else {
        message("Not gridding glass: there is no glass beyond the buffered area")
        glass_gridded_df<-NULL
    }
    if (!is.null(glass_gridded_df)) {
        # calculate properties of glass tiles
        glass_gridded_df<-glass_gridded_df %>% as.data.frame %>% 
            dplyr::rename(polygon = geometry) %>% 
            dplyr::mutate(
                polygon_centroid = st_centroid(polygon),
                polygon_area = st_area(polygon),
                SampleID = voronoi_obj_all[[i]][["metadata"]][["SampleID"]] %>% unique(),
                SampleFOV = voronoi_obj_all[[i]][["metadata"]][["SampleFOV"]]%>% unique(),
                FOV = "FOV",
                region = "glass"
            ) 
        glass_gridded_ls[[i]] <- glass_gridded_df
        
    }
    
    
    if (length(glass_gridded_ls) == 0) {
  print("The list is empty")
} else {
  print("The list is not empty")
}
    
    
} 



cellgeoms_fov <- list()    

for (i in 1:length(voronoi_obj_all)) {
    # Extract the metadata
    metadata <- voronoi_obj_all[[i]]$metadata %>%
        dplyr::select(cellID, tileID, SampleFOV, SampleID, starts_with("polygon"), region)
    
    # Check if glass_gridded_ls[[i]] is empty before combining
    if (length(glass_gridded_ls) > 0) {
        cellgeoms_fov[[i]] <- metadata %>% bind_rows(glass_gridded_ls[[i]])
    } else {
        cellgeoms_fov[[i]] <- metadata
    }
}
```




```{r}
res <- list()
fov1 <- list()
for (i in 1:length(voronoi_obj_all)){
#i=1
        message("Step 3: Build gcmat")
    # remove bg tx from tissue regions but not from the glass region
    # construct two tfiles - one for cell polygons and one for glass
    
    #cellgeoms_fov[[i]] %>% filter(region=="glass")
    #txfile[[i]]$tileID <- NULL
    #txfile[[i]]$polygon <- NULL
    
  # txfile_new <-  txfile[[i]]
   #cellgeoms_fov_new <- cellgeoms_fov[[i]]
   #txfile_tissue_glass_new <- txfile_tissue_glass
   #rm(list=ls()[! ls() %in% c("txfile_new","cellgeoms_fov_new", "txfile_tissue_glass_new")])
    
      txfile_tissue_glass<-txfile[[i]]%>% 
        st_sf %>% 
        st_join(
            cellgeoms_fov[[i]] %>%  
                st_set_geometry("polygon") %>% 
                dplyr::select(cellID, tileID, polygon, region),
            .predicate = st_intersects
        ) %>%
        as.data.table
    #txfile_tissue_glass$region <- txfile_tissue_glass$region %>% replace_na('glass')
    #txfile_tissue_glass$cellID <- txfile_tissue_glass$cell
        #txfile_tissue_glass$tileID <- txfile_tissue_glass$cellt
    
txfile_tissue_glass<-txfile_tissue_glass[!is.na(tileID)]
    txfile_tissue_glass<-txfile_tissue_glass[(cell != 0 & region == "tissue") | (region == "glass")] # select only fg tx
    
    
    #txfile_tissue_glass<-txfile_tissue_glass[!is.na(tileID)]
    # remove tx not belonging to any cells
    # message("#tx not belonging to any voronoi polygon: ", data.table(txfile_tissue_glass)[is.na(tileID), .N])
    
    message("Building gcmat")
    # generate gcmat!
    gcmat_tissue_glass<-spatula::tx_to_counts(
        txfile_tissue_glass$gene,
        txfile_tissue_glass$tileID,
        remove_bg = FALSE
    )
    message("Done")
    
    fov1$metadata<-cellgeoms_fov[[i]] %>% 
        subset(tileID %in% colnames(gcmat_tissue_glass)) %>% 
        identity
    
    
    # fov1$txfile<-txfile_tissue_glass
    fov1$counts_raw<-gcmat_tissue_glass[, as.character(fov1$metadata$tileID)]
    
    fov1$metadata<-fov1$metadata %>% dplyr::mutate(nCounts = colSums(fov1$counts_raw))
    
    fov1$metadata$FOV <- 'FOV'
    ###problem###
    #fov1$glass<-glass_gridded_ls[[i]]
    #uses updated function with try_catch to give NULL if < 5 cells
    res[[i]]<-diffuse_adjmat(fov1$metadata, fov1$counts_raw,
        lambda = 0.03, k = 10, verbose = FALSE, cut_connections = TRUE,
        prune_connections = TRUE, prune_communities = TRUE,
        community_size_thresh = 4, pxsize = 1, dist_threshold_cells = 50)
    
    
}



names(res) <- names(cellgeoms_fov)


```

```{r}


res <- list()
fov1 <- list()

for (i in 1:length(voronoi_obj_all)) {
    message("Step 3: Build gcmat")
    
    # Remove background transcripts from tissue regions but not from glass
    txfile_tissue_glass <- txfile[[i]] %>%
        st_sf() %>%
        st_join(
            cellgeoms_fov[[i]] %>%
                st_set_geometry("polygon") %>%
                dplyr::select(cellID, tileID, polygon, region),
            .predicate = st_intersects
        ) %>%
        as.data.table()

    txfile_tissue_glass <- txfile_tissue_glass[!is.na(tileID)]
    txfile_tissue_glass <- txfile_tissue_glass[
        (cell != 0 & region == "tissue") | (region == "glass")
    ]

    message("Building gcmat")
    gcmat_tissue_glass <- spatula::tx_to_counts(
        txfile_tissue_glass$gene,
        txfile_tissue_glass$tileID,
        remove_bg = FALSE
    )
    message("Done")
    
    fov1$metadata <- cellgeoms_fov[[i]] %>%
        subset(tileID %in% colnames(gcmat_tissue_glass)) %>%
        identity()
    
    # Skip processing if the number of rows in metadata is < 5
    if (nrow(fov1$metadata) < 5) {
        message("Skipping FOV ", i, " due to insufficient metadata rows (< 5).")
        next
    }

    # Processing the sample if nrow >= 5
    fov1$counts_raw <- gcmat_tissue_glass[, as.character(fov1$metadata$tileID)]
    fov1$metadata <- fov1$metadata %>%
        dplyr::mutate(nCounts = colSums(fov1$counts_raw))
    fov1$metadata$FOV <- 'FOV'
    
    res[[i]] <- diffuse_adjmat(
        fov1$metadata, fov1$counts_raw,
        lambda = 0.03, k = 10, verbose = FALSE, cut_connections = TRUE,
        prune_connections = TRUE, prune_communities = TRUE,
        community_size_thresh = 4, pxsize = 1, dist_threshold_cells = 50
    )
}





```


```{r}

#remove  empty dfs
res_f <- res[!sapply(res, is.null)]
res_f2 <- res_f[!sapply(res_f, function(x) nrow(x$metadata) == 0)]

for (i in 1:length(res_f2)){
    names(res_f2)[[i]] <-  res_f2[[i]]$metadata$SampleFOV %>% unique() 
      }


```




```{r}


library(furrr)

plots_before_pruning<-future_map(res[1:30], function(obj){
    
    graph_adj1<-igraph::graph_from_adjacency_matrix(obj$neighbors > 0, mode = "undirected")
    tidy_graph1<-tidygraph::as_tbl_graph(graph_adj1)
    polygon_coords<-st_coordinates(obj$metadata$polygon_centroid)
    p<-ggraph::ggraph(tidy_graph1, layout = polygon_coords) +
        ggraph::geom_node_point(shape = 0, alpha = 1, color = "black") +
        ggraph::geom_edge_link(edge_width = 0.1) + 
        geom_sf(data = obj$metadata, aes(geometry = polygon_centroid, color = region), fill = NA, alpha = 1) + 
        scale_color_manual(values = c("darkorange", "darkblue")) + 
        ggtitle(paste0("Adj graph before community sz pruning: ", obj$collapsemetadata$SampleFOV[1])) +
        #theme(plot.title = element_textbox_simple()) +  
      theme_ArchR()+
        NULL
}, .options = furrr_options(seed = TRUE))

  
  
  plots_after_pruning<-future_map(res[1:30], function(obj){
    
    graph_adj1<-igraph::graph_from_adjacency_matrix(obj$neighbors_collapsed[[2]] > 0, mode = "undirected")
    tidy_graph1<-tidygraph::as_tbl_graph(graph_adj1)
    polygon_coords<-st_coordinates(obj$metadata$polygon_centroid)
    p<-ggraph::ggraph(tidy_graph1, layout = polygon_coords) +
        ggraph::geom_node_point(shape = 0, alpha = 1, color = "black") +
        ggraph::geom_edge_link(edge_width = 0.1) + 
        geom_sf(data = obj$metadata, aes(geometry = polygon_centroid, color = region), fill = NA, alpha = 1) +
        scale_color_manual(values = c("darkorange", "darkblue")) + 
        ggtitle(paste0("Adj graph after community sz pruning: ", obj$metadata$SampleFOV[1])) +
        #theme(plot.title = element_textbox_simple()) +
      theme_ArchR()+
        NULL
}, .options = furrr_options(seed = TRUE))

  
  #plots_after_pruning[[i]]
  
  #for (i in 1:length(plots_after_pruning)){
  #  print(plots_after_pruning[[i]])
  #  print(plots_before_pruning[[i]])
  #}
  
  
   for (i in 1:30){
    print(plots_after_pruning[[i]])
    print(plots_before_pruning[[i]])
  }

  
  x=1
  


```






```{r}

rm(list=ls()[! ls() %in% c("res_f2","geoms_f", "geoms_f_f", "txfile", "meta_list", "meta_list_small", "xenium_chrissy")])
gc()

save.image("/rds/projects/c/croftap-mapjagx1/analysis/niches_int/repeat_clean/Step2_Clustering.Rmd.RData")
```




```{r}
source("/rds/projects/c/croftap-mapjagx1/analysis/niches_int/repeat_clean/functions.txt")

library(harmony)
library(furrr)
library(purrr)

#combine all daTA TOGETHER
obj<-list()
seurat_obj_counts <- list()
seurat_obj_raw <- list()
obj_list_k <- list()


res = res_f2

for (k in 1:10){
    obj_final<-list()
          for (i in 1:length(res)){
            collapsed_counts<-res[[i]]$counts_raw %*% t(res[[i]]$neighbors_collapsed[[k]]) 
            colnames(collapsed_counts)<-colnames(res[[i]]$counts_raw)
            obj$counts<-collapsed_counts
            obj$metadata<-res[[i]]$metadata
            obj$counts_raw<-res[[i]]$counts_raw
            obj_final[[i]] <-obj
          }
    
                obj_list_k[[k]] <- obj_final
}



```

```{r}
#Generate a list of Seurat objects from raw and collapsed counts
obj_list_k_seurat_counts <- list()
seurats_k_merged <- list()
for (k in 1:length(obj_list_k)){
  obj_list_k_seurat_counts_tmp <- list()
  for (i in 1:length(res)){
    obj_list_k_seurat_counts_tmp[[i]] <- CreateSeuratObject(counts=obj_list_k[[k]][[i]]$counts)
                      }
    obj_list_k_seurat_counts[[k]] <- obj_list_k_seurat_counts_tmp
    seurats_k_merged[[k]] <- merge(x = obj_list_k_seurat_counts[[k]][[1]], y = obj_list_k_seurat_counts[[k]][-1])

}


#crete suerat obj of all cells for each k (raw counts)
obj_list_k_seurat_RAWcounts <- list()
seurats_k_merged_raw <- list()
for (k in 1:length(obj_list_k)){
  obj_list_k_seurat_RAWcounts_tmp <- list()
  for (i in 1:length(res)){
    obj_list_k_seurat_RAWcounts_tmp[[i]] <- CreateSeuratObject(counts=obj_list_k[[k]][[i]]$counts_raw)
                      }
    obj_list_k_seurat_RAWcounts[[k]] <- obj_list_k_seurat_RAWcounts_tmp
    seurats_k_merged_raw[[k]] <- merge(x = obj_list_k_seurat_RAWcounts[[k]][[1]], y = obj_list_k_seurat_RAWcounts[[k]][-1])
}



#get metadata together
meta_data_all_k_list <- list()
for (k in 1:length(obj_list_k)){
meta_data_all <- list()
for (i in 1:length(res)){
  meta_data_all[[i]] <- obj_list_k[[k]][[i]]$metadata
   meta_data_all[[i]]$polygon <- NULL
}
  meta_data_all <- rbindlist(meta_data_all, use.names=TRUE)
meta_data_all_k_list[[k]] <- meta_data_all
}

```


```{r}

obj_ksweep <- list()
for (i in 1:length(meta_data_all_k_list)){
tmp_list <- list()
tmp_list$metadata <- meta_data_all_k_list[[i]]
tmp_list$counts <- seurats_k_merged[[i]]@assays[["RNA"]]@counts
tmp_list$counts_raw <- seurats_k_merged_raw[[i]]@assays[["RNA"]]@counts
obj_ksweep[[i]] <- tmp_list
}


```

```{r}


#need adjusting just take k=10

library(rlang)
library(harmony)
gc()
objH <- list()
for (i in 1:length(obj_ksweep)){

obj_ksweep[[i]]$metadata<-obj_ksweep[[i]]$metadata %>% 
    as.data.frame %>%
    dplyr::mutate(cellID = tileID)
objH[[i]] <- dimred_and_cluster(
  obj_ksweep[[i]], 
  do_harmony = TRUE, 
  wts = "wts",
  vars_use = c("SampleFOV", "SampleID"),
  resolution_clustering = c(0.1),
  theta = c(0,0),
  sigma = 0.2, 
  max.iter.harmony = 12,
  max.iter.cluster = 40,
  do_QC = TRUE,
  do_cluster_after = TRUE,
  do_cluster_before = TRUE,
  return_object = TRUE,
  do_umap_after = TRUE,
  do_umap_before = TRUE
)
}


```

```{r}


for (i in 1:10){
print(cbind(objH[[i]]$Humap$embedding %>% as.data.frame(),objH[[i]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=SampleFOV))+ geom_point(shape=".")+theme_ArchR()+NoLegend()+ggtitle(paste("itteration", i)))


print(cbind(objH[[i]]$Humap$embedding %>% as.data.frame(),objH[[i]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=SampleID))+ geom_point(shape=".")+theme_ArchR()+NoLegend()+ggtitle(paste("itteration", i)))

}



```




```{r}


objH[[10]]$Humap$clusters <- RunModularityClustering(
      SNN = objH[[10]]$Humap$fgraph, 
      resolution = 0.2,
      print.output = FALSE
    )


clusters <- objH[[10]]$Humap$clusters %>% as.data.frame()
colnames(clusters) <- "clust1"
clusters$clust1 %>% table()
objH[[10]]$metadata$clust1 <- as.character(clusters$clust1)


cbind(objH[[10]]$Humap$embedding %>% as.data.frame(),objH[[10]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=clust1))+ geom_point(shape=".")+theme_ArchR()+NoLegend()


```


```{r}


# Check that you get obvouse tissue morphology
unique(objH[[10]]$metadata$SampleID)


objH[[10]]$metadata %>% as.data.frame() %>% filter(SampleID == "D2518_24") %>%  ggplot() + 
  geom_sf(aes(geometry = polygon_centroid, color=clust1), alpha = 0.7) + 
            ggtitle("Tissue regions")+facet_wrap("SampleID", ncol=4) +xlim(900,3600)+ylim(5500,8200)





```






```{r}


# visulise on geoms cell object

geoms_final <- geoms_f[names(geoms_f) %in% names(res)]

meta_new <- objH[[10]]$metadata
meta_list <- list()
for (i in 1:length(unique(objH[[10]]$metadata$SampleFOV))){
  
  meta_list[[i]] <- meta_new %>% filter(SampleFOV == unique(objH[[10]]$metadata$SampleFOV)[[i]])
  
}

names(meta_list) <- names(res)



for (i in 1:length(geoms_final)){
   index <- match(geoms_final[[i]]$cell, meta_list[[i]]$cellID)
geoms_final[[i]]$clust_new <- meta_list[[i]]$clust1[index]
}

cols <- ArchR::paletteDiscrete(meta_new$clust1) %>% as.data.frame()


geoms_final[["D2518_36_Z3"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = clust_new), alpha = 0.7,
              color = "black")+ scale_fill_manual(values = cols$.) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())#+xlim(2100,2300)+ylim(4370,4550)

```
```{r}


# Choose optimal resolution
# Ofetn have to subset out a niche and subcluster

objH_final <- objH[[10]]

rm(list=ls()[! ls() %in% c("res_f2","geoms_final", "txfile", "meta_list", "xenium_chrissy", "objH_final")])

save.image("/rds/projects/c/croftap-mapjagx1/analysis/niches_int/repeat_clean/Step2_Clustering.Rmd.RData")
```


# Final cluster annotation


```{r}
library(spatula)
library(data.table)
library(sfdct)
library(sf)



process_voronoi_data <- function(txfile, meta_list) {
  voronoi_obj_all <- list()  # Initialize empty list to store results
  
  for (i in 1:length(txfile)) {
    # Generate bounding box for Voronoi tessellation
    bbox_tx <- st_rectangle(min(txfile[[i]]$x), max(txfile[[i]]$x), min(txfile[[i]]$y), max(txfile[[i]]$y))
    message("Tessellating FOV ", i)
    
    # Generate Voronoi cells for each FOV based on cell centroids
    cells_voronoi <- pts_to_voronoi(meta_list[[i]]$cell_center_x, meta_list[[i]]$cell_center_y, bbox = bbox_tx)
    
    # Add cell IDs to Voronoi polygons
    cells_voronoi <- cells_voronoi %>% 
      as.data.frame() %>% 
      cbind(meta_list[[i]] %>% dplyr::select(cellID, cell_centroid, SampleFOV, SampleID)) %>% 
      dplyr::rename(polygon = geometry) %>% 
      dplyr::mutate(tileID = cellID)
    
    # Apply alpha shape to Voronoi cells
    cells_voronoi_alpha_shaped <- alpha_shape_cells(cells_voronoi)
    
    # Generate gene count matrix (gcmat) based on Voronoi polygons
    tfile_voronoi <- txfile[[i]][cell != 0] %>% 
      st_sf() %>% 
      st_join(
        cells_voronoi_alpha_shaped %>% 
          st_sf() %>% 
          st_set_geometry("polygon") %>% 
          dplyr::select(cellID, tileID, polygon),
        .predicate = st_intersects
      ) 
    
    # Remove tx (transcripts) not belonging to any cells
    tfile_voronoi <- data.table(tfile_voronoi)[!is.na(tileID)]
    
    # Generate gene count matrix (gcmat)
    gcmat_voronoi <- spatula::tx_to_counts(
      tfile_voronoi$gene, 
      tfile_voronoi$cellID, 
      remove_bg = TRUE
    )
    
    # Prepare Voronoi object with metadata and counts
    voronoi_obj <- list()
    voronoi_obj$metadata <- cells_voronoi_alpha_shaped %>% as.data.frame()
    voronoi_obj$counts <- gcmat_voronoi[, as.character(voronoi_obj$metadata$cellID)]
    
    # Add region information to the metadata
    voronoi_obj$metadata <- voronoi_obj$metadata %>% 
      dplyr::mutate(region = "tissue")
    
    # Store the Voronoi object in the results list
    voronoi_obj_all[[i]] <- voronoi_obj
  }
  
  return(voronoi_obj_all)  # Return the list of Voronoi objects
}





library(rapportools)
library(purrr)
library(dplyr)
library(sf)
library(spatstat.geom)

process_voronoi_and_glass <- function(voronoi_obj_all, txfile) {
  
  # Initialize lists for glass gridded data and cell geometries
  glass_gridded_ls <- list()
  cellgeoms_fov <- list()
  
  # Loop through each Voronoi object
  for (i in 1:length(voronoi_obj_all)) {
    obj <- voronoi_obj_all[[i]]
    tfile <- txfile[[i]]
    
    # Step 1: Define tissue and glass regions
    tissue_region <- obj$metadata %>% 
      st_sf() %>% 
      st_set_geometry("polygon") %>% 
      st_union() %>% 
      st_cast("MULTIPOLYGON") %>% 
      st_buffer(30)
    
    fov_region <- st_bbox(c(xmin = min(tfile$x), xmax = max(tfile$x), ymin = min(tfile$y), ymax = max(tfile$y))) %>% 
      st_as_sfc()
    
    glass_region <- st_difference(fov_region, tissue_region) %>% 
      as.data.frame() %>% 
      dplyr::mutate(region = "glass")
    
    # Check if glass region exists
    res <- !(spatstat.geom::is.empty(glass_region$geometry))
    obj$tissue_region <- tissue_region
    obj$glass_region <- glass_region
    obj$glass_exists <- res
    
    # Step 2: Gridding the glass region if it exists
    if (obj$glass_exists) {
      message("Gridding glass for sample ", i, ":")
      ncounts_glass <- IQR(colSums(obj$counts)) * 2.5
      glass_gridded_df <- grid_glass(
        glass_region_input = obj$glass_region, 
        txfile_input = tfile, 
        cellgeoms_input = obj$metadata, 
        by_counts = TRUE, 
        ncounts = ncounts_glass
      )
    } else {
      message("Not gridding glass: there is no glass beyond the buffered area for sample ", i)
      glass_gridded_df <- NULL
    }
    
    # Step 3: Process the glass gridded data if it's not null
    if (!is.null(glass_gridded_df)) {
      glass_gridded_df <- glass_gridded_df %>% 
        as.data.frame() %>% 
        dplyr::rename(polygon = geometry) %>% 
        dplyr::mutate(
          polygon_centroid = st_centroid(polygon),
          polygon_area = st_area(polygon),
          SampleID = obj$metadata[["SampleID"]] %>% unique(),
          SampleFOV = obj$metadata[["SampleFOV"]] %>% unique(),
          FOV = "FOV",
          region = "glass"
        )
      
      # Append to the glass gridded list
      glass_gridded_ls[[i]] <- glass_gridded_df
    }
    
    # Step 4: Combine the gridded glass data with the metadata for each sample
    metadata <- obj$metadata %>%
      dplyr::select(cellID, tileID, SampleFOV, SampleID, starts_with("polygon"), region)
    
    if (length(glass_gridded_ls) > 0) {
      cellgeoms_fov[[i]] <- bind_rows(metadata, glass_gridded_ls[[i]])
    } else {
      cellgeoms_fov[[i]] <- metadata
    }
  }
  
  # Return both the glass gridded list and the cell geometries list
  return(list(cellgeoms_fov = cellgeoms_fov, glass_gridded_ls = glass_gridded_ls))
}






library(rapportools)
library(purrr)
library(spatula)
library(dplyr)
library(sf)
library(spatstat.geom)
library(data.table)
library(igraph)

process_voronoi_and_gcmat <- function(voronoi_obj_all, txfile, cellgeoms_fov, lambda = 0.03, k = 10) {
  res <- list()  # Initialize an empty list to store results
  
  # Loop over the Voronoi objects
  for (i in 1:length(voronoi_obj_all)) {
    message("Step 3: Build gcmat")
    
    # Remove background transcripts from tissue regions but not from glass
    txfile_tissue_glass <- txfile[[i]] %>%
      st_sf() %>%
      st_join(
        cellgeoms_fov[[i]] %>%
          st_set_geometry("polygon") %>%
          dplyr::select(cellID, tileID, polygon, region),
        .predicate = st_intersects
      ) %>%
      as.data.table()

    txfile_tissue_glass <- txfile_tissue_glass[!is.na(tileID)]
    txfile_tissue_glass <- txfile_tissue_glass[
      (cell != 0 & region == "tissue") | (region == "glass")
    ]
    
    # Build gcmat
    message("Building gcmat")
    gcmat_tissue_glass <- spatula::tx_to_counts(
      txfile_tissue_glass$gene,
      txfile_tissue_glass$tileID,
      remove_bg = FALSE
    )
    message("Done")
    
    # Set metadata for the current FOV
    fov1 <- list()
    fov1$metadata <- cellgeoms_fov[[i]] %>%
      subset(tileID %in% colnames(gcmat_tissue_glass)) %>%
      identity()
    
    # Skip processing if the number of rows in metadata is < 5
    if (nrow(fov1$metadata) < 5) {
      message("Skipping FOV ", i, " due to insufficient metadata rows (< 5).")
      next
    }
    
    # Process the sample if nrow >= 5
    fov1$counts_raw <- gcmat_tissue_glass[, as.character(fov1$metadata$tileID)]
    fov1$metadata <- fov1$metadata %>%
      dplyr::mutate(nCounts = colSums(fov1$counts_raw))
    fov1$metadata$FOV <- 'FOV'
    
    # Run diffusion-based adjacency matrix computation
    res[[i]] <- diffuse_adjmat(
      fov1$metadata, fov1$counts_raw,
      lambda = lambda, k = k, verbose = FALSE, cut_connections = TRUE,
      prune_connections = TRUE, prune_communities = TRUE,
      community_size_thresh = 4, pxsize = 1, dist_threshold_cells = 50
    )
  }
  
  # Clean up results: remove NULL and empty metadata entries
  names(res) <- names(cellgeoms_fov)
  res_f <- res[!sapply(res, is.null)]
  res_f2 <- res_f[!sapply(res_f, function(x) nrow(x$metadata) == 0)]
  
  # Update names of results with SampleFOV
  for (i in 1:length(res_f2)) {
    names(res_f2)[[i]] <- res_f2[[i]]$metadata$SampleFOV %>% unique()
  }
  
  return(res_f2)
}






process_voronoi_and_seurat <- function(res_f2, max_k = 10) {
  
  # Initialize lists to store results
  obj_list_k <- list()
  obj_list_k_seurat_counts <- list()
  seurats_k_merged <- list()
  obj_list_k_seurat_RAWcounts <- list()
  seurats_k_merged_raw <- list()
  meta_data_all_k_list <- list()
  obj_ksweep <- list()

  # Loop over the range of k (1 to max_k)
  for (k in 1:max_k) {
    
    # Step 1: Generate collapsed counts for each FOV (i)
    obj_final <- list()
    for (i in 1:length(res_f2)) {
      collapsed_counts <- res_f2[[i]]$counts_raw %*% t(res_f2[[i]]$neighbors_collapsed[[k]]) 
      colnames(collapsed_counts) <- colnames(res_f2[[i]]$counts_raw)
      
      # Store counts, metadata, and raw counts in the obj list
      obj <- list()
      obj$counts <- collapsed_counts
      obj$metadata <- res_f2[[i]]$metadata
      obj$counts_raw <- res_f2[[i]]$counts_raw
      
      obj_final[[i]] <- obj
    }
    
    # Add obj_final for each k to obj_list_k
    obj_list_k[[k]] <- obj_final
    
    # Step 2: Generate Seurat objects for collapsed counts
    obj_list_k_seurat_counts_tmp <- list()
    for (i in 1:length(res_f2)) {
      obj_list_k_seurat_counts_tmp[[i]] <- CreateSeuratObject(counts=obj_list_k[[k]][[i]]$counts)
    }
    obj_list_k_seurat_counts[[k]] <- obj_list_k_seurat_counts_tmp
    seurats_k_merged[[k]] <- merge(x = obj_list_k_seurat_counts[[k]][[1]], y = obj_list_k_seurat_counts[[k]][-1])
    
    # Step 3: Generate Seurat objects for raw counts
    obj_list_k_seurat_RAWcounts_tmp <- list()
    for (i in 1:length(res_f2)) {
      obj_list_k_seurat_RAWcounts_tmp[[i]] <- CreateSeuratObject(counts=obj_list_k[[k]][[i]]$counts_raw)
    }
    obj_list_k_seurat_RAWcounts[[k]] <- obj_list_k_seurat_RAWcounts_tmp
    seurats_k_merged_raw[[k]] <- merge(x = obj_list_k_seurat_RAWcounts[[k]][[1]], y = obj_list_k_seurat_RAWcounts[[k]][-1])
    
    # Step 4: Combine metadata for each k
    meta_data_all <- list()
    for (i in 1:length(res_f2)) {
      meta_data_all[[i]] <- obj_list_k[[k]][[i]]$metadata
      meta_data_all[[i]]$polygon <- NULL  # Remove the polygon column
    }
    meta_data_all <- rbindlist(meta_data_all, use.names = TRUE)
    meta_data_all_k_list[[k]] <- meta_data_all
  }

  # Step 5: Generate obj_ksweep with metadata and counts
  for (i in 1:length(meta_data_all_k_list)) {
    tmp_list <- list()
    tmp_list$metadata <- meta_data_all_k_list[[i]]
    tmp_list$counts <- seurats_k_merged[[i]]@assays[["RNA"]]@counts
    tmp_list$counts_raw <- seurats_k_merged_raw[[i]]@assays[["RNA"]]@counts
    obj_ksweep[[i]] <- tmp_list
  }
  
  # Return only obj_ksweep
  return(obj_ksweep)
}




```


```{r}
# Next need tx counts fro each fov and cell (this uses baysor output)

#read in required data
all_csv <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern=".csv")
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs/", tx_csv, sep="")
rm(cell_stat_csvs)
rm(all_csv)


samplefov_list <- gsub(".csv", "",tx_csv)
tfiles_list <- samples

library(Matrix)

sample_list <- c("D2518_11", "D2518_24","D2518_28", "D2518_36",  "D2518_40","D2519_04", "D2519_15", "D2519_18")



#length(tfiles_list) prev
txfile <- list()
txfile <- lapply(seq_along(tfiles_list), function(i) {
  fread(tfiles_list[[i]])[, 
    .(x, y, z_location, gene, cell, SampleID = donor, FOV = fov_name)
  ][, 
    SampleFOV := paste(SampleID, FOV, sep = "_")
  ][, 
    `:=`(x = x * 1, y = y * 1) # Multiplication by 1 is redundant unless coercing
  ] %>% 
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE) %>% 
    setDT()
})


names <- tfiles_list
names <- gsub("/rds/projects/c/croftap-mapjagx1/analysis/analysis_initial_setup/baysor_outs", "", names)
names <- gsub(".csv", "", names)
names <- gsub("/", "", names)
names(txfile) <- names


library(dplyr)
#Annotated Seurat Xenoum obj to obtain meta.data list
meta_list <- list()
for (i in 1:length(unique(xenium_chrissy@meta.data$orig.ident))){
   meta_list[[i]] <-   xenium_chrissy@meta.data %>% filter(orig.ident==unique(xenium_chrissy@meta.data$orig.ident)[[i]])
   geoms_f_f2 <- geoms_f_f[geoms_f_f$cell %in% rownames(meta_list[[i]]),]
  meta_list[[i]]$cell_centroid <- st_centroid(geoms_f_f2$geometry)
   meta_list[[i]]$cellID <- rownames(meta_list[[i]])
   }

names(meta_list) <- unique(xenium_chrissy@meta.data$orig.ident)


#make sure cells that appear in all files are the same (in this case a subset of data has been taken)
txfile <- txfile[names(txfile) %in% names(meta_list)]

for (i in 1:length(txfile)){
txfile[[i]] <- txfile[[i]][txfile[[i]]$cell %in% meta_list[[i]]$cellID,]
geoms_f[[i]] <- geoms_f[[i]][geoms_f[[i]]$cell %in% txfile[[i]]$cell,]

}

#prepare column names for downstream functions

for (i in 1:length(txfile)){
  txfile[[i]]$z <- txfile[[i]]$z_location
        meta_list[[i]]$cell_center_x <- meta_list[[i]]$Pos_X
    meta_list[[i]]$cell_center_y <- meta_list[[i]]$Pos_Y
        meta_list[[i]]$bbox_tx <- meta_list[[i]]$cell_centroid
                meta_list[[i]]$SampleFOV <- meta_list[[i]]$orig.ident
                    meta_list[[i]]$SampleID <- meta_list[[i]]$donor
                              meta_list[[i]]$geometry    <-   geoms_f[[i]]$geometry
}


ver <- process_voronoi_data(txfile, meta_list)
glass_gridded_ls <- process_voronoi_and_glass(ver, txfile)
res <- process_voronoi_and_gcmat(ver, txfile, glass_gridded_ls$cellgeoms_fov)
obj_ksweep <- process_voronoi_and_seurat(res)
```


```{r}
# Harmonisation and processing

library(rlang)
library(harmony)
gc()
objH <- list()
for (i in 1:length(obj_ksweep)){

obj_ksweep[[i]]$metadata<-obj_ksweep[[i]]$metadata %>% 
    as.data.frame %>%
    dplyr::mutate(cellID = tileID)
objH[[i]] <- dimred_and_cluster(
  obj_ksweep[[i]], 
  do_harmony = TRUE, 
  wts = "wts",
  vars_use = c("SampleFOV", "SampleID"),
  resolution_clustering = c(0.1),
  theta = c(0,0),
  sigma = 0.2, 
  max.iter.harmony = 12,
  max.iter.cluster = 40,
  do_QC = TRUE,
  do_cluster_after = TRUE,
  do_cluster_before = TRUE,
  return_object = TRUE,
  do_umap_after = TRUE,
  do_umap_before = TRUE
)
}


#plotting
for (i in 1:10){
print(cbind(objH[[i]]$Humap$embedding %>% as.data.frame(),objH[[i]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=SampleFOV))+ geom_point(shape=".")+theme_ArchR()+NoLegend()+ggtitle(paste("itteration", i)))
```


```{r}
# re clustering  
objH[[10]]$Humap$clusters <- RunModularityClustering(
      SNN = objH[[10]]$Humap$fgraph, 
      resolution = 0.2,
      print.output = FALSE
    )

# Further plotting
clusters <- objH[[10]]$Humap$clusters %>% as.data.frame()
colnames(clusters) <- "clust1"
clusters$clust1 %>% table()
objH[[10]]$metadata$clust1 <- as.character(clusters$clust1)


cbind(objH[[10]]$Humap$embedding %>% as.data.frame(),objH[[10]]$metadata) %>% ggplot(aes(x=V1, y=V2, color=clust1))+ geom_point(shape=".")+theme_ArchR()+NoLegend()
  
  
  

```
```{r}


# Check that you get obvouse tissue morphology
unique(objH[[10]]$metadata$SampleID)


objH[[10]]$metadata %>% as.data.frame() %>% filter(SampleID == "D2518_24") %>%  ggplot() + 
  geom_sf(aes(geometry = polygon_centroid, color=clust1), alpha = 0.7) + 
            ggtitle("Tissue regions")+facet_wrap("SampleID", ncol=4) +xlim(900,3600)+ylim(5500,8200)





```






```{r}


# visulise on geoms cell object

geoms_final <- geoms_f[names(geoms_f) %in% names(res)]

meta_new <- objH[[10]]$metadata
meta_list <- list()
for (i in 1:length(unique(objH[[10]]$metadata$SampleFOV))){
  
  meta_list[[i]] <- meta_new %>% filter(SampleFOV == unique(objH[[10]]$metadata$SampleFOV)[[i]])
  
}

names(meta_list) <- names(res)



for (i in 1:length(geoms_final)){
   index <- match(geoms_final[[i]]$cell, meta_list[[i]]$cellID)
geoms_final[[i]]$clust_new <- meta_list[[i]]$clust1[index]
}

cols <- ArchR::paletteDiscrete(meta_new$clust1) %>% as.data.frame()


geoms_final[["D2518_36_Z3"]] %>% as.data.frame %>% 
    ggplot()+
     geom_sf(aes(geometry = geometry, fill = clust_new), alpha = 0.7,
              color = "black")+ scale_fill_manual(values = cols$.) +theme(panel.background = element_rect(fill='white'))+
  theme(axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())#+xlim(2100,2300)+ylim(4370,4550)

```



                        
