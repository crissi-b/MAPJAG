
#rm(list=ls()[! ls() %in% c("xenium_chrissy","sce_all", "geoms_f")])

source("/rds/projects/c/croftap-celldive01/xenium/functions_roopha/TissueSegFunctions.R")
source("/rds/projects/c/croftap-celldive01/xenium/functions_roopha/workflows.R")
source("/rds/projects/c/croftap-celldive01/xenium/functions_roopha/VizFunctions.R")


plot_dim_red <- function(dim_red_embeddings, clusters = NULL, metadata, cell_id_colname = "cellID", 
                         color_by, plot_title = NULL, dim_red_type, legend_posn = "right", 
                         shape_points = ".", size_points = 0.1, alpha_value = 1, point_type = "size", 
                         legend_title_input = NULL, plot_labels = FALSE, shorten_labels = FALSE) {
    # Check if the color_by column exists in metadata
    if (!color_by %in% colnames(metadata)) {
        stop(paste("The column", color_by, "is not found in the metadata"))
    }
    
    if (is.null(plot_title)) 
        plot_title <- dim_red_type
    if (is.null(clusters)) 
        clusters <- rep(1, nrow(dim_red_embeddings)) %>% as.data.table()
    
    dim_red_embeddings <- dim_red_embeddings[, 1:2] %>% as.data.table()
    setnames(dim_red_embeddings, c("V1", "V2"))
    
    # Rename columns in metadata to avoid duplication
    metadata <- copy(metadata)
    common_cols <- intersect(names(metadata), names(dim_red_embeddings))
    if (length(common_cols) > 0) {
        setnames(metadata, old = common_cols, new = paste0(common_cols, "_meta"))
    }
    
    plt_df <- cbind(metadata, dim_red_embeddings)
    plt_df[, clusters := clusters]
    plt_df[, color_col := get(color_by)]
    
    plt_df <- plt_df[, `:=`(x_mid = median(V1), y_mid = median(V2), label_col = color_col), by = color_col][sample(.N, .N * 1)]
    cols_plot <- generate_colors_tableau(plt_df$color_col)
    
    if (is.null(legend_title_input)) 
        legend_title_input <- color_by
    
    if (plot_labels & shorten_labels) {
        plt_df[, label_col := gsub("(.*?):.*", "\\1", color_col)]
    }
    
    p <- ggplot(data = plt_df, aes(V1, V2, color = color_col, alpha = alpha_value)) + 
        {
            if (point_type == "size") 
                geom_point(size = size_points, alpha = alpha_value)
        } + 
        {
            if (point_type == "shape") 
                geom_point(shape = shape_points, alpha = alpha_value)
        } + 
        ggtitle(plot_title) + 
        scale_color_manual(values = cols_plot) + 
        {
            if (plot_labels) 
                geom_label_repel(data = unique(plt_df[, .(x_mid, y_mid, label_col, color_col)]), 
                                 aes(x = x_mid, y = y_mid, label = label_col, fill = color_col), 
                                 color = "black", size = 5, 
                                 min.segment.length = 0, max.overlaps = Inf)
        } + 
        {
            if (plot_labels) 
                scale_fill_manual(values = cols_plot)
        } + 
        {
            if (plot_labels) 
                guides(fill = guide_legend(title = legend_title_input, 
                                           override.aes = list(color = NA)))
        } + 
        {
            if (!plot_labels) 
                guides(color = guide_legend(title = legend_title_input, 
                                            override.aes = list(shape = 16, alpha = 1, size = 5)))
        } + 
        xlab(paste0(dim_red_type, "_1")) + 
        ylab(paste0(dim_red_type, "_2")) + 
        theme(plot.title = element_text(margin = margin(10, 0, 10, 0)), 
              legend.position = legend_posn, 
              legend.text = element_text(margin = margin(5, 5, 5, 5))) + 
        NULL
    
    return(p)
}

# Dummy function to generate colors (as the original code seems to refer to an undefined function)
generate_colors_tableau <- function(labels) {
    colors <- rainbow(length(unique(labels)))
    names(colors) <- unique(labels)
    return(colors)
}



#read in required data
all_csv <- dir(path="/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs", pattern=".csv")
cell_stat_csvs <- dir(path="/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs", pattern="stats.csv")
tx_csv <- all_csv[!all_csv %in% cell_stat_csvs]
samples<-paste("/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/", tx_csv, sep="")
rm(cell_stat_csvs)
rm(all_csv)


samplefov_list <- gsub(".csv", "",tx_csv)
tfiles_list <- samples

library(Matrix)

sample_list <- c("D2518_11", "D2518_24","D2518_28", "D2518_36",  "D2518_40","D2519_04", "D2519_15", "D2519_18")

xenium_chrissy$Pos_Y <- sce_all$Pos_Y
xenium_chrissy$Pos_X <- sce_all$Pos_X
geoms_f_f <- geoms_f

library(spatula)
library(sfdct)
library(sf)
meta_list <- list()
for (i in 1:length(unique(xenium_chrissy@meta.data$orig.ident))){
   meta_list[[i]] <-   xenium_chrissy@meta.data %>% filter(orig.ident==unique(xenium_chrissy@meta.data$orig.ident)[[i]])
   geoms_f_f[[i]] <- geoms_f_f[[i]][geoms_f_f[[i]]$cell %in% rownames(meta_list[[i]]),]
  meta_list[[i]]$cell_centroid <- st_centroid(geoms_f_f[[i]]$geometry)
   meta_list[[i]]$cellID <- rownames(meta_list[[i]])
   }

library(data.table)
metadata <- rbindlist(meta_list)
names(meta_list) <- (unique(xenium_chrissy@meta.data$orig.ident))

#generate txfile

fov1<-list()
txfile <- list()
for (i in 1:length(tfiles_list)){
  #samplefov_name<-gsub(".*/(D[0-9]+_[0-9]+)_.*\\.csv$", "\\1", fov_data_list[[1]]$tfile_path)
  txfile[[i]]<-fread(tfiles_list[[i]])[, .(x, y, z_location, gene, cell, SampleID = donor, FOV = fov_name)
        ][, `:=` (SampleFOV = paste(SampleID, FOV, sep = "_"))
        ][, `:=` (x = x * 1, y = y * 1)
        ] %>% 
        sf::st_as_sf(coords = c("x", "y"), remove = FALSE) %>% 
        setDT()
  
}


names <- tfiles_list
names <- gsub("/rds/projects/c/croftap-celldive01/xenium/analysis/baysor_outs/", "", names)
names <- gsub(".csv", "", names)
names(txfile) <- names


for (i in 1:length(txfile)){
      message("Step1: tessallate/tile and build gcmat")
  txfile[[i]]$z <- txfile[[i]]$z_location
        meta_list[[i]]$cell_center_x <- meta_list[[i]]$Pos_X
    meta_list[[i]]$cell_center_y <- meta_list[[i]]$Pos_Y
        meta_list[[i]]$bbox_tx <- meta_list[[i]]$cell_centroid
                meta_list[[i]]$SampleFOV <- meta_list[[i]]$orig.ident
                    meta_list[[i]]$SampleID <- meta_list[[i]]$donor
                              meta_list[[i]]$geometry    <-   geoms_f_f[[i]]$geometry
}
meta_list[[43]] <- NULL
txfile <- txfile[names(txfile)%in%names(meta_list)]

meta_list[[85]] <- NULL
txfile <- txfile[names(txfile)%in%names(meta_list)]


txfile <- Filter(function(x) nrow(x) > 1, txfile)
meta_list <- Filter(function(x) nrow(x) > 1, meta_list)
txfile <- txfile[names(txfile)%in%names(meta_list)]

#remove 132

source("~/pkgs.R")
library(spatula)
library(sfdct)
library(sf)

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
            txfile_input = txfile[[i]] , 
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
}    

glass_gridded_ls[[426]] <- "end"   

save.image("/rds/projects/c/croftap-mapjagb6/xenium2/niche_roopha/analysis_full.RData")

x=1

names(voronoi_obj_all) <- names(meta_list)
names(glass_gridded_ls) <- names(meta_list)


meta_list <- Filter(function(x) nrow(x) > 3, meta_list)
txfile <- txfile[names(txfile)%in%names(meta_list)]
voronoi_obj_all <- voronoi_obj_all[names(voronoi_obj_all)%in%names(meta_list)]
glass_gridded_ls <- glass_gridded_ls[names(glass_gridded_ls)%in%names(meta_list)]


cellgeoms_fov <- list()    
for (i in 1:length(voronoi_obj_all)){
cellgeoms_fov[[i]]<-voronoi_obj_all[[i]]$metadata %>%  dplyr::select(cellID, tileID,SampleFOV, SampleID, starts_with("polygon"), region) %>% 
        bind_rows(glass_gridded_ls[[i]])
}


res <- list()
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
    fov1$glass<-glass_gridded_ls[[i]]
    #uses updated function with try_catch to give NULL if < 5 cells
    res[[i]]<-diffuse_adjmat(fov1$metadata, fov1$counts_raw,
        lambda = 0.03, k = 10, verbose = FALSE, cut_connections = TRUE,
        prune_connections = TRUE, prune_communities = TRUE,
        community_size_thresh = 4, pxsize = 1, dist_threshold_cells = 50)
    
    
}

save.image("/rds/projects/c/croftap-mapjagb6/xenium2/niche_roopha/analysis_full.RData")

names(res) <- names(cellgeoms_fov)


library(furrr)

plots_before_pruning<-future_map(res[1:5], function(obj){
    
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

  
  
  plots_after_pruning<-future_map(res[1:5], function(obj){
    
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

  
  plots_after_pruning[[1]]
  
  for (i in 1:length(plots_after_pruning)){
    print(plots_after_pruning[[i]])
    print(plots_before_pruning[[i]])
  }
  

#roophas
library(purrr)
.tmp <- map(res_safe, ~ rownames(.x$counts_raw))
length(.tmp)


gene_panel<-furrr::future_map(res, function(.x){
    library(Matrix)
    res<-rownames(.x$counts_raw)
    return(res)
    }, .options = furrr_options(seed = TRUE)
) 

gene_panel2<-reduce(.tmp[1:2], intersect)
gene_panel2 %>% unique() %>% length()

intersect(rownames(res_safe[[1]]$counts_raw),rownames(res_safe[[2]]$counts_raw)) %>% length()

library(furrr)



genes<-res %>% map('counts') %>% map(rownames) %>% reduce(union)
counts<-res[1] %>%
    map('counts') %>%
    map(function(X) {
        genes_missing <- setdiff(genes, rownames(X))
        filler <- Matrix(nrow = length(genes_missing), ncol = ncol(X), data = 0)
        rownames(filler) <- genes_missing
        colnames(filler) <- colnames(X)
        Matrix::rbind2(X, filler)[genes, ]
    }) %>%
    reduce(Matrix::cbind2)

#send roopha code and rds of res
 map(res, ~ is.null(.x$counts_raw)) %>% reduce(`+`)

obj_ksweep<-map(1:3, function(k){

    obj_final<-map(res[1:27], function(path_obj_fov){
        tryCatch({
            obj_samplefov<-path_obj_fov
            obj<-list()
            common_genes<-intersect(gene_panel2, rownames(res$counts_raw))
            
            collapsed_counts<-obj_samplefov$counts_raw[common_genes, ] %*% t(obj_samplefov$neighbors_collapsed[[k]])
            colnames(collapsed_counts)<-colnames(obj_samplefov$counts_raw)
            
            #genes_missing <- setdiff(genes, rownames(X))
        #filler <- Matrix(nrow = length(genes_missing), ncol = ncol(X), data = 0)
        #rownames(filler) <- genes_missing
        #colnames(filler) <- colnames(X)
        #obj_samplefov$counts_raw <- Matrix::rbind2(obj_samplefov$counts_raw, filler)[genes, ]
  
            
            obj$counts<-collapsed_counts
            obj$metadata<-obj_samplefov$metadata
            obj$counts_raw<-obj_samplefov$counts_raw[gene_panel2, ]
            # print("here")
            return(obj)
            
        }, error = function(e){
            message(path_obj_fov)
            message(e)
            return(NULL)
        })
        
    })
})


intersect(res_safe[[3]]$counts_raw %>% rownames, gene_panel2)
res_safe[[3]]$counts_raw[intersect(gene_panel2, rownames(res_safe[[3]]$counts_raw)), 1:5]


    # print(obj_final[[1]]$metadata %>% head)
    obj<-list()
    obj$metadata<-map(obj_final, ~ .x$metadata) %>%  bind_rows()
    obj$counts<-map(obj_final, ~.x$counts) %>% reduce(cbind)
    obj$counts_raw<-map(obj_final, ~.x$counts_raw) %>% reduce(cbind)
    #saveRDS(obj, paste0("../dataFinal/cache/tissueSegmentation/voronoiObj/diffused_mat_sample_k_", k, ".RDS"))
    return(obj)
}) 




res_safe <- res
res <- res_safe
res = res[-which(sapply(res, is.null))]


#now make a list of serat objs for each res at each k for raw and normal counts
#combine for each k
#then extract counts/metadat for each one
obj_list_k <- list()
  obj<-list()
obj_list_k_seurat_counts_tmp_k1 <- list()
obj_list_k_seurat_RAWcounts_tmp_k1 <- list()
for (i in 1:length(res)){
          collapsed_counts<-res[[i]]$counts_raw %*% t(res[[i]]$neighbors_collapsed[[1]]) 
            colnames(collapsed_counts)<-colnames(res[[i]]$counts_raw)
            obj$counts<-collapsed_counts
            obj_list_k[[i]] <- obj
    obj_list_k_seurat_counts_tmp_k1[[i]] <- CreateSeuratObject(counts=obj_list_k[[i]]$counts %>% as.matrix(), min.cells = 0, min.features = 0)
        obj_list_k_seurat_RAWcounts_tmp_k1[[i]] <- CreateSeuratObject(counts=obj_list_k[[i]]$counts_raw, min.cells = 0, min.features = 0)
                      }


obj_list_k_seurat_counts <- list()

for (k in 1:length(obj_list_k)){
  obj_list_k_seurat_counts_tmp <- list()
  for (i in 1:length(res)){
    obj_list_k_seurat_counts_tmp <- CreateSeuratObject(counts=obj_list_k[[k]][[i]]$counts)
    
                      }
    obj_list_k_seurat_counts[[k]] <- obj_list_k_seurat_counts_tmp
}


obj_list_k_seurat_RAWcounts <- list()

for (k in 1:length(obj_list_k)){
  obj_list_k_seurat_counts_tmp <- list()
  for (i in 1:length(res)){
    obj_list_k_seurat_counts_tmp <- CreateSeuratObject(counts=obj_list_k[[k]][[i]]$counts_raw%>% as.data.frame())
            
          }
    obj_list_k_seurat_RAWcounts[[k]] <- obj_list_k_seurat_counts_tmp
}

obj_tmp <- list()
    obj_tmp$metadata<-map(obj_final, ~ res[[i]]$metadata) %>%  bind_rows()
    obj_tmp$counts<-map(obj_final, ~ res[[i]]$counts) %>% reduce(cbind)
    obj_tmp$counts_raw<-map(obj_final, ~ res[[i]]$counts_raw) %>% reduce(cbind)
        obj_ksweep[[k]] <-  obj_tmp   
}

library(furrr)
library(dplyr)
library(Matrix)

obj_ksweep <- list()

for (k in 1:10) {
  obj_final <- list()
  
  for (i in seq_along(res)) {
    path_obj_fov <- res[[i]]
    tryCatch({
      obj_samplefov <- path_obj_fov
      obj <- list()
      collapsed_counts <- obj_samplefov$counts_raw[gene_panel2, ] %*% t(obj_samplefov$neighbors_collapsed[[k]])
      colnames(collapsed_counts) <- colnames(obj_samplefov$counts_raw)
      
      obj$counts <- collapsed_counts
      obj$metadata <- obj_samplefov$metadata
      obj$counts_raw <- obj_samplefov$counts_raw[gene_panel2, ]
      
      obj_final <- c(obj_final, list(obj))
    }, error = function(e) {
      message(path_obj_fov)
      message(e)
    })
  }

  obj_final <- obj_final[!sapply(obj_final, is.null)]
  
  obj <- list()
  obj$metadata <- bind_rows(obj_final[[1]]$metadata, obj_final[[2]]$metadata)
  obj$counts <- Reduce(cbind, obj_final[[1]]$counts, obj_final[[2]]$counts)
  obj$counts_raw <- Reduce(cbind, obj_final[[1]]$counts_raw, obj_final[[2]]$counts_raw)
  
  obj_ksweep[[k]] <- obj
}

objH <- list()
for (i in 1:2){

obj_ksweep[[i]]$metadata<-obj_ksweep[[i]]$metadata %>% 
    as.data.frame %>%
    dplyr::mutate(cellID = tileID)
objH[[i]] <- dimred_and_cluster(
  obj_ksweep[[i]], 
  do_harmony = TRUE, 
  wts = "wts",
  vars_use = c("SampleFOV"), #should also include sample here as well if across multiple smaples
  resolution_clustering = c(0.1, 0.5, 0.7, 1, 2),
  theta = c(0),
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




find_common_genes <- function(matrices) {
    if (length(matrices[[1]][["counts"]]) == 0) return(character(0))
    
    # Start with the gene list of the first matrix
    common_genes <- rownames(matrices[[1]][["counts"]])
    
    # Iterate over the rest of the matrices
    for (i in 2:length(matrices)) {
        common_genes <- intersect(common_genes, rownames(matrices[[i]][["counts"]]))
        if (length(common_genes) == 0) break  # No common genes left
    }
    
    return(common_genes)
}

gene_panel <- find_common_genes(res)


library(furrr)
library(purrr)
obj_ksweep<-map(1:10, function(k){

    obj_final<-furrr::future_map(as.list(res), function(path_obj_fov){
        tryCatch({
          
            obj_samplefov<-path_obj_fov
            obj<-list()
            collapsed_counts<-obj_samplefov$counts_raw %*% t(obj_samplefov$neighbors_collapsed[[k]]) #change to k
            colnames(collapsed_counts)<-colnames(obj_samplefov$counts_raw)
            
            obj$counts<-collapsed_counts
            obj$metadata<-obj_samplefov$metadata
            obj$counts_raw<-obj_samplefov$counts_raw
            # print("here")
            return(obj)
            
        }, error = function(e){
            message(path_obj_fov)
            message(e)
            return(NULL)
        })
        
    }, .options = furrr_options(seed = TRUE))
    # print(obj_final[[1]]$metadata %>% head)
    obj<-list()
    obj$metadata<-map(obj_final, ~ .x$metadata) %>%  bind_rows()
    obj$counts<-map(obj_final, ~.x$counts) %>% reduce(cbind)
    obj$counts_raw<-map(obj_final, ~.x$counts_raw) %>% reduce(cbind)
    #saveRDS(obj, paste0("../dataFinal/cache/tissueSegmentation/voronoiObj/diffused_mat_sample_k_", k, ".RDS"))
    return(obj)
}) 

save.image("/rds/projects/c/croftap-mapjagb6/xenium2/niche_roopha/analysis_full.RData")



library(furrr)
library(purrr)
library(dplyr)
library(Matrix)

obj_ksweep <- map(1:10, function(k) {
    obj_final <- furrr::future_map(as.list(res), function(path_obj_fov) {
        tryCatch({
            obj_samplefov <- path_obj_fov
            obj <- list()
            collapsed_counts <- obj_samplefov$counts_raw %*% t(obj_samplefov$neighbors_collapsed[[k]]) #change to k
            colnames(collapsed_counts) <- colnames(obj_samplefov$counts_raw)
            
            obj$counts <- collapsed_counts
            obj$metadata <- obj_samplefov$metadata
            obj$counts_raw <- obj_samplefov$counts_raw
            return(obj)
            
        }, error = function(e) {
            message(path_obj_fov)
            message(e)
            return(NULL)
        })
        
    }, .options = furrr_options(seed = TRUE))
    
    # Filter out NULL values
    obj_final <- obj_final[!sapply(obj_final, is.null)]
    
    # Determine the union of all gene names
    all_genes <- obj_final %>%
        map(~ rownames(.x$counts)) %>%
        reduce(union)
    
    # Fill in missing genes with zeros for each matrix
    obj_final <- map(obj_final, function(x) {
        missing_genes <- setdiff(all_genes, rownames(x$counts))
        if (length(missing_genes) > 0) {
            missing_counts <- matrix(0, nrow = length(missing_genes), ncol = ncol(x$counts))
            colnames(missing_counts) <- colnames(x$counts)
            rownames(missing_counts) <- missing_genes
            x$counts <- rbind(x$counts, missing_counts)
            x$counts_raw <- rbind(x$counts_raw, missing_counts)
            x$metadata <- rbind(x$metadata, missing_counts)
        }
        return(x)
    })
    
    # Combine matrices
    obj <- list()
    obj$metadata <- map(obj_final, ~ .x$metadata) %>% bind_rows()
    obj$counts <- map(obj_final, ~ .x$counts) %>% reduce(cbind)
    obj$counts_raw <- map(obj_final, ~ .x$counts_raw) %>% reduce(cbind)
    
    return(obj)
})


save.image("/rds/projects/c/croftap-mapjagb6/xenium2/niche_roopha/analysis_full.RData")


#seems to be working but stalls?
library(furrr)
library(purrr)
library(dplyr)
library(Matrix)

reorder_fill_genes_list <- function(matrices, all_genes) {
  lapply(matrices, function(mat) {
    mat <- as.matrix(mat)
    
    mat_full <- matrix(0, nrow = length(all_genes), ncol = ncol(mat))
    rownames(mat_full) <- all_genes
    
    if (is.null(rownames(mat))) {
      stop("Matrix rownames are NULL")
    }
    
    common_genes <- intersect(rownames(mat), all_genes)
    common_indices <- match(common_genes, all_genes)
    
    #message("Dimensions of mat: ", paste(dim(mat), collapse = "x"))
    #message("Dimensions of mat_full: ", paste(dim(mat_full), collapse = "x"))
    #message("Number of common genes: ", length(common_genes))
    #message("First few common genes: ", paste(head(common_genes), collapse = ", "))
    #message("Length of common_indices: ", length(common_indices))
    #message("Length of common_genes: ", length(common_genes))
    #message("Dimensions before assignment: mat = ", paste(dim(mat[common_genes, , drop = FALSE]), collapse = "x"))
    #message("Dimensions before assignment: mat_full[common_indices, ] = ", paste(dim(mat_full[common_indices, ]), collapse = "x"))
    
    mat_common_subset <- mat[common_genes, , drop = FALSE]
    mat_full_subset <- mat_full[common_indices, ]
    #message("Dimensions of mat_common_subset: ", paste(dim(mat_common_subset), collapse = "x"))
    #message("Dimensions of mat_full_subset: ", paste(dim(mat_full_subset), collapse = "x"))
    
    if (nrow(mat_common_subset) == length(common_indices) && ncol(mat_common_subset) == ncol(mat_full_subset)) {
      mat_full[common_indices, ] <- mat_common_subset
      message("Assignment completed successfully.")
    } else {
      message("Mismatch in dimensions between mat_common_subset and mat_full_subset.")
    }
    
    return(mat_full)
  })
}

pad_matrix_list <- function(matrices) {
  max_rows <- max(sapply(matrices, nrow))
  max_cols <- max(sapply(matrices, ncol))
  
  padded_matrices <- lapply(matrices, function(mat) {
    if (nrow(mat) < max_rows) {
      pad_rows <- max_rows - nrow(mat)
      pad_matrix <- matrix(0, nrow = pad_rows, ncol = ncol(mat))
      mat <- rbind(mat, pad_matrix)
    }
    if (ncol(mat) < max_cols) {
      pad_cols <- max_cols - ncol(mat)
      pad_matrix <- matrix(0, nrow = nrow(mat), ncol = pad_cols)
      mat <- cbind(mat, pad_matrix)
    }
    return(mat)
  })
  
  return(padded_matrices)
}

obj_ksweep <- map(1:10, function(k) {
  message("Processing k = ", k)
  obj_final <- furrr::future_map(as.list(res), function(path_obj_fov) {
    tryCatch({
      obj_samplefov <- path_obj_fov
      obj <- list()
      collapsed_counts <- obj_samplefov$counts_raw %*% t(obj_samplefov$neighbors_collapsed[[k]])
      colnames(collapsed_counts) <- colnames(obj_samplefov$counts_raw)
      
      obj$counts <- collapsed_counts
      obj$metadata <- obj_samplefov$metadata
      obj$counts_raw <- obj_samplefov$counts_raw
      return(obj)
      
    }, error = function(e) {
      message(path_obj_fov)
      message(e)
      return(NULL)
    })
  }, .options = furrr_options(seed = TRUE))
  
  obj_final <- obj_final[!sapply(obj_final, is.null)]
  
  if (length(obj_final) == 0) {
    stop("No valid objects in obj_final after filtering NULLs.")
  }
  
  all_genes <- map(obj_final, ~ rownames(.x$counts)) %>% reduce(union)
  counts_matrices <- map(obj_final, ~ .x$counts)
  
  counts_matrices <- pad_matrix_list(counts_matrices)
  
  row_counts <- sapply(counts_matrices, nrow)
  col_counts <- sapply(counts_matrices, ncol)
  if (length(unique(row_counts)) != 1 || length(unique(col_counts)) != 1) {
    stop("Inconsistent dimensions in counts matrices after padding.")
  }
  
  obj <- list()
  obj$metadata <- map(obj_final, ~ .x$metadata) %>% bind_rows()
  
  # Combine counts matrices row-wise first to ensure consistency
  combined_counts_rowwise <- reorder_fill_genes_list(counts_matrices, all_genes) %>% reduce(rbind)
  combined_counts <- combined_counts_rowwise %>% reduce(cbind)
  
  obj$counts <- combined_counts
  obj$counts_raw <- map(obj_final, ~ .x$counts_raw) %>% reduce(cbind)
  
  return(obj)
})



library(purrr)
library(dplyr)
library(Matrix)

reorder_fill_genes_list <- function(matrices, all_genes) {
  lapply(seq_along(matrices), function(i) {
    mat <- as.matrix(matrices[[i]])
    
    mat_full <- Matrix(0, nrow = length(all_genes), ncol = ncol(mat), sparse = TRUE)
    rownames(mat_full) <- all_genes

    if (is.null(rownames(mat))) {
      stop("Matrix rownames are NULL")
    }

    common_genes <- intersect(rownames(mat), all_genes)
    common_indices <- match(common_genes, all_genes)

    message("Matrix ", i, " - Dimensions of mat: ", paste(dim(mat), collapse = "x"))
    message("Matrix ", i, " - Dimensions of mat_full: ", paste(dim(mat_full), collapse = "x"))
    message("Matrix ", i, " - Number of common genes: ", length(common_genes))
    message("Matrix ", i, " - First few common genes: ", paste(head(common_genes), collapse = ", "))
    message("Matrix ", i, " - Length of common_indices: ", length(common_indices))
    message("Matrix ", i, " - Length of common_genes: ", length(common_genes))
    message("Matrix ", i, " - Dimensions before assignment: mat = ", paste(dim(mat[common_genes, , drop = FALSE]), collapse = "x"))
    message("Matrix ", i, " - Dimensions before assignment: mat_full[common_indices, ] = ", paste(dim(mat_full[common_indices, ]), collapse = "x"))

    mat_common_subset <- mat[common_genes, , drop = FALSE]
    mat_full_subset <- mat_full[common_indices, ]
    message("Matrix ", i, " - Dimensions of mat_common_subset: ", paste(dim(mat_common_subset), collapse = "x"))
    message("Matrix ", i, " - Dimensions of mat_full_subset: ", paste(dim(mat_full_subset), collapse = "x"))

    if (nrow(mat_common_subset) == length(common_indices) && ncol(mat_common_subset) == ncol(mat_full_subset)) {
      mat_full[common_indices, ] <- mat_common_subset
      message("Matrix ", i, " - Assignment completed successfully.")
    } else {
      message("Matrix ", i, " - Mismatch in dimensions between mat_common_subset and mat_full_subset.")
    }

    return(mat_full)
  })
}

pad_matrix_list <- function(matrices) {
  max_rows <- max(sapply(matrices, nrow))
  max_cols <- max(sapply(matrices, ncol))

  padded_matrices <- lapply(matrices, function(mat) {
    mat <- Matrix(mat, sparse = TRUE)
    if (nrow(mat) < max_rows) {
      pad_rows <- max_rows - nrow(mat)
      pad_matrix <- Matrix(0, nrow = pad_rows, ncol = ncol(mat), sparse = TRUE)
      mat <- rbind(mat, pad_matrix)
    }
    if (ncol(mat) < max_cols) {
      pad_cols <- max_cols - ncol(mat)
      pad_matrix <- Matrix(0, nrow = nrow(mat), ncol = pad_cols, sparse = TRUE)
      mat <- cbind(mat, pad_matrix)
    }
    return(mat)
  })

  return(padded_matrices)
}

# Main processing function
obj_ksweep <- map(1:10, function(k) {
  message("Processing k = ", k)
  obj_final <- map(as.list(res), function(path_obj_fov) {
    tryCatch({
      obj_samplefov <- path_obj_fov
      obj <- list()
      collapsed_counts <- obj_samplefov$counts_raw %*% t(obj_samplefov$neighbors_collapsed[[k]])
      colnames(collapsed_counts) <- colnames(obj_samplefov$counts_raw)
      
      obj$counts <- collapsed_counts
      obj$metadata <- obj_samplefov$metadata
      obj$counts_raw <- obj_samplefov$counts_raw
      return(obj)
      
    }, error = function(e) {
      message(path_obj_fov)
      message(e)
      return(NULL)
    })
  })

  message("Completed future_map for k = ", k)
  
  # Remove NULL values from obj_final
  obj_final <- discard(obj_final, is.null)

  if (length(obj_final) == 0) {
    stop("No valid objects in obj_final after filtering NULLs.")
  }

  all_genes <- map(obj_final, ~ rownames(.x$counts)) %>% reduce(union)
  counts_matrices <- map(obj_final, ~ .x$counts)

  counts_matrices <- pad_matrix_list(counts_matrices)

  row_counts <- sapply(counts_matrices, nrow)
  col_counts <- sapply(counts_matrices, ncol)
  if (length(unique(row_counts)) != 1 || length(unique(col_counts)) != 1) {
    stop("Inconsistent dimensions in counts matrices after padding.")
  }

  obj <- list()
  obj$metadata <- map(obj_final, ~ .x$metadata) %>% bind_rows()
  
  # Combine counts matrices row-wise first to ensure consistency
  combined_counts_rowwise <- reorder_fill_genes_list(counts_matrices, all_genes) %>% reduce(rbind)
  combined_counts <- combined_counts_rowwise %>% reduce(cbind)
  
  obj$counts <- combined_counts
  obj$counts_raw <- map(obj_final, ~ .x$counts_raw) %>% reduce(cbind)
  message("Finished Processing k = ", k)
  return(obj)
})


save.image("/rds/projects/c/croftap-mapjagb6/xenium2/niche_roopha/analysis_full.RData")


