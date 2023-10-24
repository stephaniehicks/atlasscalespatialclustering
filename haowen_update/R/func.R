init <- function(){
  suppressMessages(library(spatialLIBD))
  suppressMessages(library(Seurat))
  suppressMessages(library(ggpubr))
  suppressMessages(library(tidyverse))
  suppressMessages(library(PRECAST))
  suppressMessages(library(igraph))
  suppressMessages(library(psych))
  suppressMessages(library(ggnewscale))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(circlize))
  suppressMessages(library(SpotClean))
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(reshape2))
  suppressMessages(library(graphlayouts))
  message("Initialize done!")
}

myPreprocess <- function(data_mat, scale_coef = 10000){
  data_mat <- data_mat/matrix(rep(colSums(data_mat), nrow(data_mat)), 
                              nrow = nrow(data_mat),
                              byrow = F) * scale_coef
  data_mat
}


swap <- function(vec){
  c(vec[2],vec[1])
}


# Convert to seurat list
new_seu_ls <- function(sel_assay="logcounts"){
  seu_ls <- list()
  #hvg_ls <- c()
  #sel_assay <- "logcounts" # counts logcounts
  
  for(i in unique(spe$sample_id)){
    idx <- spe$sample_id == i
    seu_ls[[i]] <- CreateSeuratObject(counts = spe@assays@data@listData[[sel_assay]][,idx],
                                      project = paste0("spe_",i),
                                      meta.data = data.frame(barcode = spe@colData@rownames[idx],
                                                             layer = spe@colData[idx,c("layer_guess_reordered_short")],
                                                             spatialLIBD = spe@colData[idx,c("spatialLIBD")],
                                                             coord_x = spe@int_colData@listData[["spatialCoords"]][idx,1],
                                                             coord_y = spe@int_colData@listData[["spatialCoords"]][idx,2],
                                                             row = spe$array_row[idx],
                                                             col = spe$array_col[idx]))
    seu_ls[[i]] <- FindVariableFeatures(seu_ls[[i]], verbose = F)
    #hvg_ls <- unique(c(hvg_ls,VariableFeatures(seurat_ls[[i]])))
  }
  seu_ls
}


new_spseu_ls <- function(){
  # Convert to seurat list
  seu_ls <- list()
  #hvg_ls <- c()
  
  for(i in sets_tokeep){
    seu_ls[[i]] <- CreateSeuratObject(counts = sp_obj_ls[[paste0("sp_",i)]]@assays@data@listData[["decont"]],
                                      project = paste0("spe_",i),
                                      meta.data = sp_obj_ls[[paste0("sp_",i)]]@metadata[["slide"]])
    seu_ls[[i]] <- FindVariableFeatures(seu_ls[[i]], verbose = F)
    #hvg_ls <- unique(c(hvg_ls,VariableFeatures(seurat_sp_ls[[i]])))
  }
  
  for(i in sets_tokeep){
    seu_ls[[i]]@meta.data[["coord_x"]] = seu_ls[[i]]@meta.data[["imagecol"]]
    seu_ls[[i]]@meta.data[["coord_y"]] = seu_ls[[i]]@meta.data[["imagerow"]]
    seu_ls[[i]]@meta.data[["layer"]]   = seu_ls[[i]]@meta.data[["layer_guess"]]
  }
  seu_ls
}


##### Create PRECASTObjec #####
PRECAST_test <- function(seu_ls, k=7){
  tic <- Sys.time()
  PRECAST_obj <- CreatePRECASTObject(seu_ls, project = "SpatialLIBD", gene.number = 2000, selectGenesMethod = "SPARK-X",
                                     premin.spots = 20, premin.features = 20, postmin.spots = 1, 
                                     postmin.features = 10, verbose = F)
  
  PRECAST_obj <- AddAdjList(PRECAST_obj, platform = "Visium")
  PRECAST_obj <- AddParSetting(PRECAST_obj, Sigma_equal = FALSE, verbose = F, int.model = NULL)
  
  PRECAST_obj <- PRECAST(PRECAST_obj, K = k)
  PRECAST_obj <- SelectModel(PRECAST_obj)
  seuInt <- IntegrateSpaData(PRECAST_obj, species = "Human")
  toc <- Sys.time()
  message(toc - tic)
  seuInt
}

draw_PRECAST <- function(seuInt_obj, k=7){
  if(k <= 20){
    cols_cluster <- chooseColors(palettes_name = "Classic 20", n_colors = k, plot_colors = TRUE)
  }else{
    cols_cluster <- chooseColors(palettes_name = "hue n", n_colors = k, plot_colors = TRUE)
  }
  
  pList <- SpaPlot(seuInt_obj, item = "cluster", batch = NULL, point_size = 1, cols = cols_cluster, combine = FALSE,
                   nrow.legend = 7)
  drawFigs(pList, layout.dim = c(2, 2), common.legend = TRUE, legend.position = "right", align = "hv")
}

##### Multi-level + Seurat NN Clustering Results #####
draw_slide_graph <- function(meta_data, edge_df=NULL, threshold=NULL, col_sel, pal = NULL){
  g <- ggplot()
  if(!is.null(edge_df)){
    edge_df[["x"]] <- meta_data[edge_df$from,'coord_x']
    edge_df[["y"]] <- meta_data[edge_df$from,'coord_y']
    edge_df[["xend"]] <- meta_data[edge_df$to,'coord_x']
    edge_df[["yend"]] <- meta_data[edge_df$to,'coord_y']
    if(is.null(threshold)){
      g <- g + geom_segment(mapping = aes(x=x,y=y,xend=xend,yend=yend, 
                                          color = weight),
                            size=abs(edge_df$weight),
                            data = edge_df) + 
        scale_colour_gradientn(colours = rev(viridis::rocket(n=10)))
    }else{ 
      edge_df[["filtered"]] <- as.numeric(edge_df$weight < threshold)
      g <- g + geom_segment(mapping = aes(x=x,y=y,xend=xend,yend=yend, 
                                          color = weight<threshold),
                            size=abs(edge_df$filtered+0.1),
                            data = edge_df) + 
        scale_color_manual(values = c("grey70","red"))
    }
     g <- g + new_scale_colour()
  }
  g <- g + geom_point(mapping = aes_string(x = "coord_x", y = "coord_y", color=col_sel), 
                      data = meta_data)
  if(is.numeric(meta_data[[col_sel]])){
    require(viridis)
    g <- g + scale_colour_gradientn(colours = viridis::viridis(n=10))
  }else{
    if(is.null(pal))    g <- g + scale_color_discrete()
    else g <- g + scale_colour_manual(values = pal)
  }
    g <- g + theme_classic() +coord_flip()
  
  return(g)
}



draw_edge_dstr <-  function(edge_df, col_sel, sample_name = ""){
  g <- ggplot(edge_df,aes_string(x=col_sel)) +
    geom_density(color="darkblue", fill="lightblue") + 
    scale_y_continuous(trans = "pseudo_log") +
    labs(title = paste("Log density of", col_sel, sample_name)) +
    theme_classic()
  return(g)
}

stage_1 <- function(seu_ls, cor_threshold = 0.2, nn = 12, cl_resolution = 10, 
                    top_pcs = 30, cl_min=5, preprocess = T, hvg = 2000, cor_slot = "PC", edge_smoothing = F, use_glmpca = F){
  tic <- Sys.time()
  for(i in names(seu_ls)){
    #seu_ls[[i]] <- NormalizeData(seu_ls[[i]], normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
    if(preprocess){
      stopifnot("HVG# Exceeded" = hvg <= nrow(seu_ls[[i]]))
      seu_ls[[i]] <- FindVariableFeatures(seu_ls[[i]], selection.method = "vst", nfeatures = hvg,verbose = F)
      seu_ls[[i]] <- ScaleData(seu_ls[[i]], features = rownames(seu_ls[[i]]),verbose = F)
      
    }else{
      VariableFeatures(object = seu_ls[[i]]) <- row.names(seu_ls[[i]])
      hvg <- nrow(seu_ls[[i]])
    }
    require(scater)
    require(scry)
    
    seu_ls[[i]] <- RunPCA(seu_ls[[i]], features = VariableFeatures(object = seu_ls[[i]]),verbose = F)
    
    if(use_glmpca == T){
      mat <- as.matrix(seu_ls[[i]]@assays$RNA@counts[VariableFeatures(object = seu_ls[[i]]),])
      mat <- nullResiduals(mat, type="deviance")
      res1 <- suppressWarnings(calculatePCA(mat,ncomponents = top_pcs, scale = TRUE, BSPARAM = BiocSingular::RandomParam()))
      seu_ls[[i]]@reductions[["pca"]]@cell.embeddings <- res1 %>% as.matrix()
      #seu_ls[[i]]@reductions[["pca"]]@feature.loadings <- res1$loadings %>% as.matrix()
    } 
    
    stopifnot("Not supported Cor Calc Method" = cor_slot %in% c("PC","HVG"))
    
    if(!edge_smoothing){
      if(cor_slot == "PC"){
        cor_mat <- cor(t(as.matrix(seu_ls[[i]]@reductions[["pca"]]@cell.embeddings[,1:top_pcs])),method = "pearson")
      } else if(cor_slot == "HVG") cor_mat <- cor(seu_ls[[i]]@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_ls[[i]]),],method = "pearson")
      
      dist_mat <- matrix(0, nrow = ncol(seu_ls[[i]]), ncol = ncol(cor_mat), 
                         dimnames = list(colnames(cor_mat), row.names(cor_mat)))
      
      for(j in seq_len(nrow(dist_mat))){
        dist_mat[,j] <- sqrt((seu_ls[[i]]@meta.data[["coord_x"]] - seu_ls[[i]]@meta.data[["coord_x"]][j])^2 + 
                               (seu_ls[[i]]@meta.data[["coord_y"]] - seu_ls[[i]]@meta.data[["coord_y"]][j])^2)
        not_nn_vec <- sort(dist_mat[,j])[(nn+2):ncol(dist_mat)] %>% names()
        cor_mat[j,not_nn_vec] <- 0
        cor_mat[j,j] <- 0
      }
    }else{
      cor_mat <- matrix(0, nrow = ncol(seu_ls[[i]]), ncol = ncol(seu_ls[[i]]), 
                        dimnames = list(colnames(seu_ls[[i]]), colnames(seu_ls[[i]])))
      dist_mat <- matrix(0, nrow = nrow(cor_mat), ncol = ncol(cor_mat), 
                         dimnames = list(colnames(cor_mat), row.names(cor_mat)))
      nn_mat <- matrix(0, nrow = nrow(cor_mat), ncol = ncol(cor_mat), 
                         dimnames = list(colnames(cor_mat), row.names(cor_mat)))
      # Find NN
      for(j in seq_len(nrow(dist_mat))){
        dist_mat[,j] <- sqrt((seu_ls[[i]]@meta.data[["coord_x"]] - seu_ls[[i]]@meta.data[["coord_x"]][j])^2 + 
                               (seu_ls[[i]]@meta.data[["coord_y"]] - seu_ls[[i]]@meta.data[["coord_y"]][j])^2)
        nn_vec <- sort(dist_mat[,j])[1:(nn+1)] %>% names()
        nn_mat[j, nn_vec] <- 1
      }
      # Calc Smoothed Edge Weight
      for(j in seq_len(nrow(dist_mat)-1)){
        for(k in (j+1):nrow(dist_mat)){
          if(nn_mat[j,k] != 1) next
          nn_vec1 <- colnames(nn_mat)[nn_mat[j,]==1]
          nn_vec2 <- colnames(nn_mat)[nn_mat[k,]==1]
          
          if(FALSE){
            ggplot(seu_ls[[i]]@meta.data[union(which(nn_mat[j,] == 1),
                                               which(nn_mat[k,] == 1)),],
                   aes(x = coord_x, y = coord_y))+
              geom_point() + 
              theme_classic()
          }
          
          nn_common <- intersect(nn_vec1, nn_vec2)
          nn_vec1 <- c(colnames(nn_mat)[j],nn_vec1[!nn_vec1 %in% nn_common])
          nn_vec2 <- c(colnames(nn_mat)[k],nn_vec2[!nn_vec2 %in% nn_common])
          
          if(cor_slot == "PC"){
              cor_val <- cor(x = colMeans(matrix(seu_ls[[i]]@reductions[["pca"]]@cell.embeddings[nn_vec1,1:top_pcs], ncol = top_pcs)),
                             y = colMeans(matrix(seu_ls[[i]]@reductions[["pca"]]@cell.embeddings[nn_vec2,1:top_pcs], ncol = top_pcs)),
                             method = "pearson")
          }else if(cor_slot == "HVG") {
            cor_mat <- cor(seu_ls[[i]]@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_ls[[i]]),],method = "pearson")
            cor_val <- cor(x = rowMeans(matrix(seu_ls[[i]]@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_ls[[i]]),nn_vec1], nrow = hvg)),
                           y = rowMeans(matrix(seu_ls[[i]]@assays[["RNA"]]@scale.data[VariableFeatures(object = seu_ls[[i]]),nn_vec2], nrow = hvg)),
                           method = "pearson")
          }
          cor_mat[j,k] <- cor_val
          cor_mat[k,j] <- cor_val
        }
      }
    }
    
    
    
    g <- graph.adjacency(cor_mat, mode = "directed", 
                         weighted = TRUE, diag = TRUE)
    
    seu_ls[[i]]@misc[["raw_edges"]] <- as_data_frame(g,"edges")
    seu_ls[[i]]@misc[["raw_graph_plot_label"]] <- 
      draw_slide_graph(seu_ls[[i]]@meta.data,seu_ls[[i]]@misc[["raw_edges"]],
                       cor_threshold,"layer")
    
    # Primary Clustering
    cor_mat[cor_mat < cor_threshold] <- 0
    g <- graph.adjacency(cor_mat, mode = "directed", 
                         weighted = TRUE, diag = TRUE)
    g <- as.undirected(g,mode = "mutual")
    
    seu_ls[[i]]@misc[["edges"]] <- as_data_frame(g,"edges")
    seu_ls[[i]]@misc[["graph_plot_label"]] <- 
      draw_slide_graph(seu_ls[[i]]@meta.data,seu_ls[[i]]@misc[["edges"]],
                       cor_threshold,"layer")
    
    louvain_res <- cluster_louvain(g, resolution = cl_resolution)
    seu_ls[[i]]@meta.data[["cluster"]] <- as.factor(louvain_res$membership)
    
    message(paste(i,"cluster number:",length(levels(as.factor(louvain_res$membership)))))
    seu_ls[[i]]@misc[["graph_plot_cluster"]] <- 
      draw_slide_graph(seu_ls[[i]]@meta.data,seu_ls[[i]]@misc[["edges"]],
                       cor_threshold,"cluster") + theme(legend.position = "none")
    
    # Cluster merging
    seu_ls[[i]]@meta.data[["merged_cluster"]] <- seu_ls[[i]]@meta.data[["cluster"]]
    while(sum(table(seu_ls[[i]]@meta.data[["merged_cluster"]]) < cl_min) > 0){
      sml_cl_idx <- names(table(seu_ls[[i]]@meta.data[["merged_cluster"]]))[table(seu_ls[[i]]@meta.data[["merged_cluster"]]) < cl_min]
      
      for(j in sml_cl_idx){
        node_ls <- colnames(seu_ls[[i]])[seu_ls[[i]]@meta.data[["merged_cluster"]]==j]
        nn_ls <- c()
        for(node in node_ls){
          nn_ls <- append(nn_ls,sort(dist_mat[,node])[1:(nn+1)] %>% names())
        }
        nn_ls <- nn_ls[!nn_ls %in% node_ls]
        nn_cl <- seu_ls[[i]]@meta.data[["merged_cluster"]][colnames(seu_ls[[i]]) %in% nn_ls]
        seu_ls[[i]]@meta.data[["merged_cluster"]][colnames(seu_ls[[i]]) %in% node_ls] <- names(sort(table(nn_cl),decreasing=TRUE)[1])
      }
      seu_ls[[i]]@meta.data[["merged_cluster"]] <- droplevels(seu_ls[[i]]@meta.data[["merged_cluster"]])
      message(paste(i,"merged cluster number:",length(levels(seu_ls[[i]]@meta.data[["merged_cluster"]]))))
    }
    
    
    seu_ls[[i]]@misc[["graph_plot_cluster_merged"]] <- 
      draw_slide_graph(seu_ls[[i]]@meta.data,seu_ls[[i]]@misc[["edges"]],
                       cor_threshold,"merged_cluster") + theme(legend.position = "none")
    
  }
  toc <- Sys.time()
  message("Elasped Time(sec):")
  message(toc - tic)
  seu_ls
}

mat_cor <- function(mat_x, mat_y){
  cor_mat <- matrix(0,nrow = nrow(mat_x),ncol = nrow(mat_y))
  for(i in seq_len(nrow(mat_x))){
    cor_mat[i,] <- sapply(seq_len(nrow(mat_y)),
                          function(x){
                            cor(x=mat_x[i,],y=mat_y[x,],method = "pearson")
                          })
  }
  row.names(cor_mat) <- row.names(mat_x)
  colnames(cor_mat) <- row.names(mat_y)
  cor_mat
}

louvain_w_cor <- function(cor_mat_, nn_=10, res_ = 1){
  for(j in seq_len(nrow(cor_mat_))){
    not_nn_vec <- sort(cor_mat_[,j],decreasing = T)[(nn_+2):ncol(cor_mat_)] %>% names()
    cor_mat_[j,not_nn_vec] <- 0
    cor_mat_[j,j] <- 0
  }
  g <- graph.adjacency(cor_mat_, mode = "directed", 
                       weighted = TRUE, diag = TRUE)
  g <- as.undirected(g,mode = "mutual")
  
  louvain_res <- cluster_louvain(g, resolution = res_)
  louvain_res
}


stage_2 <- function(seu_ls, top_pcs = 30, nn_2=10, cl_key = "merged_cluster",
                    rtn_seurat = F,method="louvain", hly_cor = 0.9, hvg = 2000, cor_slot = "PC", use_glmpca = F, resolution = 1,
                    rare_ct="none"){
  tic <- Sys.time()
  hvg_union <- c()
  gene_intersect <- row.names(seu_ls[[names(seu_ls)[1]]])
  cl_num <- c()
  
  stopifnot("Not supported Cor Calc Method" = cor_slot %in% c("PC","HVG"))
  
  for(i in names(seu_ls)){
    hvg_union <- union(hvg_union,VariableFeatures(object = seu_ls[[i]]))
    gene_intersect <- intersect(gene_intersect, row.names(seu_ls[[i]]))
    cl_num <- append(cl_num,
                     c(length(levels(droplevels(as.factor(seu_ls[[i]]@meta.data[[cl_key]]))))))
  }
  
  hvg_union <- intersect(hvg_union, gene_intersect)
  
  cl_expr <- matrix(0, nrow = sum(cl_num),ncol = length(hvg_union))
  cl_df <- data.frame(sample=rep(names(seu_ls),cl_num),cluster=0)
  
  for(i in names(seu_ls)){
    cl_df[["cluster"]][cl_df$sample==i] <- levels(droplevels(as.factor(seu_ls[[i]]@meta.data[[cl_key]])))
  }
  row.names(cl_expr) <- paste0("X",cl_df$sample,"_",cl_df$cluster)
  colnames(cl_expr) <- hvg_union
  
  for(i in seq_len(nrow(cl_df))){
    idx <- which(seu_ls[[cl_df$sample[i]]]@meta.data[[cl_key]] == cl_df$cluster[i])
    if(length(idx)>1) cl_expr[i,] <-  seu_ls[[cl_df$sample[i]]]@assays[["RNA"]]@data[hvg_union,idx] %>% rowSums(.)/length(idx)
    else cl_expr[i,] <-  seu_ls[[cl_df$sample[i]]]@assays[["RNA"]]@data[hvg_union,idx]
  }
  
  cl_expr_obj <- CreateSeuratObject(t(cl_expr),verbose = F)
  cl_expr_obj <- FindVariableFeatures(cl_expr_obj, selection.method = "vst", nfeatures = hvg,verbose = F)
  cl_expr_obj <- ScaleData(cl_expr_obj, features = rownames(cl_expr_obj),verbose = F)
  cl_expr_obj <- RunPCA(cl_expr_obj, features = VariableFeatures(object = cl_expr_obj),verbose = F)
  cl_expr_obj <- RunUMAP(cl_expr_obj, dims = 1:top_pcs,verbose = F)
  if(use_glmpca){
    mat <- as.matrix(cl_expr_obj@assays$RNA@counts[VariableFeatures(object = cl_expr_obj),])
    mat <- nullResiduals(mat, type="deviance")
    res1 <- calculatePCA(mat,ncomponents = top_pcs, scale = TRUE, BSPARAM = BiocSingular::RandomParam())
    cl_expr_obj@reductions[["pca"]]@cell.embeddings <- res1 %>% as.matrix()
    #cl_expr_obj@reductions[["pca"]]@feature.loadings <- gp_res$loadings %>% as.matrix()
  }
  
  if(method == "louvain"){
    if(cor_slot == "PC"){
        cor_mat <- cor(t(as.matrix(cl_expr_obj@reductions[["pca"]]@cell.embeddings[,1:top_pcs])),method = "pearson")
    }else if(cor_slot == "HVG") cor_mat <- cor(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj),],method = "pearson")
    
    louvain_res <- louvain_w_cor(cor_mat)
    #table(louvain_res$membership)
    cl_df[["louvain"]] <- as.character(louvain_res$membership)
    sml_cl_idx <- names(table(cl_df[["louvain"]]))[table(cl_df[["louvain"]]) < 2]
    cl_df[["louvain"]][cl_df[["louvain"]] %in% sml_cl_idx] <- NA
  }else if(method == "MNN"){
    g <- make_empty_graph(directed = F)
    node_df <- data.frame(ID=seq_along(colnames(cl_expr_obj)),
                          cl = colnames(cl_expr_obj),
                          sample = cl_expr_obj$orig.ident,
                          row.names = colnames(cl_expr_obj))
    g <- add_vertices(g, nrow(node_df))
    V(g)$cl <- node_df$cl
    V(g)$sample <- node_df$sample
    
    # Enable MNN within sample
    if(rare_ct == "m"){
      mnn_seq <- seq_along(levels(cl_expr_obj@meta.data[["orig.ident"]]))
    }else{
      mnn_seq <- seq_along(levels(cl_expr_obj@meta.data[["orig.ident"]]))[-1]
    }  
    
    for(i in mnn_seq){
      if(rare_ct == "m"){
        mnn_seq2 <- 1:i
      }else{
        mnn_seq2 <- 1:(i-1)
      }
      for(j in mnn_seq2){
        idx_i <- cl_expr_obj@meta.data[["orig.ident"]] == levels(cl_expr_obj@meta.data[["orig.ident"]])[i]
       
        idx_j <- cl_expr_obj@meta.data[["orig.ident"]] == levels(cl_expr_obj@meta.data[["orig.ident"]])[j]
        
        if(cor_slot == "PC"){
          mat_i <- as.matrix(cl_expr_obj@reductions[["pca"]]@cell.embeddings[idx_i,1:top_pcs])
          mat_j <- as.matrix(cl_expr_obj@reductions[["pca"]]@cell.embeddings[idx_j,1:top_pcs])
        }else if(cor_slot == "HVG"){
          mat_i <- as.matrix(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj), idx_i]) %>% t()
          mat_j <- as.matrix(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj), idx_j]) %>% t()
        }
        
        
        cor_mat <- mat_cor(mat_i,mat_j)
        hly_cor_mat <- cor_mat > hly_cor
        i_nn_mat <- matrix(0,nrow = nrow(mat_i),ncol = nrow(mat_j))
        j_nn_mat <- matrix(0,nrow = nrow(mat_i),ncol = nrow(mat_j))
        for(k in seq_len(nrow(mat_i))){
          i_nn_mat[k,cor_mat[k,] %>% order(decreasing = T) %>% .[1:nn_2]] <- T
        }
        for(k in seq_len(nrow(mat_j))){
          j_nn_mat[cor_mat[,k] %>% order(decreasing = T) %>% .[1:nn_2],k] <- T
        }
        graph_mat <- hly_cor_mat | (i_nn_mat & j_nn_mat)
        # Add edges
        row_id <- node_df[row.names(mat_i),1]
        col_id <- node_df[row.names(mat_j),1]
        edge_vec <- c()
        for(k in seq_len(nrow(graph_mat))){
          tmp_vec <- col_id[which(graph_mat[k,])]
          if(length(tmp_vec)>0){
            edge_vec <- c(edge_vec,
                          paste(row_id[k],tmp_vec) %>% str_split(pattern = " ") %>% unlist() %>% as.numeric())
          }
        }
        g <- add.edges(g,edges = edge_vec)
        
      }
    }
    #layout <- layout_with_kk(g)
    #f <- factor(node_df$sample)
    #vcols <- chooseColors(palettes_name = "Blink 23", 
    #                             n_colors = length(levels(f)), plot_colors = F)
    #cols_cluster <- vcols[f]
    #plot(g, layout = layout,vertex.size=3, vertex.color=cols_cluster, vertex.label.cex=0.1)
    #legend("topleft", legend = levels(f), pch = 16, col = vcols, bty = "n")
    louvain_res <- cluster_louvain(g, resolution = resolution)
    #table(louvain_res$membership)
    cl_df[["louvain"]] <- as.character(louvain_res$membership)
    sml_cl_idx <- names(table(cl_df[["louvain"]]))[table(cl_df[["louvain"]]) < 2]
    cl_df[["louvain"]][cl_df[["louvain"]] %in% sml_cl_idx] <- NA
    if(rare_ct == "a"){
      if(cor_slot == "PC"){
        cor_mat <- cor(t(as.matrix(cl_expr_obj@reductions[["pca"]]@cell.embeddings[is.na(cl_df[["louvain"]]),1:top_pcs])),method = "pearson")
      }else if(cor_slot == "HVG") cor_mat <- cor(cl_expr_obj@assays[["RNA"]]@scale.data[VariableFeatures(object = cl_expr_obj),is.na(cl_df[["louvain"]])],method = "pearson")
      
      louvain_res <- louvain_w_cor(cor_mat, nn_=20,res_ = 1)
      #table(louvain_res$membership)
      cl_df[["louvain_2"]] <- NA
      cl_df[["louvain_2"]][is.na(cl_df[["louvain"]])] <- as.character(louvain_res$membership)
      cl_df[["louvain"]] <- apply(cl_df,1,
                                  FUN = function(x){
                                    ifelse(is.na(x[["louvain"]]),
                                           paste0("R2_",x[["louvain_2"]]),
                                           paste0("R1_",x[["louvain"]]))
                                  })
    
      }
    
    
  }
  
  
  
  
  cl_df[["umap_1"]] <- cl_expr_obj@reductions[["umap"]]@cell.embeddings[,1]
  cl_df[["umap_2"]] <- cl_expr_obj@reductions[["umap"]]@cell.embeddings[,2]
  cl_df[["layer"]] <- apply(cl_df,1,
                            FUN = function(x){
                              idx <- seu_ls[[x[["sample"]] ]]@meta.data[[cl_key]] == x[["cluster"]]
                              layer_vec <- seu_ls[[x[["sample"]] ]]@meta.data[["layer"]][idx]
                              table(layer_vec) %>% sort(decreasing = T) %>% names() %>% .[1]
                              
                            })
  toc <- Sys.time()
  message("Elasped Time(sec):")
  message(toc - tic)
  if(rtn_seurat){
    if(method == "MNN") list(seurat_obj = cl_expr_obj, cl_df = cl_df, g=g)
    else list(seurat_obj = cl_expr_obj, cl_df = cl_df)
  }else{
    cl_df
  }
}

assign_label <- function(seu_ls, cl_df,nn,cor_threshold = 0.2,cl_key = "merged_cluster"){
  # Assign secondary clustering label
  for(i in names(seu_ls)){
    seu_ls[[i]]@meta.data[[paste0("sec_cluster_",nn)]] <- apply(seu_ls[[i]]@meta.data,1,
                                                                FUN = function(x){
                                                                  cl_df$louvain[cl_df$sample==i & cl_df$cluster==x[[cl_key]]]
                                                                })
    seu_ls[[i]]@misc[[paste0("graph_plot_cluster_sec_",nn)]] <- 
      draw_slide_graph(seu_ls[[i]]@meta.data,seu_ls[[i]]@misc[["edges"]],
                       cor_threshold,paste0("sec_cluster_",nn)) +
      labs(title = paste("sample:", i))
  }
  seu_ls
}

plotConfusionMatrix <- function(x,y,col_title = ""){
  require(ComplexHeatmap)
  na_index <- (is.na(x) | is.na(y))
  x <- x[!na_index]
  y <- y[!na_index]
  u_x <- unique(x)
  u_y <- unique(y)
  hm_mat <- matrix(0, nrow = length(u_x), ncol = length(u_y))
  row.names(hm_mat) <- u_x
  colnames(hm_mat) <- u_y
  for(i in seq_len(length(u_x))){
    for(j in seq_len(length(u_y))){
      hm_mat[i,j] = sum(x==u_x[i] & y ==u_y[j])
    }
  }
  hm_mat <- hm_mat/matrix(rep(rowSums(hm_mat),ncol(hm_mat)), 
                          byrow = F,nrow = nrow(hm_mat), ncol = ncol(hm_mat))
  Heatmap(hm_mat,cluster_columns = F, cluster_rows = F,
          col = colorRamp2(seq(0, 1,length.out=5), viridis::viridis(5)),
          rect_gp = gpar(col = "white", lwd = 1),column_title = col_title)
}

