library(tidyverse)
library(SpatialExperiment)
library(STexampleData)
library(scater)             # log-transformation
library(here)
library(Seurat)

# Helper functions for four Possible Spatial Patterns
# Functions that creates expression ---------------------------------------
# * Smooth & Circular -------------------------------------------------------
smth_circ_fun <- function(x, y){ (x-0.5)^2 + (y-0.5)^2 }

# * Smooth & Linear -------------------------------------------------------
smth_lnr_fun <- function(x, y) { 0.5-x + 0.5-y }


# * Layered & Circular -------------------------------------------------------
# Calcualte center of the graph
lyr_circ_fun <- function(x,y){
  # ret <- NA_real_
  dplyr::case_when(
    (0.5-x)^2 + (0.5-y)^2 < (0.125)^2 ~ 40,
    (0.5-x)^2 + (0.5-y)^2 < (0.25)^2 ~ 30,
    (0.5-x)^2 + (0.5-y)^2 < (0.375)^2 ~ 20,
    (0.5-x)^2 + (0.5-y)^2 < (1)^2 ~ 10,
    TRUE ~ NA_real_
  )
  # return(ret)
}

# * Layered & Linear -------------------------------------------------------
lyr_lnr_fun <- function(x,y){
  # ret <- NA_real_
  dplyr::case_when(
    1.5 - x - y <= 0 ~ 40,
    1 - x - y < 0 ~ 30,
    0.5 - x - y < 0 ~ 20,
    0 - x - y < 0 ~ 10,
    TRUE ~ NA_real_
  )
}

# Corresponding Functions to Create Covariates ----------------------------
# * Layered Pattern -------------------------------------------------------
# Caliberate

lnr_group_fun <- function(x,y){
  dplyr::case_when(
    1.5 - x - y <= 0 ~ "Group 1",
    1 - x - y < 0 ~ "Group 2",
    0.5 - x - y < 0 ~ "Group 3",
    0 - x - y < 0 ~  "Group 4",
    TRUE ~ NA_character_
  )
}

# # * Circular Pattern -------------------------------------------------------
circ_group_fun <- function(x,y){
  dplyr::case_when(
    (0.5-x)^2 + (0.5-y)^2 < (0.125)^2 ~ "Group 1",
    (0.5-x)^2 + (0.5-y)^2 < (0.25)^2 ~ "Group 2",
    (0.5-x)^2 + (0.5-y)^2 < (0.375)^2 ~ "Group 3",
    (0.5-x)^2 + (0.5-y)^2 < (1)^2 ~ "Group 4",
    TRUE ~ NA_character_
  )
}

semicirc_group_fun <- function(x,y){
  dplyr::case_when(
    (0.5-x)^2 + y^2 < (0.2)^2 ~ "Group 1",
    (0.5-x)^2 + y^2 < (0.4)^2 ~ "Group 2",
    (0.5-x)^2 + y^2 < (0.6)^2 ~ "Group 3",
    (0.5-x)^2 + y^2 < 1.25 ~ "Group 4",
    TRUE ~ NA_character_
  )
}

qurtcirc_group_fun <- function(x,y){
  dplyr::case_when(
    x^2 + y^2 < (0.3)^2 ~ "Group 1",
    x^2 + y^2 < (0.6)^2 ~ "Group 2",
    x^2 + y^2 < (0.8)^2 ~ "Group 3",
    x^2 + y^2 < 2 ~ "Group 4",
    TRUE ~ NA_character_
  )
}

norm <- function(x){
  (x - min(x))/ (max(x) - min(x))
}

### Helper function ends

# Sim Parameters ----------------------------------------------------------
it <- 1

tmp_spe <- Visium_humanDLPFC()# Visium Template
n_gene <- 200
n_group <- 4
n_spots <- ncol(tmp_spe)

seu_ls <- list()

set.seed(it)
# Ground truth ranking
rowData_df <- data.frame(
  gene_idx  = 1:(n_gene*n_group), 
  mu_shift = runif(n = n_gene*n_group, min = 2, max = 4),# Different mean expression level
  var_scale = runif(n = n_gene*n_group, min = 1, max = 4) # Different effect size
) |> 
  mutate(gene_name = paste0("gene_", gene_idx),
         marker=paste0("Group ",ceiling(gene_idx/200))) |> 
  column_to_rownames("gene_name")



# Choose one of the spatial pattern
for (pattern in c("lnr","circ","semicirc","qurtcirc")) {
  print(pattern)
  str_func <- switch(pattern,
                     "lnr"=lnr_group_fun,
                     "circ"=circ_group_fun,
                     "semicirc"=semicirc_group_fun,
                     "qurtcirc"=qurtcirc_group_fun) 
  set.seed(it)
  spa_str_df <- spatialCoords(tmp_spe) |> 
    data.frame() |> 
    mutate(
      x = norm(pxl_col_in_fullres),
      y = norm(pxl_row_in_fullres),
      z = str_func(
        x, y
      )
      #std_z = scale(z) # standardized
    )
  
  stopifnot(all(!is.na(spa_str_df$std_z)))
  # Check region
  #ggplot(spa_str_df,aes(x=pxl_col_in_fullres,y=pxl_row_in_fullres, color=as.factor(z))) + 
  #  geom_point() + 
  #  theme_classic()
  # Simulated raw counts following Poisson noise
  # gene by spot
  gene_count_mat <- 
    map2(.x = rowData_df$mu_shift,
         #.y = rowData_df$var_scale,
         .y = rowData_df$marker,
         .f = function(shift, markers){
           #eta_vec <- spa_str_df$std_z*scale + shift
           eta_vec <- (spa_str_df$z == markers)*shift
           stopifnot(length(eta_vec) == n_spots)
           ret_vec <- rnorm(
             n = n_spots,
             mean = eta_vec,
             sd = 0.2 # Extremely small noise
           )
           return(data.frame(`gene_` = ret_vec))
         }) |> 
    list_cbind() |> 
    t()
  
  
  rownames(gene_count_mat) <- rownames(rowData_df)
  colnames(gene_count_mat) <- rownames(spa_str_df)
  
  seu_obj <- CreateSeuratObject(counts = gene_count_mat,
                                project = "tmp_seurat",
                                meta.data = spa_str_df)
  seu_obj <- ScaleData(seu_obj, features = row.names(seu_obj), do.scale = F, do.center = F)
  seu_ls[[pattern]] <- seu_obj
}











# Create spe
sim_spe <- SpatialExperiment(
  assay = list(logcounts = gene_count_mat), # Convert to gene by spots
  colData = spa_str_df, 
  rowData = rowData_df,
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres")
)

# Create Seurat
pattern_vec <- c("lnr","circ","semicirc","qurtcirc")
DoHeatmap(seu_ls[[pattern_vec[1]]], features = row.names(seu_obj)[c(1:10,201:210,401:410,601:610)],group.by = "z") + NoLegend()
ggplot(seu_ls[["lnr"]]@meta.data,aes(x=pxl_col_in_fullres,y=pxl_row_in_fullres, color=as.factor(z))) + 
  geom_point() + 
  theme_classic()

for(pattern in pattern_vec){
  seu_ls[[pattern]]@meta.data[["layer"]] = seu_ls[[pattern]]@meta.data[["z"]]
  seu_ls[[pattern]]@meta.data[["coord_x"]] = seu_ls[[pattern]]@meta.data[["x"]]*1000
  seu_ls[[pattern]]@meta.data[["coord_y"]] = seu_ls[[pattern]]@meta.data[["y"]]*1000
  seu_ls[[pattern]]@meta.data[["row"]] = seu_ls[[pattern]]@meta.data[["pxl_row_in_fullres"]]
  seu_ls[[pattern]]@meta.data[["col"]] = seu_ls[[pattern]]@meta.data[["pxl_col_in_fullres"]]
}

seu_ls <- stage_1(seu_ls,preprocess = F)
rtn_ls <- stage_2(seu_ls,cl_key = "merged_cluster",rtn_seurat = T,nn_2 = 1,method = "MNN")
seu_ls <- assign_label(seu_ls,rtn_ls[["cl_df"]], "MNN", cl_key = "merged_cluster")
