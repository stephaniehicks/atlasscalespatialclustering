library(SpotClean)
library(ggplot2)
library(Matrix)
library(dplyr)
library(viridis)
library(cowplot)
library(reshape2)
library(SoupX)
library(celda)
library(Seurat)
library(scran)
library(scuttle)
library(spatialLIBD)

my_fetch_data <- function (type = c("sce", "sce_layer", "modeling_results", "sce_example", 
                                    "spe", "spatialDLPFC_Visium", "spatialDLPFC_Visium_pseudobulk", 
                                    "spatialDLPFC_Visium_modeling_results", "spatialDLPFC_Visium_SPG", 
                                    "spatialDLPFC_snRNAseq", "Visium_SPG_AD_Visium_wholegenome_spe", 
                                    "Visium_SPG_AD_Visium_targeted_spe", "Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe", 
                                    "Visium_SPG_AD_Visium_wholegenome_modeling_results"), destdir = tempdir(), 
                           eh = ExperimentHub::ExperimentHub(), bfc = BiocFileCache::BiocFileCache()) 
{
  sce <- sce_layer <- modeling_results <- sce_sub <- spe <- NULL
  type <- match.arg(type)
  stopifnot(methods::is(eh, "ExperimentHub"))
  if (type == "spe") {
    spe <- sce_to_spe(fetch_data("sce", destdir = destdir, 
                                 eh = eh))
    return(spe)
  }
  if (type == "sce") {
    if (!T) {
      warning(paste("Your system might not have enough memory available.", 
                    "Try with a machine that has more memory", "or use the 'sce_example'."))
    }
    tag <- "Human_Pilot_DLPFC_Visium_spatialLIBD"
    hub_title <- "Human_Pilot_DLPFC_Visium_spatialLIBD_spot_level_SCE"
    file_name <- "Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"
    url <- "https://www.dropbox.com/s/f4wcvtdq428y73p/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata?dl=1"
  }
  else if (type == "sce_layer") {
    tag <- "Human_Pilot_DLPFC_Visium_spatialLIBD"
    hub_title <- "Human_Pilot_DLPFC_Visium_spatialLIBD_layer_level_SCE"
    file_name <- "Human_DLPFC_Visium_processedData_sce_scran_sce_layer_spatialLIBD.Rdata"
    url <- "https://www.dropbox.com/s/bg8xwysh2vnjwvg/Human_DLPFC_Visium_processedData_sce_scran_sce_layer_spatialLIBD.Rdata?dl=1"
  }
  else if (type == "modeling_results") {
    tag <- "Human_Pilot_DLPFC_Visium_spatialLIBD"
    hub_title <- "Human_Pilot_DLPFC_Visium_spatialLIBD_modeling_results"
    file_name <- "Human_DLPFC_Visium_modeling_results.Rdata"
    url <- "https://www.dropbox.com/s/se6rrgb9yhm5gfh/Human_DLPFC_Visium_modeling_results.Rdata?dl=1"
  }
  else if (type == "sce_example") {
    tag <- "Human_Pilot_DLPFC_Visium_spatialLIBD"
    hub_title <- "Human_DLPFC_Visium_sce_example"
    file_name <- "sce_sub_for_vignette.Rdata"
    url <- "https://www.dropbox.com/s/5ra9o8ku9iyyf70/sce_sub_for_vignette.Rdata?dl=1"
  }
  else if (type == "spatialDLPFC_Visium") {
    if (!T) {
      warning(paste("Your system might not have enough memory available (7GB).", 
                    "Try with a machine that has more memory."))
    }
    tag <- "spatialDLPFC_Visium_VisiumSPG_snRNAseq_spatialLIBD"
    hub_title <- "spatialDLPFC_Visium_spe"
    file_name <- "spe_filtered_final_with_clusters_and_deconvolution_results.rds"
    url <- "https://www.dropbox.com/s/y2ifv5v8g68papf/spe_filtered_final_with_clusters_and_deconvolution_results.rds?dl=1"
  }
  else if (type == "spatialDLPFC_Visium_pseudobulk") {
    tag <- "spatialDLPFC_Visium_VisiumSPG_snRNAseq_spatialLIBD"
    hub_title <- "spatialDLPFC_Visium_pseudobulk_spe"
    file_name <- "sce_pseudo_BayesSpace_k09.rds"
    url <- "https://www.dropbox.com/s/pbti4strsfk1m55/sce_pseudo_BayesSpace_k09.rds?dl=1"
  }
  else if (type == "spatialDLPFC_Visium_modeling_results") {
    tag <- "spatialDLPFC_Visium_VisiumSPG_snRNAseq_spatialLIBD"
    hub_title <- type
    file_name <- "modeling_results_BayesSpace_k09.Rdata"
    url <- "https://www.dropbox.com/s/srkb2ife75px2yz/modeling_results_BayesSpace_k09.Rdata?dl=1"
  }
  else if (type == "spatialDLPFC_Visium_SPG") {
    tag <- "spatialDLPFC_Visium_VisiumSPG_snRNAseq_spatialLIBD"
    hub_title <- "spatialDLPFC_Visium_SPG_spe"
    file_name <- "spe.rds"
    url <- "https://www.dropbox.com/s/nbf13dna9ibqfaa/spe.rds?dl=1"
  }
  else if (type == "spatialDLPFC_snRNAseq") {
    tag <- "spatialDLPFC_Visium_VisiumSPG_snRNAseq_spatialLIBD"
    hub_title <- type
    file_name <- "sce_DLPFC_annotated.zip"
    url <- "https://www.dropbox.com/s/5919zt00vm1ht8e/sce_DLPFC_annotated.zip?dl=1"
  }
  else if (type == "Visium_SPG_AD_Visium_wholegenome_spe") {
    tag <- "Visium_SPG_AD_Alzheimer_Disease_ITC_spatialLIBD"
    hub_title <- type
    file_name <- "Visium_SPG_AD_spe_wholegenome.Rdata"
    url <- "https://www.dropbox.com/s/ng036m63grykdm6/Visium_SPG_AD_spe_wholegenome.Rdata?dl=1"
  }
  else if (type == "Visium_SPG_AD_Visium_targeted_spe") {
    tag <- "Visium_SPG_AD_Alzheimer_Disease_ITC_spatialLIBD"
    hub_title <- type
    file_name <- "Visium_SPG_AD_spe_targeted.Rdata"
    url <- "https://www.dropbox.com/s/kda9160awc2h8jq/Visium_SPG_AD_spe_targeted.Rdata?dl=1"
  }
  else if (type == "Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe") {
    tag <- "Visium_SPG_AD_Alzheimer_Disease_ITC_spatialLIBD"
    hub_title <- type
    file_name <- "sce_pseudo_pathology_wholegenome.rds"
    url <- "https://www.dropbox.com/s/p8foxj6t6inb8uf/sce_pseudo_pathology_wholegenome.rds?dl=1"
  }
  else if (type == "Visium_SPG_AD_Visium_wholegenome_modeling_results") {
    tag <- "Visium_SPG_AD_Alzheimer_Disease_ITC_spatialLIBD"
    hub_title <- type
    file_name <- "Visium_IF_AD_modeling_results.Rdata"
    url <- "https://www.dropbox.com/s/5plupu8bj5m0kfh/Visium_IF_AD_modeling_results.Rdata?dl=1"
  }
  file_path <- file.path(destdir, file_name)
  if (!file.exists(file_path)) {
    q <- AnnotationHub::query(eh, pattern = c(tag, hub_title))
    if (length(q) == 1) {
      res <- q[[1]]
      if (type %in% c("sce", "sce_example")) {
        res <- .update_sce(res)
      }
      else if (type == "sce_layer") {
        res <- .update_sce_layer(res)
      }
      return(res)
    }
    else {
      file_path <- BiocFileCache::bfcrpath(bfc, url)
    }
  }
  message(Sys.time(), " loading file ", file_path)
  if (grepl(".Rdata", file_path)) {
    load(file_path, verbose = FALSE)
    if (type == "sce") {
      return(spatialLIBD:::.update_sce(sce))
    }
    else if (type == "sce_layer") {
      return(spatialLIBD:::.update_sce_layer(sce_layer))
    }
    else if (type == "modeling_results" || type == "spatialDLPFC_Visium_modeling_results" || 
             type == "Visium_SPG_AD_Visium_wholegenome_modeling_results") {
      return(modeling_results)
    }
    else if (type == "sce_example") {
      return(spatialLIBD:::.update_sce(sce_sub))
    }
    else if (type == "Visium_SPG_AD_Visium_wholegenome_spe" || 
             type == "Visium_SPG_AD_Visium_targeted_spe") {
      return(spe)
    }
  }
  else if (grepl(".rds", file_path)) {
    return(readRDS(file_path))
  }
  else {
    file_path
  }
}



myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

ehub <- ExperimentHub::ExperimentHub()
sce_LIBD <- my_fetch_data(type = "sce", eh = ehub,destdir = "E:/share_hdd/research/sce/")

sce_layer <- my_fetch_data(type = "sce_layer", eh = ehub)
modeling_results <- my_fetch_data("modeling_results", eh = ehub)

sample_id <- unique(sce_LIBD@colData@listData[["sample_name"]])
sample_names <- paste0("sp_",sample_id)



raw_dir <- "E:/share_hdd/research/spatialLIBD"
raw_files <- list.files(raw_dir, "h5", full.names = T)
sample_id <- gsub(".*LIBD/(.*)_raw.*","\\1",raw_files)

# Load raw matrix data
raw_matrix <- list()
for(i in seq_along(raw_files)){
  raw_matrix[[sample_names[i]]] <- Read10X_h5(raw_files[i])
}
names(raw_matrix) <- sample_names

# Load raw matrix data
#raw_matrix <- list()
#for(i in seq_along(sample_names)){
#  idx = sce_LIBD@colData@listData[["sample_name"]] == sample_id[i]
#  raw_matrix[[sample_names[i]]] <- sce_LIBD@assays@data@listData[["counts"]][,idx]
#}
#names(raw_matrix) <- sample_names

# Generate per-sample slides
sce_slide <- sce_LIBD@colData@listData
slide_list <- list()
for(dset in seq_along(sample_names)){
  slide_list[[sample_names[dset]]] <- sce_slide %>% data.frame %>% filter(sample_name==sample_id[dset]) %>% arrange()
}


#ref_slide <- slide_list_impute$`151507`
# df contains all spots
ref_slide <- readRDS("SpotClean_ref_slide.RDS")

slide_list_impute <- list()
for(dset in sample_names){
  slide_impute <- slide_list[[dset]] %>% select(barcode, tissue, row, col, imagerow, imagecol) %>% arrange(row, col)
  bg_slide <- ref_slide %>% select(barcode, tissue, row, col, imagerow, imagecol) %>%
    filter(!barcode%in%slide_impute$barcode) %>% 
    arrange(row, col)
  bg_slide$tissue <- 0
  
  # estimate average row distance and column distance
  lm_row <- lm(imagerow~row, data=slide_impute)
  row_dist <- coef(lm_row)[2]
  lm_col <- lm(imagecol~col, data=slide_impute)
  col_dist <- coef(lm_col)[2]
  
  # impute imagerow and imagecol of background spots
  bg_slide$imagerow <- predict(lm_row, data.frame(row=bg_slide$row))
  bg_slide$imagecol <- predict(lm_col, data.frame(col=bg_slide$col))
  
  slide_list_impute[[dset]] <- rbind(slide_impute, bg_slide)
  
  slide_list_impute[[dset]]$sum_umi <- colSums(raw_matrix[[dset]][,as.character(slide_list_impute[[dset]]$barcode)])
  slide_list_impute[[dset]]$width <- slide_list[[dset]]$width[1]
  slide_list_impute[[dset]]$height <- slide_list[[dset]]$height[1]
}

# Add layers to imputed slides
for(dset in sample_names){
  slide_list_impute[[dset]] <- merge(
    slide_list_impute[[dset]],
    slide_list[[dset]][,c("barcode","layer_guess")],
    by="barcode", all.x=TRUE)
  
  # Change barcodes to chr
  slide_list[[dset]]$barcode <- as.character(slide_list[[dset]]$barcode)
  slide_list_impute[[dset]]$barcode <- as.character(slide_list_impute[[dset]]$barcode)
}

# build slide object
#load("E:/share_hdd/research/spatialLIBD/12brains.RData")
#ref_slide <- slide_list_impute$`151507`

sets_tokeep <- c("sp_151507","sp_151508","sp_151669","sp_151670","sp_151673","sp_151674")


slide_obj_list <- list()
for(dset in sets_tokeep){
  slide_obj_list[[dset]] <- createSlide(raw_matrix[[dset]],slide_list_impute[[dset]])
}

for(dset in sets_tokeep){
  gene_tokeep <- keepHighGene(raw_matrix[[dset]])
  print(range(rowMeans(raw_matrix[[dset]][gene_tokeep,])))
}


decont_obj_list <- list()
for(dset in names(slide_obj_list)){
  decont_obj_list[[dset]] <- spotclean(slide_obj_list[[dset]])
}

saveRDS(decont_obj_list,"SpotCleaned_obj_ls.RDS")




