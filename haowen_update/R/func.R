init <- function(){
  suppressMessages(library(tidyverse))
  suppressMessages(library(igraph))
  suppressMessages(library(psych))
  suppressMessages(library(reshape2))
  suppressMessages(library(ComplexHeatmap))
  suppressMessages(library(circlize))
  suppressMessages(library(clusterSim))
  message("Initialize done!")
}

myPreprocess <- function(data_mat, scale_coef = 10000){
  data_mat <- data_mat/matrix(rep(colSums(data_mat), nrow(data_mat)), 
                              nrow = nrow(data_mat),
                              byrow = F) * scale_coef
  data_mat
}


myCorGraph <- function(data_mat, threshold = 0.5){
  
}

swap <- function(vec){
  c(vec[2],vec[1])
}

