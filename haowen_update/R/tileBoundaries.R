# From https://github.com/SydneyBioX/SpatialUtils/blob/main/R/tileBoundaries.R

#moleculesAssay = "detected"
#boundariesAssay = "tiles"
#tile_width = 80
#returnBoundariesOnly = FALSE

tileBoundaries = function(me, 
                          moleculesAssay = "detected",
                          boundariesAssay = "tiles",
                          tile_width = 0.1,
                          returnBoundariesOnly = FALSE) {
  
  require(dplyr)
  require(deldir)
  
  bds_list = list()
  # I would like to have something like samples(me) give the sample names
  # of the molecules, instead of
  samples = names(MoleculeExperiment::molecules(me, assayName = moleculesAssay)[[moleculesAssay]])
  for (sample in samples){
    
    # enabling something like 
    # me_i = me[[sample]]
    message(sample)
    # would be very nice
    # to subset the moleculeExperiment object to just that sample
    
    molecules_sub = MoleculeExperiment::molecules(me, flatten = TRUE, assayName = moleculesAssay) |>
      dplyr::filter(sample_id == sample)
    
    # select boundary seed points given grid structure
    # it would be nice if I could use 
    # molecules(me_i, flatten = TRUE)
    xr = range(molecules_sub$x_location)
    yr = range(molecules_sub$y_location)
    by = tile_width
    
    bds_points = expand.grid(x_location = seq(from = xr[1] - by/2, to = xr[2] + by, by = by),
                             y_location = seq(from = yr[1] - by/2, to = yr[2] + by, by = by))
    
    # Calculate Voronoi Tesselation and tiles
    tesselation <- deldir(bds_points$x_location, bds_points$y_location)
    tiles <- tile.list(tesselation)
    
    # the order of cells doesnt work if we have cell_id corresponding
    # to a numeric or integer value
    bds_i = do.call(rbind, lapply(tiles, function(x) data.frame(x_location = c(x$x, x$x[1]),
                                                                y_location = c(x$y, x$y[1]),
                                                                cell_id = paste0("cell_", x$ptNum))))
    bds_i$sample_id = sample
    
    bds_list[[sample]] <- bds_i
    
  }
  
  bds = do.call(rbind, bds_list)
  
  boundaries_ls <- dataframeToMEList(bds,
                                     dfType = "boundaries",
                                     assayName = "tiles",
                                     sampleCol = "sample_id",
                                     factorCol = "cell_id",
                                     xCol = "vertex_x",
                                     yCol = "vertex_y")
  
  if (returnBoundariesOnly) return(boundaries_ls)
  # str(boundaries_ls, max.level = 2)
  
  boundaries(me, "tiles") <- boundaries_ls
  
  return(me)
}
