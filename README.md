# Atlas-scale Spatial Clustering

This repository is used to share code for the atlas-scale spatial clustering CZI project. 

## Notes 

This is a [google doc linking to notes](https://docs.google.com/document/d/1KPqFb7MLW5MPU0u8n1NvsvgSispunT-df23IfOTkBqg/edit?usp=sharing) that we find about possible approaches to try for atlas-scale spatial clustering. A similar set of notes that Stephanie created when developing mbkmeans can be found [here](https://github.com/stephaniehicks/benchmark-hdf5-clustering/blob/2020-05-07-freeze/notes/NOTES.md). 

## Algorithms 

Specifically, we compare the following algorithms:  

- mini-batch _k_-means (non-spatial clustering)
    - available using the `mbkmeans()` function from the `mbkmeans` [Bioconductor](https://www.bioconductor.org/packages/mbkmeans) package. Can use data stored (i) in memory OR (ii) data stored on-disk (e.g. HDF5 files). 
- TBA



## Contributors

- Haowen Zhou
- Shila Ghazanfar
- Stephanie Hicks
- Boyi Guo
