## Spatial Clustering

### Jun 06 Update

#### Methods related

##### Benchmark:

Benchmarking atlas-level data integration in single-cell genomics: https://www.nature.com/articles/s41592-021-01336-8

Benchmarking cell-type clustering methods for spatially resolved transcriptomics data: https://academic.oup.com/bib/article/24/1/bbac475/6835380

Benchmarking Computational Integration Methods for Spatial Transcriptomics Data: https://www.biorxiv.org/content/10.1101/2021.08.27.457741v2.full



##### Spatial Clustering Methods

**SpaGCN:** https://www.nature.com/articles/s41592-021-01255-8

<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-021-01255-8/MediaObjects/41592_2021_1255_Fig1_HTML.png" alt="Fig. 1" style="zoom: 25%;" />

Based on RGB channel, $$z_v = \frac{{r_v \times V_r + g_v \times V_g + b_v \times V_b}}{{V_r + V_g + V_b}}, V_i\rightarrow\text{variance of }i_v$$, $z_v^*=\frac{z_v-\mu_z}{\sigma_z}\times \max \left( {\sigma _x,\sigma _y} \right) \times s$

Distance $d\left( {u,v} \right) = \sqrt {(x_u - x_v)^2 + (y_u - y_v)^2 + (z_u^ \ast - z_v^ \ast )^2} .$

Edge weight: $w\left( {u,v} \right) = {{{\mathrm{exp}}}}\left( { - \frac{{d\left( {u,v} \right)^2}}{{2l^2}}} \right).$

GCN: $f\left( {{{{{X}}}},\,{{{{A}}}}} \right) = \delta \left( {A_{n\times n}X_{n\times 50}B_{50\times 50}} \right), $ 

$\text{A: Graph Adjacency Matrix, B: Filter Parameters, X: Top 50 PCs Embedding Matrix}$

Use Louvain to initialize cluster centroid 

**BayesSpace:** https://www.nature.com/articles/s41587-021-00935-2#Sec10

<img src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41587-021-00935-2/MediaObjects/41587_2021_935_Fig1_HTML.png" alt="figure 1"  />

Gibbs sampling for updating most of the parameters, Metropolis–Hastings algorithm updating clusters

##### Data Integration Methods - ST Data

*scVI (TBF)*

*MOFA+ (TBF)*

**~~Seurat v4~~ (WNN for multi-modality)**: https://www.sciencedirect.com/science/article/pii/S0092867421005833

Steps for Weighted Nearest Neighbor Analysis:

- (1) Constructing independent *k*-nearest neighbor (KNN) graphs for both modalities. 
- (2) Performing within and across-modality prediction 
  - With-in: $\hat{r}= \frac{\sum_{i\in \{\text{r KNN}\}} r_i}{k}, \hat{p}= \frac{\sum_{i\in \{\text{p KNN}\}} p_i}{k}$
  - Cross-modality:  $\hat{r}= \frac{\sum_{i\in \{\text{p KNN}\}} p_i}{k}, \hat{p}= \frac{\sum_{i\in \{\text{r KNN}\}} p_i}{k}$
- (3) Calculating cell-specific modality weights. 
  - $\theta_w(i,j) = w_r(i)\theta_r(i,j)+w_p(i)\theta_p(i,j)$
- (4) Calculating a WNN graph.

**Seurat v3** (Anchors)

![Figure thumbnail fx1](https://www.cell.com/cms/attachment/a61e553c-d662-414a-be80-a043cbdf03c3/fx1.jpg)

- Feature selection for integrated analysis of multiple datasets
- Identification of anchor correspondences between two datasets
- Anchor scoring & Weighting
- Data integration for reference assembly

*SNF (TBF)*

*CIMLR (TBF)*

**PASTE:** https://www.nature.com/articles/s41592-022-01459-6

<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-022-01459-6/MediaObjects/41592_2022_1459_Fig1_HTML.png" alt="Fig. 1" style="zoom: 33%;" />

$(X_{p\times n},D_{n\times n},g_{1\times n})$, X: expression matrix, D: Distance Matrix, g: cell weight vector

$F({{\Pi }}\,;\,X,D,X^{\prime} ,D^{\prime} ,c,\alpha )=(1-\alpha )\mathop{\sum}\limits_{i,j}c({x}_{\cdot i},{x}_{\cdot j}^{\prime}){\pi }_{ij}+\alpha \mathop{\sum}\limits_{i,j,k,l}{({d}_{ik}-{d}_{jl}^{\prime})}^{2}{\pi }_{ij}{\pi }_{kl}.$ c is a function that measures a nonnegative cost between the expression profiles of two spots over all genes.

Gromov–Wasserstein optimal transport problem

PASTE v2 (TBF)

**Ideal for replicates**

**PRECAST:** https://www.nature.com/articles/s41467-023-35947-w

<img src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41467-023-35947-w/MediaObjects/41467_2023_35947_Fig1_HTML.png" alt="figure 1" style="zoom: 80%;" />



To promote spatial smoothness in the space of cluster labels, we assume each latent class label, $y_{ri}$, is interconnected with the class labels of its neighborhoods via a discrete hidden Markov random field (HMRF).

##### Data Integration Methods - ST Data

*Scanorama (TBF)*

*fastMNN (TBF)*

**Harmony**

![Fig. 1](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41592-019-0619-0/MediaObjects/41592_2019_619_Fig1_HTML.png)

### Possible Approach

#### Two-level Clustering for data integration

##### Step 0: Preprocessing

Data representation:

Gene Expression (Scaled or raw counts?): $X_{n\times g}$

Euclidean Distance Matrix (symmetric): $D_{n\times n}$

Cell/Spot Similarity Matrix (symmetric): $M_{n\times n},m_{i,j}=f(x_i,x_j)$

Primary Weighted Undirected Graph Matrix: $G_{n\times n} = g(D,M)$

Function $g$ needs to be sensitive to $D$, therefore, for each spot, its nearest neighbors in the network should be physically near the spot. Thus the weighted graph can preserve the spatial information in later clustering steps

Filter out edges with low weight: $G'=G[G>v]$

##### Step 1: Build KNN from $G'$:

For each node, find its k neighbors with highest weights, remove connection to any other nodes. This step is designed to reduce the complexity of graphs and improve time & memory efficiency.

##### Step 2: Graph-based clustering with high resolution

Use Louvain to conduct graph-based clustering, increase the resolution to get about $m$ (50-100) small clusters (called niche?). 

We have the averaged gene expression matrix for each niche in each sample: $Y_{i}\in \mathbb{R}^{m\times g}$

Then, create an abstract graph for niches in each sample

##### Step 3: Graph integration

This step can be done directly using Seurat v3 data integration method, as the input will be the averaged expression data of each niche.

Or, calculate the similarity between two niches $(k,l)$ from two different samples $(i,j)$: $s_{i,j,k,l}=h(y_{ik},y_{jl})$ to 'connect' these two niches. Another threshold filter will also be applied to reduce the complexity of the integrated graph.

##### Step 4: Top-level clustering for cell-typing

Another round of Louvain clustering to assign cluster labels for each niche and for each cell within.





