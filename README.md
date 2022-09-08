SCellBow
====

  This is a long description

   -   [Installation](#desktop-installation)
   -   [Tutorial](#vignette-tutorial)
       -  [Class Initialisation](#class-initialisation)
       -  [Pre-processing](#pre-processing)
       -  [Loading data](#setting-source-and-target-data)
       -  [Clustering](#clustering)
       -  [Visualizing](#visualizing-clusters)
       -  [Phenotype Algebra](#phenotype-algebra)
       -  [Differential gene analysis](#find-cluster-specific-differentially-expressed-genes)



Installation
===============

The developer version of the python package can be installed with the following commands:

```
pip install -i https://test.pypi.org/simple/ SCellBOW==0.0.1
```



Vignette tutorial
------------------
This vignette uses a small data set from the 10X website (3K PBMC dataset [here](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) ) to demonstrate a standard pipeline. This vignette can be used as a tutorial as well.


Class Initialisation
------------------
```python
scb = SCellBOW(read_path, write_path)
```


Pre-processing
--------------
dropClust performs pre-processing to remove poor quality cells and genes. dropClust is also equipped to mitigate batch-effects that may be present. The user does not need to provide any information regarding the source of the batch for individual transcriptomes. However, the batch-effect removal step is optional.

Cells are filtered based on the total UMI count in a cell specified by parameter `min_count`.  Poor quality genes are removed based on the minimum number of cells `min_count` with expressions above a given threshold `min_count`. 

``` python
adata_source = scb.preprocessing(path, min_genes, min_cells, target_sum, n_top_genes, max_value)
adata_target = scb.preprocessing(path, min_genes, min_cells, target_sum, n_top_genes, max_value)
```

Setting source and target data
-----------------------

dropClust loads UMI count expression data from three input files. The files follow the same structure as the datasets available from the 10X website, i.e.:

-   count matrix file in sparse format
-   transcriptome identifiers as a TSV file and
-   gene identifiers as a TSV file

``` python
scb.set_source_data(adata_source)
scb.set_target_data(adata_target)
```


Gene selection based on PCA
---------------------------
Another gene selection is performed to reduce the number of dimensions. PCA is used to identify genes affecting major components. 

``` python

# Find PCA top 200 genes. This may take some time.
sce<-RankPCAGenes(sce)

```


Clustering
------------------

### Fine tuning the clustering process

By default best-fit, Louvain based clusters are returned. However, the user can tune the parameters to produce the desired number of clusters. The un-sampled transcriptomes are assigned cluster identifiers from among those identifiers produced from fine-tuning clustering. The post-hoc assignment can be controlled by setting the confidence value `conf`. High `conf` values will assign cluster identifiers to only those transcriptomes sharing a majority of common nearest neighbours. 


``` python
# When `method = hclust`
# Adjust Minimum cluster size with argument minClusterSize (default = 20)
# Adjust tree cut with argument level deepSplit (default = 3), higher value produces more clusters.
sce<-Cluster(sce, method = "default", conf = 0.8)
```

Visualizing clusters
--------------------

Compute 2D embeddings for samples followed by post-hoc clustering.

``` python

sce<-PlotEmbedding(sce, embedding = "umap", spread = 10, min_dist = 0.1)

plot_data = data.frame("Y1" = reducedDim(sce,"umap")[,1], Y2 = reducedDim(sce, "umap")[,2], color = sce$ClusterIDs)

ScatterPlot(plot_data,title = "Clusters")
```

![](doc/vignette_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Find cluster specific Differentially Expressed genes
----------------------------------------------------

``` python

DE_genes_all = FindMarkers(sce, selected_clusters=NA, lfc_th = 1, q_th =0.001, nDE=30)

write.csv(DE_genes_all$genes, 
          file = file.path(tempdir(),"ct_genes.csv"),
          quote = FALSE)

```



## Visualizing clusters

Compute 2D embeddings for samples followed by post-hoc clustering.

``` python
ScatterPlot(dc.corr, title = "Clusters")
```

![Batch corrected dropClust based
Clustering.](doc/batchCorrection_files/figure-gfm/unnamed-chunk-5-1.png)

## Optional Batch correction

Users can use `fastmnn` method for batchcorrection. Specific arguments of fastmnn can also be passed through the `Correction` module.

``` python
merged_data.fastmnn<-Merge(all.objects,use.de.genes = FALSE)
set.seed(1)
mnn.corr <-  Correction(merged_data.fastmnn,  method="fastmnn", d = 10)
mnn.corr = Cluster(mnn.corr,method = "kmeans",centers = 3)
ScatterPlot(mnn.corr, title = "Clusters")
```

