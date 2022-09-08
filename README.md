SCellBow
====

  This is a long description

   -   [Installation](#installation)
   -   [Tutorial](#vignette-tutorial)
       -  [Class Initialisation](#class-initialisation)
       -  [Pre-processing](#pre-processing)
       -  [Loading data](#setting-source-and-target-data)
       -  [Model Creation](#model-creation)
       -  [Clustering and Visualisation](#clustering-and-visualisation)
       -  [Phenotype Algebra](#phenotype-algebra)
       -  [Differential gene analysis](#find-cluster-specific-differentially-expressed-genes)



Installation
=============

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
-------------------------------

dropClust loads UMI count expression data from three input files. The files follow the same structure as the datasets available from the 10X website, i.e.:

-   count matrix file in sparse format
-   transcriptome identifiers as a TSV file and
-   gene identifiers as a TSV file

``` python
scb.set_source_data(adata_source)
scb.set_target_data(adata_target)
```

Model Creation
---------------
```python
scb.SCellBOW_source(epochs, vec_size)
```


Clustering and Visualisation
--------------------------

### Fine tuning the clustering process

By default best-fit, Louvain based clusters are returned. However, the user can tune the parameters to produce the desired number of clusters. The un-sampled transcriptomes are assigned cluster identifiers from among those identifiers produced from fine-tuning clustering. The post-hoc assignment can be controlled by setting the confidence value `conf`. High `conf` values will assign cluster identifiers to only those transcriptomes sharing a majority of common nearest neighbours. 


``` python
scb.SCellBOW_test(self, svd_solver, n_neighbors, n_pcs, resolution)
```

Phenotype Algebra
----------------------------------------------------


Find cluster specific Differentially Expressed genes
----------------------------------------------------

``` python

????????????

```


