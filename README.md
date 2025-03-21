# Glioblastoma Spatial Transcriptomics
Practicing using Seurat on a [glioblastoma dataset](https://www.10xgenomics.com/datasets/gene-and-protein-expression-library-of-human-glioblastoma-cytassist-ffpe-2-standard) to look at spatial changes in gene expression. 

## Steps
1. `Quality control`; removing spots with low gene detection and a high percentage of mitochondrial reads
2. `Normalisation` using SCTransform
3. `UMAP` visualisation clustering
4. `DE analysis` between different clusters

## Output files
Output file contains UMAP and spatial representations of spot clusters, as well as a spatial expression plot of one DE gene (CPLX2; Complexin 2) between clusters 3 & 7.
