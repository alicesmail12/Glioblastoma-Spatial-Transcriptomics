# Packages
library('Seurat')
library('hdf5r')
library('tidyverse')
library('extrafont')
options(future.globals.maxSize = 8000 * 1024^2)

# GET DATA #####################################################################
# Create seurat image
seuratObj <- Load10X_Spatial(
  data.dir='C:/Users/c4064951/OneDrive - Newcastle University/Documents/Bioinformatics/Coding-Club/Spatial-Files',
  filename = "CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.h5",
  assay = "Spatial"
)

# Access different seurat attributes
seuratObj@assays$Spatial
seuratObj@meta.data$nCount_Spatial
seuratObj@meta.data$nFeature_Spatial

# PRE-PROCESSING ###############################################################
# Violin plot
VlnPlot(seuratObj, features = "nFeature_Spatial", ncol = 1, pt.size = 0) + 
  theme(text=element_text(family='Roboto'), axis.text.x=element_text(angle=0, hjust=0.5))+
  guides(fill='none')
VlnPlot(seuratObj, features = "nCount_Spatial", ncol = 1, pt.size = 0) + 
  theme(text=element_text(family='Roboto'), axis.text.x=element_text(angle=0, hjust=0.5))+
  guides(fill='none')

# Spatial plot
SpatialFeaturePlot(seuratObj, features = "nCount_Spatial") + theme(legend.position = "right") +
  theme(text=element_text(family='Roboto'), axis.text.x=element_text(angle=0, hjust=0.5))+
  scale_fill_gradientn(colours = c('#90bdcf','#75a450','#cfc963','#ba5346'), 
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) 

# UMI COUNTS & GENES ###########################################################
#' 3 metrics can be used to assess quality of spots in a spatial dataset
#' UMI Counts: number of transcripts per spot
#' Detected Genes: a low number of detected genes indicates low-quality spots
#' Mitochondrial Gene Expression: a high percentage of mitochondrial transcripts
#' can indicate stress or dying cells

# Get mito genes
counts <- seuratObj[["Spatial"]]$counts
MitoGenes <- grep("^MT-", rownames(counts), value = TRUE)

# Mitochondrial percentage
seuratObj[["mito_percent"]] <- (colSums(counts[MitoGenes, ]) / seuratObj[["nCount_Spatial"]]) * 100

# Plot QC histogram
ggplot(seuratObj@meta.data, aes(x=mito_percent)) +
  geom_histogram() +
  ggtitle("Histogram of Mitochondrial Percentage") +
  xlab("Mitochondrial Percentage") +
  ylab("Count")+
  theme_classic()+
  theme(text=element_text(family='Roboto'))

# Filter out spots that have more than 15% mitochondrial reads
seuratObj <- subset(seuratObj, mito_percent<15)

# Filter out spots that have less than 500 genes detected
seuratObj <- subset(seuratObj, nFeature_Spatial>500)

# Spatial plot
SpatialFeaturePlot(seuratObj, features = "mito_percent") + theme(legend.position = "right") +
  theme(text=element_text(family='Roboto'), axis.text.x=element_text(angle=0, hjust=0.5))+
  scale_fill_gradientn(colours = c('#90bdcf','#75a450','#cfc963','#ba5346'), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) 

# Filter out mitochondrial genes
seuratObj <- seuratObj[!grepl("^MT-", rownames(seuratObj)), ]

# NORMALISATION ################################################################
#seuratObj <- SCTransform(seuratObj, assay = "Spatial", verbose = TRUE)

# PCA/UMAP #####################################################################
# These are now standard steps in the Seurat workflow for visualization and clustering
seuratObj <- RunPCA(seuratObj, verbose = FALSE)
seuratObj <- RunUMAP(seuratObj, dims = 1:30, verbose = FALSE)

# Clustering
seuratObj <- FindNeighbors(seuratObj, dims = 1:30, verbose = FALSE)
seuratObj <- FindClusters(seuratObj, verbose = FALSE)

# Plot
palette <- list(colorRampPalette(colors = c('#ba5346','#cfc963','#75a450','#90bdcf'))(max(as.numeric(seuratObj@meta.data$seurat_clusters))+1))
DimPlot(seuratObj, label = TRUE)+ theme(legend.position = "right") +
  theme_classic()+
  theme(text=element_text(family='Roboto', size=14), axis.text.x=element_text(angle=0, hjust=0.5))+
  scale_colour_manual(values=palette[[1]]) +
  labs(x='UMAP1', y='UMAP2')+
  guides(colour='none')

# Where are clusters located
SpatialDimPlot(seuratObj, label = FALSE, label.size = 3)+
  theme_bw()+
  theme(text=element_text(family='Roboto'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())+
  scale_fill_manual(values=palette[[1]])+
  guides(fill='none')

# differential expression between cluster 7 and cluster 3
#de_markers <- FindMarkers(seuratObj, ident.1 = 7, ident.2 = 3)

# plot top marker
SpatialFeaturePlot(object = seuratObj, features = rownames(de_markers)[1:1],
                   alpha = c(1, 1), ncol = 1)+
  theme(text=element_text(family='Roboto'), legend.position = "right")+
  scale_fill_gradientn(colours = c('#90bdcf','#75a450','#cfc963','#ba5346'), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) 











