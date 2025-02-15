### Stephano project
library(umap)
library(glue)
library(ggplot2)
library(Seurat)

# Load annotations - the file is in the github repository

setwd("/path/to/data")

PrepareData <- function(annotations, seurat_10x, genes) {
  data <- annotations
  
  fetched_data <- FetchData(object = seurat, vars = genes)
  
  for (gene in genes) {
    data[[gene]] <- fetched_data[[gene]]
    data[[paste0(gene, "_above_zero")]] <- data[[gene]] > 0 # This is for counting expressing cells
  }
  
  data["one"] <- 1 # This is for counting number of cells
  
  
  data
}

SubsampleNormalizeAge <- function(data) {
  # Subsampling to make young and aged the same amount
  young_cells <- dim(data[data$age == "Young",])[1]
  data.young <- data[data$age == "Young",]
  data.aged <- data[data$age == "Aged",]
  data.aged <- data.aged[sample(nrow(data.aged), young_cells),]
  
  data <- rbind.data.frame(data.young, data.aged)
}

GenerateDotplot <- function(data, condition, genes) {
  result_df = data.frame()
  for (gene in genes) {
    aggregation <- aggregate(
      data[,gene], 
      list(data[,condition]), 
      sum
    )
    
    groupCounts <- aggregate(
      data[,"one"], 
      list(data[,condition]),
      sum
    )
    
    groupAboveZeroCounts <- aggregate(
      data[, paste0(gene, "_above_zero")], 
      list(data[,condition]),
      sum
    )
    
    
    # aggregation$mean_value <- aggregation$x # Sum
    # aggregation$mean_value <- aggregation$x / groupAboveZeroCounts$x # Across expressed cells in condition
    # aggregation$mean_value <- aggregation$x / length(rownames(data)) # Across all CP cells
    aggregation$mean_value <- aggregation$x / groupCounts$x # Across all cells in condition
    
    aggregation$percentage <- groupAboveZeroCounts$x / groupCounts$x # Across condition
    # aggregation$percentage <- groupAboveZeroCounts$x / length(rownames(data)) # Across all CP cells
    
    aggregation["gene"] = gene # This is for naming the X axis
    
    result_df <- rbind.data.frame(result_df, aggregation)
  }
  
  ggplot(result_df, aes(x=gene, y=Group.1, colour=mean_value, size=percentage*100)) + 
    geom_point(binaxis='y', stackdir='center') +
    scale_colour_gradient(
      low = "#e6e6e6", high = "#ff6a00"
    ) + 
    theme_classic() +
    theme(
      axis.text = element_text(sizecr = 12), 
      axis.text.x = element_text(face = "italic"),
    ) +
    scale_size(range = c(4,18)) +
    labs(size = "Percent expressed\nper cell type", colour = "Mean expression\nper cell type", y=NULL, x=NULL) +
    guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
    scale_y_discrete(labels=c("Endothelial", "Epithelial", "Hematopoietic", "Mesenchymal", "Neuron associated"))
}

### Cells across celltypes
# annotations <- read.csv("annotations.csv") # Don't need this - we already have it in memory from annotations.RData file
data <- PrepareData(
  annotations, 
  seurat, 
  c('Cyp46a1', 'Ch25h', 'Cyp7a1', 'Cyp11a1', 'Cyp27a1', 'Tnfrsf1a', 'Tnfrsf1b', 'Nr1h2', 'Nr1h3', 'Cdh1', 'Sp1')
)
data <- SubsampleNormalizeAge(data)

ggplot(data, aes(y=age)) + geom_bar() # Make sure we equalized populations
number_of_cells <- nrow(data)

for (gene in c('Cyp46a1', 'Ch25h', 'Cyp7a1', 'Cyp11a1', 'Cyp27a1', 'Tnfrsf1a', 'Tnfrsf1b', 'Nr1h2', 'Nr1h3', 'Cdh1', 'Sp1')) {
  p <- GenerateDotplot(data, "cellType", gene)
  ggsave(glue("{gene}_{condition}.svg"), p, width=5, height=5)
}

### Nr1h couple Tnfrsf couple
p <- GenerateDotplot(data, "cellType", c('Nr1h2', "Nr1h3"))
ggsave(glue("Nr1h2 + Nr1h3 dotplot.svg"), p, width=6, height=5)

p <- GenerateDotplot(data, "cellType", c('Tnfrsf1a', 'Tnfrsf1b'))
ggsave(glue("Tnfrsf1a + Tnfrsf1b dotplot.svg"), p, width=6, height=5)

### Only epithelial cells (DO THIS AFTER GENERATING UMAPS, DATA IS NEEDED FOR ALL CELLTYPES)
# annotations <- read.csv("annotations.csv") # Don't need this - we already have it in memory from annotations.RData file
epithelial_data <- PrepareData(
  annotations, 
  seurat, 
  c('Cyp46a1', 'Ch25h', 'Cyp7a1', 'Cyp11a1', 'Cyp27a1', 'Tnfrsf1a', 'Tnfrsf1b', 'Nr1h2', 'Nr1h3', 'Cdh1', 'Sp1')
)
epithelial_data <- epithelial_data[epithelial_data$cellType == "Epithelial cell",] # Filter only epithelial cells
epithelial_data <- SubsampleNormalizeAge(epithelial_data)

ggplot(epithelial_data, aes(y=age)) + geom_bar() # Make sure we equalized populations
number_of_cells <- nrow(epithelial_data)

p <- GenerateDotplot(epithelial_data, "cellType", c('Cyp46a1', 'Ch25h', 'Cyp7a1', 'Cyp11a1', 'Cyp27a1'))
p = p + coord_flip() + theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(face = "plain")) + scale_y_discrete(labels=c("Epithelial"))
ggsave(glue("5_Cyp_genes_epithelial_cells.svg"), p, width=4, height=5)

p <- GenerateDotplot(epithelial_data, "cellType", c('Nr1h2', "Nr1h3"))
p = p + coord_flip() + theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(face = "plain")) + scale_y_discrete(labels=c("Epithelial"))
ggsave(glue("2_Nr1h_genes_epithelial_cells.svg"), p, width=4, height=5)

p <- GenerateDotplot(epithelial_data, "cellType", c('Tnfrsf1a', 'Tnfrsf1b'))
p = p + coord_flip() + theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(face = "plain")) + scale_y_discrete(labels=c("Epithelial"))
ggsave(glue("2_Tnfrsf_genes_epithelial_cells.svg"), p, width=4, height=5)


### MTX seurat. We already loaded this part from the state.RData file
# annotations <- read.csv("annotations.csv")
data_10x = Read10X(
  "/path/to/10x_data_folder", # Should contain matrix.mtx, genes.tsv, and barcodes.tsv files. 
  # You can find the files in the download section of https://singlecell.broadinstitute.org/single_cell/study/SCP1366/choroid-plexus-nucleus-atlas#study-download.
  # Rename the files [sc_rna_data.mtx, sc_rna_genes.tsv, sc_rna_barcodes.tsv] to the names above
  gene.column = 1
)
seurat <- CreateSeuratObject(counts = data_10x)

# Set Idents: age, ventricle, cellType
seurat@meta.data$condition_ventricle = annotations$ventricle
seurat@meta.data$condition_age = annotations$age
seurat@meta.data$condition_cellType = annotations$cellType
seurat@meta.data$condition_cellTypeSubset = annotations$cellTypeSubset

seurat <- SetIdent(object = seurat, value = "condition_vetricle") # It takes the value from seurat@meta.data$cluster_name
seurat <- SetIdent(object = seurat, value = "condition_age")
seurat <- SetIdent(object = seurat, value = "condition_cellType")
seurat <- SetIdent(object = seurat, value = "condition_cellTypeSubset")


# Filter out Embryos, and make aged same number as adults
subsetted_seurat <- seurat
seurat_young_subset = subset(subsetted_seurat, subset = condition_age == "Young")
seurat_aged_subset = subset(subsetted_seurat, subset = condition_age == "Aged")
num_of_young_cells <- ncol(seurat_young_subset)
seurat_aged_subset <- seurat_aged_subset[, sample(colnames(seurat_aged_subset), size=num_of_young_cells, replace=F)]
subsetted_seurat <- merge(seurat_aged_subset, seurat_young_subset)

# OR just use the same cells in data variable from before
subsetted_seurat <- subset(seurat, cells = data$cell)

# Make sure the age amounts are equalized
seurat_young_subset
seurat_aged_subset

# Downsample (don't do this for the finalized umaps)
# subsetted_seurat = subset(subsetted_seurat, downsample = 1000)

# Normalize + Scale
#subsetted_seurat = NormalizeData(subsetted_seurat) # Data is already normalized
subsetted_seurat <- FindVariableFeatures(subsetted_seurat)
subsetted_seurat = ScaleData(subsetted_seurat)

# Remove X and Y chromosome genes
var_features <- VariableFeatures(object = subsetted_seurat)
xy_genes_to_remove <- intersect(VariableFeatures(object = subsetted_seurat), c("Uba1","Rbmx","Zfx","Kdm5c","Hnrnpdl","Gapdh","Sycp3","Zfy1","Uba1y","Gm28588","Kdm5d","Eif2s3y","Gm29650","Uty","Ddx3y","Usp9y","Zfy2","Gm10256","Gm21693","Gm21704","Gm28919","Gm29554","Gm21454","Gm20831","Gm28595","Gm21865","Gm29433","Gm29276","Gm21739","Gm28887","Gm29221","Gm29206","Gm29209","Gm29060","Gm29380","Gm28612","Gm29082","Gm29644","Gm29566","Gm28664","Gm28930", "Xist", "Tsix"))
var_features <- var_features[-which(var_features %in% xy_genes_to_remove)]

# UMAP
subsetted_seurat <- RunUMAP(
  subsetted_seurat, 
  features = var_features,
  reduction.name = "umap_reduction_scaled",
  slot = "scale.data" # Using scaled data instead of just "data" by default
)

colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")

# UMAP by cell type, ventricle, age
p <- DimPlot(subsetted_seurat, group.by = "condition_cellType", reduction="umap_reduction_scaled", 
        cols = colors) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) + 
  ggtitle('UMAP by cell type')
ggsave(glue("UMAP_by_cell_type.svg"), p, width=9, height=7)

p <- DimPlot(subsetted_seurat, group.by = "condition_ventricle", reduction="umap_reduction_scaled", 
        cols = colors) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) + 
  ggtitle('UMAP by ventricle')
ggsave(glue("UMAP_by_ventricle.svg"), p, width=7.7, height=7)

p <- DimPlot(subsetted_seurat, group.by = "condition_age", reduction="umap_reduction_scaled", 
        cols = colors) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) + 
  ggtitle('UMAP by age')
ggsave(glue("UMAP_by_age.svg"), p, width=7.8, height=7)

# UMAP by each genes expression
for (gene in c('Cyp46a1', 'Ch25h', 'Cyp7a1', 'Cyp11a1', 'Cyp27a1', 'Tnfrsf1a', 'Tnfrsf1b', 'Nr1h2', 'Nr1h3', 'Cdh1', 'Sp1')) {
  p <- FeaturePlot(
    subsetted_seurat, 
    features = gene, 
    reduction="umap_reduction_scaled",
    cols = c("#e6e6e6", "#ff6a00"),
    pt.size = 1,
    order = TRUE
  ) + 
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  ggsave(glue("UMAP_{gene}.svg"), p, width=8, height=7)
}

# Scatter correlation plot
p <- FeatureScatter(
  subsetted_seurat,
  feature1 = "Cyp46a1",
  feature2 = "Sp1",
  group.by = "condition_cellType",
  cols = colors,
  slot = "data"
) + 
  ggtitle('Cyp46a1 x Sp1') + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(glue("Scatter_correlation_cyp46a1_x_sp1.svg"), p, width=7, height=5)

# Distribution of expressing cells across cellTypes
expressing_data = p$data[p$data$Cyp46a1 > 0 & p$data$Sp1 > 0,]
table(expressing_data$colors)

# Correlation test
cor.test(expressing_data$Cyp46a1, 
         expressing_data$Sp1,
         method = "spearman")



