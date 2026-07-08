## -----------------------------------------------------------------------------
library(SpatialExperiment)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(scran)
library(CATALYST)
library(kohonen)
library(patchwork)
library(pheatmap)
library(dplyr)
library(cytomapper)
library(grid)
library(paletteer)
library(RColorBrewer)
library(openxlsx)
library(scater)
library(Seurat)
library(tidyverse)
library(scuttle)

## -----------------------------------------------------------------------------
# Set working directory
working_dir <- "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/all_IMC"
setwd(working_dir)
getwd()

# Create output directory
out_dir <- file.path(working_dir, "out/cell_type_SP")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## -----------------------------------------------------------------------------
# Set variables
id= "image_ID"
group_by= "slide_ID"
group_by_2= "condition"
group_by_3= "treatment_arm"
group_by_4= "patient_ID"
batch_correction= "seurat"
umap_correction= "UMAP_seurat"
k <- 80
phenotyped_spe <- "phenograph"
seed <- 2026

# Set .png
open_png <- function(file, width = 10, height = 8, res = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename = file, width = width, height = height,
                 units = "in", res = res, type = "cairo")
}
close_device <- function() try(grDevices::dev.off(), silent = TRUE)

## -----------------------------------------------------------------------------
# Set seed
set.seed(seed)

## -----------------------------------------------------------------------------
# Reload the relevant objects
input_dir = file.path("/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/ephemeral/IMC_output_folder_11_06_2026/06a_Cell_Labelling/spe_clusters.rds")
spe <- readRDS(input_dir)

cur_cells <- sample(seq_len(ncol(spe)), 2000)

## -----------------------------------------------------------------------------
# Set color palettes
cmap_colors_blue <- paletteer_c("grDevices::Oslo", 30)
pal <- as.character(paletteer::paletteer_c("ggthemes::Red-Gold", 30))

# Set colors for clusters
color_vectors <- list()
colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#3100AD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#969696", "#9ECAE1", "#A1D99B", "#FDBF6F", "#CAB2D6",
  "#FB9A99", "#B2DF8A", "#FFFF99", "#6A3D9A", "#B15928",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"
)

u <- unique(spe$pg_clusters_harmony)
pg_clusters_harmony_colors <- setNames(
  colors[seq_along(u)],
  u
)
color_vectors$pg_clusters_harmony <- pg_clusters_harmony_colors
metadata(spe)$color_vectors <- color_vectors

# Set color palette for patient IDs
palette_100 <- c(
  # blues / teals
  "#5B8FA8", "#8BB0C2", # controls 
  "#7EB5A6", "#6FAE9E", "#9BD0C0", "#5FA593","#3D6B5E", "#2F5C50", # MICA
  
  # warm oranges / browns
  "#C4956A", "#8F653E", "#D4A97A", "#E0B68B","#A67A52", # PM08
  
  # reds / salmon / rust
  "#C17B6E", "#8F4D3E", "#D28B7D", "#A95C4E", "#B56355",
  "#A85446", "#C97C6C", "#D4957F", "#E0A08D", "#C68A77",
  "#A05A4A", "#D9A08C",
  "#C88C78",
  
  # purples / lavenders
  "#8B7CB3", "#9A8BC0", "#7C6AA8", "#6E5B9C", "#A08DB8",
  "#B19CC6", "#8F7FB0", "#6B5E9E", "#5E5090", "#4F4282",
  "#B8A9CC", "#C7B8DA", "#A996C2", "#7A6EA8", "#695C96",
  "#C5B8D6", "#D3C7E2",
  
  # greens (T/NK-like tones)
  "#7E9E6E", "#8BAA7A", "#6F8F60", "#5F7F52", "#92B080",
  "#A0BC8E", "#7B996A", "#5E7E50", "#4F6E43", "#3F5E37",
  "#B5C99A", "#A8C08C", "#9AB67E", "#6D8C5A", "#5A784A",
  
  # yellows / immune neutrals
  "#B8A06E", "#C4AD7A", "#A9925F", "#9E8B5A", "#8C7A4A",
  "#D4C07A", "#7C8A7A", "#6E7C6C",
  
  # greys / structural / fallback
  "#9E9E8A", "#8F8F7C", "#A7A795", "#7F7F6D", "#7A9E9A",
  "#6F8F8B", "#8AA6A3", "#A0A0A0", "#8C8C8C", "#B0B0B0"
)
metadata(spe)$color_vectors$patient_id <-palette_100

## -----------------------------------------------------------------------------
# Plot clusters on UMAP
plot <- dittoDimPlot(
  spe,
  var = "pg_clusters_harmony",
  reduction.use = "UMAP_seurat",
  size = 0.5,
  do.label = TRUE
) +
  ggtitle("Phenograph clusters on UMAP") +
  scale_color_manual(values = metadata(spe)$color_vectors$pg_clusters_harmony)
plot

ggsave(
  file.path(out_dir, "UMAP_seurat_pg_clusters_harmony.png"),
  plot = plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

## -----------------------------------------------------------------------------
# Heatmap
open_png(file.path(
  out_dir,
  file.path("heatmap_harmony_", k, ".png")
), width = 10, height = 10, res = 300)

ht <- dittoHeatmap(
  spe,
  genes = rownames(spe)[rowData(spe)$use_channel],
  assay = "exprs", scale = "none",
  heatmap.colors = cmap_colors_blue,
  annot.by = c("pg_clusters_harmony"),
  main = paste0(
    "Phenograph clusters on heatmap k=",
    k,
    " Correction: ",
    batch_correction
  ),
  annot.colors = c(color_vectors$pg_clusters_harmony)
)

close_device()

ggsave(
  paste0(file.path, "heatmap_harmony_", k, ".png"),
  plot = plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

## -----------------------------------------------------------------------------
# Loop through each marker and create a violin plot
for (feature in rownames(spe)) {
  print(feature)
  plot <- plotExpression(
    spe,
    features = feature,  # Use the current feature
    x = "pg_clusters_harmony",
    exprs_values = "exprs",  # Use the appropriate expression values slot
    colour_by = "pg_clusters_harmony"
  )
  
  # Add title and theme settings
  plot <- plot +
    ggtitle(feature) +  # Set the title to the current feature
    theme(
      plot.title = element_text(size = 20),  # Adjust the title size
      axis.title.x = element_text(size = 15),  # Adjust the x-axis title size
      axis.title.y = element_text(size = 15),  # Adjust the y-axis title size
      # Adjust the x-axis text size and angle
      axis.text.x = element_text(size = 12, angle = 90),
      axis.text.y = element_text(size = 12)  # Adjust the y-axis text size
    ) +
    scale_color_manual(values = metadata(spe)$color_vectors$pg_clusters_harmony)
  
  # Save plot
  ggsave(
    file.path(out_dir, paste0("cluster_violinplot_", feature, ".png")),
    plot = plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
}

## -----------------------------------------------------------------------------
# Plot all violin plots in single figures
plot <- plotExpression(
  spe,
  features = rownames(spe)[rowData(spe)$use_channel],
  x = "pg_clusters_harmony",
  exprs_values = "exprs",
  colour_by = "pg_clusters_harmony"
) +
  scale_color_manual(
    values = metadata(spe)$color_vectors$pg_clusters_harmony
  ) +
  theme_classic()
plot

ggsave(
  file.path(out_dir, "cluster_violinplot_all.png"),
  plot = plot,
  width = 12,
  height = 30,
  units = "in",
  dpi = 300
)


## -----------------------------------------------------------------------------
# Unscaled dotplot
plot <- dittoDotPlot(
  spe,
  vars = rownames(spe)[rowData(spe)$use_channel],
  group.by = "pg_clusters_harmony",
  scale = FALSE,
  min.color = "#F1F1F1",
  max.color = "#2C5985",
  main = "Marker expression across clusters - unscaled",
  legend.color.title = "Colors title",
  legend.size.title = "Legend"
)
plot

ggsave(
  file.path(out_dir, "dotplot_cluster_all_noscale_blues.png"),
  plot = plot,
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

# Reds
plot <- dittoDotPlot(
  spe,
  vars = rownames(spe)[rowData(spe)$use_channel],
  group.by = "pg_clusters_harmony",
  scale = FALSE
) +
  scale_color_gradientn(colors = pal)
plot

ggsave(
  file.path(out_dir, "dotplot_cluster_all_noscale_reds.png"),
  plot = plot,
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

## -----------------------------------------------------------------------------
# Scaled dot plot
plot <- dittoDotPlot(
  spe,
  vars = rownames(spe)[rowData(spe)$use_channel],
  group.by = "pg_clusters_harmony",
  scale = TRUE,
  min.color = "#F1F1F1",
  max.color = "#2C5985",
  main = "Marker expression across clusters - scaled",
  legend.color.title = "Colors title",
  legend.size.title = "Legend"
)

ggsave(
  file.path(out_dir, "dotplot_cluster_all_scale_blues.png"),
  plot = plot,
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

# Reds
plot <- dittoDotPlot(
  spe,
  vars = rownames(spe)[rowData(spe)$use_channel],
  group.by = "pg_clusters_harmony",
  scale = TRUE
) +
  scale_color_gradientn(colors = pal)
plot

ggsave(
  file.path(out_dir, "dotplot_cluster_all_scale_reds.png"),
  plot = plot,
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

## -----------------------------------------------------------------------------
### Dotplot scaled and percent expression
# Ensure expression matrix and clusters are available
expr_matrix <- assays(spe)$exprs
clusters <- colData(spe)$pg_clusters_harmony
markers <- rownames(spe)[rowData(spe)$use_channel]

# Convert to binary (1 = expressed, 0 = not expressed)
# expression of the marker needs to be greater than 0.25
# to be considered positive
expr_binary <- expr_matrix > 0.25

# Create dataframe
expr_df <- as.data.frame(t(expr_binary))
expr_df$cluster <- clusters

# Calculate percent expressing per marker per cluster
percent_expr <- expr_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), ~ mean(.x) * 100))

# Convert wide to long format
percent_long <- pivot_longer(
  percent_expr,
  cols = -cluster,
  names_to = "marker",
  values_to = "percent_cells"
)

# Ensure marker names are correct
colnames(percent_expr)[-1] <- markers

# Compute mean expression per marker per cluster
expr_long <- as.data.frame(t(expr_matrix))  # Convert to dataframe
expr_long$cluster <- clusters

# Compute mean expression
mean_expr <- expr_long %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean))

# Reshape to long format
mean_expr_long <- pivot_longer(
  mean_expr,
  cols = -cluster,
  names_to = "marker",
  values_to = "mean_expression"
)

# Merge percent expression and mean expression
dotplot_data <- left_join(
  percent_long,
  mean_expr_long,
  by = c("cluster", "marker")
)

# Normalize percent_cells for dot size (0-1 scale)
dotplot_data <- dotplot_data %>%
  group_by(marker) %>%
  mutate(
    size_scaled = scales::rescale(percent_cells, to = c(0, 1))
  )

# Plot
plot <- ggplot(
  dotplot_data,
  aes(x = marker, y = cluster, size = size_scaled, color = mean_expression)
) +
  geom_point() +
  scale_size(range = c(0.1, 8)) +  # Adjust dot sizes
  scale_color_gradientn(colors = pal) +  # Use the Oslo color scale
  theme_classic() +
  labs(
    x = "",
    y = "",
    size = "Percent Expressing",
    color = "Mean Expression"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save plot
ggsave(file.path(out_dir,"dotplot_cluster_percent_expression.png"),plot = plot, width = 10, height = 8, units = "in", dpi = 300)

## -----------------------------------------------------------------------------
# Heatmap
# Define which variable to use for the heatmap's second annotation
group_by_heatmap <- "patient_id"

metadata(spe)$color_vectors$patient_id <- palette_100

# Make sure cluster IDs sort in the order you expect (avoids "1","10","11"... alphabetic sort issue)
spe$pg_clusters_harmony <- factor(
  spe$pg_clusters_harmony,
  levels = as.character(sort(as.numeric(unique(spe$pg_clusters_harmony))))
)
cluster_levels <- levels(spe$pg_clusters_harmony)

# Pull cluster colors from metadata, in matching order
cluster_colors_ordered <- metadata(spe)$color_vectors$pg_clusters_harmony[cluster_levels]

# Pull patient_id colors, sized and ordered to match actual levels present
group_levels <- sort(unique(spe[[group_by_heatmap]]))
group_colors_named <- setNames(
  metadata(spe)$color_vectors$patient_id[seq_along(unique(spe[[group_by_heatmap]]))],
  group_levels
)
group_colors_ordered <- group_colors_named[group_levels]

# Heatmap
plot <- dittoHeatmap(
  spe,
  genes = rownames(spe)[rowData(spe)$use_channel],
  assay = "exprs",
  scale = "none",
  heatmap.colors = cmap_colors_blue,
  annot.by = c("pg_clusters_harmony", group_by_heatmap),
  annot.colors = c(cluster_colors_ordered, group_colors_ordered)
)

ggsave(file.path(out_dir, "heatmap_clusters_short.png"), plot = plot,
       width = 10, height = 10, units = "in", dpi = 300)

## -----------------------------------------------------------------------------
# Check number of cells in each cluster

# NUMBER IN: pg_clusters_corrected
cluster_counts <- as.data.frame(table(spe$pg_clusters_corrected))
colnames(cluster_counts) <- c("Cluster", "Cell_count")
cluster_counts <- cluster_counts[order(cluster_counts$Cell_count), ]
cluster_counts

# NUMBER IN: pg_clusters_harmony
cluster_counts <- as.data.frame(table(spe$pg_clusters_harmony))
colnames(cluster_counts) <- c("Cluster", "Cell_count")
cluster_counts <- cluster_counts[order(cluster_counts$Cell_count), ]
cluster_counts

## -----------------------------------------------------------------------------
#### CELL TYPE MAPPING ###
celltype_mapping <- list(
  "1"  = "Neutrophils (CD66a+MPO+)",
  "2"  = "Basal cells (KI67 int)",
  "3"  = "Endothelial cells",
  "4"  = "Basal cells (KI67 neg)",
  "5"  = "Fibroblasts",
  "6"  = "Unknown 1",
  "7"  = "Basal cells (KI67 hi)",
  "8"  = "CD8+ T cells (epithelial)",
  "9"  = "Smooth muscle cells",
  "10" = "Airway macrophages",
  "11" = "Blood endothelial cells",
  "12" = "Unknown - Pi16hi",
  "13" = "Smooth muscle cells",
  "14" = "CD11b+ cells",
  "15" = "Smooth muscle cells",
  "16" = "CD4+ T cells",
  "17" = "Mast cells",
  "18" = "Goblet cells",
  "19" = "CD4+ T cells (epithelial)",
  "20" = "GATA3+PDGFRb+",
  "21" = "Airway macrophages, (Tryptase hi)",
  "22" = "Unknown 2 - KRT5+",
  "23" = "Unknown 3",
  "24" = "Epithelial cells (HLA-DRhi)",
  "25" = "CD8+ T cells",
  "26" = "Secretory epithelial cells (ECP hi)",
  "27" = "Dendritic cells",
  "28" = "B cells",
  "29" = "CD4+ T cells",
  "30" = "Goblet cells (CD11c hi)",
  "31" = "Neutrophils/Eosinophils (CD66a+)"
)

## -----------------------------------------------------------------------------
# Convert the mappings to a named vector (required for recode)
celltype_mapping <- unlist(celltype_mapping)

# Apply the mappings using recode
celltype <- recode(spe$pg_clusters_harmony, !!!celltype_mapping)
# Assign the cell types to the object
spe$celltype <- celltype

# Set cluster colors
cluster_colors <- metadata(spe)$color_vectors$pg_clusters_harmony

celltype_colors <- setNames(
  cluster_colors[names(celltype_mapping)],
  unname(celltype_mapping)
)

metadata(spe)$color_vectors$cluster_celltype <- celltype_colors

## -----------------------------------------------------------------------------
# Check number of cells in each cluster

# NUMBER IN: celltype
cluster_counts <- as.data.frame(table(spe$celltype))
colnames(cluster_counts) <- c("Cluster", "Cell_count")
cluster_counts <- cluster_counts[order(cluster_counts$Cell_count), ]
cluster_counts

## -----------------------------------------------------------------------------
# We can save the generated data objects
# for further downstream processing and analysis.
data_path = "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/all_IMC/data"
#saveRDS(spe, file.path(data_path, "spe_celltyped.rds"))

# Read saved spe
spe <- readRDS(file.path(data_path, "spe_celltyped.rds"))
