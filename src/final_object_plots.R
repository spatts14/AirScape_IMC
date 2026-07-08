
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
# We can save the generated data objects
# for further downstream processing and analysis.
data_path = "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/all_IMC/data"

# Read saved spe
spe <- readRDS(file.path(data_path, "spe_celltyped.rds"))

## -----------------------------------------------------------------------------
# Plot umap with cell type names
plot <- dittoDimPlot(
  spe,
  var = "celltype",
  reduction.use = "UMAP_seurat",
  size = 0.5,
  do.label = TRUE,
  labels.size = 3
) +
  scale_color_manual(values = metadata(spe)$color_vectors$cluster_celltype)

# Save plot
ggsave(file.path(out_dir,"celltype_umap.png"),plot = plot, width = 12, height = 8, units = "in", dpi = 300)
ggsave(file.path(out_dir,"celltype_umap.pdf"),plot = plot, width = 12, height = 8, units = "in", dpi = 300)


## -----------------------------------------------------------------------------
# Heatmap
group_by_heatmap <- "patient_id"
metadata(spe)$color_vectors$patient_id <- palette_100

# Dedupe while preserving cluster-number order
celltype_order <- unique(unname(celltype_mapping[as.character(sort(as.numeric(names(celltype_mapping))))]))
spe$celltype <- factor(spe$celltype, levels = celltype_order)
cluster_levels <- levels(spe$celltype)

# Celltype colors, keyed by unique labels (first color wins for shared labels)
cluster_colors_ordered <- metadata(spe)$color_vectors$cluster_celltype[cluster_levels]

# Patient_id colors
group_levels <- sort(unique(spe[[group_by_heatmap]]))
group_colors_named <- setNames(
  metadata(spe)$color_vectors$patient_id[seq_along(unique(spe[[group_by_heatmap]]))],
  group_levels
)
group_colors_ordered <- group_colors_named[group_levels]

# Heatmap: order cells by celltype, no column clustering/dendrogram at all
plot <- dittoHeatmap(
  spe,
  genes = rownames(spe)[rowData(spe)$use_channel],
  assay = "exprs",
  scale = "none",
  heatmap.colors = cmap_colors_blue,
  annot.by = c("celltype", group_by_heatmap),
  annot.colors = c(cluster_colors_ordered, group_colors_ordered),
  cluster_cols = FALSE,
  show_colnames = FALSE   # cell IDs on x-axis are meaningless at this scale anyway
)
plot
ggsave(file.path(out_dir, "heatmap_expression_celltype_patient_.png"), plot = plot,
       width = 10, height = 15, units = "in", dpi = 300)

## -----------------------------------------------------------------------------
### Dotplot scaled and percent expression
# Ensure expression matrix and clusters are available
expr_matrix <- assays(spe)$exprs
clusters <- colData(spe)$celltype
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
ggsave(file.path(out_dir,"dotplot_cluster_percent_expression_celltype.png"),
       plot = plot, width = 12, height = 10, units = "in", dpi = 300)

## -----------------------------------------------------------------------------
# Plot cluster composition...
group_list <- c(group_by, group_by_2, group_by_3, group_by_4)

for (view in group_list){
  print(paste0("Plotting ", view, " for compositional analysis..."))
  # Cell populations composition
  plot <- dittoBarPlot(
    spe,
    var = "celltype",
    group.by = view
  ) +
    scale_fill_manual(values = metadata(spe)$color_vectors$cluster_celltype)
  
  ggsave(
    file.path(out_dir, paste0("Composition_", view, "_cell_type.png")),
    plot = plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
  )
}


## -----------------------------------------------------------------------------
### Visualize cells
# Ensure factors are consistent
spe$celltype <- factor(spe$celltype)
celltype_colors <- metadata(spe)$color_vectors$cluster_celltype

# Create output directory if needed
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Iterate over all samples
for (sample in unique(spe$sample_id)) {
  
  idx <- spe$sample_id == sample
  
  df <- data.frame(
    Pos_X = spatialCoords(spe)[idx, "Pos_X"],
    Pos_Y = spatialCoords(spe)[idx, "Pos_Y"],
    celltype = droplevels(spe$celltype[idx])
  )
  
  plot <- ggplot(df, aes(Pos_X, Pos_Y, color = celltype)) +
    geom_point(size = 0.4) +
    coord_fixed() +
    scale_y_reverse() +
    scale_color_manual(
      values = celltype_colors,
      drop = FALSE
    ) +
    theme_void() +
    ggtitle(sample)
  
  ggsave(
    filename = file.path(out_dir, paste0("celltype_", sample, "_visualize.png")),
    plot = plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
  )
}


## -----------------------------------------------------------------------------
# Plot specific cell type
cell_to_show <- "Airway macrophages"
spe$celltype <- factor(spe$celltype)
celltype_colors <- metadata(spe)$color_vectors$cluster_celltype

# Safe folder name
safe_name <- gsub("[^A-Za-z0-9]+", "_", cell_to_show)
safe_name <- gsub("^_|_$", "", safe_name)

out_dir_celltype <- file.path(out_dir, safe_name)
dir.create(out_dir_celltype, recursive = TRUE, showWarnings = FALSE)

# Iterate over samples
for (sample in unique(spe$sample_id)) {
  
  idx <- spe$sample_id == sample
  
  df <- data.frame(
    Pos_X = spatialCoords(spe)[idx, "Pos_X"],
    Pos_Y = spatialCoords(spe)[idx, "Pos_Y"],
    celltype = droplevels(spe$celltype[idx])
  )
  
  # Color mapping: highlight one cell type, rest grey
  levs <- levels(spe$celltype)
  
  color_map <- setNames(rep("grey80", length(levs)), levs)
  
  if (cell_to_show %in% levs) {
    color_map[cell_to_show] <- celltype_colors[cell_to_show]
  } else {
    warning(sprintf("'%s' not found in celltype levels", cell_to_show))
  }
  
  # Plot
  plot <- ggplot(df, aes(Pos_X, Pos_Y, color = celltype)) +
    geom_point(size = 0.8) +
    coord_fixed() +
    scale_y_reverse() +
    scale_color_manual(values = color_map, drop = FALSE) +
    theme_void() +
    ggtitle(sample)
  
  # Save
  ggsave(
    filename = file.path(
      out_dir_celltype,
      paste0("highlight_", safe_name, "_", sample, ".png")
    ),
    plot = plot,
    width = 15,
    height = 10,
    units = "in",
    dpi = 300
  )
}

## -----------------------------------------------------------------------------
# Ensure factor is set once
spe$celltype <- factor(spe$celltype)

# Color mapping
celltype_colors <- metadata(spe)$color_vectors$cluster_celltype

# All cell types to iterate over
all_celltypes <- levels(spe$celltype)

for (cell_to_show in all_celltypes) {
  
  # Safe folder name
  safe_name <- gsub("[^A-Za-z0-9]+", "_", cell_to_show)
  safe_name <- gsub("^_|_$", "", safe_name)
  
  out_dir_celltype <- file.path(out_dir, safe_name)
  dir.create(out_dir_celltype, recursive = TRUE, showWarnings = FALSE)
  
  message("Processing: ", cell_to_show)
  
  # Iterate over samples
  for (sample in unique(spe$sample_id)) {
    
    idx <- spe$sample_id == sample
    
    df <- data.frame(
      Pos_X = spatialCoords(spe)[idx, "Pos_X"],
      Pos_Y = spatialCoords(spe)[idx, "Pos_Y"],
      celltype = droplevels(spe$celltype[idx])
    )
    
    # Base color map: all grey
    levs <- levels(spe$celltype)
    color_map <- setNames(rep("grey80", length(levs)), levs)
    
    # Highlight current cell type if it exists
    if (cell_to_show %in% levs) {
      color_map[cell_to_show] <- celltype_colors[cell_to_show]
    } else {
      warning(sprintf("'%s' not found in celltype levels", cell_to_show))
    }
    
    # Plot
    p <- ggplot(df, aes(Pos_X, Pos_Y, color = celltype)) +
      geom_point(size = 0.4) +
      coord_fixed() +
      scale_y_reverse() +
      scale_color_manual(values = color_map, drop = FALSE) +
      theme_void() +
      ggtitle(paste(sample, "-", cell_to_show))
    
    # Save
    ggsave(
      filename = file.path(
        out_dir_celltype,
        paste0("highlight_", safe_name, "_", sample, ".png")
      ),
      plot = p,
      width = 15,
      height = 10,
      units = "in",
      dpi = 300
    )
  }
}

# -----------------------------------------------------------------------------
target_cells <- c("CD4+ T cells", "CD8+ T cells", "B cells")

# Base color map: all grey
levs <- levels(spe$celltype)
color_map <- setNames(rep("grey80", length(levs)), levs)

# Assign distinct colors to multiple selected cell types
selected_present <- intersect(target_cells, levs)

missing <- setdiff(target_cells, levs)
if (length(missing) > 0) {
  warning("Missing cell types: ", paste(missing, collapse = ", "))
}

# Use predefined palette or fallback colors
highlight_colors <- celltype_colors[selected_present]
names(highlight_colors) <- selected_present

color_map[selected_present] <- highlight_colors

p <- ggplot(df, aes(Pos_X, Pos_Y, color = celltype)) +
  geom_point(size = 0.4) +
  coord_fixed() +
  scale_y_reverse() +
  scale_color_manual(values = color_map, drop = FALSE) +
  theme_void() +
  ggtitle(paste(sample, paste(target_cells, collapse = ", ")))
p

