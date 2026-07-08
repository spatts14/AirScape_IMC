## -----------------------------------------------------------------------------
library(batchelor)
library(scater)
library(cowplot)
library(dittoSeq)
library(viridis)
library(harmony)
library(BiocSingular)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(paletteer)

## -----------------------------------------------------------------------------
working_dir <- "/Volumes/phenotypingsputumasthmaticsaurorawellcomea1/live/Sara_Patti/all_IMC"
setwd(working_dir)
getwd()

# Create output directory
out_dir <- file.path(working_dir, "out/dimensionality_reduction")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## -----------------------------------------------------------------------------
# Reload the relevant objects
input_dir = file.path(working_dir, "data/spe_celltyped.rds")
spe <- readRDS(input_dir)


## -----------------------------------------------------------------------------
# Look at col data
cd <- colData(spe)

cd[] <- lapply(cd, function(x) {
  if (is.character(x)) {
    toupper(x)
  } else if (is.factor(x)) {
    factor(toupper(as.character(x)))
  } else {
    x
  }
})

colData(spe) <- cd
#colnames(colData(spe)) <- toupper(colnames(colData(spe)))
df <- as.data.frame(colData(spe))


# RENAME META
cd <- colData(spe)

cd$condition[cd$condition == "HEALTHY"] <- "PROXIMAL CONTROL"
cd$condition[cd$condition == "CONTROL"]  <- "DISTAL CONTROL"

colData(spe) <- cd


# Set color palette
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

# Create color vector container if it doesn't exist
if (is.null(metadata(spe)$color_vectors)) {
  metadata(spe)$color_vectors <- list()
}

metadata(spe)$color_vectors$condition <- c(
  "IPF"  = "#6A7FB5",
  "DISTAL CONTROL" = "#B07D4A",
  "COPD" = "#7EB0B8",
  "PROXIMAL CONTROL" = "#A67B8A"
)

metadata(spe)$color_vectors$diagnosis <- c(
  "IPF"         = "#6A7FB5",
  "LUNG_CANCER" = "#B07D4A",
  "COPD"        = "#7EB0B8",
  "HEALTHY"     = "#8EA882",
  "NO_CRD"      = "#A67B8A"
)

metadata(spe)$color_vectors$study <- c(
  "RBH"        = "#6A7FB5",
  "PM08"       = "#B07D4A",
  "REJUVENAIR" = "#7EB0B8",
  "MICA_III"   = "#A67B8A"
)

metadata(spe)$color_vectors$treatment_arm <- c(
  "SHAM"      = "#6A7FB5",
  "TREATMENT" = "#B07D4A"
)

metadata(spe)$color_vectors$timepoint <- c(
  "V1" = "#4A6699",
  "V2"  = "#B07D4A",
  "V3" = "#7EB0B8"
)

metadata(spe)$color_vectors$gender <- c(
  "MALE"   = "#6A7FB5",
  "FEMALE" = "#B07D4A"
)

metadata(spe)$color_vectors$lung_location <- c(
  "PROXIMAL" = "#4A6699",
  "DISTAL"   = "#B07D4A"
)

metadata(spe)$color_vectors$biopsy_type <- c(
  "CRYOBIOPSY"            = "#6A7FB5",
  "TRANSBRONCHIAL" = "#B07D4A",
  "RESECTION"             = "#7EB0B8",
  "ENDOBRONCHIAL"  = "#A67B8A"
)

metadata(spe)$color_vectors$smoker <- c(
  "NEVER"          = "#8EA882",
  "EX_SMOKER"      = "#C4A85A",
  "CURRENT_SMOKER" = "#A05A4A",
  "UNKNOWN"        = "#8A8A85"
)



# Plot UMAP
group_by = "patient_id"

# unique(colData(spe)$batch)

# CORRECTED UMAP
p1 <- dittoDimPlot(spe, var = group_by, reduction.use = "UMAP_seurat", size = 0.1) +
  scale_color_manual(values = metadata(spe)$color_vectors[[group_by]]) +
  ggtitle(paste0(group_by, " on UMAP"))
p1

# Save plot
filename <- paste0(group_by, "_UMAP_corrected.png")
ggsave(
  paste(out_dir, filename, sep = "/"),
  plot = p1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 900
)



# UNCORRECTED TSNE
p1 <- dittoDimPlot(spe, var = group_by, reduction.use = "TSNE", size = 0.2) +
  scale_color_manual(values = metadata(spe)$color_vectors[[group_by]]) +
  ggtitle(paste0(group_by, " on TSNE"))
p1

# Save plot
filename <- paste0(group_by, "_TSNE_uncorrected.png")
ggsave(
  paste(out_dir, filename, sep = "/"),
  plot = p1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 900
)

# UNCORRECTED UMAP
p1 <- dittoDimPlot(spe, var = group_by, reduction.use = "UMAP", size = 0.2) +
  scale_color_manual(values = metadata(spe)$color_vectors[[group_by]]) +
  ggtitle(paste0(group_by, " on UMAP"))
p1

# Save plot
filename <- paste0(group_by, "_UMAP_uncorrected.png")
ggsave(
  paste(out_dir, filename, sep = "/"),
  plot = p1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 900
)



# DROP NA FROM TREATMENT ARM
spe <- spe[, !is.na(spe$treatment_arm)]

group_by= "timepoint"

p1 <- dittoDimPlot(spe, var = group_by, reduction.use = "UMAP_seurat", size = 0.1) +
  scale_color_manual(values = metadata(spe)$color_vectors[[group_by]]) +
  ggtitle(paste0(group_by, " on UMAP"))
p1

# Save plot
filename <- paste0(group_by, "_UMAP_corrected.png")
ggsave(
  paste(out_dir, filename, sep = "/"),
  plot = p1,
  width = 8,
  height = 6,
  units = "in",
  dpi = 900
)


