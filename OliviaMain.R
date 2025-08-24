if (!requireNamespace("BigVAR", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("wbnicholson/BigVAR/BigVAR", build_vignettes = FALSE, force = TRUE)
}

if (!requireNamespace("scCATCH", quietly = TRUE)) 
  install.packages("scCATCH")

install.packages("remotes")
remotes::install_github("njetzel/ngc", force = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

if (!requireNamespace("msigdbr", quietly = TRUE)) 
  install.packages("msigdbr")

if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
  BiocManager::install("clusterProfiler")

if (!requireNamespace("slingshot", quietly = TRUE)) 
  BiocManager::install("slingshot")

BiocManager::install("org.Hs.eg.db", force = TRUE)

library(BigVAR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cluster)
library(msigdbr)
library(clusterProfiler)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(stringr)
library(scCATCH)
library(rTensor)
library(NGC)
library(igraph)

counts <- Read10X_h5("/Users/oliviaxu/Desktop/GSM8552950_PA11-2/filtered_feature_bc_matrix.h5")

positions <- read.csv(
  "/Users/oliviaxu/Desktop/GSM8552950_PA11-2/spatial/tissue_positions.csv",
  header = FALSE
)
colnames(positions) <- c("barcode", "in_tissue", "array_row", "array_col",
                         "pxl_row_in_fullres", "pxl_col_in_fullres")
rownames(positions) <- positions$barcode

setwd("/Users/oliviaxu/Desktop/wavelet_project/4band_2reg")
wavelet_mat <- read.csv("wavelet_4band_2reg_view3_for_r.csv", header = TRUE, row.names = 1)

setwd("/Users/oliviaxu/Desktop/wavelet_project")
spot_ids <- readLines("spots.csv")
rownames(wavelet_mat) <- spot_ids
wavelet_mat <- t(as.matrix(wavelet_mat))
colnames(wavelet_mat) <- gsub("\\.", "-", colnames(wavelet_mat))

hvg_genes <- readLines("hvg.txt")
valid_genes <- intersect(hvg_genes, rownames(counts))
counts_filtered <- counts[valid_genes, spot_ids]

seu <- CreateSeuratObject(counts = counts_filtered)
seu <- AddMetaData(seu, metadata = positions[colnames(seu), ])
seu[["wavelet"]] <- CreateAssayObject(counts = wavelet_mat)

DefaultAssay(seu) <- "wavelet"
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = rownames(seu), npcs = 50)

pca_var <- seu[["pca"]]@stdev^2
var_exp <- cumsum(pca_var) / sum(pca_var)
num_pcs <- which(var_exp >= 0.9)[1]
cat(paste("Using", num_pcs, "PCs for clustering\n"))

seu <- FindNeighbors(seu, dims = 1:num_pcs)
seu <- FindClusters(seu, resolution = 0.7, algorithm = 4)
seu <- RunUMAP(seu, dims = 1:num_pcs)

DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 1) +
  labs(title = "Wavelet-Derived Clusters on UMAP") +
  theme_minimal()

library(colorspace)
levs <- levels(Idents(seu))
pal <- colorspace::rainbow_hcl(length(levs), c = 80, l = 70)
names(pal) <- levs

p_umap <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 1.2,
  alpha = 0.5
) +
  scale_color_manual(values = pal, limits = names(pal)) +
  labs(title = "Low", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.line = element_line(size = 1.0, colour = "black"),
    axis.ticks = element_line(size = 1.0, colour = "black"),
    axis.ticks.length = unit(6, "pt"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.box.margin = margin(0, 8, 0, 8),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    override.aes = list(size = 5)
  ))

print(p_umap)
#ggsave("~/Desktop/umap_without_dwt.png", p_umap, width = 7, height = 5.5, dpi = 600)


seu$pxl_row_in_fullres <- as.numeric(seu$pxl_row_in_fullres)
seu$pxl_col_in_fullres <- as.numeric(seu$pxl_col_in_fullres)

x_min <- min(seu$pxl_col_in_fullres)
x_max <- max(seu$pxl_col_in_fullres)

p_spatial <- ggplot(seu@meta.data, aes(
  x = pxl_col_in_fullres,
  y = pxl_row_in_fullres,
  color = seurat_clusters
)) +
  geom_point(size = 1.8, shape = 16) +
  scale_color_manual(values = pal, limits = names(pal), drop = TRUE) +
  labs(title = "High 3", x = "X (pixels)", y = "Y (pixels)", color = "Cluster") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0)),
    limits = c(x_min - 200, x_max + 200)
  ) +
  scale_y_reverse(expand = expansion(mult = c(0, 0))) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 30),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.line = element_line(size = 1.0, colour = "black"),
    axis.ticks = element_line(size = 1.0, colour = "black"),
    axis.ticks.length = unit(6, "pt"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    override.aes = list(size = 5)
  ))

print(p_spatial)

ggsave(
"~/Desktop/spatial_4band_2reg_high3.png",
p_spatial,
width = 7.5, height = 5.5, dpi = 300, bg = "white"
)

DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, assay = "RNA")

markers <- FindAllMarkers(
  seu,
  assay = "RNA",
  layer = "data",
  group.by = "seurat_clusters",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30) %>%
  ungroup()

print(top_markers)

cluster_genes <- split(top_markers$gene, top_markers$cluster)

muraro_pancreas <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  filter(str_starts(gs_name, "MURARO_PANCREAS")) %>%
  dplyr::select(gs_name, gene_symbol)

c7_all <- msigdbr(species = "Homo sapiens", category = "C7")

keywords <- c(
  "NEUTROPHIL","MACROPHAGE","MONOCYTE","DENDRITIC","NK","MAST","EOSINOPHIL",
  "BASOPHIL","GRANULOCYTE","ILC","MDSC","T_CELL","B_CELL","CD4","CD8","TREG",
  "CYTOTOXIC","HELPER","ANTIBODY","ANTIGEN_PRESENTING","MHC","IFN","IL1","IL2",
  "IL4","IL6","IL10","IL17","IL21","TH1","TH2","TH17","NAIVE","EFFECTOR",
  "ACTIVATED","EXHAUSTED","SUPPRESSOR","IMMUNOSUPPRESSIVE","INFLAMMATORY",
  "CYTOKINE","IMMUNE_RESPONSE"
)

max_per_keyword <- 15

final_c7_filtered <- purrr::map_dfr(keywords, function(kw) {
  hits <- c7_all %>%
    dplyr::filter(grepl(kw, gs_name, ignore.case = TRUE, fixed = TRUE))
  top_hits <- head(unique(hits$gs_name), max_per_keyword)
  dplyr::filter(hits, gs_name %in% top_hits)
})

c7_term2gene <- final_c7_filtered %>%
  dplyr::select(gs_name, gene_symbol)

stromal_sets <- tibble::tribble(
  ~gs_name,       ~gene_symbol,
  "F07_FIBRO_NRP2", "LUM",
  "F07_FIBRO_NRP2", "COL1A1",
  "F07_FIBRO_NRP2", "DCN",
  "F07_FIBRO_NRP2", "PDGFRA",
  "F07_FIBRO_NRP2", "NGFR",
  "F07_FIBRO_NRP2", "AUTS2",
  "F07_FIBRO_NRP2", "RGMA",
  "F07_FIBRO_NRP2", "LAMA2",
  "F07_FIBRO_NRP2", "CD34",
  "F07_FIBRO_NRP2", "TMEM100",
  "F07_FIBRO_NRP2", "SOX9",
  "F07_FIBRO_NRP2", "SOX5",
  "F07_FIBRO_NRP2", "GSK3B",
  "F07_FIBRO_NRP2", "APOD",
  "F07_FIBRO_NRP2", "ABCA10",
  "F07_FIBRO_NRP2", "FMO2",
  "CAF_GENERAL", "MMP11",
  "CAF_GENERAL", "SDC1",
  "CAF_GENERAL", "COL11A1",
  "CAF_GENERAL", "DPT",
  "CAF_GENERAL", "CFD",
  "CAF_GENERAL", "PI16",
  "CAF_GENERAL", "FAP",
  "CAF_GENERAL", "MME",
  "CAF_GENERAL", "TMEM158",
  "CAF_GENERAL", "NDRG1",
  "CAF_GENERAL", "HSPA1A",
  "CAF_GENERAL", "MKI67",
  "CAF_GENERAL", "MALAT1",
  "CAF_GENERAL", "COLEC11",
  "CAF_GENERAL", "NRP2",
  "CAF_GENERAL", "IGFBP5"
)

neural_sets <- tibble::tribble(
  ~gs_name,       ~gene_symbol,
  "SCHWANN_ALL", "ABCA8",
  "SCHWANN_ALL", "NEGR1",
  "SCHWANN_ALL", "AHR",
  "SCHWANN_ALL", "MPZ",
  "SCHWANN_ALL", "MBP",
  "SCHWANN_ALL", "NDRG1",
  "SCHWANN_ALL", "TGFBI",
  "SCHWANN_ALL", "FN1",
  "SCHWANN_ALL", "COL12A1",
  "SCHWANN_ALL", "CCN3",
  "SCHWANN_ALL", "FGF5",
  "SCHWANN_ALL", "MDK",
  "SCHWANN_ALL", "BTC",
  "SCHWANN_ALL", "SERPINA3",
  "SCHWANN_ALL", "JUN",
  "SCHWANN_ALL", "JUNB",
  "SCHWANN_ALL", "FOS",
  "SCHWANN_ALL", "FOSB",
  "SCHWANN_ALL", "GFAP",
  "SCHWANN_ALL", "CCL2"
)

combined_term2gene <- bind_rows(
  muraro_pancreas,
  c7_term2gene,
  stromal_sets,
  neural_sets
)

combined_term2name <- combined_term2gene %>%
  distinct(gs_name) %>%
  mutate(description = gs_name)


universe_genes <- rownames(seu)

all_results <- lapply(names(cluster_genes), function(cl) {
  genes <- cluster_genes[[cl]]
  res <- tryCatch(
    enricher(
      gene = genes,
      TERM2GENE = combined_term2gene,
      TERM2NAME = combined_term2name,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      universe = universe_genes
    ),
    error = function(e) NULL
  )
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    as.data.frame(res) %>% mutate(cluster = cl)
  } else NULL
})

combined_results <- bind_rows(all_results)

top_hits <- combined_results %>%
  filter(p.adjust < 0.2) %>%
  group_by(cluster) %>%
  slice_min(order_by = p.adjust, n = 1) %>%
  ungroup() %>%
  dplyr::select(cluster, Description, p.adjust) %>%
  arrange(as.integer(cluster))

ggplot(top_hits, aes(x = factor(cluster), y = -log10(p.adjust), fill = Description)) +
  geom_col() +
  labs(x = "Cluster", y = "-log10(adjusted p-value)",
       title = "Top Enriched Pathway per Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(top_hits, n = Inf)

seu_immune <- seu
sce <- as.SingleCellExperiment(seu_immune)

reducedDims(sce)$PCA <- Embeddings(seu_immune, "pca")[, 1:10]

sce <- slingshot(
  sce,
  clusterLabels = 'seurat_clusters',
  reducedDim = 'PCA',
  start.clus = names(sort(table(seu_immune$seurat_clusters), decreasing = TRUE))[1]
)

cl <- factor(colData(sce)$seurat_clusters)
levs <- levels(cl)
pal <- pal[levs]
cols <- pal[as.character(cl)]

#png("~/Desktop/slingshot_pseudotime_4band_2reg_high3.png", width = 3000, height = 2000, res = 300, bg = "white")

par(mar = c(5, 5, 4, 14), xpd = NA)

cols_alpha <- adjustcolor(cols, alpha.f = 0.5)
pal_alpha <- adjustcolor(pal, alpha.f = 0.5)

plot(
  reducedDims(sce)$PCA[, 1:2],
  col = cols_alpha,
  pch = 16,
  cex = 1,
  xlab = "PC1",
  ylab = "PC2",
  main = "High 3",
  cex.main = 2.2, font.main = 2,
  cex.lab = 1.8, font.lab = 2,
  cex.axis = 1.4, font.axis = 2,
  bty = "l", las = 1
)

lines(SlingshotDataSet(sce), lwd = 2, col = "black")

usr <- par("usr")
xleg <- usr[2] + 0.02 * diff(usr[1:2])
yleg <- mean(usr[3:4])

legend(
  x = xleg, y = yleg,
  xjust = 0, yjust = 0.5,
  legend = c("Lineage curve", paste("Cluster", levs)),
  col = c("black", pal_alpha),
  lty = c(1, rep(NA, length(levs))),
  lwd = c(2, rep(NA, length(levs))),
  pch = c(NA, rep(16, length(levs))),
  bty = "n",
  cex = 1.67,
  pt.cex = 2.6,
  y.intersp = 1,
  x.intersp = 0.5,
  text.font = 2
)

#dev.off()

pt_mat <- slingPseudotime(sce)

allowed_cols <- c("Celltype_TFH CD4 T Cell", "Celltype_Pancreas Acinar Cell")

top10_markers <- top_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()
marker_genes <- intersect(unique(top10_markers$gene), rownames(seu_immune))

cnt <- GetAssayData(seu_immune, assay = "RNA", slot = "counts")

setwd("/Users/oliviaxu/Desktop/wavelet_project/without_dwt")
labels <- read.csv("without_dwt_annotations.csv", header = FALSE, stringsAsFactors = FALSE)[[1]]

for (i in seq_len(ncol(pt_mat))) {
  lin <- as.character(i)
  
  pt <- pt_mat[colnames(seu_immune), i]
  valid_cells <- names(pt)[!is.na(pt)]
  ord <- order(pt[valid_cells], decreasing = FALSE)
  cells_all <- valid_cells[ord]
  
  cnt_all <- cnt[marker_genes, cells_all, drop = FALSE]
  M_lognorm_all <- t(as.matrix(log1p(cnt_all)))
  #write.csv(M_lognorm_all, paste0("spot_by_gene_lognorm_", lin, ".csv"), row.names = FALSE)
  
  lab_all <- as.integer(as.character(seu_immune$seurat_clusters[cells_all]))
  present_ids <- sort(unique(lab_all))
  idx <- if (any(present_ids == 0L)) present_ids + 1L else present_ids
  
  M_spot_by_ct_all <- do.call(cbind, lapply(present_ids, function(id) as.integer(lab_all == id)))
  colnames(M_spot_by_ct_all) <- labels[idx]
  rownames(M_spot_by_ct_all) <- cells_all
  #write.csv(M_spot_by_ct_all, paste0("spot_by_celltype_", lin, ".csv"), row.names = FALSE)
  }
}

