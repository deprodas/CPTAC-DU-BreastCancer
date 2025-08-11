setwd("C:\\03. BRCA ZHH\\02. Data analysis\\02. scRNA-seq (PBMC)")

# Author : Depro Das, Department of Biochemistry and Molecular Biology, University of Dhaka, Bangladesh 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(SingleCellExperiment)
library(Seurat)
library(scater)
library(biomaRt) 
library(openxlsx)
library(tidyverse) 
library(patchwork) 
library(SingleR) 
library(celldex) 
library(UCell) 
library(dittoSeq) 
library(ComplexHeatmap) 
library(circlize) 
library(ggplotify) 
library(presto) 
library(msigdbr) 
library(fgsea) 
library(ggplot2) 
library(viridis)
library(RColorBrewer) 
library(tibble) 

# ── Data input and integration ──────────────────────────────────────────────── 

# RDS files are SCE objects 

rds_files <- list.files(pattern = "\\.rds$")


# Map ensembl ids to gene symbols in SCE object 

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

convert_ensembl_to_symbols <- function(sce, mart) {
  ensembl_ids <- rownames(sce) 
  ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids) 
  gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                    filters = "ensembl_gene_id",
                    values = ensembl_ids_clean,
                    mart = mart) 
  
  gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  keep_idx <- ensembl_ids_clean %in% gene_map$ensembl_gene_id 
  sce <- sce[keep_idx, ]
  ensembl_ids_clean <- ensembl_ids_clean[keep_idx]
  gene_symbols <- gene_map$hgnc_symbol[match(ensembl_ids_clean, gene_map$ensembl_gene_id)]
  rownames(sce) <- gene_symbols
  
  sce <- sce[!duplicated(rownames(sce)), ] 
  return(sce)
}


# Process all SCEs (map genes and convert to seurat) 

seurat_list <- lapply(rds_files, function(file) {
  message("Processing: ", file)
  sce <- readRDS(file)
  sce <- convert_ensembl_to_symbols(sce, mart)
  
  sample_name <- sub("(_.*|\\.rds$)", "", file)
  sce$sample <- sample_name 
  
  seu <- as.Seurat(sce, counts = "counts", data = NULL)
  seu$sample <- sample_name
  
  return(seu)
})


# Steps before integration (both cca and rpca)

seurat_list <- lapply(seurat_list, function(seu) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000) 
  seu <- ScaleData(seu)
  seu <- RunPCA(seu) 
  return(seu)
})


# Perform integration

features <- SelectIntegrationFeatures(seurat_list)

anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                  reduction = "cca", 
                                  anchor.features = features) 

integrated <- IntegrateData(anchorset = anchors) 

head(rownames(integrated)) 


# Steps after integration  

DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

# Additionally insert sct assay 

integrated <- SCTransform(integrated, assay = "originalexp", new.assay.name = "SCT", verbose = FALSE)

# Plot to check 

p1 <- DimPlot(integrated, reduction = "umap", group.by = "intrinsic.subtype")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

# ── Observe metadata ────────────────────────────────────────────────────────── 

# Modify metadata names  

integrated <- readRDS("PBMC integrated (Mangiola).RDS")
colnames(integrated@meta.data) 

integrated@meta.data <- integrated@meta.data %>%
  dplyr::rename(Stage = STAGE.WHEN.BIOPSY.TAKEN..EBC.vs..MBC., 
                Treatment.response = Treatment.response.at.time.sample.taken..progressing..responding..stable.disease., 
                Menopausal.status = menopausal.status..premenopausal..perimenopausal..postmenopausal., 
                Neoadjuvant.therapy = neoadjuvant.therapy..specified., 
                Radiotherapy = radiotherapy..y.n., 
                Endocrine.therapy = endocrine.therapy..y.n., 
                Her2.therapy = her2.targeted.therapy..y.n., 
                Adjuvant.chemotherapy = adjuvant.chemotherapy..y.n., 
                Nr.of.previous.chemo = nr.of.previous.chemotherapies.for.metastatic.disease, 
                Sample.OMBC = AT.time.of.INDEX.sample.OMBC..1..de.novo..2..recurrent..3..residual..or.NO, 
                Comment.of.OMBC = Comment.of.OMBC.or.not, 
                CNS.other = CNS..other., 
                Other.specify = other..specify., 
                Alive.or.dead = alive.or.dead.at.last.follow.up, 
                Current.treatment = current.treatment.at.last.follow.up) 

# ── Cell annotation ───────────────────────────────────────────────────────────  

# SingleR annotation 

DefaultAssay(integrated) <- "originalexp" 

hpca.se <- HumanPrimaryCellAtlasData() 
count_matrix <- GetAssayData(integrated, slot = "data") 

singler_res <- SingleR(test = count_matrix, ref = hpca.se, labels = hpca.se$label.main)
integrated$singleR.annotation <- singler_res$labels

singler_p1 <- DimPlot(integrated, group.by = "singleR.annotation", reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE)
singler_p1 
ggsave("01. Singler annotation umap.pdf", plot = singler_p1, width = 9.6, height = 6, units = c("in")) 

table(integrated$singleR.annotation) 


# Marker-based annottaion (PBMC data - all are immune cells found in blood) 

DefaultAssay(integrated) <- "integrated" 
Idents(object = integrated) <- "seurat_clusters"

all_markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(all_markers, file = "Result cluster-wise markers (PBMC).csv", row.names = FALSE)

# Plot the markers 

DefaultAssay(integrated) <- "SCT"
Idents(object = integrated) <- "seurat_clusters" 

# Cluster 0, 5, 8, 12 - macrophage / monocyte 

monocyte <- c("CD68", "CD14", "LYZ", "FCGR3A", "S100A8", "S100A9", "VCAN", "MS4A7")
DotPlot(integrated, features = monocyte) 

# Cluster 13, 14 - dendritic cells 

DCs <- c("CD1C", "CLEC10A", "FCER1A", "CLEC9A", "XCR1", "LILRA4", "IL3RA", "IRF8", "TCF4", "CD207", "LAMP3")
DotPlot(integrated, features = DCs) 

# Cluster 2, 3, 6, 10, 11 - T cells  

t.cell <- c("CD45", "CD69", "CD25", "CD3D", "CD28", "CD3E", "CD3G", "CD4", "CD8A", "IL7R", "CCR7", "CD122", "CD45RA", "CD27", "CD127", "CD95") 
DotPlot(integrated, features = t.cell)

# Cluster 1, 4, 15 - NK cells 

nk.cell <- c("NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A", "GZMB", "PRF1")
DotPlot(integrated, features = nk.cell) 

# Cluster 7 - Bcells 

b.cell <- c("MS4A1", "CD19", "CD79A", "CD79B")
DotPlot(integrated, features = b.cell) 

# Cluster 16 - Proliferating plasmablasts / B cells 

plasmablast <- c("MKI67", "IGHA1", "IGHG1", "TOP2A") 
DotPlot(integrated, features = plasmablast) 

# Cluster 9 - Erythroblasts 

erythroblast <- c("HBB", "HBA1", "HBA2", "AHSP", "ALAS2", "HBM", "HBD", "HEMGN", "TRIM58", "TENT5C", "SLC25A37", "CA1", "SELENBP1", "BNIP3L")
DotPlot(integrated, features = erythroblast) 

# Cluster 17 - Platelets

platelet <- c("PPBP", "PF4", "ITGA2B", "GP9")
DotPlot(integrated, features = platelet) 

# Cluster 18 - Monocyte progenitors 

mo_prog <- c("SOX4", "GATA2", "CD34")
DotPlot(integrated, features = mo_prog) 

# Classify the cells - final classification 

integrated$cell_type <- ifelse(integrated$seurat_clusters %in% c(0, 5, 8, 12), "Monocyte/Macrophage",
                        ifelse(integrated$seurat_clusters %in% c(13, 14), "Dendritic cells",
                        ifelse(integrated$seurat_clusters %in% c(2, 3, 6, 10, 11), "T cells",
                        ifelse(integrated$seurat_clusters %in% c(1, 4, 15), "NK cells",
                        ifelse(integrated$seurat_clusters == 7, "B cells",
                        ifelse(integrated$seurat_clusters == 16, "Proliferating plasmablasts",
                        ifelse(integrated$seurat_clusters == 9, "Erythroblasts",
                        ifelse(integrated$seurat_clusters == 17, "Platelets",
                        ifelse(integrated$seurat_clusters == 18, "Monocyte progenitor", "Unknown")))))))))

# ── Plot everything ───────────────────────────────────────────────────────────  

# Plot metadata 

colnames(integrated@meta.data) 

sn_metavars <- c("type", "Stage.at.diagnosis", "intrinsic.subtype", "Treatment.response", "Menopausal.status", "Nr.of.previous.chemo", "Alive.or.dead", "Current.treatment", "integrated_snn_res.0.5", "seurat_clusters", "cell_type") 
cell_colors <- c("Monocyte/Macrophage" = "#2CA02C", "Dendritic cells" = "#9ACD32", "T cells" = "#FF69B4", "NK cells" = "#ff0090", "B cells" = "#6495ED", "Proliferating plasmablasts" = "#7600bc", "Erythroblasts" = "#FFB300", "Platelets" = "#D62728", "Monocyte progenitor" = "#800000")

umap_plots <- lapply(sn_metavars, function(group) {
  if (group == "cell_type") {
    DimPlot(integrated, group.by = group, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE, cols = cell_colors)
  } else {
    DimPlot(integrated, group.by = group, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE)
  }
})

p1.all <- wrap_plots(umap_plots, ncol = 3) 
ggsave("02. Metadata umap.pdf", plot = p1.all, width = 24, height = 18, units = c("in")) 


# Plot marker genes 

DefaultAssay(integrated) <- "SCT"
Idents(object = integrated) <- "cell_type"

final_marker <- c("CD14", "LYZ", "S100A9", "FCER1A", "IRF8", "TCF4", "IL7R", "CCR7", "CD3D", "NKG7", "GNLY", "KLRD1", "MS4A1", "CD79A", "CD79B", "MKI67", "IGHA1", "TOP2A", "HBB", "HBA1", "CA1", "PPBP", "PF4", "GP9", "SOX4", "GATA2", "CD34") 
cell_order <- c("Monocyte/Macrophage", "Dendritic cells", "T cells", "NK cells", "B cells", "Proliferating plasmablasts", "Erythroblasts", "Platelets", "Monocyte progenitor")

integrated@active.ident <- factor(integrated@active.ident, levels = cell_order)

p.mark_cell <- DotPlot(integrated, features = final_marker) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p.mark_cell 
ggsave("03. Markers cell-wise.pdf", plot = p.mark_cell, width = 10, height = 3, units = c("in")) 


# ── Explore immune-related genes ──────────────────────────────────────────────  

# Plot immune-related bulk-identified genes 

immune_genes <- read.xlsx("Selected immune related genes.xlsx") %>% pull(1) 

d.imm.genes <- DotPlot(integrated, features = immune_genes) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
d.imm.genes
ggsave("04. Immune-related genes dot (SCT).pdf", plot = d.imm.genes, width = 9, height = 3, units = c("in")) 


# Umap of immune-related genes  

DefaultAssay(integrated) <- "originalexp"

f.imm.genes <- FeaturePlot(integrated, 
                           features = immune_genes, 
                           min.cutoff = 0.3,
                           max.cutoff = 0.7,
                           raster = TRUE) 
ggsave("04. Immune-related genes feature (RNA).pdf", plot = f.imm.genes, width = 13, height = 24, units = c("in")) 


# ── Explore claudin-low subtype ───────────────────────────────────────────────   

# Ucell enrichment of the subtypes  

gene_dir <- "C:\\03. BRCA ZHH\\02. Data analysis\\01. RNA-seq (CPTAC, Cell 2020)"
subtype_genesets <- openxlsx::read.xlsx(file.path(gene_dir, "NMF breast subtype geneset.xlsx"))

set.seed(123) 
subtype_signatures_list <- lapply(subtype_genesets, as.character) 
integrated <- AddModuleScore_UCell(integrated, features = subtype_signatures_list)

feature_names <- c("Her2_UCell", "LumA_UCell", "LumB_UCell", "Basal_UCell", "Claudin_low_UCell") 
f.subtype_ucell <- FeaturePlot(integrated, 
                               reduction = "umap", 
                               features = feature_names, 
                               ncol = 3, 
                               min.cutoff = 0.03,
                               max.cutoff = 0.12,
                               raster = TRUE, 
                               order = TRUE) 
f.subtype_ucell 
ggsave(f.subtype_ucell, file = "05. Subtypes ucell feature.pdf", width = 11, height = 6, units = "in")

# Classify the data based on subtype enrichment  

# Category 1 

ucell_scores <- integrated@meta.data[, feature_names]
integrated$UCell_type <- colnames(ucell_scores)[apply(ucell_scores, 1, which.max)]
integrated$UCell_type <- gsub("_UCell", "", integrated$UCell_type) 

# Category 2 

ucell_category <- integrated@meta.data[, "UCell_type", drop = FALSE]
ucell_category <- ucell_category %>% mutate(claudin_class = ifelse(ucell_category$UCell_type == "Claudin_low", "Claudin_low", "Others")) 
integrated <- AddMetaData(integrated, metadata = ucell_category$claudin_class, col.name = "claudin_class")

# Plot the categories 

ucell_color <- c(Claudin_low = "#B4E50D", Basal = "#640D5F", Her2 = "#FF2DD1", LumA = "#4ED7F1", LumB = "#FF0B55")
claud_color <- c(Claudin_low = "#B4E50D", Others = "#640D5F")

d.ucell <- DimPlot(integrated, group.by = "UCell_type", cols = ucell_color, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) + 
  DimPlot(integrated, group.by = "claudin_class", cols = claud_color, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE)
d.ucell 
ggsave("06. Subtypes ucell umap.pdf", plot = d.ucell, width = 8.35, height = 3, units = "in")


# Calculate proportion - bar plots 

new_barplot <- function(feature, group) {
  dittoBarPlot(integrated, feature, group.by = group) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

compares <- list(c("UCell_type", "cell_type"), c("cell_type", "sample"), c("claudin_class", "sample"), c("cell_type", "claudin_class"))
barplots <- lapply(compares, function(x) new_barplot(x[1], x[2]))

p.b.all <- wrap_plots(barplots, ncol = 3) 
p.b.all 
ggsave("07. Barplot metadata.pdf", plot = p.b.all, width = 14, height = 9, units = "in")


# Correlation between immune-related gene expression and claudin-low scores 

immune_counts <- GetAssayData(integrated, slot = "data")[immune_genes, ] %>% as.data.frame() %>% t() %>% as.data.frame()
df_iu <- cbind(ucell_scores, immune_counts) 
genes <- setdiff(colnames(df_iu), feature_names)

cor_matrix <- sapply(feature_names, function(subtype) {
  sapply(genes, function(gene) {
    cor(df_iu[[subtype]], df_iu[[gene]], method = "pearson", use = "complete.obs")
  })
})

cor_matrix <- t(cor_matrix)
cor_matrix_clean <- cor_matrix[, colSums(is.na(cor_matrix)) == 0] 

set.seed(123) 
p.cor.1 <- Heatmap(cor_matrix_clean, 
                   col = colorRamp2(c(0, 0.3, 0.6), c("slateblue1", "white", "red")),
                   rect_gp = gpar(col = "white", lwd = 1), 
                   border = TRUE) 
p.cor.1
p.cor.1 <- as.ggplot(p.cor.1) 
ggsave("08. Immune & subtype correlation-1.pdf", plot = p.cor.1, width = 11, height = 2.9, units = "in")


# ── Macrophage states ─────────────────────────────────────────────────────────    

# Collect macrophage signatures: https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00094-1

tam_dir <- "C:\\03. BRCA ZHH\\01. Data search\\Signatures\\TAMs"
tam_genesets <- openxlsx::read.xlsx(file.path(tam_dir, "TAMs signatures.xlsx"))

# Run ucell 

set.seed(123) 
signatures_list_tam <- lapply(tam_genesets, as.character) 

integrated <- AddModuleScore_UCell(integrated, features = signatures_list_tam)
feature_names_tams <- paste0(names(signatures_list_tam), "_UCell")

f.tam_ucell <- FeaturePlot(integrated,
                           reduction = "umap",
                           features = feature_names_tams,
                           ncol = 3,
                           min.cutoff = 0.1,
                           max.cutoff = 0.5,
                           raster = TRUE,
                           order = TRUE)
f.tam_ucell 
ggsave(f.tam_ucell, file = "09. TAMs ucell feature.pdf", width = 10, height = 13, units = "in")


# Correlation between claudin-low and tams 

feature_tams_selected <- c("Claudin_low_UCell", "IFN_TAM_UCell", "Inflam_TAM_UCell", "LA_TAM_UCell", "Angio_TAM_UCell", "Reg_TAM_UCell", "Prolif_TAM_UCell", "Kupffer_RTM_like_UCell", "Alveolar_RTM_like_UCell", "Microglia_RTM_like_UCell", "MT_RTM_like_UCell", "Interstitial_RTM_like_UCell", "Classical_TIMs_UCell", "Nonclassical_Monocytes_UCell") 
meta_tam <- FetchData(integrated, vars = feature_tams_selected)

cor_mat = cor(meta_tam) 

# Plot corrplot 

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("forestgreen", "white", "deeppink"))

cell_fun_diag = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
  if(i == j) {
    grid.text(nm[i], x = x, y = y)
  } else if(i > j) {
    grid.circle(x = x, y = y, r = abs(cor_mat[i, j])/2 * min(unit.c(width, height)), gp = gpar(fill = col_fun(cor_mat[i, j]), col = NA))
  } else {
    grid.text(sprintf("%.1f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
} 

cell_fun_nodiag = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
  if(i > j) {
    grid.circle(x = x, y = y, r = abs(cor_mat[i, j])/2 * min(unit.c(width, height)), gp = gpar(fill = col_fun(cor_mat[i, j]), col = NA))
  } else if(i < j) {
    grid.text(sprintf("%.1f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
}

h.corr.1 <- Heatmap(cor_mat, 
                    col = col_fun, 
                    cell_fun = cell_fun_nodiag, 
                    name = "correlation", 
                    rect_gp = gpar(type = "none"), 
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE, 
                    show_row_names = TRUE, 
                    show_column_names = TRUE)
h.corr.1 
h.corr.1 <- as.ggplot(h.corr.1)
ggsave(h.corr.1, file = "10. TAMs & claudin correlation corrplot.pdf", width = 8, height = 7.2, units = "in")

# Plot heatmap 

h.corr.2 <- Heatmap(cor_mat, 
                    col = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "red")),
                    rect_gp = gpar(col = "white", lwd = 1), 
                    border = TRUE)
h.corr.2 
h.corr.2 <- as.ggplot(h.corr.2)
ggsave(h.corr.2, file = "11. TAMs & claudin correlation heatmap.pdf", width = 7.5, height = 6.75, units = "in")


# Save the object 

saveRDS(integrated, file = "PBMC integrated (Mangiola).RDS")


# ── Gene set enrichment ───────────────────────────────────────────────────────   

# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html 

# Fast wilcoxon test - presto 

claudin.gr.marker_prs <- wilcoxauc(integrated, group_by = 'claudin_class', assay = "data", seurat_assay = "originalexp") 

# Create ranked files 

claudin_genes <- claudin.gr.marker_prs %>% dplyr::filter(group == "Claudin_low") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc) 
others_genes <- claudin.gr.marker_prs %>% dplyr::filter(group == "Others") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc) 

claudin_ranks <- deframe(claudin_genes)
others_ranks <- deframe(others_genes)


# Fgsea - multiple groups (the more the rank file the more the gsea result) : group-wise top pathways 

top_fgsea_pathways <- function(ranks_list, category, subcategory = NULL, species = "Homo sapiens", nperm = 1000, top_n = 10) { 
  geneset_df <- msigdbr(species = species, category = category, subcategory = subcategory)
  fgsea_sets <- geneset_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  fgsea_results <- purrr::imap(ranks_list, function(ranks, label) { 
    fgsea(pathways = fgsea_sets, stats = ranks, nperm = nperm) %>% as_tibble() %>% arrange(desc(NES)) %>% mutate(group = label)
  })
  
  combined_results <- bind_rows(fgsea_results) 
  
  top_pathways <- combined_results %>%
    group_by(group) %>%
    slice_max(order_by = NES, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  return(top_pathways)
}


# Before runing the function it is important to make sure that labels naming is accurate 

collections <- msigdbr_collections() 
claudin.gr.marker_prs %>% dplyr::count(group)  

# Ranked files can be multiple, not only limited for two groups 

ranked_list <- list(Claudin_low = claudin_ranks, Others = others_ranks)

# Run fgsea  

top_h <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "H") 
top_h$database <- "Hallmark" 

top_gobp <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C5", subcategory = "GO:BP") 
top_gobp$database <- "GO_BP" 

top_gocc <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C5", subcategory = "GO:CC") 
top_gocc$database <- "GO_CC" 

top_gomf <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C5", subcategory = "GO:MF") 
top_gomf$database <- "GO_MF" 

top_rect <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C2", subcategory = "CP:REACTOME") 
top_rect$database <- "REACTOME" 

top_wiki <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C2", subcategory = "CP:WIKIPATHWAYS") 
top_wiki$database <- "WIKIPATHWAYS"

top_onco <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C6") 
top_onco$database <- "Oncogenic" 

top_immune <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C7", subcategory = "IMMUNESIGDB") 
top_immune$database <- "IMMUNESIGDB" 

top_cpbio <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C2", subcategory = "CP:BIOCARTA") 
top_cpbio$database <- "BIOCARTA" 

top_keggleg <- top_fgsea_pathways(ranks_list = ranked_list, top_n = 10, category = "C2", subcategory = "CP:KEGG_LEGACY") 
top_keggleg$database <- "KEGG_LEGACY" 

top_all <- bind_rows(top_h, top_gobp, top_gocc, top_gomf, top_rect, top_wiki, top_onco, top_immune, top_cpbio, top_keggleg)

# Plot bubble 

bub_all <- ggplot(top_all, aes(x = group, y = pathway, color = NES, size = -log10(pval))) +
  geom_point(alpha = 0.7) + 
  scale_color_viridis(option = "C") + 
  facet_wrap(~database, scales = "free", ncol = 2) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
bub_all
ggsave(bub_all, file = "12. Claudin classified pathways.pdf", width = 14, height = 28, units = "in")


# Plot heatmap 

top.new_h <- top_h %>%
  group_by(group) %>%
  slice_max(order_by = NES, n = 5, with_ties = FALSE) %>%
  ungroup() 

top.new_gobp <- top_gobp %>%
  group_by(group) %>%
  slice_max(order_by = NES, n = 5, with_ties = FALSE) %>%
  ungroup() 

top.new_rect <- top_rect %>%
  group_by(group) %>%
  slice_max(order_by = NES, n = 5, with_ties = FALSE) %>%
  ungroup() 

top.new_wiki <- top_wiki %>%
  group_by(group) %>%
  slice_max(order_by = NES, n = 5, with_ties = FALSE) %>%
  ungroup() 

top_fix <- bind_rows(top.new_h, top.new_gobp, top.new_rect, top.new_wiki) 

claudin_path <- top_fix %>% 
  filter(group == "Claudin_low") %>%
  group_by(pathway, database) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = database, values_from = NES, values_fill = NA) %>%
  column_to_rownames("pathway") %>%
  as.data.frame() 

others_path <- top_fix %>% 
  filter(group == "Others") %>%
  group_by(pathway, database) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = database, values_from = NES, values_fill = NA) %>%
  column_to_rownames("pathway") %>%
  as.data.frame() 

# Prepare annotation data 

claudin_ann_meta <- top_fix %>%
  filter(group == "Claudin_low") %>% 
  column_to_rownames(var = "pathway") %>%
  select(size, database) 

others_ann_meta <- top_fix %>%
  filter(group == "Others") %>% 
  column_to_rownames(var = "pathway") %>%
  select(size, database) 

all_databases <- unique(top_fix$database)
db_colors <- setNames(RColorBrewer::brewer.pal(length(all_databases), "Set1"), all_databases)

# Top pathways - claudin-low 

col_ha1_claudin <- rowAnnotation(Database = claudin_ann_meta$database, 
                                 col = list(Database = db_colors),
                                 annotation_name_side = "top", 
                                 border = TRUE) 

col_ha2_claudin <- rowAnnotation(Gene_count = anno_barplot(claudin_ann_meta$size, gp = gpar(fill = "yellow2"), bar_width = 0.8, border = FALSE), 
                                 annotation_name_side = "top") 

set.seed(123)
ht_claudin <- Heatmap(as.matrix(claudin_path), 
                      col = colorRampPalette(brewer.pal(9, "PiYG"))(10),
                      right_annotation = c(col_ha1_claudin, col_ha2_claudin),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      na_col = "white", 
                      border = TRUE, 
                      column_title_gp = gpar(fontsize = 10, fontface = "bold")) 
ht_claudin 
ht_claudin <- as.ggplot(ht_claudin)
ggsave(ht_claudin, file = "13. Claudin pathways heatmap.pdf", width = 5.75, height = 8, units = "in")


# Top pathways - others 

col_ha1_others <- rowAnnotation(Database = others_ann_meta$database, 
                                col = list(Database = db_colors),
                                annotation_name_side = "top", 
                                border = TRUE) 

col_ha2_others <- rowAnnotation(Gene_count = anno_barplot(others_ann_meta$size, gp = gpar(fill = "yellow2"), bar_width = 0.8, border = FALSE), 
                                annotation_name_side = "top") 

set.seed(123)
ht_others <- Heatmap(as.matrix(others_path), 
                     col = colorRampPalette(brewer.pal(9, "PiYG"))(10),
                     right_annotation = c(col_ha1_claudin, col_ha2_claudin),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     na_col = "white", 
                     border = TRUE, 
                     column_title_gp = gpar(fontsize = 10, fontface = "bold")) 
ht_others 
ht_others <- as.ggplot(ht_others)
ggsave(ht_others, file = "13. Others pathways heatmap.pdf", width = 5.75, height = 8, units = "in")


