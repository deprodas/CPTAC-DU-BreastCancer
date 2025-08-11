setwd("C:\\03. BRCA ZHH\\02. Data analysis\\01. RNA-seq (CPTAC, Cell 2020)")

# Author : Depro Das, Department of Biochemistry and Molecular Biology, University of Dhaka, Bangladesh 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(openxlsx) 
library(tidyverse) 
library(dplyr) 
library(tibble) 
library(tidyr) 
library(ggplot2) 
library(ggpubr) 
library(ggplotify)
library(ggrepel) 
library(patchwork) 
library(decoupleR) 
library(rlang) 
library(ComplexHeatmap) 
library(circlize) 
library(RColorBrewer) 
library(viridis) 
library(limma) 
library(edgeR) 
library(fgsea) 
library(msigdbr) 
library(readr) 
library(stringr) 
library(tools)
library(clusterProfiler) 
library(tidyestimate) 
library(grid) 
library(purrr) 

# ── Data load and manipulation ──────────────────────────────────────────────── 

# Count countdata  

countdata <- read.table(file = "data_mrna_seq_fpkm.txt", sep = "\t", header = TRUE)
countdata[ , -1] <- lapply(countdata[ , -1], function(x) as.numeric(as.character(x))) 
countdata$Hugo_Symbol <- make.names(countdata$Hugo_Symbol, unique = TRUE) 
countdata <- countdata %>% remove_rownames() %>% column_to_rownames(var = "Hugo_Symbol")

# Metadata 

metadata.pt <- read.table(file = "data_clinical_patient.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata.sm <- read.table(file = "data_clinical_sample.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- merge(metadata.pt, metadata.sm, by = "PATIENT_ID")
metadata <- metadata %>% column_to_rownames(var = "PATIENT_ID") 

# Subset count and metadata for matching 

metadata <- metadata %>% subset(rownames(metadata) %in% colnames(countdata)) 
countdata <- countdata %>% dplyr::select(rownames(metadata))
all(colnames(countdata) %in% rownames(metadata))
all(rownames(metadata) %in% colnames(countdata))

# Convert NA values to zero (0) in count matrix 

any(is.na(countdata)) 
all(is.numeric(countdata)) 
countdata[is.na(countdata)] <- 0 
any(is.na(countdata)) 
all(is.numeric(countdata)) 

# ── Claudin group classification ────────────────────────────────────────────── 

# Custom gene set 

genesets <- openxlsx::read.xlsx("Claudin-low geneset.xlsx")
genesets <- genesets %>% dplyr::select(Gene_symbol) %>% dplyr::rename(claudin.low = "Gene_symbol")
genesets <- genesets %>% pivot_longer(c(claudin.low), names_to = "source", values_to = "target")
genesets <- genesets %>% mutate(weight = 1) %>% filter(target!=" ") %>% filter(target!="")

# Prepare count and metadata 

counts <- countdata %>% as.matrix() 
is.na(counts)
design <- metadata

# Run GSEA 

res_decoupler <- decoupleR::decouple(mat = counts, 
                                     network = genesets, 
                                     .source ='source', 
                                     .target ='target',
                                     minsize = 0) 

res_ssgsea <- decoupleR::run_gsva(mat = counts, 
                                  network = genesets, 
                                  .source ='source', 
                                  .target ='target', 
                                  minsize = 2L, 
                                  method = c("ssgsea")) # "gsva", "plage", "ssgsea", "zscore"


# Each metric tried (ssGSEA and ULM provided the same results in terms of percentage of basal-like and high-grade proportions)  

# Claudin-low group classification (greater equal then the 3rd quartile) 

res_decoupler %>% dplyr::count(statistic)

res_long <- res_decoupler %>% 
  dplyr::filter(statistic %in% "ulm") %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(claudin.low.enrichment = "claudin.low") 

res_long <- res_long %>% mutate(claudin.groups = ifelse(claudin.low.enrichment >= quantile(claudin.low.enrichment, 0.75, na.rm = TRUE), "claudin.low", "Others"))
meta.res <- cbind(res_long, metadata) 


# Tumor grade classification (low : I, II, high : III, IV) 

colnames(meta.res) 
meta.res %>% dplyr::count(TUMOR_STAGE)

meta.res$our.tumor.class <- dplyr::case_when(meta.res$TUMOR_STAGE %in% c("Stage IA", "Stage IIA", "Stage IIB") ~ "Low_Operable",
                                             meta.res$TUMOR_STAGE %in% c( "Stage III", "Stage IIIA", "STAGE IIIB", "STAGE IIIB", "STAGE IIIC") ~ "High_Inoperable",
                                             TRUE ~ "Unknown") 

# ── Check proportion ────────────────────────────────────────────────────────── 

# Function for pie chart proportion 

plot_pie <- function(data, column, group_col = "claudin.groups") {
  column_sym <- sym(column)
  group_col_sym <- sym(group_col) 
  
  plot_data <- data %>%
    group_by(!!group_col_sym, !!column_sym) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(!!group_col_sym) %>%
    mutate(percentage = count / sum(count) * 100) 
  
  ggplot(plot_data, aes(x = "", y = percentage, fill = !!column_sym)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    facet_wrap(vars(!!group_col_sym)) +
    labs(title = paste("Distribution of", column, "by", group_col), fill = column, y = "Percentage", x = NULL) +
    theme_void() +
    theme(strip.text = element_text(face = "bold"))
}

# Function for bar chart proportion 

plot_bar <- function(data, column, group_col = "claudin.groups", position = c("fill", "stack")) {
  position <- match.arg(position)
  column_sym <- sym(column)
  group_col_sym <- sym(group_col)
  
  plot_data <- data %>%
    count(!!group_col_sym, !!column_sym) %>%
    group_by(!!group_col_sym) %>%
    mutate(prop = n / sum(n), percent_label = paste0(round(prop * 100, 1), "%")) %>%
    ungroup() 
  
  contingency_tbl <- table(data[[group_col]], data[[column]])
  if (any(chisq.test(contingency_tbl)$expected < 5)) {
    test_result <- fisher.test(contingency_tbl)
    test_name <- "Fisher's exact"
  } else {
    test_result <- chisq.test(contingency_tbl)
    test_name <- "Chi-squared"
  }
  pval <- signif(test_result$p.value, 3)
  p_label <- paste0(test_name, " p = ", pval)
  
  p <- ggplot(plot_data, aes(x = !!group_col_sym, y = n, fill = !!column_sym)) +
    geom_bar(stat = "identity", position = position)
  if (position == "fill") {
    p <- p + 
      geom_text(aes(label = percent_label, y = prop / 2), position = position_stack(vjust = 0.5), size = 3, color = "white") +
      scale_y_continuous(labels = scales::percent_format())
    y_label <- "Proportion"
  } else {
    p <- p + 
      geom_text(aes(label = percent_label, y = n / 2), position = position_stack(vjust = 0.5), size = 3, color = "white")
    y_label <- "Count"
  } 
  p <- p +
    labs(title = paste("Distribution of", column, "by", group_col), subtitle = p_label, fill = column, x = group_col, y = y_label) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Proportion plots 

class(meta.res)
colnames(meta.res) 

plot_pie(meta.res, column = "PAM50", group_col = "claudin.groups")
plot_pie(meta.res, column = "TUMOR_STAGE", group_col = "claudin.groups") 
plot_pie(meta.res, column = "our.tumor.class", group_col = "claudin.groups") 

p1 <- plot_bar(meta.res, column = "NMF_CLUSTER", group_col = "claudin.groups") 
p2 <- plot_bar(meta.res, column = "TUMOR_STAGE", group_col = "claudin.groups") 
p3 <- plot_bar(meta.res, column = "our.tumor.class", group_col = "claudin.groups") 
p4 <- plot_bar(meta.res, column = "CANCER_TYPE_DETAILED", group_col = "claudin.groups") 
p5 <- plot_bar(meta.res, column = "PAM50", group_col = "claudin.groups", position = "fill")

p.all <- p1 + p2 + p3 + p4 + p5
p.all

dir.create("01. Classification", showWarnings = FALSE)
ggsave(filename = "01. Classification/01. Proportions (claudin groups).pdf", plot = p.all, width = 12, height = 9, units = "in")


# ── Enrichment gradient with metadata ───────────────────────────────────────── 

meta.ord <- meta.res %>% rownames_to_column(var = "Patients")
colnames(meta.ord)

meta.ord %>% dplyr::count(NMF_CLUSTER)
meta.ord %>% dplyr::count(PAM50)
meta.ord %>% dplyr::count(TUMOR_STAGE)
meta.ord %>% dplyr::count(our.tumor.class)

# Color palettes  

gradient_colors <- list(CIBERSORT_ABSOLUTE_SCORE = c("white", "skyblue", "navy"),
                        ESTIMATE_IMMUNE_SCORE = c("white", "lightgreen", "darkgreen"),
                        STEMNESS_SCORE = c("white", "orange", "darkred"),
                        ESTIMATE_STROMAL_SCORE = c("white", "orchid", "purple4"))

category_colors <- list(NMF_CLUSTER = c("Basal-I" = "tomato", "HER2-I" = "gold", "LumA-I" = "deepskyblue", "LumB-I" = "#D50B8B"),
                        PAM50 = c("Basal" = "forestgreen", "Her2" = "orange", "LumA" = "red", "LumB" = "purple", "Normal-like" = "grey70"),
                        our.tumor.class = c("High_Inoperable" = "darkblue", "Low_Operable" = "darkred", "Unknown" = "grey70"),
                        TUMOR_STAGE = c("Stage IA" = "#F8766D", "Stage IIA" = "#7CAE00", "Stage IIB" = "#00BFC4", "Stage III" = "#C77CFF", "Stage IIIA" = "#FF2DF1", "Stage IIIB" = "#4ED7F1", "Stage IIIC" = "#FA812F"))

# Functions for plots 

plot_numeric <- function(varname, colors) {
  ggplot(meta.ord_ordered, aes(x = Patients, y = 1, fill = .data[[varname]])) +
    geom_tile() +
    scale_fill_gradientn(colors = colors) +
    theme_void() +
    theme(legend.position = "bottom", plot.margin = margin(0, 5, 0, 5)) +
    labs(fill = varname)
}

plot_categorical <- function(varname, colors) {
  ggplot(meta.ord_ordered, aes(x = Patients, y = 1, fill = .data[[varname]])) +
    geom_tile() +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(legend.position = "bottom", plot.margin = margin(0, 5, 0, 5)) +
    labs(fill = varname)
}

# Order by claudin-low enrichment

meta.ord_ordered <- meta.ord %>%
  arrange(desc(claudin.low.enrichment)) %>%
  mutate(Patients = factor(Patients, levels = Patients))

# Create the barplot on which the other annotation plots will be added 

p_claudin <- ggplot(meta.ord_ordered, aes(x = Patients, y = claudin.low.enrichment, fill = claudin.groups)) +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = quantile(meta.ord_ordered$claudin.low.enrichment, 0.75, na.rm = TRUE), linetype = "dashed", color = "red", size = 0.8) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(5, 5, 2, 5)) +
  ylab("Claudin-low\nEnrichment")

# Create and combine the annotation plots

num_vars <- names(gradient_colors)
cat_vars <- names(category_colors)

num_plots <- mapply(plot_numeric, varname = num_vars, colors = gradient_colors, SIMPLIFY = FALSE) 
cat_plots <- mapply(plot_categorical, varname = cat_vars, colors = category_colors, SIMPLIFY = FALSE) 
all_plots <- list(p_claudin) %>% append(cat_plots) %>% append(num_plots)

# Fix the heights (main plot = 4x)

rel_heights <- c(4, rep(1, length(cat_plots) + length(num_plots)))
combined_plot <- wrap_plots(all_plots, ncol = 1, heights = rel_heights)

# Final version of the plot 

f.plot <- print(combined_plot)
ggsave(filename = "01. Classification/02. Gradient barplot (claudin enrichment).pdf", plot = f.plot, width = 12, height = 9, units = "in")


# ── PCA plot generation ─────────────────────────────────────────────────────── 

# Use log2-normalized data 

pca_res <- prcomp(t(counts), scale. = TRUE)
pca_df <- data.frame(pca_res$x[, 1:2])
pca_df$patients <- rownames(pca_df)
pca_df <- cbind(pca_df, meta.res) 

plot.pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = claudin.groups)) +
  stat_ellipse(aes(color = claudin.groups), type = "norm", level = 0.95, size = 0.75) +  # Add ellipses
  geom_point(size = 3, pch = 21) +
  xlab(paste0("PC1: ", round(100 * summary(pca_res)$importance[2, 1], 1), "%")) +
  ylab(paste0("PC2: ", round(100 * summary(pca_res)$importance[2, 2], 1), "%")) +
  theme_bw()
plot.pca
ggsave(filename = "01. Classification/03. PCA (claudin groups).pdf", plot = plot.pca, width = 5.5, height = 4, units = c("in")) 


# ── Limma analysis ────────────────────────────────────────────────────────────  

dir.create("02. Limma", showWarnings = FALSE)
limma_dir <- "02. Limma"

# Create new group for limma condition 

meta.res <- meta.res %>% dplyr::mutate(limma.condition = paste(claudin.groups, our.tumor.class, sep = "_"))
meta.res %>% dplyr::count(limma.condition) 
all(colnames(countdata) == rownames(meta.res)) 


# Prepare count matrix for fpkm counts to avoild any errors due to the negative values 

# Limitation : As raw counts could not be retrieved from the fpkm values, limma analysis was limited to genes with non-zero fpkm values 

# Be careful : Avoid log-normalization more than once, as it can distort the expression data and impact downstream analyses 

counts.lim <- countdata %>% as.matrix() 
counts.lim[counts.lim < 0] <- NA 

# Run limma 

run_limmaFPKM <- function(group_col_name, group1, group2, fpkm_matrix, coldata, outname_prefix = "Result_Limma") {
  group_col <- coldata[[group_col_name]] 
  selected <- group_col %in% c(group1, group2)
  sub_fpkm <- fpkm_matrix[, selected]
  sub_meta <- coldata[selected, ]
  
  log_fpkm <- sub_fpkm # use log2(sub_fpkm + 1) instead of only sub_fpkm if the fpkm values are not log2 transformation
  
  group <- factor(sub_meta[[group_col_name]], levels = c(group1, group2))
  design <- model.matrix(~ group)
  fit <- lmFit(log_fpkm, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, adjust = "fdr") 
  
  out_filename <- paste0(outname_prefix, "_", gsub("[^a-zA-Z0-9]", "_", group1), "_vs_", gsub("[^a-zA-Z0-9]", "_", group2), ".xlsx")
  write.xlsx(res, file.path(limma_dir, out_filename), row.names = T)
  return(res)
} 

# group1 = control condition, group2 = experimental condition 

meta.res %>% dplyr::count(claudin.groups)
meta.res %>% dplyr::count(our.tumor.class) 
meta.res %>% dplyr::count(limma.condition) 

lim.res1 <- run_limmaFPKM(group_col_name = "claudin.groups", 
                          group1 = "Others", 
                          group2 = "claudin.low",  
                          fpkm_matrix = counts.lim, 
                          coldata = meta.res, 
                          outname_prefix = "Result_Limma") 

lim.res2 <- run_limmaFPKM(group_col_name = "our.tumor.class", 
                          group1 = "Low_Operable", 
                          group2 = "High_Inoperable",  
                          fpkm_matrix = counts.lim, 
                          coldata = meta.res, 
                          outname_prefix = "Result_Limma") 

lim.res3 <- run_limmaFPKM(group_col_name = "limma.condition", 
                          group1 = "claudin.low_Low_Operable", 
                          group2 = "claudin.low_High_Inoperable",  
                          fpkm_matrix = counts.lim, 
                          coldata = meta.res, 
                          outname_prefix = "Result_Limma") 

lim.res4 <- run_limmaFPKM(group_col_name = "limma.condition", 
                          group1 = "Others_Low_Operable", 
                          group2 = "Others_High_Inoperable",  
                          fpkm_matrix = counts.lim, 
                          coldata = meta.res, 
                          outname_prefix = "Result_Limma") 

# Visualize the results from each limma comparison 

comparisons.limma <- list(list(name = "lim.res1", data = lim.res1, group1 = "Others", group2 = "claudin.low"),
                          list(name = "lim.res2", data = lim.res2, group1 = "Low_Operable", group2 = "High_Inoperable"),
                          list(name = "lim.res3", data = lim.res3, group1 = "claudin.low_Low_Operable", group2 = "claudin.low_High_Inoperable"),
                          list(name = "lim.res4", data = lim.res4, group1 = "Others_Low_Operable", group2 = "Others_High_Inoperable"))

# Volcano plotting loop

for (comp in comparisons.limma) {
  limmadata <- comp$data %>% rownames_to_column(var = "SYMBOL") %>% na.omit()
  limmadata$nlog10 <- -log10(limmadata$P.Value) 
  limmadata$expression <- ifelse(limmadata$P.Value < 0.8912509 & abs(limmadata$logFC) >= 1, 
                                 ifelse(limmadata$logFC > 1, 'Up-regulated', 'Down-regulated'), 'Stable') 
  print(dplyr::count(limmadata, expression))
  
  top <- 10
  top_limma.genes <- bind_rows(limmadata %>%
                                 filter(expression == 'Up-regulated') %>%
                                 arrange(P.Value, desc(abs(logFC))) %>%
                                 head(top),
                               limmadata %>%
                                 filter(expression == 'Down-regulated') %>%
                                 arrange(P.Value, desc(abs(logFC))) %>%
                                 head(top)) 
  
  p_vol <- ggplot(data = limmadata, aes(x = logFC, y = nlog10, fill = expression)) +
    geom_vline(xintercept = c(-1 , 1), lty = 2, col = "black", lwd = 0.5) +
    geom_hline(yintercept = 0.05, lty = 2, col = "black", lwd = 0.5) +
    geom_point(size = 1.5, alpha = 1, pch = 21, stroke = 0.0001) +
    scale_fill_manual(values = c("#FD841F", "#D6D5A8", "#38E54D")) +
    theme_bw() + 
    xlab("logFC") +
    ylab("-log10(P-value)") +
    ggrepel::geom_text_repel(aes(label = SYMBOL), data = top_limma.genes, size = 2, force = 10, fontface = "italic", max.overlaps = Inf)
  
  filename <- paste0(limma_dir, "/Volcano_", comp$group1, "_vs_", comp$group2, ".pdf")
  ggsave(filename = filename, plot = p_vol, width = 5, height = 4, units = "in")
}


# ── Run fast gene set enrichment ──────────────────────────────────────────────

dir.create("03. GSEA", showWarnings = FALSE)
gsea_dir <- "03. GSEA"

# Get MSigDB genesets 

gene.categories <- list(Hallmark = msigdbr(species = "Homo sapiens", category = "H"),
                        # KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"),
                        GOBP = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
                        Reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
                        Wiki = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")) 

geneset_list <- lapply(gene.categories, function(df) split(df$gene_symbol, df$gs_name)) 


# Create a wrapper function around fgsea 

fgsea_wrapper <- function(limma_file, genesets, output_dir = "03. GSEA") {
  df <- openxlsx::read.xlsx(limma_file) 
  colnames(df) <- make.names(colnames(df)) 
  
  df <- df %>% dplyr::rename(SYMBOL = 1) %>% dplyr::filter(!is.na(SYMBOL), !is.na(logFC))
  
  ranked_genes <- df$logFC
  names(ranked_genes) <- df$SYMBOL
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (set_name in names(genesets)) {
    fgsea_res <- fgsea::fgseaMultilevel(pathways = genesets[[set_name]],
                                        stats = ranked_genes,
                                        minSize = 15,
                                        maxSize = 500) 
    
    fgsea_df <- fgsea_res %>% dplyr::arrange(padj)
    out_name <- paste0("fgsea_", tools::file_path_sans_ext(basename(limma_file)), "_", set_name, ".xlsx") 
    openxlsx::write.xlsx(fgsea_df, file = file.path(output_dir, out_name), asTable = TRUE)
  }
}

# List all xlsx limma result files

limma_files <- list.files(path = limma_dir, pattern = "\\.xlsx$", full.names = TRUE)
limma_files 

# Loop through files and run fgsea

for (file in limma_files) {
  fgsea_wrapper(limma_file = file, genesets = geneset_list, output_dir = gsea_dir) 
}


# ── Plot gene set enrichment ────────────────────────────────────────────────── 

# List all xlsx gsea result files

gsea_files <- list.files(path = gsea_dir, pattern = "\\.xlsx$", full.names = TRUE)
gsea_files 


# Find the top up- and down-regulated pathways in each condition 

all_pathways_list <- list()

for (file in gsea_files) {
  df <- openxlsx::read.xlsx(file)
  if (!"NES" %in% colnames(df)) next
  
  top_df <- df %>% arrange(desc(NES)) %>% slice_head(n = 10) %>% mutate(direction = "Top", type = "upregulated") 
  bottom_df <- df %>% arrange(NES) %>% slice_head(n = 10) %>% mutate(direction = "Bottom", type = "downregulated") 
  combined_df <- bind_rows(top_df, bottom_df)
  
  filename <- basename(file)
  database <- str_extract(filename, "(?<=_)[^_]+(?=\\.xlsx$)")
  condition <- str_replace(filename, paste0("^fgsea_(.*?)_", database, "\\.xlsx$"), "\\1")
  
  combined_df <- combined_df %>% mutate(database = database, condition = condition)
  all_pathways_list[[file]] <- combined_df
}

pathways_combined <- bind_rows(all_pathways_list)
pathways_combined %>% dplyr::count(condition) 


# Create a complex heatmap 

# Get all the conditions

conditions <- unique(pathways_combined$condition)

# Colors for the databases 

all_databases <- unique(pathways_combined$database)
db_colors <- setNames(RColorBrewer::brewer.pal(length(all_databases), "Set1"), all_databases)


# Loop over conditions to - (i) convert to pathways to long format, (ii) create annotation metadata, (iii) create complex heatmap 

for (cond in conditions) {
  for (reg_type in c("upregulated", "downregulated")) {
    
    df <- pathways_combined %>% 
      filter(type == reg_type, condition == cond) %>%
      group_by(pathway, database) %>%
      summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = database, values_from = NES, values_fill = NA) %>%
      column_to_rownames("pathway") %>%
      as.data.frame() 
    
    if (nrow(df) == 0) next
    
    ann_meta <- pathways_combined %>%
      filter(type == reg_type, condition == cond) %>%
      column_to_rownames(var = "pathway") %>%
      select(size, database) 
    ann_meta <- ann_meta[rownames(df), , drop = FALSE] 
    
    col_ha1 <- HeatmapAnnotation(Gene_count = anno_barplot(ann_meta$size, gp = gpar(fill = "yellow2"), bar_width = 0.8, border = FALSE), 
                                 annotation_name_side = "left") 
    col_ha2 <- HeatmapAnnotation(Database = ann_meta$database, 
                                 col = list(Database = db_colors),
                                 annotation_name_side = "left", 
                                 border = TRUE) 
    set.seed(123)
    ht <- Heatmap(t(df), 
                  # col = viridisLite::plasma(10), 
                  col = colorRampPalette(brewer.pal(9, "PiYG"))(10),
                  top_annotation = c(col_ha1, col_ha2),
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  na_col = "white", 
                  border = TRUE, 
                  column_title = paste(reg_type, "pathways -", cond),
                  column_title_gp = gpar(fontsize = 10, fontface = "bold")) 
    
    safe_cond <- gsub("[^A-Za-z0-9_\\-]", "_", cond)
    pdf_file <- file.path("03. GSEA", paste0("heatmap_", reg_type, "_", safe_cond, ".pdf"))
    pdf(pdf_file, width = 14, height = 4.5)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    dev.off() 
  }
}


# ── mRNAsi stemness analysis ────────────────────────────────────────────────── 

stem_genes <- openxlsx::read.xlsx("mRNAsi stemness geneset.xlsx")
stem_genes <- stem_genes %>% dplyr::select(TCGA_gene_id) %>% dplyr::rename(stemness_genes = "TCGA_gene_id")
stem_genes <- stem_genes %>% pivot_longer(c(stemness_genes), names_to = "source", values_to = "target")
stem_genes <- stem_genes %>% mutate(weight = 1) 
stem_genes <- stem_genes %>% filter(target!=" ")
stem_genes <- stem_genes %>% filter(target!="")
stem_genes <- unique(stem_genes)

# Run GSEA 

res_decoupler.stem <- decoupleR::decouple(mat = counts, 
                                          network = stem_genes, 
                                          .source ='source', 
                                          .target ='target',
                                          minsize = 0) 

res_ssgsea.stem <- decoupleR::run_gsva(mat = counts, 
                                       network = stem_genes, 
                                       .source ='source', 
                                       .target ='target', 
                                       minsize = 2L, 
                                       method = c("ssgsea")) # "gsva", "plage", "ssgsea", "zscore"

res_long.stem <- res_decoupler.stem %>% 
  dplyr::filter(statistic %in% "ulm") %>% 
  column_to_rownames('condition') %>% 
  dplyr::select(score, p_value) %>% 
  rename(stemness_score_depro = "score", p_value_stemness = "p_value")

meta.res <- cbind(res_long.stem, meta.res)


# ── ESTIMATE analysis ─────────────────────────────────────────────────────────  

dir.create("04. ESTIMATE", showWarnings = FALSE)
estimate_dir <- "04. ESTIMATE"

# Run ESTIMATE 

estim.scores <- counts |> 
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix = T)  

estim.scores <- estim.scores %>% column_to_rownames(var = "sample") 
colnames(estim.scores) <- paste0(colnames(estim.scores), "_depro")

meta.res <- cbind(estim.scores, meta.res)

# Data modification 

estim.df <- meta.res %>% dplyr::select(stromal_depro, immune_depro, estimate_depro, purity_depro, stemness_score_depro, claudin.groups, our.tumor.class, PAM50, limma.condition)
estim.df[, 1:5] <- scale(estim.df[, 1:5]) 

estim.long <- estim.df %>%
  pivot_longer(cols = c(stromal_depro, immune_depro, estimate_depro, purity_depro, stemness_score_depro), names_to = "estimate_parameters", values_to = "normalized_scores")


# Barplots 

colnames(estim.long)  

p1_est <- ggplot(estim.long, aes(x = estimate_parameters, y = normalized_scores, fill = estimate_parameters)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") + 
  facet_wrap(~ claudin.groups) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p1_est 
ggsave(plot = p1_est, filename = file.path(estimate_dir, "01. Estimate barplot-1.pdf"), width = 5, height = 3, units = c("in")) 


p2_est <- ggplot(estim.long, aes(x = claudin.groups, y = normalized_scores, fill = claudin.groups)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") + 
  facet_wrap(~ estimate_parameters) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p2_est 
ggsave(plot = p2_est, filename = file.path(estimate_dir, "01. Estimate barplot-2.pdf"), width = 5, height = 4, units = c("in")) 


p3_est <- ggplot(estim.long, aes(x = claudin.groups, y = normalized_scores, fill = estimate_parameters)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") + 
  facet_wrap(~ PAM50) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p3_est 
ggsave(plot = p3_est, filename = file.path(estimate_dir, "01. Estimate barplot-3.pdf"), width = 7, height = 4, units = c("in")) 


p4_est <- ggplot(estim.long, aes(x = claudin.groups, y = normalized_scores, fill = PAM50)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") + 
  facet_wrap(~ estimate_parameters) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p4_est 
ggsave(plot = p4_est, filename = file.path(estimate_dir, "01. Estimate barplot-4.pdf"), width = 7.5, height = 4, units = c("in")) 


# ── Subtype enrichment ──────────────────────────────────────────────────────── 

dir.create("05. Subtypes", showWarnings = FALSE) 
subtype_dir <- "05. Subtypes" 

subtype_genes <- openxlsx::read.xlsx("NMF breast subtype geneset.xlsx")
subtype_genes <- subtype_genes %>% pivot_longer(c(NMF_Her2_I, NMF_Basal_I, NMF_LumB_I, NMF_LumA_I), names_to = "source", values_to = "target")
subtype_genes <- subtype_genes %>% mutate(weight = 1) 
subtype_genes <- subtype_genes %>% filter(target!=" ")
subtype_genes <- subtype_genes %>% filter(target!="")
subtype_genes <- unique(subtype_genes)

# Run GSEA 

res_decoupler.subtype <- decoupleR::decouple(mat = counts, 
                                             network = subtype_genes, 
                                             .source ='source', 
                                             .target ='target',
                                             minsize = 0) 

res_ssgsea.subtype <- decoupleR::run_gsva(mat = counts, 
                                          network = subtype_genes, 
                                          .source ='source', 
                                          .target ='target', 
                                          minsize = 2L, 
                                          method = c("ssgsea")) # "gsva", "plage", "ssgsea", "zscore"

res_long.subtype <- res_decoupler.subtype %>% 
  dplyr::filter(statistic %in% "ulm") %>% 
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>% 
  column_to_rownames('source') %>% 
  t() %>% 
  as.data.frame() 

meta.res <- cbind(res_long.subtype, meta.res)

# Save all the metadata 

dir.create("Result sheets (CPTAC, Cell 2020)", showWarnings = FALSE) 
sheet_dir <- "Result sheets (CPTAC, Cell 2020)" 
save.meta <- meta.res %>% rownames_to_column(var = "Patients")
write.xlsx(save.meta, file = file.path(sheet_dir, "Result metadata.xlsx")) 

plot.heat.sub_df <- meta.res 
plot.heat.sub_df %>% dplyr::count(TUMOR_STAGE)

# Heatmap 

col_ha1 <- HeatmapAnnotation(stemness = anno_barplot(plot.heat.sub_df$stemness_score_depro, gp = gpar(fill = "yellow2"), bar_width = 0.8, border = T), 
                             annotation_name_side = "left") 

col_ha2 <- HeatmapAnnotation(claudin.groups = plot.heat.sub_df$claudin.groups, 
                             col = list(claudin.groups = c("claudin.low" = "#E41A1C", "Others" = "#377EB8")), 
                             annotation_name_side = "left", 
                             show_annotation_name = TRUE) 

col_ha3 <- HeatmapAnnotation(PAM50 = plot.heat.sub_df$PAM50, 
                             col = list(PAM50 = c("Basal" = "forestgreen", "Her2" = "orange", "LumA" = "red", "LumB" = "purple", "Normal-like" = "grey70")), 
                             annotation_name_side = "left", 
                             show_annotation_name = TRUE) 

col_ha4 <- HeatmapAnnotation(our.tumor.class = plot.heat.sub_df$our.tumor.class, 
                             col = list(our.tumor.class = c("High_Inoperable" = "darkblue", "Low_Operable" = "darkred", "Unknown" = "grey70")), 
                             annotation_name_side = "left", 
                             show_annotation_name = TRUE) 

col_ha5 <- HeatmapAnnotation(TUMOR_STAGE = plot.heat.sub_df$TUMOR_STAGE, 
                             col = list(TUMOR_STAGE = c("Stage IA" = "#F8766D", "Stage IIA" = "#7CAE00", "Stage IIB" = "#00BFC4", "Stage III" = "#C77CFF", "Stage IIIA" = "#FF2DF1", "Stage IIIB" = "#4ED7F1", "Stage IIIC" = "#FA812F", "NA" = "#BDBDBD")), 
                             annotation_name_side = "left", 
                             show_annotation_name = TRUE) 

row_sep_df <- data.frame(row_sep = ifelse(grepl("^NMF_", colnames(plot.heat.sub_df)[1:9]), "NMF", "ESTIMATE"), row.names = colnames(plot.heat.sub_df)[1:9])

set.seed(123) 
subtype.ht <- Heatmap(as.matrix(t(scale(plot.heat.sub_df[, 1:9]))), 
                      col = colorRamp2(breaks = c(-2, 0, 2), colors = c("darkblue", "#FFFFFF", "#FF0000")),
                      top_annotation = c(col_ha1, col_ha2, col_ha3, col_ha4, col_ha5),
                      cluster_rows = T,
                      cluster_columns = T, 
                      column_split = plot.heat.sub_df$claudin.groups,
                      row_split = row_sep_df$row_sep,
                      na_col = "white", 
                      border = T, 
                      show_column_names = T, 
                      column_names_gp = gpar(fontsize = 5), 
                      column_title_gp = gpar(fontsize = 10, fontface = "bold")) 
subtype.ht
subtype.ht <- as.ggplot(subtype.ht)
ggsave(filename = file.path(subtype_dir, "01. Subtype enrichment heatmap.pdf"), plot = subtype.ht, width = 11, height = 5)


# ── Correlation plots ───────────────────────────────────────────────────────── 

nmf_var <- c("NMF_Her2_I", "NMF_LumA_I", "NMF_LumB_I", "NMF_Basal_I", "immune_depro", "stromal_depro", "estimate_depro", "purity_depro", "stemness_score_depro")

plot.corr.sub_df <- meta.res %>% dplyr::select(NMF_Her2_I, NMF_LumA_I, NMF_LumB_I, NMF_Basal_I, immune_depro, stromal_depro, estimate_depro, purity_depro, stemness_score_depro, claudin.low.enrichment, claudin.groups, PAM50)
plot.corr.sub_df_scaled <- plot.corr.sub_df %>% mutate(across(where(is.numeric), scale))

# Correlation 1 

nmf.corr_plots.1 <- map(nmf_var, function(yvar) {
  ggplot(plot.corr.sub_df, aes(x = claudin.low.enrichment, y = .data[[yvar]])) + 
    geom_point(position = position_jitter(width = 0.5), color = "black", aes(fill = PAM50), pch = 21, size = 2) + 
    stat_smooth(method = "lm", formula = y ~ x) + 
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) + 
    facet_wrap(~ claudin.groups, nrow = 1, scales = "free") + 
    theme_bw() 
})

nmf.corr.1_g <- wrap_plots(nmf.corr_plots.1, ncol = 1)
nmf.corr.1_g 
ggsave(filename = file.path(subtype_dir, "02. Claudin correlation-1.pdf"), plot = nmf.corr.1_g, width = 6, height = 22)

# Correlation 2 

nmf.corr_plots.2 <- map(nmf_var, function(yvar) {
  ggplot(plot.corr.sub_df, aes(x = claudin.low.enrichment, y = .data[[yvar]])) + 
    geom_point(position = position_jitter(width = 0.5), color = "black", aes(fill = claudin.groups), pch = 21, size = 2) + 
    stat_smooth(method = "lm", formula = y ~ x) + 
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) + 
    theme_bw() 
})

nmf.corr.2_g <- wrap_plots(nmf.corr_plots.2, ncol = 3)
nmf.corr.2_g 
ggsave(filename = file.path(subtype_dir, "02. Claudin correlation-2.pdf"), plot = nmf.corr.2_g, width = 12, height = 7.25)

# Correlation 3 

nmf.corr_plots.3 <- map(nmf_var, function(yvar) {
  ggplot(plot.corr.sub_df, aes(x = immune_depro, y = .data[[yvar]])) + 
    geom_point(position = position_jitter(width = 0.5), color = "black", aes(fill = claudin.groups), pch = 21, size = 2) + 
    stat_smooth(method = "lm", formula = y ~ x) + 
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) + 
    facet_wrap(~ claudin.groups, nrow = 1, scales = "free") +
    theme_bw() 
})

nmf.corr.3_g <- wrap_plots(nmf.corr_plots.3, ncol = 1)
nmf.corr.3_g 
ggsave(filename = file.path(subtype_dir, "02. Claudin correlation-3.pdf"), plot = nmf.corr.3_g, width = 6, height = 22)

# Correlation 4 

nmf.corr_plots.4 <- map(nmf_var, function(yvar) {
  ggplot(plot.corr.sub_df, aes(x = NMF_Basal_I, y = .data[[yvar]])) + 
    geom_point(position = position_jitter(width = 0.5), color = "black", aes(fill = claudin.groups), pch = 21, size = 2) + 
    stat_smooth(method = "lm", formula = y ~ x) + 
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) + 
    facet_wrap(~ claudin.groups, nrow = 1, scales = "free") +
    theme_bw() 
})

nmf.corr.4_g <- wrap_plots(nmf.corr_plots.4, ncol = 1)
nmf.corr.4_g 
ggsave(filename = file.path(subtype_dir, "02. Claudin correlation-4.pdf"), plot = nmf.corr.4_g, width = 6, height = 22)


# ── Proportion plots ──────────────────────────────────────────────────────────  

plot.prop.sub_df <- meta.res 
plot.prop.sub_df <- plot.prop.sub_df %>% 
  dplyr::select(immune_depro, stromal_depro, estimate_depro, purity_depro, stemness_score_depro, claudin.groups) %>% 
  pivot_longer(cols = c(immune_depro, stromal_depro, estimate_depro, purity_depro, stemness_score_depro), names_to = "estimate_parameters", values_to = "normalized_scores")

plot.prop.sub_df_scaled <- plot.prop.sub_df %>% mutate(across(where(is.numeric), scale))

# Violin plot statistical comparison 

vio.prop.sub.p1 <- ggplot(plot.prop.sub_df, aes(x = claudin.groups, y = normalized_scores, fill = claudin.groups)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black", fatten = 1) +  # Median line
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", label.y.npc = "top") + 
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~ estimate_parameters, nrow = 1, scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
vio.prop.sub.p1
ggsave(filename = file.path(subtype_dir, "03. Violin plot-1.pdf"), plot = vio.prop.sub.p1, width = 10, height = 3)

