setwd("C:\\03. BRCA ZHH\\02. Data analysis\\01. RNA-seq (CPTAC, Cell 2020)")

# Author : Depro Das, Department of Biochemistry and Molecular Biology, University of Dhaka, Bangladesh 

# ── Libraries ───────────────────────────────────────────────────────────────── 

library(tidyverse) 
library(tibble)
library(tidyr) 
library(openxlsx) 
library(readxl) 
library(ggplot2)
library(ggrepel) 
library(decoupleR)
library(ComplexHeatmap) 
library(circlize) 
library(RColorBrewer) 
library(ggplotify) 
library(purrr) 
library(patchwork) 
library(msigdbr) 
library(ggvenn) 

# ── Data load and prepare inputs ────────────────────────────────────────────── 

# Count data  

countdata <- read.table(file = "data_mrna_seq_fpkm.txt", sep = "\t", header = TRUE)
countdata[ , -1] <- lapply(countdata[ , -1], function(x) as.numeric(as.character(x))) 
countdata$Hugo_Symbol <- make.names(countdata$Hugo_Symbol, unique = TRUE) 
countdata <- countdata %>% remove_rownames() %>% column_to_rownames(var = "Hugo_Symbol")

# Remove null values 

countdata <- countdata %>%
  dplyr::mutate_if(~ any(is.na(.x)), ~ dplyr::if_else(is.na(.x), 0, .x)) %>% 
  as.matrix() 

# Metadata 

metadata <- openxlsx::read.xlsx("C:\\03. BRCA ZHH\\02. Data analysis\\01. RNA-seq (CPTAC, Cell 2020)\\Result sheets (CPTAC, Cell 2020)\\Result metadata.xlsx")
metadata <- metadata %>% column_to_rownames(var = "Patients")

all(colnames(countdata) %in% rownames(metadata))
all(rownames(metadata) %in% colnames(countdata))

# Limma data load and manipulation 

limma.dirc <- "C:\\03. BRCA ZHH\\02. Data analysis\\01. RNA-seq (CPTAC, Cell 2020)\\02. Limma"
limma.data <- read.xlsx(file.path(limma.dirc, "Result_Limma_Others_vs_claudin_low.xlsx"))

colnames(limma.data)[1] <- "SYMBOL" 

limma.data$nlog10 <- -log10(limma.data$P.Value) 
limma.data$expression <- ifelse(limma.data$P.Value < 0.8912509 & abs(limma.data$logFC) >= 1,
                                ifelse(limma.data$logFC > 1, 'Up-regulated', 'Down-regulated'), 'Stable') 


# ── Immune data and gene identification ─────────────────────────────────────── 

dir.create("06. Immune genes", showWarnings = FALSE)
immune_dir <- "06. Immune genes"

# Up-regulated genes 

up_genes <- limma.data %>% filter(expression == "Up-regulated")
up_genes_list <- unique(up_genes$SYMBOL) 

# Immune genes-1 (InnateDB) 

main.dirc <- "C:\\03. BRCA ZHH\\02. Data analysis\\01. RNA-seq (CPTAC, Cell 2020)"
immune_genes <- read_excel(file.path(main.dirc, "InnateDB.xls")) 
immune_genes_list <- unique(immune_genes$Gene_Symbol) 

# Immune genes-2 (MSigDB C7 immunologic signatures) 

msigdb_c7_immune <- msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB")
msigdb_genes_list <- unique(msigdb_c7_immune$gene_symbol)


# Intersect and find genes 

venn_list <- list(InnateDB = immune_genes_list, MSigDB = msigdb_genes_list, Upregulated = up_genes_list)

ven.1 <- ggvenn(venn_list, 
                # fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                stroke_size = 0.5, 
                set_name_size = 4) 
ven.1 
ggsave(filename = file.path(immune_dir, "01. Subset venn.pdf"), plot = ven.1, width = 5, height = 5, units = "in")


# Upset plot 

ups_mtx <- make_comb_mat(venn_list) 

row_ha <- rowAnnotation("Set Size" = anno_barplot(set_size(ups_mtx), gp = gpar(fill = c("#4E71FF", "#A4DD00", "#FF2DD1")), width = unit(3, "cm")))
col_ha <- HeatmapAnnotation("Intersection Size" = anno_barplot(comb_size(ups_mtx), gp = gpar(fill = "#EFC000FF"), height = unit(3, "cm")))

sub.ups.1 <- UpSet(ups_mtx, 
                   pt_size = unit(4, "mm"), 
                   top_annotation = col_ha, 
                   right_annotation = row_ha,
                   comb_order = order(comb_size(ups_mtx))) 
sub.ups.1
sub.ups.1 <- as.ggplot(sub.ups.1)
ggsave(filename = file.path(immune_dir, "02. Subset upset.pdf"), plot = sub.ups.1, width = 4.75, height = 3, units = "in")


# ── DocupleR default ────────────────────────────────────────────────────────── 

dir.create("07. DecoupleR", showWarnings = FALSE)
decup_dir <- "07. DecoupleR"

# Prepare collectri regulatory network 

net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE) 

# Run ulm 

sample_acts <- decoupleR::run_ulm(mat = countdata, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 5) 

meta.filt <- metadata %>% select(claudin.groups, claudin.low.enrichment, PAM50, limma.condition) %>% rownames_to_column(var = "Patients")
merged_df <- sample_acts %>% left_join(meta.filt, by = c("condition" = "Patients"))

# Top tfs just to narrow down the tf matrix 

n_tfs <- 50
top_tfs <- merged_df %>%
  group_by(source) %>%
  summarise(std = sd(score, na.rm = TRUE)) %>%
  arrange(desc(std)) %>%
  slice_head(n = n_tfs) %>%
  pull(source)

subset_df <- merged_df %>% filter(source %in% top_tfs)

median_scores <- subset_df %>%
  group_by(source, claudin.groups) %>%
  summarise(median_score = mean(score, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = claudin.groups, values_from = median_score)

heatmap_matrix <- median_scores %>%
  column_to_rownames(var = "source") %>% 
  scale() %>% 
  as.matrix() 

set.seed(123) 
ht.decup.1 <- Heatmap(t(heatmap_matrix), 
                      col = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "red")),
                      rect_gp = gpar(col = "white", lwd = 1), 
                      border = TRUE)  
ht.decup.1 
ht.decup.1 <- as.ggplot(ht.decup.1)  
ggsave(plot = ht.decup.1, filename = file.path(decup_dir, "01. TFs heatmap.pdf"), width = 14, height = 2, units = c("in")) 


# ── Docupler limma ────────────────────────────────────────────────────────────  

# Prepare limma result and extract t-values per gene

deg <- limma.data %>%
  dplyr::select(SYMBOL, logFC, t, P.Value) %>% 
  dplyr::filter(!is.na(t)) %>% 
  tibble::column_to_rownames(var = "SYMBOL") %>%
  as.matrix() 

# Run ulm on limma 

contrast_limm <- decoupleR::run_ulm(mat = deg[, 't', drop = FALSE], 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5) 

# Create target gene groups (CD36)

net_fil <- net %>% filter(target %in% "CD36") 
contrast_limm <- contrast_limm %>% mutate(CD36_TFs = ifelse(contrast_limm$source %in% net_fil$source, "CD36_TF", "None"))

# Rank the tfs 

f_contrast_limm <- contrast_limm %>% dplyr::mutate(rnk = NA)
msk <- f_contrast_limm$score > 0
f_contrast_limm[msk, 'rnk'] <- rank(-f_contrast_limm[msk, 'score'])
f_contrast_limm[!msk, 'rnk'] <- rank(-abs(f_contrast_limm[!msk, 'score']))

# Select the top tfs based on ranks 

n_tfs_limm <- 50

tfs_limm <- f_contrast_limm %>% dplyr::arrange(rnk) %>% head(n_tfs_limm) %>% dplyr::pull(source)
f_contrast_limm <- f_contrast_limm %>% filter(source %in% tfs_limm)

# Plot the tfs 

colors_limm <- brewer.pal(9, "PiYG") 
f_contrast_limm$CD36_colors <- ifelse(f_contrast_limm$CD36_TFs == "CD36_TF", "#FF0B55", "#00CAFF")

p_bar_limm <- ggplot(data = f_contrast_limm, mapping = aes(x = stats::reorder(source, score), y = score)) + 
  geom_bar(mapping = aes(fill = score), color = "black", stat = "identity") +
  geom_tile(aes(y = min(score) - 0.05 * diff(range(score))), fill = f_contrast_limm$CD36_colors, height = 0.02 * diff(range(f_contrast_limm$score))) +
  scale_fill_gradientn(colors = rev(colors_limm), limits = c(min(f_contrast_limm$score), max(f_contrast_limm$score))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  xlab("TFs") 
p_bar_limm
ggsave(plot = p_bar_limm, filename = file.path(decup_dir, "02. TFs ranked bar.pdf"), width = 12, height = 2.5, units = c("in")) 


# ── Individual tfs and their targets ────────────────────────────────────────── 

tfs_all <- c("HIF1A", "CEBPB", "SMAD3", "SP1")  

# Map basic activating and deactivating genes for specific tfs 

tf_basic_volcano <- function(tf, net, deg) {
  df.tf <- net %>%
    dplyr::filter(source == tf) %>%
    dplyr::arrange(target) %>%
    dplyr::mutate(ID = target, color = "3") %>%
    tibble::column_to_rownames('target')
  
  inter.tf <- sort(dplyr::intersect(rownames(deg), rownames(df.tf))) 
  
  df.tf <- df.tf[inter.tf, ]
  df.tf[, c('logfc', 't_value', 'p_value')] <- deg[inter.tf, ]
  
  df.tf <- df.tf %>%
    dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>%
    dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p.tf <- ggplot(data = df.tf, mapping = aes(x = logfc, y = -log10(p_value), color = color, size = abs(mor))) + 
    geom_point(size = 1.5, alpha = 1, pch = 21, stroke = 0.0001) +
    scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
    ggrepel::geom_label_repel(mapping = aes(label = ID, size = 0.5)) + 
    theme_bw() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(tf) 
  
  return(p.tf)
} 

maping_tfs_basic <- map(tfs_all, ~ tf_basic_volcano(.x, net, deg))
names(maping_tfs_basic) <- tfs_all
p_tfs_basic <- wrap_plots(maping_tfs_basic, nrow = 1)
p_tfs_basic 
ggsave(plot = p_tfs_basic, filename = file.path(decup_dir, "03. TFs basic volcano.pdf"), width = 10, height = 4, units = c("in")) 


# Map immune-related activating and deactivating genes for specific tfs 

imm_selected <- read.xlsx("Selected immune related genes.xlsx")
immune_genes <- as.character(imm_selected[[1]])

tf_immune_volcano <- function(tf, net, deg, immune_genes) {
  df.tf <- net %>%
    dplyr::filter(source == tf) %>%
    dplyr::arrange(target) %>%
    dplyr::mutate(ID = target, color = "3") %>%
    tibble::column_to_rownames('target')
  
  inter.tf <- sort(dplyr::intersect(rownames(deg), rownames(df.tf))) 
  
  df.tf <- df.tf[inter.tf, ]
  df.tf[, c('logfc', 't_value', 'p_value')] <- deg[inter.tf, ]
  
  df.tf <- df.tf %>%
    dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>%
    dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>%
    dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  df.tf$label <- ifelse(df.tf$ID %in% immune_genes, df.tf$ID, NA)
  
  p.tf <- ggplot(data = df.tf, mapping = aes(x = logfc, y = -log10(p_value), color = color, size = abs(mor))) + 
    geom_point(size = 1.5, alpha = 1, pch = 21, stroke = 0.0001) +
    scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
    ggrepel::geom_label_repel(mapping = aes(label = label), size = 2, na.rm = TRUE, max.overlaps = 100) + 
    theme_bw() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    ggtitle(tf) 
  
  return(p.tf)
} 

maping_tfs_immune <- purrr::map(tfs_all, ~ tf_immune_volcano(.x, net, deg, immune_genes)) 
names(maping_tfs_immune) <- tfs_all
p_tfs_immune <- wrap_plots(maping_tfs_immune, nrow = 1)
p_tfs_immune 
ggsave(plot = p_tfs_immune, filename = file.path(decup_dir, "04. TFs immune volcano.pdf"), width = 10, height = 4, units = c("in")) 

