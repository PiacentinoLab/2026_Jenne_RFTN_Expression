# Import packages for DE-normalization, line plots and heatmaps
library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

############################## Data Prep ##############################

#Read in counts matrix and ensure sample columns appropriate
data <- read.delim("Leung2019_featureCounts_WT_hiNCC_final.txt", header=TRUE,
                   comment.char="#")

# Define sample info including Cell Type, Day, and Replicate
sample_info <- data.frame(
  Sample = c("ESC_Day0_rep1", "ESC_Day0_rep2", "PreNCC_Day3_rep1", 
             "PreNCC_Day3_rep2", "NCC_Day5_rep1", "NCC_Day5_rep2"),
  Cell_Type = c("ESC", "ESC", "PreNCC", "PreNCC", "NCC", "NCC"),
  Day = c(0, 0, 3, 3, 5, 5),
  Replicate = c(1, 2, 1, 2, 1, 2)
)

sample_cols <- sample_info$Sample[sample_info$Sample %in% colnames(data)]
if (length(sample_cols) == 0) stop("No sample columns found. Check 
                                   sample_info$Sample vs colnames(data).")

# remove rows with missing Geneid or with all-sample NAs, save removed rows
missing_id <- is.na(data$Geneid) | trimws(as.character(data$Geneid)) == ""
all_na_samples <- apply(data[, sample_cols, drop = FALSE], 1, 
                        function(x) all(is.na(x)))
remove_mask <- missing_id | all_na_samples
data_clean <- data[!remove_mask, ] # Keep only valid rows

############################## Normalization ##############################

# Prepare count matrix for DESeq2
# Extract count columns and set gene IDs as rownames
count_matrix <- data_clean[, sample_cols]
rownames(count_matrix) <- data_clean$Geneid

# Convert to integer matrix (DESeq2 requires integer counts)
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer"

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ Cell_Type  # or ~ Day, depending on your analysis
)

# Perform VST normalization
vst_data <- vst(dds, blind = TRUE)  # blind=TRUE for exploratory analysis

# Extract VST-normalized values as data frame
# Extract VST-normalized values as data frame
vst_matrix <- assay(vst_data)
vst_df <- as.data.frame(vst_matrix)

# Add Geneid and GeneName directly from data_clean
vst_df$Geneid <- data_clean$Geneid
vst_df$GeneName <- data_clean$GeneName

# Reorder columns
vst_df <- vst_df[, c("Geneid", "GeneName", sample_cols)]

# Check result
cat("VST normalization complete!\n")
cat("Dimensions:", nrow(vst_df), "genes x", length(sample_cols), "samples\n")
cat("Columns:", colnames(vst_df), "\n")
head(vst_df)


############################## Line Plots ##############################
# Define genes of interest
genes_of_interest <- c("NANOG", "PAX3", "TFAP2B", "RFTN1", "RFTN2")

# Filter VST data for genes of interest
vst_genes <- vst_df %>%
  filter(GeneName %in% genes_of_interest)

# Check which genes were found
cat("Genes found:", unique(vst_genes$GeneName), "\n")
cat("Genes missing:", setdiff(genes_of_interest, vst_genes$GeneName), "\n")

# Reshape data for plotting (long format)
vst_long <- vst_genes %>%
  select(GeneName, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), 
               names_to = "Sample", 
               values_to = "VST_Expression") %>%
  left_join(sample_info, by = "Sample")

# Create line plot
ggplot(vst_long, aes(x = Day, y = VST_Expression, color = GeneName, group = GeneName, fill = GeneName)) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, color = NA) +  # error bands
  stat_summary(fun = mean, geom = "line", size = 1.2) +  # line connects means
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +  # matching fill colors
  scale_x_continuous(breaks = c(0, 3, 5), labels = c("ESC", "Pre-NCC", "NCC")) +
  labs(title = "Gene Expression Over Differentiation",
       x = "Day",
       y = "VST-Normalized Expression",
       color = "Gene",
       fill = "Gene") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))


############################## Heatmap ##############################
library(viridis)

# Define genes of interest
genes_of_interest <- c("NANOG", 
                       "PAX3", "PAX7", "SOX10", "SNAI2", 
                       "TFAP2B", "RFTN1", "RFTN2")

# Filter VST data for genes of interest
vst_genes <- vst_df %>%
  filter(GeneName %in% genes_of_interest)


# Heatmap with individual replicates
heatmap_all <- as.data.frame(vst_genes)
rownames(heatmap_all) <- heatmap_all$GeneName
heatmap_all <- as.matrix(heatmap_all[, sample_cols])

# Create annotation for samples
col_annotation <- data.frame(
  Cell_Type = sample_info$Cell_Type,
  Day = as.factor(sample_info$Day),
  row.names = sample_info$Sample
)

# Create heatmap
pheatmap(heatmap_all,
         #scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = col_annotation,
         color = viridis(100),
         #color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         #color = colorRampPalette(c("blue", "white", "red"))(10),
         main = "Gene Expression Over Differentiation",
         fontsize = 12)

