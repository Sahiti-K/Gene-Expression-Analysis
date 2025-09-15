# ------------------------------------------------------------
# Differential Gene Expression Analysis
# ------------------------------------------------------------

# Read the file, splitting on any whitespace
data <- read.table(
  "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt",
  header = TRUE,
  sep = "",            # empty string = split on any white space
  stringsAsFactors = FALSE
)

# Verify the data structure
str(data)
# Check the first few rows
head(data)

# Calculate -log10(p-value) for volcano plot y-axis
data$logP <- -log10(data$pvalue)

# Classify genes based on fold change and p-value
data$status <- "Not Significant"  # default
data$status[data$log2FoldChange > 1 & data$pvalue < 0.01] <- "Up"
data$status[data$log2FoldChange < -1 & data$pvalue < 0.01] <- "Down"

# 3. Save updated table
write.csv(data, "classified_gene_expression.csv", row.names = FALSE)

# Save only significant genes
sig_data <- subset(data, status != "Not Sig")
write.csv(sig_data, "significant_genes.csv", row.names = FALSE)

# Optional: Preview first few rows
head(data)

library(ggplot2)

# 3. Volcano plot
volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = logP, color = status)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(p-value)")

ggsave("volcano_plot.png", volcano_plot, width = 7, height = 5)

# Filter for significant genes
up_genes <- subset(data, log2FoldChange > 1 & pvalue < 0.01)
down_genes <- subset(data, log2FoldChange < -1 & pvalue < 0.01)

# Top 5 up-regulated genes (largest positive log2FC)
top5_up <- head(up_genes[order(-up_genes$log2FoldChange), ], 5)

# Top 5 down-regulated genes (most negative log2FC)
top5_down <- head(down_genes[order(down_genes$log2FoldChange), ], 5)

# Create a summary table
summary_table <- data.frame(
  Upregulated_Gene = top5_up$Gene,
  Up_Log2FC = round(top5_up$log2FoldChange, 3),
  Up_Pvalue = signif(top5_up$pvalue, 3),
  Downregulated_Gene = top5_down$Gene,
  Down_Log2FC = round(top5_down$log2FoldChange, 3),
  Down_Pvalue = signif(top5_down$pvalue, 3)
)

# View the table
print(summary_table)

# Save to file
write.csv(summary_table, "summary_table.csv", row.names = FALSE)
