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
data$status <- "Neutral"  # default - not significant
data$status[data$log2FoldChange > 1 & data$pvalue < 0.01] <- "Up"
data$status[data$log2FoldChange < -1 & data$pvalue < 0.01] <- "Down"

# Save only significant genes
sig_data <- subset(data, status != "Neutral")
write.csv(sig_data, "significant_genes.csv", row.names = FALSE)

# Optional: Preview first few rows
head(data)

# Colors
cols <- c("Down"="blue", "Neutral"="grey", "Up"="red")

# Open a PDF device
pdf("volcano_plot.pdf", width = 7, height = 5)  # size in inches

# Volcano plot
plot(
  data$log2FoldChange, data$logP,
  col = adjustcolor(cols[data$status]),
  pch = 16, cex = 0.8,
  xlab = "Log2 Fold Change", ylab = "-log10(p-value)",
  main = "Volcano Plot"
)
abline(v=c(-1,1), lty=3)
abline(h=-log10(0.01), lty=3)
legend("topright", legend=names(cols), col=cols, pch=16,
       cex=0.7, pt.cex=0.7, bty="o")

# Close the device
dev.off()

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




