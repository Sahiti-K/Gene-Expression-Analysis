# Gene Expression Analysis: R Project 

**Project: DumbDeseq — RNA Expression Software**  
*Processed RNA-seq dataset from quantitated gene expression data.*  

## Overview 
This repository provides a reproducible workflow for **differential gene expression analysis** using R. It processes expression data, classifies genes into **upregulated**, **downregulated**, or **not significant**, generates a **volcano plot**, and summarizes the **top regulated genes**.  

## Data Source
Dataset used here is available from:  
[Stephen Turner’s gist: `results.txt`](https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt)  

This dataset contains gene expression results with `log2FoldChange` and `pvalue` columns for differential expression analysis.  

## Workflow

| Steps | Tasks | Details |
|------|------|---------|
| 1 | **Data Import** | Inspect structure and preview rows |
| 2 | **Preprocessing** | Compute `-log10(p-value)` for improved visualization of significance; classify genes as **Up**, **Down**, or **Not Significant** based on thresholds |
| 3 | **Export Results** | Save the full dataset with classification; save filtered table of significantly regulated genes |
| 4 | **Visualization** | Create a volcano plot using `ggplot2` |
| 5 | **Summary Table** | Extract the top 5 Upregulated and the top 5 Downregulated genes |
| 6 | **Biological Interpretation** | Annotate top genes using functional summaries from GeneCards |  

## Outputs  
- `classified_gene_expression.csv` — full dataset with classification  
- `significant_genes.csv` — filtered significant genes  
- `volcano_plot.png` — visual summary of differential expression  
- `summary_table.csv` — top regulated genes with fold change and p-values  
-  HTML report generated via R Markdown  

## How to Run  
Make sure you have R and the required packages installed:

```r
install.packages(c("ggplot2", "dplyr", "readr", "knitr"))
```
Then knit the R Markdown file to generate HTML and PDF reports:
```r
rmarkdown::render("gene_expression_analysis.Rmd")
```

## Repository Structure

| Folder/File | Description |
|-------------|-------------|
| `scripts/`  | R Markdown and helper scripts |
| `results/`  | Output tables and plots |
| `report/`   | HTML |
| `README.md` | Project overview |
| `LICENSE`   | Usage terms |


## References
Stelzer, G., et al. (2016). The GeneCards Suite: From Gene Data Mining to Disease Genome Sequence Analyses. Current Protocols in Bioinformatics, 54(1).

## License
This project is licensed under the MIT License.
