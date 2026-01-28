library(dplyr)
setwd("/Users/sarah/Duke Bio_Ea Dropbox/Sarah Marion/Lethals Project/Empircal Data Analysis R/Relatedness")
## In DCC, merged.vcf.gz of all sequenced lethal lines dv SNPs (plus Illumina sequenced L14, L230, L310)
## run vcftools --gzvcf merged.vcf.gz --relatedness2 --out relatedness_data

relatedness_data <- read.table("relatedness_data.txt", header = TRUE)
relatedness_data$Relatedness_fix <- ifelse(relatedness_data$INDV1 == relatedness_data$INDV2, NA, relatedness_data$RELATEDNESS_PHI)
relatedness_clean <- select(relatedness_data, c("INDV1","INDV2","Relatedness_fix"))
LL_remove <- c("L139","L284","L83","L230","L114","L207","L247","L278")
relatedness_extra_rmvd <- relatedness_clean %>%  filter(!(INDV1 %in% LL_remove | INDV2 %in% LL_remove))

# Reshape data into a matrix
relatedness_matrix <- reshape(relatedness_extra_rmvd, idvar = "INDV1", timevar = "INDV2", direction = "wide")

# Remove the 'INDV1' column (we only want the relatedness values)
relatedness_matrix <- relatedness_matrix[,-1]

# Replace the column names to remove the prefix "relatedness."
colnames(relatedness_matrix) <- gsub("Relatedness_fix.","", colnames(relatedness_matrix))
rownames(relatedness_matrix) <- colnames(relatedness_matrix)
relatedness_matrix <- as.matrix(relatedness_matrix)

# Manually define the order of the samples
sample_order <- c("L19","L256","L59","L96","L14","L121","L257","L49","L310","L34","L170","L30","L122",
                  "L173","L302","L72","L95","L10","L12","L212","L271","L274","L288","L297","L301","L303","L158","L75",
                  "L282","L98","L127","L208","L224")

# Manually name rows by the gene they mapped to and columns by year sampmle was collected
row_names_gene <- c("drosha","drosha","drosha","drosha","drosha","CG33155","CG33155","CG33155","CG33155","shot","shot","shot","Nipped-A",
                    "Nipped-A","Nipped-A","Nipped-A","Nipped-A","Nipped-A","l(2)gl","l(2)gl","l(2)gl","l(2)gl","l(2)gl","l(2)gl","l(2)gl","l(2)gl","CG13185","CG13185",
                    "Ca-alpha1D","Ca-alpha1D","Fgop2","1psc","kr-h1")

col_names_year <- c("2018","2018","2018","2018","2018","2018","2018","2018","2021","2018","2021_SM","2018","2018",
                    "2021","2021","2020","2018","2018","2018","2018","2018","2021","2021","2021","2021","2021","2021","2020",
                    "2021","2018","2018","2018","2018")

row_names_geneyear <- c("drosha, 2018","drosha, 2018","drosha, 2018","drosha, 2018","drosha, 2018","CG33155, 2018","CG33155, 2018","CG33155, 2018","CG33155, 2021","shot, 2018","shot, 2021","shot, 2018","Nipped-A, 2018",
                        "Nipped-A, 2021","Nipped-A, 2021","Nipped-A, 2020","Nipped-A, 2018","Nipped-A, 2018","l(2)gl, 2018","l(2)gl, 2018","l(2)gl, 2018","l(2)gl, 2021","l(2)gl, 2021","l(2)gl, 2021","l(2)gl, 2021","l(2)gl, 2021","CG13185, 2021","CG13185, 2020",
                        "Ca-alpha1D, 2021","Ca-alpha1D, 2018","Fgop2, 2018","1psc, 2018","kr-h1, 2018")
order_labels_new <- data.frame(Sample= sample_order,
                               gene_year = row_names_geneyear)
order_labels <- data.frame(Sample = sample_order, 
                           Gene = row_names_gene, 
                           Year = col_names_year)
relatedness_matrix <- relatedness_matrix[sample_order, sample_order]
print(rownames(relatedness_matrix))

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(viridis)
#col_fun = colorRamp2(c(0, 0.25, 0.5), viridis(3))
#col_fun <- colorRamp2(c(0.2, 0.35, 0.5), c("lightgray", "yellow", "red"))
col_fun <- colorRamp2(c(0.2, 0.4, 0.5), rev(hcl.colors(3, "YlGnBu")))
# Set the upper triangle to NA before passing it to Heatmap
#relatedness_matrix[lower.tri(relatedness_matrix)] <- NA

# Plot the heatmap with the modified matrix
Heatmap(relatedness_matrix, col = col_fun,
        row_order = order_labels_new$Sample, 
        row_labels = order_labels_new$gene_year[match(order_labels_new$Sample, sample_order)],  # Ensure alignment
        column_order = order_labels_new$Sample,
        #column_labels = order_labels_new$gene_year[match(order_labels_new$Sample, sample_order)],
        rect_gp = gpar(col = "gray", lwd = 1),
        na_col = "white",  
        heatmap_legend_param = list(title = "Kinship Coefficient"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x, y, width, height, gp = gpar(col = "white", fill = fill))
        })
#Plot without column labels
Heatmap(relatedness_matrix, col = col_fun,
        row_order = order_labels_new$Sample, 
        row_labels = order_labels_new$gene_year[match(order_labels_new$Sample, sample_order)],  # Ensure alignment
        column_order = order_labels_new$Sample,
        rect_gp = gpar(col = "gray", lwd = 1),
        na_col = "white",  
        heatmap_legend_param = list(title = "Kinship Coefficient"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x, y, width, height, gp = gpar(col = "white", fill = fill))
        },
        show_column_names = FALSE) 







# Install and load the necessary library
install.packages("pheatmap")
library(pheatmap)

pheatmap(relatedness_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         main = "Relatedness Heatmap")
