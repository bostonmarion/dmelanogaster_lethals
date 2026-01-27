if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
install.packages("dplyr")
library(GenomicRanges)
library(dplyr)

setwd("/Users/sarah/Duke Bio_Ea Dropbox/Sarah Marion/Lethals Project/Empircal Data Analysis R")


##TEST DATA

rr_test <- read.csv("RR_Data/rr_test_data.csv")
rr_test_gr <- makeGRangesFromDataFrame(rr_test,
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE,
                                  seqinfo = NULL,
                                  seqnames.field = "arm",
                                  start.field="start",
                                  end.field = "end",
                                  starts.in.df.are.0based = FALSE)
rr_test_gr

df_test <- read.csv("DF_mapping_data/DF_test_data.csv")
df_test_gr <- makeGRangesFromDataFrame(df_test,
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE,
                                  seqinfo = NULL,
                                  seqnames.field = "arm",
                                  start.field="start",
                                  end.field ="end",
                                  starts.in.df.are.0based = TRUE)
df_test_gr

##############______________________________________________________________________________

##REAL DATA
#Read in Comeron RR data and combine into one RR data set differentiated by ARM
RR_2L <- read.table("RR_Data/Comeron_100kb_chr2L_ref6.txt", sep = '')
RR_2L <- RR_2L %>% mutate(arm = "L") %>% mutate(end=V1+99999) 
RR_2L <- dplyr::rename(RR_2L, start = V1, rr = V2)

RR_2R <- read.table("RR_Data/Comeron_100kb_chr2R_ref6.txt", sep = '')
RR_2R <- RR_2R %>% mutate(arm = "R") %>% mutate(end=V1+99999) %>% dplyr::rename(start=V1, rr=V2)

RR <- rbind(RR_2L,RR_2R)

rr_gr <- makeGRangesFromDataFrame(RR,
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE,
                                  seqinfo = NULL,
                                  seqnames.field = "arm",
                                  start.field="start",
                                  end.field = "end",
                                  starts.in.df.are.0based = FALSE)
rr_gr

#Read in DF Data
df <- read.csv("DF_mapping_data/DF_R_Oct_11_2024_FINAL.csv")

df_gr <- makeGRangesFromDataFrame(df,
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE,
                                  seqinfo = NULL,
                                  seqnames.field = "arm",
                                  start.field="start",
                                  end.field ="end",
                                  starts.in.df.are.0based = TRUE)
df_gr

# Find Overlaps

hits <- findOverlaps(df_gr, rr_gr)

# Step 4: Calculate percent overlap and weighted recombination rate
# Compute overlapping regions and their widths
overlapping_DF <- pintersect(df_gr[queryHits(hits)], rr_gr[subjectHits(hits)])
overlap_widths <- width(overlapping_DF)
total_widths <- width(df_gr[queryHits(hits)])

# Calculate percent overlap and weighted recombination rates
percentOverlap <- overlap_widths / total_widths
weighted_rr <- mcols(rr_gr[subjectHits(hits)])$rr * percentOverlap

# Step 5: Create a data frame to summarize results
summary_df <- data.frame(
  df_id = mcols(df_gr[queryHits(hits)])$df_id,  # Assuming df_id is in the original df
  rr_weighted = weighted_rr
)

# Calculate the average recombination rate by df_id
avg_rr <- summary_df %>% 
  group_by(df_id) %>% 
  summarize(rr_avg = sum(rr_weighted, na.rm = TRUE), .groups = 'drop')

# Step 6: Merge average recombination rate with DF data frame
DF_data <- merge(avg_rr, df, by = "df_id", all.y = TRUE)  # Use all.y to keep all df rows
DF_data <- DF_data %>% 
  mutate(map_freq = n_mapped / n_crossed) %>% 
  mutate(map_freq_weightbysize = map_freq / size)
#------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------
##Read in and calculate gene density by 100,000 bp section with GenomicRanges (Gene Density bed files from FlyBase to GRanges)

gd_2L <- as.data.frame(read.table("GD_Data/Gene span-2L-1..23513712.bed",header = FALSE, sep="",stringsAsFactors=FALSE))
gd_2L$V1[gd_2L$V1 == '2L'] <- 'L'
head(gd_2L)
gd_2R <- as.data.frame(read.table("GD_Data/Gene span-2R-1..25286936.bed",header = FALSE, sep="",stringsAsFactors=FALSE))
gd_2R$V1[gd_2R$V1 == '2R'] <- 'R'
gd <- rbind(gd_2L, gd_2R)

gd_gr <- makeGRangesFromDataFrame(gd,
                                  keep.extra.columns = TRUE,
                                  ignore.strand = TRUE,
                                  seqinfo = NULL,
                                  seqnames.field = "V1",
                                  start.field="V2",
                                  end.field ="V3",
                                  starts.in.df.are.0based = TRUE)
gd_gr



# Find Overlaps between Deficiencies and gene ranges
hits_gd <- findOverlaps(df_gr, gd_gr)

# For each gene that overlaps with a deficiency, calculate the intersection length (proportion in the region)
overlapping_genes_DF <- pintersect(df_gr[queryHits(hits_gd)], gd_gr[subjectHits(hits_gd)])
overlap_lengths <- width(overlapping_genes_DF)


# Split the overlapping gene regions by deficiency region
overlaps_by_def <- split(overlapping_genes_DF, queryHits(hits_gd))

# Apply the `reduce` function to merge overlapping gene regions within each deficiency
reduced_overlaps_by_def <- lapply(overlaps_by_def, reduce)

# Calculate the total length of reduced overlapping regions for each deficiency
gene_overlap_lengths_by_def <- lapply(reduced_overlaps_by_def, function(gr) sum(width(gr)))

# Convert the list to a numeric vector (sums of overlap lengths for each deficiency)
sum_gene_overlap_per_def <- unlist(gene_overlap_lengths_by_def)

# Get the width of each deficiency region (ensure it matches the order in the split list)
def_lengths <- width(df_gr)[unique(queryHits(hits_gd))]

# Calculate the proportion of base pairs covered by genes (gene density) for each deficiency region
gene_density_per_def <- sum_gene_overlap_per_def / def_lengths


## Merge this with the other data set to finally get Deficiencies with average recombination rate, and gene density
mcols(df_gr)$gene_density <- gene_density_per_def
df_gr

df_gd <- as.data.frame(df_gr) %>% select(df_id, gene_density)
df_gd <- df_gd[, c("df_id", "gene_density")]
DF_data <- merge(DF_data, df_gd, by="df_id") 
# Remove outlier ed50001 because outlier and because we don't really know breakpoints...
DF_data <- DF_data[DF_data$df_id != "DF(2L)ED50001", ]

####_------------------------------------------------------------------------------------
##The above code calculates average basepair in a gene but not necessarily CODING base pair. Downloaded CDS regions from flybase to get coding basepairs instead

#install.packages("BiocManager")
#BiocManager::install("rtracklayer")
library(rtracklayer)

##import gff cds files with rtracklayer so that they can be granges objects
cds_2L <- import("GD_Data/CDS-2L-1..23513712.gff3", format = "gff3")
cds_2R <- import("GD_Data/CDS-2R-1..25286936.gff3", format = "gff3")
seqlevels(cds_2L)
seqlevels(cds_2R)
seqlevels(cds_2L) <- "L"
seqlevels(cds_2R) <- "R"

# Combine the two GRanges objects
cds_gr <- c(cds_2L, cds_2R)

# Find Overlaps between Deficiencies and cds ranges
hits_cds <- findOverlaps(df_gr, cds_gr)

# For each gene that overlaps with a deficiency, calculate the intersection length (proportion in the region)
overlapping_cds_DF <- pintersect(df_gr[queryHits(hits_cds)], cds_gr[subjectHits(hits_cds)])
overlap_lengths_cds <- width(overlapping_cds_DF)


# Split the overlapping coding regions by deficiency region
overlapcds_by_def <- split(overlapping_cds_DF, queryHits(hits_cds))

# Apply the `reduce` function to merge overlapping coding regions within each deficiency
reduced_overlapcds_by_def <- lapply(overlapcds_by_def, reduce)

# Calculate the total length of reduced overlapping regions for each deficiency
cds_overlap_lengths_by_def <- lapply(reduced_overlapcds_by_def, function(gr) sum(width(gr)))

# Convert the list to a numeric vector (sums of overlap lengths for each deficiency)
sum_cds_overlap_per_def <- unlist(cds_overlap_lengths_by_def)

# Get the width of each deficiency region (ensure it matches the order in the split list)
def_lengths_cds <- width(df_gr)[unique(queryHits(hits_cds))]

# Calculate the proportion of base pairs covered by genes (gene density) for each deficiency region
cds_density_per_def <- sum_cds_overlap_per_def / def_lengths_cds


## Merge this with the other data set to finally get Deficiencies with average recombination rate, and gene density
mcols(df_gr)$coding_density <- cds_density_per_def
df_gr

df_cds <- as.data.frame(df_gr) 
df_cds <- df_cds[, c("df_id", "coding_density")]

DF_data <- merge(DF_data, df_cds, by="df_id") 


## Look at GC content relationships with lethals mpaping

install.packages("Biostrings")
library(Biostrings)

##Read in reference genome sequence for 2L and 2R
chr2L <- readDNAStringSet("sequence/2L_ref6.fas")
chr2L <- unlist(chr2L)
chr2R <- readDNAStringSet("sequence/2R_ref6.fas")
chr2R <- unlist(chr2R)
print(class(chr2R))
# divide df into df_2L and df_2R

df_2L <- df %>% filter(arm=="L")
df_2L[13, "end"] <- 23513711
print(df_2L)
df_2R <- df %>% filter(arm=="R")

# # Initialize vectors to store results
at_props_2L <- numeric()
cg_props_2L <- numeric()

# Loop through each specified range
for (i in 1:nrow(df_2L)) {
  start <- df_2L$start[i]
  end <- df_2L$end[i]
  
  # Extract the sequence region
  region_2L <- subseq(chr2L, start=start, end=end)
  
  # Count nucleotides
  counts_2L <- alphabetFrequency(region_2L, baseOnly=TRUE)
  at_count_2L <- counts_2L["A"] + counts_2L["T"]
  cg_count_2L <- counts_2L["C"] + counts_2L["G"]
  total_2L <- at_count_2L + cg_count_2L
  
  # Compute proportions
  at_props_2L <- c(at_props_2L, ifelse(total_2L > 0, at_count_2L / total_2L, NA))
  cg_props_2L <- c(cg_props_2L, ifelse(total_2L > 0, cg_count_2L / total_2L, NA))
}

# Add proportions to the dataframe
df_2L$AT_Proportion <- at_props_2L
df_2L$CG_Proportion <- cg_props_2L

## Now same thing for 2R

# # Initialize vectors to store results
at_props_2R <- numeric()
cg_props_2R <- numeric()

# Loop through each specified range
for (i in 1:nrow(df_2R)) {
  start <- df_2R$start[i]
  end <- df_2R$end[i]
  
  # Extract the sequence region
  region_2R <- subseq(chr2R, start=start, end=end)
  
  # Count nucleotides
  counts_2R <- alphabetFrequency(region_2R, baseOnly=TRUE)
  at_count_2R <- counts_2R["A"] + counts_2R["T"]
  cg_count_2R <- counts_2R["C"] + counts_2R["G"]
  total_2R <- at_count_2R + cg_count_2R
  
  # Compute proportions
  at_props_2R <- c(at_props_2R, ifelse(total_2R > 0, at_count_2R / total_2R, NA))
  cg_props_2R <- c(cg_props_2R, ifelse(total_2R > 0, cg_count_2R / total_2R, NA))
}

# Add proportions to the dataframe
df_2R$AT_Proportion <- at_props_2R
df_2R$CG_Proportion <- cg_props_2R

df_CG <- bind_rows(df_2R, df_2L)
DF_data <- merge(DF_data, df_CG, by = "df_id", all.y = TRUE)  # Use all.y to keep all df rows



##Plot it all colored by low or high recombination rate

DF_data <- DF_data  %>%
  mutate(rr_range = case_when(rr_avg < 1 ~ 'low',
                              rr_avg > 1 ~ 'high')) %>%
  mutate(map_freq_weightbysize_cdense = map_freq_weightbysize*coding_density)


# DF by ARM

DF_2R <- DF_data[(DF_data$arm.x == "R"),]
DF_2L <- DF_data[(DF_data$arm.x == "L"),]


## filter to high_rr

DF_data_rrhigh <- DF_data %>% filter(rr_range == "high")
DF_data_rrlow <- DF_data %>% filter(rr_range == "low")
DF_2L_rrhigh <- DF_2L %>% filter(rr_range == "high")
DF_2L_rrlow <- DF_2L %>% filter(rr_range == "low")
DF_2R_rrhigh <- DF_2R %>% filter(rr_range == "high")
DF_2R_rrlow <- DF_2R %>% filter(rr_range == "low")

DF_data_summary_cd <- DF_data %>%
  summarise(cd_avg = mean(coding_density, na.rm = TRUE), 
            cd_se = sd(coding_density, na.rm = TRUE),
            map_freq_mean = mean(map_freq_weightbysize, na.rm = TRUE),
            map_freq_max = max(map_freq_weightbysize, na.rm = TRUE),
            map_freq_min = min(map_freq_weightbysize, na.rm = TRUE),
            map_freq_top_quartile = quantile(map_freq_weightbysize, 0.75, na.rm = TRUE))

DF_data_highestmap <- DF_data %>% filter(map_freq_weightbysize > 4.182493e-08)


#ggplot recombination rate by lethals mapped/crosses complete/DF size
library(ggplot2)

###  positive correlation between rr and coding density 
rr_gd_dotplot <- ggplot(data=DF_data, aes(x=rr_avg, y=coding_density)) + 
  geom_point(aes(color="red")) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Mean recombination rate in deficiency (cM/Mb)") +
  ylab("gene density (Proportion of DF Coding)") +
  ggtitle("Gene Density by Recombination Rate per Deficiency") +
  theme(plot.title = element_text(hjust = 0.5)) 
#theme_light()

rr_gd_dotplot

cor.test(DF_data$rr_avg,DF_data$coding_density, method ="pearson")

### No correlation between lethals and rr
rr_lethals_dotplot <- ggplot(data=DF_data, aes(x=rr_avg, y=map_freq_weightbysize)) + 
  geom_point(aes(color="red")) +
  #geom_smooth(method = "lm", se = FALSE) +
  xlab("Mean recombination rate in deficiency (cM/Mb)") +
  ylab("Lethals mapped / Crosses completed / Deficiency size") +
  ggtitle("Lethal Mutation Distribution by Recombination Rate") +
  theme(plot.title = element_text(hjust = 0.5)) 
  theme_light()
  
rr_lethals_dotplot

cor.test(DF_data$rr_avg, DF_data$map_freq_weightbysize, method ="pearson")

###from df_data, only low rec (rr <1 cM/Mb)

rr_lethals_dotplot_lowrr <- ggplot(data=DF_data_rrlow, aes(x=rr_avg, y=map_freq_weightbysize)) + 
  geom_point(aes(color="red")) +
  #geom_smooth(method = "lm", se = FALSE) +
  xlab("2L Mean recombination rate in deficiency (In(2L)t removed) (cM/Mb)") +
  ylab("Lethals mapped / Crosses completed / Deficiency size") +
  ggtitle("Lethal Mutation Distribution by Recombination Rate") +
  theme(plot.title = element_text(hjust = 0.5)) 
theme_light()

rr_lethals_dotplot_lowrr

cor.test(DF_data_rrlow$rr_avg, DF_data_rrlow$map_freq_weightbysize, method ="pearson")
### Positive/almost significant correlation between lethals and coding density
rr_gd <- ggplot(data=DF_data, aes(x=coding_density, y=map_freq_weightbysize)) + 
  geom_point(aes(color="red")) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("gene density (Proportion basepairs that are coding)") +
  ylab("Lethals mapped / Crosses completed / Deficiency size") +
  ggtitle("Lethal Frequency by Gene Density") +
  theme(plot.title = element_text(hjust = 0.5)) 
#theme_light()

rr_gd

cor.test(DF_data$map_freq_weightbysize, DF_data$gene_density, method ="pearson")

### ?? correlation between lethals and GC content
CGprop_lethals <- ggplot(data=DF_data, aes(x=map_freq_weightbysize, y=CG_Proportion)) + 
  geom_point(aes(color="red")) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC content") +
  ylab("Lethals mapped / Crosses completed / Deficiency size") +
  ggtitle("Lethal Frequency by GC content") +
  theme(plot.title = element_text(hjust = 0.5)) 
#theme_light()
CGprop_lethals


cor.test(DF_data$coding_density, DF_data$CG_Proportion.x, method ="pearson")




#ggplot x-axis as chromosome arm


p_2R <- ggplot(data=DF_2R, aes(x=midpoint)) + geom_line(aes(y=map_freq_weightbysize*50000000, color="red")) + geom_line(aes(y=rr_avg))
p_2R
p_2L <- ggplot(data=DF_2L, aes(x=midpoint)) + geom_line(aes(y=map_freq_weightbysize*50000000, color="red")) + geom_line(aes(y=rr_avg))
p_2L




## Test

library(tidyr)

# Reshape data into long format
DF_2R_long <- DF_2R %>%
  pivot_longer(cols = c(map_freq_weightbysize, rr_avg, coding_density, CG_Proportion),
               names_to = "variable", values_to = "value")

p_2R <- ggplot(data = DF_2R_long, aes(x = midpoint.x, y = value)) +
  geom_line(aes(color = variable)) +
  facet_grid(variable ~ ., scales = "free_y", switch = "y", labeller = as_labeller(c(
    "coding_density" = "Coding Density",
    "map_freq_weightbysize" = "Mapping Frequency",
    "rr_avg" = "Recombination Rate"
  ))) +  
  scale_color_manual(values = c("coding_density" = "red", 
                                "map_freq_weightbysize" = "black",  
                                "rr_avg" = "blue"),
                     labels = c("Coding Density", 
                                "Mapping Frequency", 
                                "Recombination Rate")) + 
  labs(color = "Legend", x = "2R") + 
  theme_minimal() + 
  theme(strip.text.y = element_text(angle = 0)) +
  #panel.grid.major.x = element_line(color = "gray", size = 0.5),  # Default gridlines
  # panel.grid.minor.x = element_blank()) +  # Hide minor gridlines if desired
  geom_vline(xintercept = c(6000000, 8000000, 9950000, 11800000, 13900000, 21000000),  # Example x-positions
             color = "black", linetype = "dashed", size = 0.2) +  # Bold lines
  geom_rect(aes(xmin = 15391154, xmax = 20276334, ymin = -Inf, ymax = Inf),  # Define x-axis range
            fill = "gray", alpha = 0.007)
p_2R

# Reshape data into long format
DF_2L_long <- DF_2L %>%
  pivot_longer(cols = c(map_freq_weightbysize, rr_avg, coding_density),
               names_to = "variable", values_to = "value")

p_2L <- ggplot(data = DF_2L_long, aes(x = midpoint.x, y = value)) +
  geom_line(aes(color = variable)) +
  facet_grid(variable ~ ., scales = "free_y", switch = "y", labeller = as_labeller(c(
    "coding_density" = "Coding Density",
    "map_freq_weightbysize" = "Mapping Frequency",
    "rr_avg" = "Recombination Rate"
  ))) +  
  scale_color_manual(values = c("coding_density" = "red", 
                                "map_freq_weightbysize" = "black",  
                                "rr_avg" = "blue"),
                     labels = c("Coding Density", 
                                "Mapping Frequency", 
                                "Recombination Rate")) + 
  labs(color = "Legend", x = "2R") + 
  theme_minimal() + 
  theme(strip.text.y = element_text(angle = 0)) +
  #panel.grid.major.x = element_line(color = "gray", size = 0.5),  # Default gridlines
  # panel.grid.minor.x = element_blank()) +  # Hide minor gridlines if desired
  geom_vline(xintercept = c(6000000, 8000000, 9950000, 11800000, 13900000, 21000000),  # Example x-positions
             color = "black", linetype = "dashed", size = 0.2) +  # Bold lines
  geom_rect(aes(xmin = 15391154, xmax = 20276334, ymin = -Inf, ymax = Inf),  # Define x-axis range
            fill = "gray", alpha = 0.007)
p_2L


##Try whole chromosome
# Reshape data into long format


DF_wholechrom_long <- DF_data %>%
  filter(!is.na(arm.x)) %>%  # Remove rows where arm is NA
  pivot_longer(cols = c(map_freq_weightbysize, rr_avg, coding_density),  
               names_to = "variable", values_to = "value")

ggplot(data = DF_wholechrom_long, aes(x = midpoint.x, y = value)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable), size=0.5) +
  facet_grid(variable ~ arm.x, scales = "free_y", switch = "y", 
             labeller = as_labeller(c(
               "coding_density" = "Coding Density",
               "map_freq_weightbysize" = "Mapping Frequency",
               "rr_avg" = "Recombination Rate (cM/Mb)"
             ))) +  
  scale_x_continuous(labels = scales::label_comma()) +  # Format x-axis labels
  scale_color_manual(values = c("coding_density" = "salmon", 
                                "map_freq_weightbysize" = "darkseagreen4",  
                                "rr_avg" = "orchid4"),
                     labels = c("Coding Density", 
                                "Mapping Frequency", 
                                "Recombination Rate (cM/Mb)")) + 
  labs(color = "Legend", x = "Deficiency MidPoint") + 
  theme_minimal() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.y = element_text(angle = 0))


##Gene dnsity
DF_wholechrom_long <- DF_data %>%
  filter(!is.na(arm.x)) %>%  # Remove rows where arm is NA
  pivot_longer(cols = c(map_freq_weightbysize, rr_avg, gene_density),  
               names_to = "variable", values_to = "value")

ggplot(data = DF_wholechrom_long, aes(x = midpoint.x, y = value)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable), size=0.5) +
  facet_grid(variable ~ arm.x, scales = "free_y", switch = "y", 
             labeller = as_labeller(c(
               "gene_density" = "Gene Density",
               "map_freq_weightbysize" = "Mapping Frequency",
               "rr_avg" = "Recombination Rate (cM/Mb)"
             ))) +  
  scale_x_continuous(labels = scales::label_comma()) +  # Format x-axis labels
  scale_color_manual(values = c("gene_density" = "salmon", 
                                "map_freq_weightbysize" = "darkseagreen4",  
                                "rr_avg" = "orchid4"),
                     labels = c("Gene Density", 
                                "Mapping Frequency", 
                                "Recombination Rate (cM/Mb)")) + 
  labs(color = "Legend", x = "Deficiency MidPoint") + 
  theme_minimal() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.y = element_text(angle = 0))
