library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(VennDiagram)            


# Read in dataframe.
datae0040raw <- read.table("data/ps_all.txt.gz", header=TRUE, stringsAsFactors = FALSE)
# Remove columns that are no longer necessary for analysis.
datae0040 <- datae0040raw %>%
  dplyr::select(well, OTU, Kingdom, Phylum, Class, Order, Family, Genus, count, relAbundance)

# Set up a metadata file. Recommended columns include well, community type 
# community type (preAbx, postAbxv1, postAbxv2, BacteOnly, BacteSwap, blank).
# media, replicate, Bacte1, Bacte2. Export as TSV.
# Use left_join(datae0040, metadata, by=c("well")) to add the metadata to your 16S data.
# Format your wells as A1, A2...H11, H12. You can use unique(datae0040$well) to check.

# Set file path for the metadata table
appendCol_path <- "C:/Users/aparr/e0040-bottom-up-communities-analyses/data/append_col_e0040.tsv"

# Import metadata table
appendCol <- read.table(appendCol_path, header = TRUE)

# Join the metadata table to the original data frame
datae0040meta <- left_join(datae0040, appendCol, by=c("well"))
write.table(datae0040meta, "datae0040meta.txt", row.names = FALSE, quote = FALSE)

# Set file path for the updated replicates column
appendRep2_path <- "C:/Users/aparr/e0040-bottom-up-communities-analyses/data/replicates2.tsv"

# Import metadata table
appendRep2 <- read.table(appendRep2_path, header = TRUE)

# Join the metadata table to the original data frame
datae0040meta <- left_join(datae0040meta, appendRep2, by=c("well"))
write.table(datae0040meta, "datae0040meta.txt", row.names = FALSE, quote = FALSE)

# Filter the data frame to include only rows with relAbundance greater than 0.1%
datae0040meta <- datae0040meta %>% filter(relAbundance > 0.001)

# e0040.A - XEA b-u
# Run the OTUs through dada 2 to get the ASVs
# Relative abundance of top 15 families


# Top 20 most abundant species
# and calculate the cumulative sum of the relative abundance represented by the species.
# pre-antibiotic
output_preAbx_specTop = datae0040meta %>% filter(communityType=="preAbx") %>%
  top_n(20, relAbundance) %>% arrange(desc(relAbundance)) %>%
  dplyr::select(count, communityType, relAbundance, media, Family) %>%
  mutate(cumulative_abundance=cumsum(relAbundance))
write.table(output_preAbx_specTop, "e0040A_preAbx_asv20.txt", row.names = FALSE, quote = FALSE)


#post-antibiotic V1
output_postAbxV1_specTop = datae0040meta %>% filter(communityType=="postAbxV1") %>%
  top_n(20, relAbundance) %>% arrange(desc(relAbundance)) %>%
  dplyr::select(count, communityType, relAbundance, media, Family) %>%
  mutate(cumulative_abundance=cumsum(relAbundance))
write.table(output_postAbxV1_specTop, "e0040A_postAbxV1_asv20.txt", row.names = FALSE, quote = FALSE)


# post-antibiotic V2
output_postAbxV2_specTop = datae0040meta %>% filter(communityType=="postAbxV2") %>%
  top_n(20, relAbundance) %>% arrange(desc(relAbundance)) %>%
  dplyr::select(count, communityType, relAbundance, media, Family) %>%
  mutate(cumulative_abundance=cumsum(relAbundance))
write.table(output_postAbxV2_specTop, "e0040A_postAbxV2_asv20.txt", row.names = FALSE, quote = FALSE)


# Raw alpha diversity and alpha diversity w/ limit of detection**
# Get alpha diversity by well
alpha_diversity_e0040 <- datae0040meta %>%
  group_by(well) %>%
  summarize(alpha_diversity_e0040 = sum(count > 0))

# Join alpha_diversity to the original table
e0040_alpha_diversity_table <- datae0040meta <- datae0040meta %>%
  left_join(alpha_diversity_e0040, by = "well")

# Save table to "out" folder
write.table(e0040_alpha_diversity_table, file = "a_diversity_by_well.txt", quote = FALSE, row.names = FALSE)

# Arrange communities by alpha diversity in descending order
dec_alpha_diversity <- datae0040meta %>%
  arrange(desc(alpha_diversity_e0040)) 


# Filter sorted alpha diversity by community type
# Pre-Abx and save table
alpha_diversity_preAbx <- dec_alpha_diversity %>%
  filter(communityType == "preAbx") %>%
  distinct(well, .keep_all = TRUE) %>%
  dplyr::select(well, relAbundance, alpha_diversity_e0040, communityType, media)
write.table(alpha_diversity_preAbx, file = "alpha_diversity_preAbx.txt", quote = FALSE, row.names = FALSE)

# Post-Abx V1 and save table
alpha_diversity_postAbxV1 <- dec_alpha_diversity %>%
  filter(communityType == "postAbxV1") %>%
  distinct(well, .keep_all = TRUE) %>%
  dplyr::select(well, relAbundance, alpha_diversity_e0040, communityType, media)
write.table(alpha_diversity_preAbx, file = "alpha_diversity_postAbxV1.txt", quote = FALSE, row.names = FALSE)

# Post-Abx V2 and save table
alpha_diversity_postAbxV2 <- dec_alpha_diversity %>%
  filter(communityType == "postAbxV2") %>%
  distinct(well, .keep_all = TRUE) %>%
  dplyr::select(well, relAbundance, alpha_diversity_e0040, communityType, media)
write.table(alpha_diversity_preAbx, file = "alpha_diversity_postAbxV2.txt", quote = FALSE, row.names = FALSE)


# Compare pre-abx in mBHI vs mGAM, YCFA, and RCM by species/rel abundance
alpha_diversity_media <- dec_alpha_diversity %>%
  filter(media %in% c("mGAM", "RCM", "YCFA") | well %in% c("A1", "B1", "C1", "D1")) %>%
  distinct(well, .keep_all = TRUE) %>%
  dplyr::select(well, relAbundance, alpha_diversity_e0040, communityType, media)
write.table(alpha_diversity_media, file = "alpha_diversity_media.txt", quote = FALSE, row.names = FALSE)


# Ensure controls look as expected; i.e. A4-D4 are blank
#A4_diversity <- datae0040meta[datae0040meta$well == "A4", "alpha_diversity_e0040"]
# Define well vector
wells <- c("A4", "B4", "C4", "D4")

# Create a new data frame with blank wells and alpha diversity
blankWellVector <- unique(datae0040meta[datae0040meta$well  %in% wells, "alpha_diversity_e0040"])

# Create a data frame with the well and a_div values
blankWellTable <- data.frame(well = wells, alpha_diversity = blankWellVector)

# Save blank well table
write.table(blankWellTable, file = "blankWellTable.txt", quote = FALSE, row.names = FALSE)






# Bar plot showing the count or relative abundance of each community type 
# Summary of rel abundance by comm
community_abundance <- aggregate(relAbundance ~ communityType, data = datae0040meta, FUN = sum)

# Create a bar plot using ggplot2
community_abundance_bar_v2 <- datae0040meta %>% ggplot() +
  geom_bar(aes(x=well, y=relAbundance, fill=factor(Family)), stat= "identity") +
  xlab("Community Type") +
  ylab("Relative Abundance") +
  ggtitle("Distribution of Community Types based on Relative Abundance")

# Save plot
save_plot("community_abundance_bar_v2.png", community_abundance_bar_v2)






# Stacked bar plots to visualize the relative abundance of different taxonomic levels within each comm
# Ensure taxonomic_levels is character vector
taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Convert taxonomic levels to character if they are not already
datae0040meta[taxonomic_levels] <- lapply(datae0040meta[taxonomic_levels], as.character)

# Calculate relative abundance for each taxonomic level within each community
community_taxa <- datae0040meta %>%
  group_by(communityType) %>%
  summarise(across(all_of(taxonomic_levels), ~ sum(relAbundance, na.rm = TRUE)))

# I need to fix this table and stackoverflow mentioned pivot_longer?
community_taxa_long <- tidyr::pivot_longer(community_taxa, 
                                           cols = -communityType,
                                           names_to = "TaxonomicLevel",
                                           values_to = "RelativeAbundance")

# Create the stacked bar plot
community_taxa_bars <- ggplot(data = community_taxa_long, aes(x = communityType, y = RelativeAbundance, fill = TaxonomicLevel)) +
  geom_bar(stat = "identity") +
  xlab("Community Type") +
  ylab("Relative Abundance") +
  ggtitle("Composition of Taxonomic Levels in Different Communities") +
  theme(legend.position = "right") +
  scale_fill_brewer(palette = "Set3")

# Save plot
save_plot("community_taxa_bars.png", community_taxa_bars)









# boxplots to compare the distribution of relative abundance values among different community types
# Create andombine the data frames for preAbx and postAbx
# Subset the data by community type
preAbx_data <- datae0040meta[datae0040meta$communityType == "preAbx", ]
postAbxV1_data <- datae0040meta[datae0040meta$communityType == "postAbxV1", ]
postAbxV2_data <- datae0040meta[datae0040meta$communityType == "postAbxV2", ]

# combined_prePost_data <- rbind(preAbx_data, postAbxV1_data, postAbxV2_data)

# Create a boxplot using ggplot2
alpha_diversity_boxPlot <- ggplot() +
  geom_boxplot(data = preAbx_data, aes(x = communityType, y = alpha_diversity_e0040), fill = "steelblue", color = "black") +
  geom_boxplot(data = postAbxV1_data, aes(x = communityType, y = alpha_diversity_e0040), fill = "steelblue", color = "black") +
  geom_boxplot(data = postAbxV2_data, aes(x = communityType, y = alpha_diversity_e0040), fill = "steelblue", color = "black") +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Distribution of Alpha Diversity by Community Type") +
  theme_bw()

# Save plot
save_plot("relAbundance_boxPlot.png", relAbundance_boxPlot)




# Plot the alpha diversity of microbial communities against the different media types used
# Create a scatter plot using ggplot2
alphaDiversity_media_scatter <- ggplot(alpha_diversity_media, aes(x = media, y = alpha_diversity_e0040)) +
  geom_point() +
  xlab("Media") +
  ylab("Alpha Diveristy") +
  ggtitle("Relative Abundance vs. Media") +
  theme_bw()

# Save plot
save_plot("alphaDiversity_media_scatter.png", alphaDiversity_media_scatter)




# Bar plot comparing the alpha diversity for each replicate group
# Create a subset of the data for the replicate groups
replicate_groups <- list(
  preAbx = c("A1", "B1", "C1", "D1"),
  postAbxV1 = c("A2", "B2", "C2", "D2"),
  postAbxV2 = c("A3", "B3", "C3", "D3"),
  mGAM = c("E1", "F1", "G1", "H1"),
  RCM = c("E2", "F2", "G2", "H2"),
  YCFA = c("E3", "F3", "G3", "H3"),
  bacteOnly = c("E4", "F4", "G4", "H4")
)

# Create a new data frame containing the replicate groups only
replicate_data <- datae0040meta[datae0040meta$well %in% unlist(replicate_groups), ]

# Create a boxplot for each replicate group
replicateReliability_boxplot <- ggplot(replicate_data, aes(x = well, y = relAbundance, fill = well)) +
  geom_boxplot() +
  xlab("Replicate Group") +
  ylab("Relative Abundance") +
  ggtitle("Comparison of Replicate Groups") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create a new factor variable for the well column based on the replicate groups
replicate_order <- unlist(replicate_groups)
replicate_data$well <- factor(replicate_data$well, levels = replicate_order)

# Create a color palette for the replicate groups
group_colors <- c("preAbx" = "steelblue",
                  "postAbxV1" = "darkorange",
                  "postAbxV2" = "darkgreen",
                  "mGAM" = "purple",
                  "RCM" = "deeppink",
                  "YCFA" = "gold",
                  "bacteOnly" = "darkcyan")

# Create a boxplot for each replicate group with customized color palette
replicateReliability_boxplot <- ggplot(replicate_data, aes(x = well, y = relAbundance, fill = well)) +
  geom_boxplot() +
  scale_fill_manual(values = group_colors) +
  xlab("Replicate Group") +
  ylab("Relative Abundance") +
  ggtitle("Comparison of Replicate Groups") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot
save_plot("replicateReliability_boxplot.png", replicateReliability_boxplot)

  



# Venn diagram to display the shared taxa between pre-antibiotic, post-antibiotic V1, and post-antibiotic V2 communities

# Extract microbial families for each community type
preAbx_families <- unique(datae0040meta[datae0040meta$communityType == "preAbx", "Family"])
postAbxV1_families <- unique(datae0040meta[datae0040meta$communityType == "postAbxV1", "Family"])
postAbxV2_families <- unique(datae0040meta[datae0040meta$communityType == "postAbxV2", "Family"])

# Prepare the input for the Venn diagram, removing NAs
taxa_venn_input <- list(
  preAbx = preAbx_families[!is.na(preAbx_families)],
  postAbxV1 = postAbxV1_families[!is.na(postAbxV1_families)],
  postAbxV2 = postAbxV2_families[!is.na(postAbxV2_families)]
)

# Create the Venn diagram
taxa_venn_plot <- venn.diagram(
  x = taxa_venn_input,
  category.names = c("Pre-Abx", "Post-Abx V1", "Post-Abx V2"),
  fill = c("steelblue", "darkorange", "darkgreen"),
  alpha = 0.7, # Controls the transparency
  cex = 1.5, # Adjusts size
  filename = NULL
)

# Add a title to the Venn diagram
title <- "Shared Families"
taxa_venn_plot <- grid.arrange(
  grobs = list(taxa_venn_plot, textGrob(title, gp = gpar(fontsize = 16, fontface = "bold"))),
  nrow = 2,
  heights = c(0.8, 0.2)
)

# Display the Venn diagram
grid.draw(taxa_venn_plot)

# Save plot
save_plot("taxa_venn_plot.png", taxa_venn_plot)


# Compare average alpha diversity of e0040.A across all conditions (w/post expected to be lowest)
#     ## Plot x-diff cond and y-alpha
#     ## Can do a boxplot…points with mean




# e0040.B - XBA b-u

# Relative abundance of top 15 families


# Top 20 most abundant species
# Bacte Swap 
output_e0040A_specTop = datae0040meta %>% filter(communityType=="bacteSwap") %>%
  top_n(20, relAbundance) %>% arrange(desc(relAbundance)) %>%
  dplyr::select(count, communityType, relAbundance, media, Family) %>%
  mutate(cumulative_abundance=cumsum(relAbundance))
write.table(output_e0040A_specTop, "e0040B_bacteSwap_asv20.txt", row.names = FALSE, quote = FALSE)

# Bacte Only
output_e0040A_specTop = datae0040meta %>% filter(communityType=="bacteOnly") %>%
  top_n(20, relAbundance) %>% arrange(desc(relAbundance)) %>%
  dplyr::select(count, communityType, relAbundance, media, Family) %>%
  mutate(cumulative_abundance=cumsum(relAbundance))
write.table(output_e0040A_specTop, "e0040B_bacteOnly_asv20.txt", row.names = FALSE, quote = FALSE)


# Alpha diversity
# BacteOnly and save table
alpha_diversity_bacteOnly <- dec_alpha_diversity %>%
  filter(communityType == "bacteOnly") %>%
  distinct(well, .keep_all = TRUE) %>%
  dplyr::select(well, relAbundance, alpha_diversity_e0040, communityType, media)
write.table(alpha_diversity_preAbx, file = "alpha_diversity_bacteOnly.txt", quote = FALSE, row.names = FALSE)

# BacteSwap and save table
alpha_diversity_bacteSwap <- dec_alpha_diversity %>%
  filter(communityType == "bacteSwap") %>%
  distinct(well, .keep_all = TRUE) %>%
  dplyr::select(well, relAbundance, alpha_diversity_e0040, communityType, media)
write.table(alpha_diversity_preAbx, file = "alpha_diversity_bacteSwap.txt", quote = FALSE, row.names = FALSE)


# Parse out all of the present bacteroides species along with their rel abundance





# Compare which bacteroides species were able to coexist by making a table of all the bacteroides swap comm and appending two col w rel abundance of each bacteroides strain that was in the well
