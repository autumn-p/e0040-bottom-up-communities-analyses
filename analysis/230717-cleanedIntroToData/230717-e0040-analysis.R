#Import libraries
library(gridExtra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(VennDiagram) 
library(DT)

# Set up file paths
outPath = "analysis/out/71723_FinalOut" # out
dataframePath = "data/ps_all.txt.gz" #raw data
appendCol_path <- "data/metadatae0040.tsv" #metadata
KCHpalette_path <- "config/KCHcolors-Silva-partial.txt"

# Import the color palette
KCHpalette <- read.table(KCHpalette_path, header = TRUE)

# Make a named list
KCHpalettevector <- KCHpalette$hex
names(KCHpalettevector) <- KCHpalette$taxashort

# Read in dataframe.
datae0040raw <- read.table(dataframePath, header=TRUE, stringsAsFactors = FALSE)
# Remove columns that are no longer necessary for analysis.
datae0040 <- datae0040raw %>%
  dplyr::select(well, OTU, Kingdom, Phylum, Class, Order, Family, Genus, count, relAbundance)

# Import metadata table
appendCol <- read.table(appendCol_path, header = TRUE)

# Join the metadata table to the original data frame
datae0040meta <- left_join(datae0040, appendCol, by=c("well"))
#write.table(datae0040meta, "datae0040meta.txt", row.names = FALSE, quote = FALSE)

# Filter the data frame to include only rows with relAbundance greater than 0.1%
datae0040meta <- datae0040meta %>% filter(relAbundance > 0.001)

# Raw alpha diversity and alpha diversity w/ limit of detection**
# Get alpha diversity by well
alpha_diversity_e0040 <- datae0040meta %>%
  group_by(well) %>%
  summarize(alpha_diversity_e0040 = sum(count > 0))

# Join alpha_diversity to the original table
e0040_alpha_diversity_table <- datae0040meta <- datae0040meta %>%
  left_join(alpha_diversity_e0040, by = "well")

# Save metadata table to "out" folder
write.table(datae0040meta, paste0(outPath, "/datae0040meta.txt"), quote = FALSE, row.names = FALSE)

# Create a stripped-down color palette containing only the families
# present in the dataset.
datae0040meta <- datae0040meta %>%
  mutate(fullSilvataxonomy=paste(Kingdom,Phylum,Class,Order,Family, sep="."))
KCHpalettee0040 <- KCHpalette %>%
  filter(taxa %in% sort(unique(datae0040meta$fullSilvataxonomy))) %>%
  mutate(taxashort=ifelse(taxashort=="", gsub(".*\\.","",taxa), taxashort))
# Make a named list
KCHpalettee0040vector <- KCHpalettee0040$hex
names(KCHpalettee0040vector) <- KCHpalettee0040$taxashort


# Bar plot showing the count or relative abundance of each community type 
# Summary of rel abundance by comm
community_abundance <- aggregate(relAbundance ~ communityType, data = datae0040meta, FUN = sum)

# Create a bar plot using ggplot2
community_abundance_bar <- datae0040meta %>% 
  filter(media== "mBHI", communityType!="bacteSwap", communityType!="blank") %>% 
  ggplot() +
  geom_bar(aes(x=replicate, y=relAbundance, fill=factor(Family)), color = "black", stat= "identity") +
  scale_fill_manual(values=KCHpalettee0040vector) +
  facet_wrap(~communityType) +
  xlab("Community Type") +
  ylab("Relative Abundance") +
  ggtitle("Distribution of Community Types based on Relative Abundance")
community_abundance_bar

# Save plot
save_plot(paste0(outPath, "/community_abundance_bar_preAbx_colored.png"), community_abundance_bar, base_width = 5, base_height = 5)




# I want all the taxonomic levels grouped by color with Kingdom at the top and Genus at the
# bottom with black lines within each color dividing the levels into individual different 
# Kingdoms; Phylums, Classes, etc.For example, the thickness of the family Bacteroidaceae will
# correspond to the sum of all relative abundances of bacteroidaceae in rows that are labelled A1
# Stacked bar plots to visualize the relative abundance of different taxonomic levels within each community

# Ensure taxonomic_levels is character vector
taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Convert taxonomic levels to character
datae0040meta[taxonomic_levels] <- lapply(datae0040meta[taxonomic_levels], as.character)

# Calculate relative abundance for each taxonomic level within each community
community_taxa <- datae0040meta %>%
  filter(well == "A1") %>%
  #group_by(communityType) %>%
  summarise(across(all_of(taxonomic_levels), sum, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_longer(cols = -communityType, names_to = "TaxonomicLevel", values_to = "value")

# Create the stacked bar plot
community_taxa_bars <- ggplot(community_taxa) +
  geom_bar(aes(x = communityType, y = value, fill = TaxonomicLevel), color = "black", stat = "identity", position = "stack") +
  xlab("Community Type") +
  ylab("Relative Abundance") +
  ggtitle("Composition of Taxonomic Levels in Different Communities") +
  theme(legend.position = "right") +
  scale_fill_brewer(palette = "Set3")

# Save plot
save_plot(paste0(outPath, "/community_taxa_bars.png"), community_taxa_bars, base_width = 5, base_height = 5)

### Trying a DIFFERENT way

# Ensure taxonomic_levels is character vector
taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Convert taxonomic levels to factors
datae0040meta[taxonomic_levels] <- lapply(datae0040meta[taxonomic_levels], as.factor)

# Check data types of columns
print(sapply(datae0040meta[taxonomic_levels], class))

# Identify unexpected values
unexpected_values <- lapply(datae0040meta[taxonomic_levels], function(x) unique(x[!is.na(x) & !is.na(as.numeric(as.character(x)))]))
print(unexpected_values)

# Calculate relative abundance for each taxonomic level within each community
community_taxa <- datae0040meta %>%
  filter(well == "A1") %>%
  summarise(across(all_of(taxonomic_levels), ~sum(as.numeric(as.character(.)), na.rm = TRUE, na.fail = FALSE)), .groups = "drop") %>%
  tidyr::pivot_longer(cols = -communityType, names_to = "TaxonomicLevel", values_to = "value")

# Create the stacked bar plot
community_taxa <- datae0040meta %>%
  filter(well == "A1") %>%
  group_by(communityType) %>%
  summarise(across(all_of(taxonomic_levels), ~sum(as.numeric(as.character(.)), na.rm = TRUE), .names = "{.col}_sum")) %>%
  tidyr::pivot_longer(cols = -communityType, names_to = "TaxonomicLevel", values_to = "value") %>%
  mutate(value = replace_na(value, 0))


# Print the plot
print(community_taxa_bars)

na_values <- community_taxa$value[is.na(community_taxa$value)]
print(na_values)

community_taxa <- community_taxa %>%
  mutate(value = as.numeric(replace_na(value, 0)))

# Create the stacked bar plot
community_taxa_bars <- ggplot(community_taxa) +
  geom_bar(aes(x = communityType, y = value, fill = TaxonomicLevel), color = "black", stat = "identity", position = "stack") +
  xlab("Community Type") +
  ylab("Relative Abundance") +
  ggtitle("Composition of Taxonomic Levels in Different Communities") +
  theme(legend.position = "right") +
  scale_fill_brewer(palette = "Set3")

# Print the plot
print(community_taxa_bars)


###Trying a different WAY


# Boxplots to compare the distribution of relative abundance values among different community types
# Create andombine the data frames for preAbx and postAbx
# Subset the data by community type
# combined_prePost_data <- rbind(preAbx_data, postAbxV1_data, postAbxV2_data)

# Create a scatterplot using ggplot2
alpha_diversity_boxPlot <- datae0040meta %>%
  filter(communityType %in% c("preAbx", "postAbxV1", "postAbxV2") & media == "mBHI") %>%
  ggplot() +
  geom_point(aes(x = communityType, y = alpha_diversity_e0040), color = "black", position = position_jitter(width = 0.2, height = 0)) +
  xlab("Community Type") +
  ylab("Alpha Diversity") +
  ggtitle("Distribution of Alpha Diversity by Community Type") +
  theme_bw()


# Save plot
save_plot(paste0(outPath, "/alphaDiversity_boxPlot.png"), alpha_diversity_boxPlot, base_width = 5, base_height = 5)


# jittering and make into a scatter plot instead of a boxplot


# Plot the alpha diversity of microbial communities against the different media types used
# Create a scatter plot using ggplot2
alphaDiversity_media_scatter <- datae0040meta %>%
  filter(media %in% c("mBHI", "mGAM", "YCFA", "RCM") & communityType == "preAbx") %>%
  distinct(well, .keep_all = TRUE) %>% ggplot() +
  geom_point(aes(x = media, y = alpha_diversity_e0040)) +
  xlab("Media") +
  ylab("Alpha Diveristy") +
  ggtitle("Relative Abundance vs. Media") +
  theme_bw()

# Save plot
save_plot(paste0(outPath, "/alphaDiversity_media_scatter.png"), alphaDiversity_media_scatter)





# Plot comparing the alpha diversity for each replicate group

filtered_data <- datae0040meta %>%
  filter(media == "mBHI" & communityType == "preAbx") %>%
  distinct(well, .keep_all = TRUE)

alpha_diversity_boxPlot <- filtered_data %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(replicate), y = alpha_diversity_e0040), fill = "steelblue", color = "black") +
  xlab("Replicate") +
  ylab("Alpha Diversity") +
  ggtitle("Alpha Diversity by Replicate") +
  theme_bw()

# Save plot
save_plot(paste0(outPath, "/replicateReliability_boxplot.png"), alpha_diversity_boxPlot)






# Venn diagram to display the shared taxa between pre-antibiotic, post-antibiotic V1, and post-antibiotic V2 communities
# Extract microbial families for each community type
community_types <- c("preAbx", "postAbxV1", "postAbxV2")
family_lists <- lapply(community_types, function(ctype) {
  families <- unique(datae0040meta[datae0040meta$communityType == ctype, "Family"])
  families <- families[!is.na(families)]
})

# Prepare the input for the Venn diagram
taxa_venn_input <- setNames(family_lists, community_types)

# Create the Venn diagram
taxa_venn_plot <- venn.diagram(
  x = taxa_venn_input,
  category.names = community_types,
  fill = c("steelblue", "darkorange", "darkgreen"),
  alpha = 0.7,
  cex = 1.5,
  filename = NULL
)

# Add a title to the Venn diagram
title <- "Shared Families"
taxa_venn_plot <- grid.arrange(
  grobs = list(taxa_venn_plot, textGrob(title, gp = gpar(fontsize = 16, fontface = "bold"))),
  nrow = 2,
  heights = c(0.8, 0.2)
)

# Display and save the Venn diagram
grid.draw(taxa_venn_plot)
save_plot(paste0(outPath, "/taxa_venn_plot2.png"), taxa_venn_plot)




# Table of the ASVs per families in pre, post v1 and post v2
# Filter data for the specified communities
community_types_list <- c("preAbx", "postAbxV1", "postAbxV2")

# Create a table with families and their counts for each recipient community
e0040_OTU_table <- datae0040meta %>%
  filter(communityType %in% community_types_list) %>%
  group_by(communityType, Family) %>%
  summarize(family_count = n_distinct(OTU)) %>%
  spread(Family, family_count, fill = 0)

# Convert the table to a format suitable for interactive display
# factor and levels sets the custom order of the factor levels based on the order of community_types_list
e0040_OTU_table_for_display <- e0040_OTU_table %>%
  mutate(community_type = factor(communityType, levels = community_types_list)) %>%
  arrange(community_type)

# Create an interactive dataTable
datatable(e0040_OTU_table_for_display,
          options = list(pageLength = 10,
                         lengthMenu = c(10, 20, 50),
                         dom = 't',
                         scrollX = TRUE),
          rownames = FALSE)

# Saving datatable to an HTML file
e0040_OTU_per_recipient_datatable <- datatable(e0040_OTU_table_for_display)
saveWidget(e0040_OTU_per_recipient_datatable, paste0(outPath, "/e0040_OTU_table.html"))



# Filter data for the specified communities
community_types_list <- c("preAbx", "postAbxV1", "postAbxV2")

# Create a table with families for each recipient community
e0040_family_table <- datae0040meta %>%
  filter(communityType %in% community_types_list & replicate == "1" & media == "mBHI") %>%
  group_by(communityType, Family) %>%
  summarize(present = any(!is.na(Family))) %>%
  spread(Family, present, fill = FALSE)

# Convert the table to a format suitable for interactive display
e0040_family_table_for_display <- e0040_family_table %>%
  mutate(community_types = factor(communityType,
                                      levels = c("preAbx", "postAbxV1", "postAbxV2"))) %>%
  arrange(community_types)  # Reorder rows based on recipient_community levels

# Select only the columns I want to display (excluding "community_types")
e0040_columns_to_display <- setdiff(names(e0040_family_table_for_display), "community_types")

# Create an interactive dataTable with background color formatting
e0040_Family_presence_per_recipient_datatable <- datatable(e0040_family_table_for_display[, e0040_columns_to_display],
                                                           options = list(pageLength = 10,
                                                                          lengthMenu = c(10, 20, 50),
                                                                          dom = 't',
                                                                          scrollX = TRUE),
                                                           rownames = FALSE) %>%
  formatStyle(e0040_columns_to_display, valueColumns = e0040_columns_to_display,
              backgroundColor = styleEqual(c(FALSE, TRUE), c("red", "lightgreen")))

# Saving datatable to an HTML file
saveWidget(e0040_Family_presence_per_recipient_datatable, paste0(outPath, "/e0040_family_presence_table.html"))




