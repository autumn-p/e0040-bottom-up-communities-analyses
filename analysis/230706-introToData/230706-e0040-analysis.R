library(tidyverse)

# Read in dataframe.
datae0040raw <- read.table("data/ps_all.txt.gz", header=TRUE, stringsAsFactors = FALSE)
# Remove columns that are no longer necessary for analysis.
datae0040 <- datae0040raw %>%
  dplyr::select(well, OTU, Kingdom, Phylum, Class, Order, Family, Genus, count, relAbundance)
