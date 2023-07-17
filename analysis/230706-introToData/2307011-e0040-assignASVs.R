# VERSION 1.1
# Script version 1.1: Change parsing of sample names
# to reflect new standard naming formats.
# Samples should be named in form e####_HKX###-[A-H][01-12].
# Also revise the species names of sequences that have a family ID but no species ID
# to reflect the family ID (i.e. Bacte.NA instead of just NA in the species name).

library(sangeranalyseR)
library(seqinr)
library(dada2)
library(phyloseq)
library(tidyverse)
library(cowplot)
library(ggfittext)


# Update version each time new output files are generated to keep track.
VERSION <- "V1.1"


# Generate individual FASTA files. ----------------------------------------

# Code below from Becca C, with modifications to increase the quality cutoff.

files <- "C:/Users/aparr/e0040-bottom-up-communities-analyses/data/ps_all.txt.gz"
lapply(files, function(x) {
  sangerReadF <- SangerRead(inputSource           = "ABIF",
                            readFeature           = "Forward Read",
                            readFileName          = x,
                            geneticCode           = GENETIC_CODE,
                            TrimmingMethod        = "M1",
                            M1TrimmingCutoff      = 0.01,
                            M2CutoffQualityScore  = NULL,
                            M2SlidingWindowSize   = NULL,
                            baseNumPerRow         = 100,
                            heightPerRow          = 200,
                            signalRatioCutoff     = 0.33,
                            showTrimmed           = TRUE)
  # apply function
  writeFasta(sangerReadF,
             outputDir         = 'fastaFinal',
             compress          = FALSE,
             compression_level = NA)
})

# Generate concatenated FASTA file. -------------------------------------------------------

# Import trimmed Sanger sequences and concatenate into a single FASTA file.

FASTAfiles <- list.files(path="fastaFinal", pattern="*.fa", full.names=TRUE, recursive=FALSE)
seqs <- sapply(FASTAfiles, read.fasta)
write.fasta(sequences=sapply(seqs,getSequence),
            names=sapply(seqs, getName),
            file.out="fastaFinal/seqsAll.fa",
            nbchar=100)

# Assign taxonomy. --------------------------------------------------------

# Import the trimmed Sanger sequences for each isolate.

seqsraw <- read.fasta("fastaFinal/seqsAll.fa")
seqs <- sapply(seqsraw, function(x) getSequence(x, as.string=TRUE))
# Convert to a data frame.
dataSeqs <- data.frame(names(seqs), unname(unlist(seqs)))
colnames(dataSeqs) <- c("sample","seq")

# Remove sequences shorter than 50 nts so that assignTaxonomy will not throw an error. 
# Replace sequences with 50 "N" so that they will remain in dataframe and can be exported
# to strain log without causing gaps.
dataSeqsTrimmed <- dataSeqs %>% 
  mutate(seqSize=nchar(seq), 
         keep=ifelse(seqSize>=50, TRUE, FALSE))
shortSeqs <- which(dataSeqsTrimmed$keep==FALSE)
seqsTrimmed <- seqs
seqsTrimmed[shortSeqs] <- strrep("n", 50)

# Use dada2/phyloseq to assign taxonomy based on the Silva database.
#ref_fasta <- "../../config/Sil138_BAhq_V4_train_set.fa.gz"
ref_fasta_custom <- "../../config/dada2-ref-custom.fa"
set.seed(0)
#taxa <- assignTaxonomy(unlist(seqs),refFasta = ref_fasta,multithread=FALSE)
# Use dada2/phyloseq to assign taxonomy based on our custom database, with lower confidence.
taxa <- assignTaxonomy(unlist(seqsTrimmed),refFasta = ref_fasta_custom,multithread=FALSE, minBoot=15)
dataTaxa <- as.data.frame(tax_table(taxa))
dataTaxa$seq <- gsub("\\..*","",rownames(dataTaxa))
dataTaxa$sample <- names(seqsTrimmed)

# Reset taxonomy of dummy sequences of 50 "N" back to NA.
dataTaxa <- dataTaxa %>% 
  mutate(Family=replace(Family, seq==strrep("N", 50), NA),
         Species=replace(Species, seq==strrep("N", 50), NA),
         Family=replace(Family, Family=="NA", NA))

# For sequences whose species is annotated as NA but have a family-level ID,
# replace the species name of NA with the five-letter family code followed by NA.
# Import five-letter family codes and bind to data.
familyCoding <- read.table("../../config/familyCodingFull.tsv", header=TRUE) %>% 
  dplyr::rename(Family=family)
dataTaxa <- dataTaxa %>% left_join(familyCoding, by="Family")
dataTaxa <- dataTaxa %>%
  mutate(Species=ifelse(!is.na(Species), Species,
                        ifelse(!is.na(Family), paste0(code, ".NA"), "NA.NA")))

# Reformat and export the taxonomy table for analysis.
dataTaxaExport <- dataTaxa %>%
  dplyr::select(sample, Kingdom, Phylum, Class, Order,
                Family, Genus, Species, seq) %>%
  mutate(Silva=paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep="-"))
write.table(dataTaxaExport, paste0("seqClassification-",VERSION,".txt"),
            row.names=FALSE)

# Separate the sample name into pieces based on the assumption that
# the samples are named in the format e####_HKX###-[A-H][01-12].
# CHANGE THESE IF SAMPLE NAME FORMATTING IS DIFFERENT.
dataTaxa <- dataTaxa %>%
  separate(sample, into=c(NA,"Elim","experiment","plateWell","blank","ElimWell"), sep="_", remove=FALSE) %>%
  separate(plateWell, into=c("plate","well"), sep="-", remove=FALSE)

# Reformat dataTaxa for export to strain log.
# CHANGE THESE IF SAMPLE NAME FORMATTING IS DIFFERENT.
dataTaxa <- dataTaxa %>% 
  mutate(well=sub("\\..*","",well),
         row=substr(well,1,1),
         col=substr(well,2,3)) %>%
  arrange(col,row)

# Export taxonomy to genus level for strain log.
taxToGenus <- dataTaxa %>% 
  select(well, experiment, Species, Phylum, Class, Order, Family, Genus, seq) %>% 
  filter(!is.na(Family))
write.table(taxToGenus, paste0("taxToGenus-",VERSION,".txt"), row.names=FALSE, quote=FALSE, sep="\t")

# Generate plate map. -----------------------------------------------------

# Change plot default appearance
theme_set(theme_cowplot())

# Import taxa color palette.
palette <- read.table("../../config/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>% 
  mutate(taxashort=gsub(".*\\.","",taxa))

# Bind taxa colors by family.
dataTaxa <- dataTaxa %>% left_join(palette %>% select(taxashort, taxa, hex) %>%
                                     dplyr::filter(taxashort %in% dataTaxa$Family) %>% 
                                     dplyr::rename(Family=taxashort) %>% 
                                     group_by(Family) %>% 
                                     dplyr::mutate(hexnum=row_number()) %>% 
                                     filter(hexnum==1) %>% 
                                     dplyr::select(!hexnum),
                                   by="Family")

# Extract color palette. Convert unknown taxa to dark gray.
taxaPalette <- dataTaxa %>% ungroup() %>%
  select(Family, taxa, hex) %>%
  unique() %>%
  mutate(hex=ifelse(is.na(hex),"#615c5c", hex)) %>% 
  arrange(taxa)
taxaPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaPaletteList) <- taxaPalette$taxa

# Plot a map of the samples, with their family assignment.
p_plateMap <- dataTaxa %>%
  mutate(code=replace(code, is.na(code), "NA"),
         species=paste0(code,"\n",substr(Species,7,10)),
         col=as.numeric(col)) %>%
  ggplot() +
  geom_tile(aes(x=col, y=fct_rev(as.factor(row)), fill=factor(taxa)), 
            color="black") +
  geom_fit_text(aes(x=col, y=fct_rev(as.factor(row)), label=species),
                padding.x=grid::unit(0.5, "mm"), padding.y=grid::unit(0.5,"mm")) +
  scale_fill_manual(values=taxaPaletteList) +
  scale_x_continuous(breaks=seq(1, 12)) +
  theme(legend.position="none") +
  theme(axis.title = element_blank()) +
  facet_wrap(~plate)
save_plot(paste0("plateMap-",VERSION,".png"),p_plateMap, ncol=1.5)
