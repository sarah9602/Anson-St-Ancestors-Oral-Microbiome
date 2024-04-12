## Script to (1) create PCA plots based on Metaphalan4 relative abundance data with COMPARATIVE dental calclus datasets;
## (2) Color PCA plots based on time period, location, and number of species
## (3) Run PERMANOVA association testing to determine significance
## (4) Create a comparative PCA bi-plot
## Script created by R. Fleskes (2024) modified from scripts provided from Velsko et al. 2019
## https://github.com/ivelsko/smoking_calculus


## Load libraries
library(devtools)
library(mixOmics)
library(phyloseq)
library(vegan)
library(ape)
library(ggplot2)
library(viridis)
library(tidyverse)
library(ggrepel)
library(pander)
library(dplyr)
library(rstatix)

## Read in species-level table from MetaPhlAn3 and Sample metadata file
setwd("~/Dropbox (Dartmouth College)/14_CODING/R-Studio Projects/DentalCalculus/DentalCalculus/Github")
OTUIN<-otu_table(read.table("comparative_metaphalan.txt", header=TRUE, row.names=1),taxa_are_rows = TRUE)
SAMPDATA<-sample_data(read.table("PCA_Metadata.txt",header=TRUE, row.names=1,sep="\t"))

## Create PCA plot
psobject<-phyloseq(OTUIN,SAMPDATA)
PRobject<-t(otu_table(psobject))+1
pca.result<-pca(PRobject,logratio='CLR')
coords_pca<-as.data.frame(pca.result$variates)
loadings_pca<-as.data.frame(pca.result$loadings)
FinalTable<-merge(coords_pca,as.data.frame(SAMPDATA),by="row.names")

## Obtain number of species
NUMSPECIES_OTUIN <- read.table("comparative_metaphalan.txt", header=TRUE, row.names=1)
row_counts <- NUMSPECIES_OTUIN %>%
  summarise(across(everything(), ~sum(. > 0)))
print(row_counts)
transposed_numberofspecies <- t(row_counts)
print(transposed_numberofspecies)
write.table(transposed_numberofspecies, file="comparative_metaphalan_numberofspecies.tsv", row.names = TRUE, col.names = TRUE, quote = FALSE, sep ="\t")

## Merge species counts info with FinalTable
SpeciesCounts <- read.table("comparative_metaphalan_numberofspecies.tsv", header=TRUE, row.names = NULL)
colnames(SpeciesCounts)[1] <- "Row.names"
colnames(SpeciesCounts)[2] <- "SpeciesCount"
FinalTable_Speces<-merge(FinalTable,as.data.frame(SpeciesCounts),by="Row.names")

## Format other Variables in FinalTable
FinalTable$TimePeriod <- factor(FinalTable$TimePeriod)
FinalTable$Population <- factor(FinalTable$Population)
FinalTable$TimePeriodLocation <- factor(FinalTable$TimePeriodLocation)
FinalTable$Location <- factor(FinalTable$Location)

## Create color and shape codes used for Plots
timeperiod_code_colors <- c("AncientNorthAmerica"="lightblue",
                            "AncientAfrica"="lightblue",
                            "HistoricAfrica"="#7AD151FF",
                            "HistoricEurope"="#7AD151FF",
                            "HistoricAsia"="#7AD151FF",
                            "HistoricNorthAmerica_AnsonStreetAncestors"="black",
                            "ModernEurope"="orange")
timeperiod_code_colors_black <- c("AncientNorthAmerica"="lightblue",
                            "AncientAfrica"="lightblue",
                            "HistoricAfrica"="#7AD151FF",
                            "HistoricEurope"="#7AD151FF",
                            "HistoricAsia"="#7AD151FF",
                            "HistoricNorthAmerica_AnsonStreetAncestors"="black",
                            "ModernEurope"="orange")
timeperiod_code_size <- c("AncientNorthAmerica"=2,
                          "AncientAfrica"=2,
                          "HistoricAfrica"=2,
                          "HistoricEurope"=2,
                          "HistoricAsia"=2,
                          "HistoricNorthAmerica_AnsonStreetAncestors"=3,
                          "ModernEurope"=2)
timeperiod_code_shapes_different <- c("AncientNorthAmerica"=0,
                            "AncientAfrica"=1,
                            "HistoricAfrica"=2,
                            "HistoricEurope"=3,
                            "HistoricAsia"=4,
                            "HistoricNorthAmerica_AnsonStreetAncestors"=12,
                            "ModernEurope"=5)

## Make Comparative PCA Plots
plot1<-ggplot(data=FinalTable, aes(X.PC1, X.PC2, color = TimePeriodLocation, shape = TimePeriodLocation, size = TimePeriodLocation))+
  geom_point()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = timeperiod_code_colors_black, name="")+  
  scale_size_manual(values = timeperiod_code_size, name="")+
  scale_shape_manual(values = timeperiod_code_shapes_different, name="")+
  xlim(-5.5,5.5)+ylim(-5.5,5.5)+
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 14), legend.position = "none")
ggsave("PCA_ComparativePlotSpecies_ALL.png", plot = plot1, width = 6, height = 6, dpi = 400)

SpeciesNumberplot<-ggplot(data=FinalTable_Speces, aes(X.PC1, X.PC2, color = SpeciesCount, shape = TimePeriodLocation, size = TimePeriodLocation))+
  geom_point()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_size_manual(values = timeperiod_code_size, name="")+
  scale_shape_manual(values = timeperiod_code_shapes_different, name="")+
  scale_color_gradientn(colours = c("cyan", "blue", "green3", "orange", "red"),
                        values = scales::rescale(c(20, 60, 100, 140, 160))) +
  xlim(-5.5,5.5)+ylim(-5.5,5.5)+
  labs(x = "PC1", y = "PC2")
ggsave("PCAComparativePlotSpecies_ALL_SpeciesCount.png", plot = SpeciesNumberplot, width = 9, height = 6, dpi = 400)


##########################################
## PERMANOVA TESTING USING VEGAN ON THE ##
## PCA AND COMPARATIVE METADATA RESULTS ##
##########################################

## Create Euclidean distance matrix using Vegan Package
diss_matrix <- vegdist(pca.result$x, method = "euclidean", na.rm = F)

## Create Variables based on Metadata
TimePeriod_variable <- FinalTable$TimePeriod
TimePeriodLocation_variable <- FinalTable$TimePeriodLocation
Location_variable <- FinalTable$Location
Population_variable <- FinalTable$Population

## Run PERMANOVA Association Testing using Adonis2
perm_result_TimePeriod <- adonis2(diss_matrix ~ TimePeriod, FinalTable, permutations = 999, by = "margin")
perm_result_TimePeriodLocation <- adonis2(diss_matrix ~ TimePeriodLocation, FinalTable, permutations = 999, by = "margin")
perm_result_Location <- adonis2(diss_matrix ~ Location, FinalTable, permutations = 999, by = "margin")
perm_result_Population <- adonis2(diss_matrix ~ Population, FinalTable, permutations = 999, by = "margin")

## Filter and Write PERMANOVA results
all_permanova_results <- bind_rows(perm_result_TimePeriod, perm_result_TimePeriodLocation, perm_result_Location, perm_result_Population) %>% rownames_to_column("ID")

all_permanova_results_filtered <- all_permanova_results %>%
  as_tibble() %>%
  drop_na(`F`) %>%
  arrange(desc(`F`)) %>%
  arrange(desc(R2))

write.csv(all_permanova_results_filtered, file = "all_permanova_timeperiodlocation.csv", quote =  F)


###################################################
## COMPARATIVE BIPLOT ANALYSIS USING PCA RESULTS ##
###################################################

# Extract explained variance
pca_expl_vars <- paste0(round(pca.result$variance * 100, 2), "%")

## Extract loadings (eigenvectors) from PC analysis
pca_loadings_df <- as.data.frame(pca.result$loadings$X) %>%
  rownames_to_column(var = "Species")

## Identify the 10 species that have highest positive and negative loadings in PC1 and PC2
pc1_pws <- pca_loadings_df %>%
  top_n(10, PC1) %>%
  arrange(desc(PC1)) %>%
  mutate(Species_Order = 1:10)

pc1_neg <- pca_loadings_df %>%
  top_n(-10, PC1) %>%
  arrange(PC1) %>%
  mutate(Species_Order = -1:-10)

pc2_pws <- pca_loadings_df %>%
  top_n(10, PC2)%>%
  arrange(desc(PC2)) %>%
  mutate(Species_Order = 1:10)

pc2_neg <- pca_loadings_df %>%
  top_n(-10, PC2) %>%
  arrange(PC2) %>%
  mutate(Species_Order = -1:-10)

## Calculate correction factor from ggfortify::autoplot
corr_lam <- pca.result$sdev[c("PC1", "PC2")] * sqrt(nrow(FinalTable))

## Create BiPlot
plot_bi <- FinalTable %>%
  mutate(PC1 = X.PC1 / corr_lam[1],
         pc2_residuals = X.PC2 / corr_lam[2]) %>%
  ggplot(.,
         aes(x = PC1, y = pc2_residuals)) +
  geom_point(aes(color = TimePeriodLocation, shape = TimePeriodLocation, size = TimePeriodLocation))+
  geom_segment(data = pc1_pws,
               aes(xend = PC1, yend = PC2),
               x = 0, y = 0, colour = "black", size = 0.5,
               arrow = arrow(length = unit(0.03, "npc"))) +
  geom_label_repel(data = pc1_pws,
                   aes(x = PC1, y = PC2, label = Species_Order),
                   size = 2.5, colour = "black", label.padding = 0.2) +
  geom_segment(data = pc1_neg,
               aes(xend = PC1, yend = PC2),
               x = 0, y = 0, colour = "grey50", size = 0.5,
               arrow = arrow(length = unit(0.03, "npc"))) +
  geom_label_repel(data = pc1_neg,
                   aes(x = PC1, y = PC2, label = Species_Order),
                   size = 2.5, colour = "black", label.padding = 0.2) +
  labs(x = "PC1", y = "PC2") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = timeperiod_code_colors, name="")+  
  scale_size_manual(values = timeperiod_code_size, name="")+
  scale_shape_manual(values = timeperiod_code_shapes_different, name="") +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  guides(
    color = guide_legend(title = NULL, ncol = 2),
    shape = guide_legend(title = NULL, ncol = 2),
    size = guide_legend(title = NULL))
plot_bi
ggsave("BiPlot_metaphalan4_COMPARATIVE.png", plot = plot_bi, width = 6, height =5.5, dpi = 400)

## Create Species Loading tables
pandoc.table(pc1_pws %>% arrange(desc(PC1)),
             digits = 5,
             keep.trailing.zeros = T,
             split.tables = Inf,
             caption = "PC1 positive loadings")
pandoc.table(pc1_neg %>% arrange(PC1),
             digits = 5,
             keep.trailing.zeros = T,
             split.tables = Inf,
             caption = "PC1 negative loadings")

