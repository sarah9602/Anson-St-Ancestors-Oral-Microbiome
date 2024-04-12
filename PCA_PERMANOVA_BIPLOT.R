## Script to (1) create PCA plots based on Metaphalan4 relative abundance data;
## (2) To color PCA plots based on dental pathology variables;
## (3) Run PERMANOVA association testing to determine significance;
## (4) Create a PCA biplot
## Script Created by R. Fleskes (2024) based on scripts provided from Velsko et al. 2022
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
library(gridExtra)
library(grid)


## Read in species-level table from MetaPhlAn3 and Sample metadata file
setwd("~/Dropbox (Dartmouth College)/14_CODING/R-Studio Projects/DentalCalculus/DentalCalculus/Github")
OTUIN<-otu_table(read.table("CHS_metaphalan.txt",header=TRUE,row.names=1),taxa_are_rows = TRUE)
SAMPDATA<-sample_data(read.table("PCA_Metadata.txt",header=TRUE,row.names=1,sep="\t"))

###############################################################
## CREATE PCA PLOTS BASED ON METAPHALAN4 RELATIVE ABUNDANCES ##
## AND COLOR THEM BY DENTAL PATHOLOGY METADATA               ##
###############################################################

## Create PCA plot
psobject<-phyloseq(OTUIN,SAMPDATA)
PRobject<-t(otu_table(psobject))+1
pca.result<-pca(PRobject,logratio='CLR')
coords_pca<-as.data.frame(pca.result$variates)
loadings_pca<-as.data.frame(pca.result$loadings)
FinalTable<-merge(coords_pca,as.data.frame(SAMPDATA),by="row.names")

## Format Pathologies Variables
FinalTable$ToothNumberPathologies <- factor(FinalTable$ToothNumberPathologies)
FinalTable$ArcadeNumberPathologies <- factor(FinalTable$ArcadeNumberPathologies)
FinalTable$ArcadeOralHealthNumberPathologies <- factor(FinalTable$ArcadeOralHealthNumberPathologies)
FinalTable$BirthPlace <- factor(FinalTable$BirthPlace)
FinalTable$Abscess <- factor(FinalTable$Abscess)
FinalTable$CariesAll <- factor(FinalTable$CariesAll)
FinalTable$CariesTooth <- factor(FinalTable$CariesTooth)

## Create base PCA
plot1<-ggplot(data=FinalTable,aes(X.PC1,X.PC2))+geom_point(alpha=0.0001)+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(-6,6)+ylim(-6,6)

plot1<-plot1+
  geom_point(data=subset(FinalTable,Population %in% c("AnsonStreet_SouthCarolina")), shape=17, colour="black", size=3) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(-6,6)+ylim(-6,6)+
  labs(x = "PC1", y = "PC2")
plot1


## Create color codes for PCAs
tooth_code_colors_2 <- c("Non-NorthAmerica"="#440154FF","NorthAmerica"="orange")
tooth_code_colors_2A <- c("NotPresent"="#440154FF","Present"="orange")
tooth_code_colors_3 <- c("0"="#440154FF","1"="#21908CFF", "2"="orange")
tooth_code_colors_6 <- c("0"="#440154FF","1"="#46337EFF", "2"="#365C8DFF", "3"="#22A884FF", "4"="#7AD151FF", "5"="orange")
tooth_code_colors_8 <- c("0"="#440154FF","1"="#46337EFF", "2"="#365C8DFF", "3"="#277F8EFF", "4"="#21908CFF", "5"="#27AD81FF", "6"="#5DC863FF", "7"="#AADC32FF", "8"="orange")

## Make PCA Plots colored by Variable
ToothVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = ToothNumberPathologies), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_3, name="")+
  labs(title = "Number of Pathologies on Teeth/Tooth with Dental Calculus")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.title.position = "plot",        
        text = element_text(size = 12))

ArcadeVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = ArcadeNumberPathologies), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_8, name = "")+
  labs(title = "Number of Pathologies in Dental Arcade") +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.title.position = "plot",
        text = element_text(size = 12))

OralHealthVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = ArcadeOralHealthNumberPathologies), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_6, name = "")+
  labs(title = "Number of Pathologies Related to Oral Health in Dental Arcade")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),        
        plot.title.position = "plot", 
        text = element_text(size = 12))

ResidencyVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = BirthPlace), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_2, name = "")+
  labs(title = "Residency During Childhood Based on Enamel Strontium Ratio")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.title.position = "plot", 
        text = element_text(size = 12))

AbscessVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = Abscess), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_2A, name = "")+
  labs(title = "Presence or Abscence of Abscesses in Dental Arcade")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.title.position = "plot", 
        text = element_text(size = 12))

CariesAllVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = CariesAll), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_2A, name = "")+
  labs(title = "Presence or Abscence of Caries in Dental Arcade")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.title.position = "plot", 
        text = element_text(size = 12))

CariesToothVariablePlot <- plot1 +
  geom_point(
    data = subset(FinalTable, Population %in% c("AnsonStreet_SouthCarolina")),
    aes(color = CariesTooth), shape = 17, size = 2) +
  scale_color_manual(values = tooth_code_colors_2A, name = "")+
  labs(title = "Presence or Abscence of Caries on Sampled Tooth")+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.title.position = "plot", 
        text = element_text(size = 12))

## Arrange generated Plots into Main Figures
SI_PCAS <- grid.arrange(ResidencyVariablePlot, ArcadeVariablePlot, CariesToothVariablePlot, 
                        ToothVariablePlot, ncol=2)
SI_PCAS_2 <- grid.arrange(AbscessVariablePlot,CariesAllVariablePlot, OralHealthVariablePlot,ncol=2)
ggsave("SI_PCAS_layout.png", plot = SI_PCAS, width = 10, height = 10, units = c("in"), dpi = 400)
ggsave("SI_PCAS_2_layout.png", plot = SI_PCAS_2, width = 10, height = 10, units = c("in"), dpi = 400)


##############################################
## PERMANOVA TESTING USING VEGAN ON THE     ##
## PCA AND DENTAL PATHOLOGY METADATARESULTS ##
##############################################

## Create Euclidean distance matrix using Vegan Package
diss_matrix <- vegdist(pca.result$x, method = "euclidean", na.rm = F)

## Create Variables based on Metadata
ArcadeOralHealth_variable <- FinalTable$ArcadeOralHealthNumberPathologies
Arcade_variable <- FinalTable$ArcadeNumberPathologies
Tooth_variable <- FinalTable$ToothNumberPathologies
Residency_variable <- FinalTable$BirthPlace
Abscess_variable <- FinalTable$Abscess
CariesAll_variable <- FinalTable$CariesAll
CariesTooth_variable <- FinalTable$CariesTooth

## Run PERMANOVA Association Testing using Adonis2
perm_result_ArcadeOralHealth <- adonis2(diss_matrix ~ ArcadeOralHealthNumberPathologies, FinalTable, permutations = 999, by = "margin")
perm_result_Arcade <- adonis2(diss_matrix ~ ArcadeNumberPathologies, FinalTable, permutations = 999, by = "margin")
perm_result_Tooth <- adonis2(diss_matrix ~ ToothNumberPathologies, FinalTable, permutations = 999, by = "margin")
perm_result_Residency <- adonis2(diss_matrix ~ BirthPlace, FinalTable, permutations = 999, by = "margin")
perm_result_Abscess <- adonis2(diss_matrix ~ Abscess, FinalTable, permutations = 999, by = "margin")
perm_result_CariesAll <- adonis2(diss_matrix ~ CariesAll, FinalTable, permutations = 999, by = "margin")
perm_result_CariesTooth <- adonis2(diss_matrix ~ CariesTooth, FinalTable, permutations = 999, by = "margin")

## Filter and Write PERMANOVA results
all_permanova_results <- bind_rows(perm_result_ArcadeOralHealth, perm_result_Arcade, perm_result_Tooth, perm_result_Residency, perm_result_Abscess, perm_result_CariesAll, perm_result_CariesTooth) %>% rownames_to_column("ID")
  
all_permanova_results_filtered <- all_permanova_results %>%
  as_tibble() %>%
  drop_na(`F`) %>%
  arrange(desc(`F`)) %>%
  arrange(desc(R2))

write.csv(all_permanova_results_filtered, file = "all_permanova_testing.csv", quote =  F)


#######################################
## BIPLOT ANALYSIS USING PCA RESULTS ##
#######################################

# Extract explained variance
pca_expl_vars <- paste0(round(pca.result$variance * 100, 2), "%")

# Extract loadings (eigenvectors) from PC analysis
pca_loadings_df <- as.data.frame(pca.result$loadings$X) %>%
  rownames_to_column(var = "Species")

# Identify the 10 species that have highest positive and negative loadings in PC1 and PC2
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

## Create Color Variables for BiPlot
color_palette_13 <- viridis(13, begin = 0, end = 1, direction = -1)
sample_code_colors_13 <- c("Banza (CHS01)"="orangered", "Lima (CHS03)"="orange", "Juba (CHS14)"="#8FD744FF", "Lisa (CHS22)"="#21908CFF", "Risu (CHS26)"="#287C8EFF", "Amina (CHS27)"="#31688EFF", "Pita (CHS29)"="#443A83FF", "Tima (CHS31)"="#481F70FF", "Isi (CHS36)"="#440154FF")
sample_code_shape_13 <- c("Banza (CHS01)"=15, "Lima (CHS03)"=16, "Juba (CHS14)"=17,"Lisa (CHS22)"=15, "Risu (CHS26)"=16, "Amina (CHS27)"=17, "Pita (CHS29)"=15, "Tima (CHS31)"=16, "Isi (CHS36)"=17)

## Update Names
FinalTable[1,1] <-"Banza (CHS01)"
FinalTable[2,1] <-"Lima (CHS03)"
FinalTable[3,1] <-"Juba (CHS14)"
FinalTable[4,1] <-"Lisa (CHS22)"
FinalTable[5,1] <-"Risu (CHS26)"
FinalTable[6,1] <-"Amina (CHS27)"
FinalTable[7,1] <-"Pita (CHS29)"
FinalTable[8,1] <-"Tima (CHS31)"
FinalTable[9,1] <-"Isi (CHS36)"

## Create BiPlot
plot_bi <- FinalTable %>%
  mutate(PC1 = X.PC1 / corr_lam[1],
         pc2_residuals = X.PC2 / corr_lam[2]) %>%
  ggplot(.,
         aes(x = PC1, y = pc2_residuals)) +
  geom_point(size = 3.5, aes(shape = factor(FinalTable$Row.names), color = factor(FinalTable$Row.names))) +
  scale_shape_manual(values = sample_code_shape_13, name = "") +
  scale_color_manual(values = sample_code_colors_13, name = "") +
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
  theme_light() +
  theme(text = element_text(size = 14), legend.position = "right")
plot_bi
ggsave("BiPlot_metaphalan4.png", plot = plot_bi, width = 12, height =4.5, dpi = 400)


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

