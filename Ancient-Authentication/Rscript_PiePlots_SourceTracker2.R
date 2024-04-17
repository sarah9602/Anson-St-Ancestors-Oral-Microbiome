### Barplots for sourcetracker2
setwd("/PATH/TO/Metaphlan/sourcetracker2.20000.s")

tab=read.csv(file="mixing_proportions.txt",sep="\t",header=TRUE)

library(reshape2)
tab1=setNames(data.frame(t(tab[,-1])),tab[,1])
tab1<-tab1[order(tab1$ModernDentalCalculus),]
tab2= tab1 %>% rownames_to_column(var="SampleID") 
tab2=melt(tab2,id.vars="SampleID")
ord <- c("ModernDentalCalculus","ModernSubgingivalPlaque", "ModernSupragingivalPlaque","ruralGut","urbanGut","skin","soil","Unknown","BLANK")
tab2$variable<-factor(tab2$variable,levels=ord)
tab2$SampleID=as.character(tab2$SampleID)
tab2$SampleID=factor(tab2$SampleID,levels=unique(tab2$SampleID))

library(ggplot2)
library(RColorBrewer)


cbPalette <- c("#A4D3A6", "#518244","#4E687C","#223FEC","#C15FDC","#5F3634","#70A3E9","#656565")
pie_charts <- ggplot(tab2, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~SampleID, nrow = 1) +  # One row for each sample
  scale_fill_manual(values=cbPalette) +
  theme_void() +  # Removes unnecessary elements
  theme(legend.position = "right")  # Adjust legend position if needed
pie_charts

plot.type <- 'pie'
envs <- as.factor(tab2$SampleID)
pdf("all_pie_charts.pdf", width = 10, height = 5)
# Loop through each environment and plot the pie charts
for (env in unique(envs)) {
  # Subset the data for the current environment
  subset_data <- tab2[tab2$SampleID == env, ]

  # Create the pie chart
  pie_chart <- ggplot(subset_data, aes(x = "", y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = env) +  # Set the title to the current environment
    scale_fill_manual(values=cbPalette) +
    theme_void() +  # Removes unnecessary elements
    theme(legend.position = "right")  # Adjust legend position if needed

  print(pie_chart)
}
dev.off()


