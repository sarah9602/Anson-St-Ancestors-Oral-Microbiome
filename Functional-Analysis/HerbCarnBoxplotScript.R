## Script created by Sarah Johnson (2024) to visualize comparative output from HUMANn3

library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(rstatix)
FIL<-"CHS-Vel-WAMP"
setwd("PATH/TO/DATA/TABLES")
otherdir <- "/PATH/TO/FEATURE/MAPS"

# Load merged, renorm, regrouped, meat consumption-filtered table
data<-read.csv(sprintf("%s.HerbCarn.tsv",FIL),header=T,sep='\t',row.names=1)

# Load Meat Consumption mapping file
featuremap<-read.csv(sprintf("%s/MeatConsumption_ECnums.txt",otherdir),header=T,sep='\t')

# Load metadata
metadata<-read.csv(sprintf("%s.metadata.txt",FIL),header=T,sep='\t')
# Specify order of populations
ord <- c("Anson Street Ancestors","Wichita","Historic Radcliffe","Modern")
metadata$Population <- factor(metadata$Population, levels=ord)
# Get counts of each population
for (i in ord){
  print(i)
  print(nrow(metadata[metadata$Population == i,]))
}

## Filter dataframe to only include functional groups present in at least N amount of samples
## N can vary
# At least 75% of samples have greater than 2 cpm of functional group
data <- data[ rowSums(data > 2) >= .75*ncol(data),]
# Transpose data and merge with metadata file
data.t <- data %>% t() %>% 
  as.data.frame() %>%
  rownames_to_column(var="SampleName") %>%
  mutate(Population= metadata$Population[match(SampleName, metadata$Sample)])
# Mutate the data
dfm <- data.t %>% 
  gather(Feature,Value,-SampleName,-Population) %>%
  mutate(Diet=featuremap$Consumption_category[match(Feature,featuremap$Name)])
## Drop those that don't align to dietary EC numbers
dfm <-dfm[rowSums(is.na(dfm)) == 0,]

# Plot the output and color by population
b=ggplot(dfm,aes(x=Population,y=Value, fill=Population)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#DA9101ff","#56b4e9ff","#D81B60","#009e73ff")) +
  labs(x="Population",y="Copies per million") +
  theme_minimal() +
  scale_x_discrete(labels=function(x) str_wrap(x,width=10)) +
  ## Optional - provide points on plot
  # geom_jitter(color="black", size=0.4, alpha=0.9) +
  facet_wrap(~ Diet, scales="free_y")

b
# Write mutated table to file
write.table(dfm,sprintf("%s.HerbCarn-Boxplots.txt",FIL),row.names=FALSE,sep='\t')
## Initiate empty df to add Kruskal Wallis statistics to
kwdf <- data.frame(row.names=c('statistic','p.value','parameter','method'))

for (diet in unique(dfm$Diet)) {
  dfx<-dfm[dfm$Diet == diet,]
  # print(diet)
  tb <-(group_by(dfx,Population) %>%
    summarise(count=n(),
              mean=mean(Value,na.rm=T),
              sd=sd(Value,na.rm=T),
              median=median(Value,na.rm=T),
              IQR=IQR(Value,na.rm=T)
    ))
  write.table(tb, file=sprintf("%s.%s.summary.txt",diet,FIL),col.names=T,row.names=F,sep='\t')
  # Kruskal-Wallis statistic
  kw<-tidy(kruskal.test(Value ~ Population, data=dfx))
  kw <- t(as.data.frame(kw))
  colnames(kw)[1]  <- diet
  kwdf<-cbind(kwdf,kw)
  # Pairwise statistics
  pww <- tidy(pairwise.wilcox.test(dfx$Value,dfx$Population,p.adjust.method="BH",exact=F))
  write.table(pww,file=sprintf("%s.%s.Pairwise.txt",diet,FIL),col.names=T,row.names=F,sep='\t')
}
kwdf <- cbind(Measure = rownames(kwdf),kwdf)
# Write kruskal-wallis statistic to file
write.table(kwdf,file=sprintf("%s-MeatConsumption.KW.txt",FIL),col.names=T,row.names=F,sep='\t')

