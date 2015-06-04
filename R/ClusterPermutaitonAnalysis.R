library(ggplot2)

## Parameters set to get the correct results.
setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10")
GOI = 9
GIC = 16
Cluster = 3
NumberOfPertubations = 10000

PertubationInformationFile = "gene_list2.permutations.txt"


## Gene info
distribution = read.table(PertubationInformationFile, sep = "\t", header =  FALSE, stringsAsFactors=FALSE) 
colnames(distribution) <- c("Name", "distance","Clusters", "GIC", "GOI")

GOIsAbove = length(which(distribution$GOI >= GOI))  
GOI_pvalue = GOIsAbove/NumberOfPertubations

GICsAbove = length(which(distribution$GIC >= GIC))  
GIC_pvalue = GICsAbove/NumberOfPertubations

ClustersAbove = length(which(distribution$Clusters >= Cluster))  
Cluster_pvalue = ClustersAbove/NumberOfPertubations

test = stack(distribution[,c("Clusters", "GIC", "GOI")])


plot <- ggplot(test, aes(x=values, colour=ind, group=ind))
plot + geom_density(fill=NA) +  scale_x_continuous(limits=c(0, 15))
geom_point(data=ClusterInfo, mapping=aes(y=Chr, x=Location , color = Source,size=ClusterInfo$total))
  
