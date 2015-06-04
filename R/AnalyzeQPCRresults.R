library(ggplot2)

setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/QPCR_values")

# Static values
replicates = 3
ClusterRowNames =  c("Ma st1-9","Mb st1-9","Mc st1-9","Ma st10","Mb st10","Mc st10","Da st1-9","Db st1-9","Dc st1-9","Da st10","Db st10","Dc st10","Dd st10","Ma st11","Mb st11","Mc st11","Ma st12","Mb st12","Mc st12","Da-2 st11","Db-2 st11","Dc st11","Dd st11","Da st12","Da-2 st12","Db st12","Db-2 st12","Dc st12")


Cluster13 = read.table("Cluster_13.tab.txt", sep = "\t", header =  TRUE,stringsAsFactors=FALSE)
#removing last lines BECAUSE THEY ARE EMPTY VALUES
Cluster13 = Cluster13 [1:(dim(Cluster13)[1]-9),]
Cluster13genes = colnames(Cluster13)
Cluster13Mean  = getMeanValues(Cluster13,replicates,ClusterRowNames )
Cluster13_Mean_Relative  = as.data.frame(t(t(Cluster13Mean)/apply(Cluster13Mean,2,max)))
Cluster13_Mean_Relative_MetaInfo = addMetaInfo(Cluster13_Mean_Relative,ClusterRowNames)
Cluster13_Mean_Relative_MetaInfo_ggplot = GGplotInfo(Cluster13_Mean_Relative,Cluster13_Mean_Relative_MetaInfo)
Cluster13_Final = Cluster13_Mean_Relative_MetaInfo_ggplot
Cluster13_Final$Cluster  = "Cluster 15"

Cluster3 = read.table("Cluster_3.tab.txt", sep = "\t", header =  TRUE,stringsAsFactors=FALSE)
Cluster3genes = colnames(Cluster3)
Cluster3Mean  = getMeanValues(Cluster3,replicates,ClusterRowNames )
Cluster3_Mean_Relative  = as.data.frame(t(t(Cluster3Mean)/apply(Cluster3Mean,2,max)))
Cluster3_Mean_Relative_MetaInfo = addMetaInfo(Cluster3_Mean_Relative,ClusterRowNames)
Cluster3_Mean_Relative_MetaInfo_ggplot = GGplotInfo(Cluster3_Mean_Relative,Cluster3_Mean_Relative_MetaInfo)
Cluster3_Final = Cluster3_Mean_Relative_MetaInfo_ggplot
Cluster3_Final$Cluster  = "Cluster 3"


BothClusters = rbind(Cluster3_Final,Cluster13_Final)
BothClusters$Cluster = factor(BothClusters$Cluster, levels = c("Cluster 3","Cluster 15"))




MS1pvalues = matrix(nrow = 2, ncol = 2)
colnames(MS1pvalues) <- c("DEX", "Mock")
rownames(MS1pvalues) <- c("Cluster 3","Cluster 13")

MS1pvalues[1,1] = wilcox.results.MS1(Cluster =Cluster3_Final,Treatment = "DEX", TimePoint = c("st1-9","st11","st12")) 
MS1pvalues[2,1] = wilcox.results.MS1(Cluster =Cluster13_Final,Treatment = "DEX", TimePoint = c("st1-9","st11","st12")) 
MS1pvalues[1,2] = wilcox.results.MS1(Cluster =Cluster3_Final,Treatment = "Mock", TimePoint = c("st1-9","st11","st12")) 
MS1pvalues[2,2] = wilcox.results.MS1(Cluster =Cluster13_Final,Treatment = "Mock", TimePoint = c("st1-9","st11","st12")) 

write.table(x = MS1pvalues, file = "Wilcox.test.results.compare.st1-10_to_st11-12",quote = FALSE)



DEXpvalues = matrix(nrow = 2, ncol = 6)
colnames(DEXpvalues) <- c("st1-9", "st10","st11","st12","st1-9","st11-12")
rownames(DEXpvalues) <- c("Cluster 3","Cluster 13")

DEXpvalues[1,1] = wilcox.results.Treatment(Cluster = Cluster3_Final,TimePoint = "st1-9")
DEXpvalues[1,2] = wilcox.results.Treatment(Cluster = Cluster3_Final,TimePoint = "st10")
DEXpvalues[1,3] = wilcox.results.Treatment(Cluster = Cluster3_Final,TimePoint = "st11")
DEXpvalues[1,4] = wilcox.results.Treatment(Cluster = Cluster3_Final,TimePoint = "st12")
DEXpvalues[1,5] = wilcox.results.Treatment(Cluster = Cluster3_Final,MS1 = "STAGE 1-9")
DEXpvalues[1,6] = wilcox.results.Treatment(Cluster = Cluster3_Final,MS1 = "STAGE 11-12")
DEXpvalues[2,1] = wilcox.results.Treatment(Cluster = Cluster13_Final,TimePoint = "st1-9")
DEXpvalues[2,2] = wilcox.results.Treatment(Cluster = Cluster13_Final,TimePoint = "st10")
DEXpvalues[2,3] = wilcox.results.Treatment(Cluster = Cluster13_Final,TimePoint = "st11")
DEXpvalues[2,4] = wilcox.results.Treatment(Cluster = Cluster13_Final,TimePoint = "st12")
DEXpvalues[2,5] = wilcox.results.Treatment(Cluster = Cluster13_Final,MS1 = "STAGE 1-9")
DEXpvalues[2,6] = wilcox.results.Treatment(Cluster = Cluster13_Final,MS1 = "STAGE 11-12")

write.table(x = DEXpvalues, file = "Wilcox.test.results.compare.Mock_to_DEX",quote = FALSE)

##cannot compute exact p-value with ties




head(Cluster3_Final)




p = ggplot(BothClusters, aes(x = factor(Treatment), fill = factor(MS1), y = values),order=factor(Treatment)) +
  geom_boxplot()

p + facet_grid(. ~ Cluster)+
  ggtitle("Overall relative expression before, during and after stage 10") + 
  theme(plot.title = element_text( face="bold"),legend.title=element_blank())+
  ylab("Relative expression")+
  xlab("Adding DEX (activate MS1) or mock (does not activate MS1)")+ theme(legend.title=element_blank())
  ggsave(filename = "Both_Cluster_Overall_expression_pattern.pdf")




p = ggplot(Cluster3_Final, aes(x = factor(Treatment), fill = factor(MS1), y = values),order=factor(Treatment)) +
geom_boxplot()
  
p + facet_grid(. ~ ind,margins=TRUE)+
  ggtitle("Expression before (st1-10) and after (st11-12) MS1 onset") + 
  theme(plot.title = element_text( face="bold"),legend.title=element_blank())+
  ylab("Relative expression")+
  xlab("Adding DEX (activate MS1) or mock (does not activate MS1)")+ theme(legend.title=element_blank())+
  ggsave(filename = "cluster3_expression_MS1_boxplot.pdf")

p = ggplot(Cluster3_Final, aes(x = factor(Treatment), fill = factor(TimePoint), y = values),order=factor(Treatment)) +
  geom_boxplot()

 p + facet_grid(. ~ ind,margins=TRUE)+
  ggtitle("Expression of genes at different flower development stages") + 
  theme(plot.title = element_text( face="bold"),legend.title=element_blank())+
  ylab("Relative expression")+
  xlab("Adding DEX (activate MS1) or Mock (does not activate MS1)")+ theme(legend.title=element_blank())
ggsave(filename = "cluster3_expression_ST1_to_13_boxplot.pdf")


p = ggplot(Cluster13_Final, aes(x = factor(Treatment), fill = factor(MS1), y = values),order=factor(Treatment)) +
  geom_boxplot()

p + facet_grid(. ~ ind,margins=TRUE)+
  ggtitle("Expression before (st1-10) and after (st11-12) MS1 onset") + 
  theme(plot.title = element_text( face="bold"),legend.title=element_blank())+
  ylab("Relative expression")+
  xlab("Adding DEX (activate MS1) or mock (does not activate MS1)")+ theme(legend.title=element_blank())
ggsave(filename = "cluster13_expression_MS1_boxplot.pdf")

p = ggplot(Cluster13_Final, aes(x = factor(Treatment), fill = factor(TimePoint), y = values),order=factor(Treatment)) +
  geom_boxplot()

p + facet_grid(. ~ ind,margins=TRUE)+
  ggtitle("Expression of genes at different flower development stages") + 
  theme(plot.title = element_text( face="bold"),legend.title=element_blank())+
  ylab("Relative expression")+
  xlab("Adding DEX (activate MS1) or Mock (does not activate MS1)")+ theme(legend.title=element_blank())
ggsave(filename = "cluster13_expression_ST1_to_13_boxplot.pdf")





Cluster13Dex = TEST[TEST$TimePoint == "st1-9", ]
wilcox.test(compare1 ~ Treatment, data=Cluster13Dex)

wilcox.test.conditions(gplot)



Cluster13Dex = TEST[TEST$TimePoint == "st10", ]
wilcox.test(relativeValue ~ Treatment, data=Cluster13Dex)

Cluster13Dex = TEST[TEST$TimePoint == "st11", ]
wilcox.test(relativeValue ~ Treatment, data=Cluster13Dex)

Cluster13Dex = TEST[TEST$TimePoint == "st12", ]
wilcox.test(relativeValue ~ Treatment, data=Cluster13Dex)

Cluster13After =Cluster13DexOrdered[ Cluster13DexOrdered$MS1 =="AFTER",]
Cluster13Afterlapply(Cluster13After,Cluster13DexOrdered[ Cluster13DexOrdered$MS1 =="AFTER",])




wilcox.results.MS1 <- function(Cluster,Treatment = "all",  replicate= "all", TimePoint = c("all")){
  testCluster = Cluster
  if(Treatment != "all"){
    testCluster = testCluster[testCluster$Treatment == Treatment , ]
  }
  if(replicate != "all"){
    testCluster = testCluster[testCluster$replicate == replicate , ]
  }
  if(TimePoint[[1]] != "all"){
    testCluster = testCluster[which(testCluster$TimePoint %in% TimePoint) , ]
  }
  wctest = wilcox.test(values ~ MS1, data=testCluster)
  return (wctest$p.value)
}

wilcox.results.Treatment <- function(Cluster, MS1 = "all",  replicate= "all", TimePoint = "all"){
  testCluster = Cluster
  if(replicate != "all"){
    testCluster = testCluster[testCluster$replicate == replicate , ]
  }
  if(TimePoint != "all"){
    testCluster = testCluster[testCluster$TimePoint == TimePoint , ]
  }
  if(MS1 != "all"){
    testCluster = testCluster[testCluster$MS1 == MS1 , ]
  }
  wctest =  wilcox.test(values ~ Treatment, data=testCluster)
  return (wctest$p.value)
}





addMetaInfo  <- function(cluster,rowNames,NoMS1expression = c("st1-9"), MS1induction = c("st10")){
  foo <- data.frame(do.call('rbind', strsplit(as.character(rowNames),' ',fixed=TRUE)))
  cluster$Treatment = "DEX"
  cluster$Treatment[grep("M",foo$X1)] = "Mock"
  cluster$Treatment = factor(cluster$Treatment,levels = c("Mock","DEX"))
  cluster$replicate = substring(foo$X1,2)
  cluster$TimePoint = foo$X2
  cluster$MS1 = "STAGE 11-12"
  for(i in 1:length(NoMS1expression)){
    cluster$MS1[cluster$TimePoint == NoMS1expression[i]] = "STAGE 1-9"
  }
  for(i in 1:length(MS1induction)){
    cluster$MS1[cluster$TimePoint == MS1induction[i]] = "STAGE 10"
  }
  
  cluster$MS1 = factor(cluster$MS1, levels = c("STAGE 1-9","STAGE 10","STAGE 11-12"))
  return(cluster)
}
  


GGplotInfo <- function(onlyGenes,metaInfo){
  stackedInfo = stack(onlyGenes)
  stackedInfo$Treatment = metaInfo$Treatment 
  stackedInfo$replicate = metaInfo$replicate 
  stackedInfo$TimePoint =metaInfo$TimePoint 
  stackedInfo$MS1 = metaInfo$MS1 
  
  
  return(stackedInfo)
  
}







getMeanValues <- function(Cluster, replicates,ClusterRowNames ){
  ClusterDistribution = matrix(nrow = dim (Cluster)[1]/replicates,ncol = dim (Cluster)[2])
  colnames(ClusterDistribution) <- colnames(Cluster) 
  rownames(ClusterDistribution) <-ClusterRowNames
  
  for(col in 1:dim (Cluster)[2]){
    for(i in seq(from = 1,to = dim (Cluster)[1],by=replicates)){
      ClusterDistribution[(i-1)/replicates+1,col] = getMean(Cluster,col,i,replicates) 
    }
  }
  ClusterMean = as.data.frame(ClusterDistribution)
  
  return(ClusterMean)
} 



getFoldChange <- function(cluster,group1 = c("st1-9","st10"), group2= c("st11","st12") , genes){
  group1Cluster = cluster[cluster$TimePoint %in% group1,] 
  group2Cluster = cluster[cluster$TimePoint %in% group2,] 
  group1ClusterMock = group1Cluster[group1Cluster$Treatment %in% "Mock",] 
  group2ClusterMock = group2Cluster[group2Cluster$Treatment %in% "Mock",] 
  group1ClusterDEX = group1Cluster[group1Cluster$Treatment %in% "DEX",] 
  group2ClusterDEX = group2Cluster[group2Cluster$Treatment %in% "DEX",] 
  genes = Cluster3genes
  
  group1GenesMock = as.matrix(group1ClusterMock[,colnames(group1ClusterMock) %in% genes ])
  group1MeansMock = colMeans(group1GenesMock)
  group1GenesDEX = as.matrix(group1ClusterDEX[,colnames(group1ClusterDEX) %in% genes ])
  group1MeansDEX = colMeans(group1GenesDEX)
  
  group2GenesMock = as.matrix(group2ClusterMock[,colnames(group2ClusterMock) %in% genes ])
  group2MeansMock = colMeans(group2GenesMock)
  
  group2GenesDEX = as.matrix(group2ClusterDEX[,colnames(group2ClusterDEX) %in% genes ])
  group2MeansDEX = colMeans(group2GenesDEX)
  

  group1FoldChange = group1MeansDEX/group1MeansMock
  changes =  as.data.frame(group1FoldChange)
  changes$group2FoldChange = group2MeansDEX/group2MeansMock
  changes$DexFoldChange =  group2MeansDEX / group1MeansDEX
  changes$MockFoldChange = group2MeansMock / group1MeansMock

  
  
  group1MedianMock = apply(group1GenesMock, 2, median)
  group1MedianDEX = apply(group1GenesDEX, 2, median)
  group2MedianMock = apply(group2GenesMock, 2, median)
  group2MedianDEX = apply(group2GenesDEX, 2, median)

  
  group1FoldChange = group1MedianDEX/group1MedianMock
  MedianChanges =  as.data.frame(group1FoldChange)
  MedianChanges$group2FoldChange = group2MedianDEX/group2MedianMock
  MedianChanges$DexFoldChange =  group2MedianDEX / group1MedianDEX
  MedianChanges$MockFoldChange = group2MedianMock / group1MedianMock
  
  
}



getMean <- function(dataFrame, column, start,replicates){
  mean = 0
  count = 0
  for(i in 0:(replicates-1)){
    if(dataFrame[start+i,column] != "MV"){
      count= count+1
      print(column)
      print(start+i)
      print(dataFrame[start+i,column])
      print(mean)
      
      mean = mean + as.numeric(dataFrame[start+i,column])
      print(mean)
    }
  }
  if(mean >0 ){
    mean = mean/count
  }
  return(mean)
}


  