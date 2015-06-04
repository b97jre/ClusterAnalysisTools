library(ggplot2)

setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/QPCR_values")


CondensationValues = read.table("Condensation measurements_JR.txt", sep = "\t", header =  TRUE)

CondensationValuesTapetal =CondensationValues[CondensationValues$CellType =="Tapetal",] 

p1 = ggplot(CondensationValuesTapetal, aes(x = factor(Stage), fill = factor(Treatment), y = nm),order=factor(Sample)) 
p1 + geom_violin(adjust=1) + facet_grid(. ~ Cluster)+
  ggtitle("Condensation of cluster") + 
  theme(plot.title = element_text( face="bold"),legend.title=element_blank())+
  ylab("Length(nm)")+
  xlab("Floral stages")+ theme(legend.title=element_blank())
ggsave(filename = "Floral_stages_distribution_Tapetal.pdf")



CondensationValuesNonTapetal =CondensationValues[CondensationValues$CellType !="Tapetal",] 

p1 = ggplot(CondensationValuesTapetal, aes(x = factor(Stage), fill = factor(Treatment), y = nm),order=factor(Sample)) 




factor()



p+ geom_violin(adjust=2)+ facet_grid(CellType ~ Cluster)


p+ geom_violin()+ facet_grid(CellType ~ Cluster)

+ stat_bin()

DEX_KS_pvalues = matrix(nrow = 2, ncol = 4)
colnames(DEX_KS_pvalues) <- c("DEX_10vs11","MOCK_10vs11","10_DEXvsMOCK","11_DEXvsMOCK" )
rownames(DEX_KS_pvalues) <- c("Cluster 3","Cluster 13")

DEX_KS_pvalues[1,1] = ks.test.Stage(info =CondensationValues, Treatment = "DEX",Celltype = "Tapetal", Cluster = 3 )
DEX_KS_pvalues[1,2] = ks.test.Stage(info =CondensationValues, Treatment = "Mock",Celltype = "Tapetal", Cluster = 3 )
DEX_KS_pvalues[1,3] = ks.test.Treatment(info =CondensationValues, Stage = 10,Celltype = "Tapetal", Cluster = 3 )
DEX_KS_pvalues[1,4] = ks.test.Treatment(info =CondensationValues, Stage = 11,Celltype = "Tapetal", Cluster = 3 )
DEX_KS_pvalues[2,1] = ks.test.Stage(info =CondensationValues, Treatment = "DEX",Celltype = "Tapetal", Cluster = 15 )
DEX_KS_pvalues[2,2] = ks.test.Stage(info =CondensationValues, Treatment = "Mock",Celltype = "Tapetal", Cluster = 15 )
DEX_KS_pvalues[2,3] = ks.test.Treatment(info =CondensationValues, Stage = 10,Celltype = "Tapetal", Cluster = 15 )
DEX_KS_pvalues[2,4] = ks.test.Treatment(info =CondensationValues, Stage = 11,Celltype = "Tapetal", Cluster = 15 )





Celltype= "Tapetal"
Stage = "10"
Cluster = "3"


ks.test.Treatment <- function(info, Stage = "all",  Celltype= "all", Cluster = "all"){
  testCluster = info
  if(Celltype != "all"){
    testCluster = testCluster[testCluster$CellType == Celltype , ]
  }
  if(Stage != "all"){
    testCluster = testCluster[testCluster$Stage == Stage , ]
  }
  if(Cluster != "all"){
    testCluster = testCluster[testCluster$Cluster == Cluster , ]
  }
  x = testCluster$nm[testCluster$Treatment == "Mock"]
  y = testCluster$nm[testCluster$Treatment == "DEX"]
  ks_results = ks.test(x, y)
  return (ks_results$p.value)
}


Treatment = "DEX"
ks.test.Stage <- function(info, Treatment = "all",  Celltype= "all", Cluster = "all"){
  testCluster = info
  if(Celltype != "all"){
    testCluster = testCluster[testCluster$CellType == Celltype , ]
  }
  if(Treatment != "all"){
    testCluster = testCluster[testCluster$Treatment == Treatment , ]
  }
  if(Cluster != "all"){
    testCluster = testCluster[testCluster$Cluster == Cluster , ]
  }
  
  x = testCluster$nm[testCluster$Stage == "10"]
  y = testCluster$nm[testCluster$Stage == "11"]
  ks_results = ks.test(x, y)
  return (ks_results$p.value)
}

ggplot(testCluster, aes(nm, colour = as.factor(Stage),linetype = as.factor(Treatment)))+ geom_density(adjust=3)
p
testCluster$Stage
factor(CondensationValues$nm)
39*39+56*56

