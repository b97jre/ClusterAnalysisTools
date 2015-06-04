library(ggplot2)
library(plyr)


setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/figure1")


trueDistribution = read.table("trueDistribution.txt")
permutedDistribution = read.table("BackogroundDistribution.txt")

colnames(trueDistribution) <- c("DataSet","Length","NrOfGOI","NrOfClusters","NrOfGenesInCluster", "NumberOfGenesOfInterestInClusters")
colnames(permutedDistribution) <- c("DataSet","Length","NrOfClusters","NrOfGenesInCluster", "NumberOfGenesOfInterestInClusters")

# removing 100000 clusters
trueDistribution = trueDistribution[trueDistribution$Length < 80000, ]
permutedDistribution = permutedDistribution[permutedDistribution$Length < 80000, ]

trueDistribution$sigLevel = 2

trueDistribution$DataSet2 = trueDistribution$DataSet

DSfactors = levels(factor(trueDistribution$DataSet))
LengthFactors = levels(factor(trueDistribution$Length))

PvalueMatrix = matrix(nrow  = length(DSfactors), ncol = length(LengthFactors))

for(i in 1:length(DSfactors)){
  for(j in 1:length(LengthFactors)){
    subDistribution =  permutedDistribution[permutedDistribution$DataSet == DSfactors[i] 
                                            & permutedDistribution$Length == LengthFactors[j],]
    nrOfPermutations = dim(subDistribution)[1]
    trueCluster = trueDistribution[trueDistribution$DataSet == DSfactors[i] 
                         & trueDistribution$Length == LengthFactors[j],]
    position = length(which(subDistribution$NrOfClusters< trueCluster$NrOfClusters))
    pvalue = 1/(nrOfPermutations+1)
    if(position !=   nrOfPermutations){
      pvalue = (nrOfPermutations-position)/(nrOfPermutations)
    }  
    PvalueMatrix[i,j] = pvalue
    if(floor(-log(pvalue, base = 10)) > 2){
      trueDistribution$sigLevel[trueDistribution$DataSet == DSfactors[i] 
                         & trueDistribution$Length == LengthFactors[j]] = floor(-log(pvalue,base = 10))
    }
  }
}



p <- ggplot(permutedDistribution, aes(x = factor(Length), y = NrOfClusters))
p + geom_boxplot(aes(fill = factor(DataSet)))+ 
  theme_bw()+
  geom_point(data = trueDistribution, size = 5, aes(factor(Length), NrOfClusters, shape = factor(sigLevel),colour = factor(DataSet),fill = factor(DataSet)))
ggsave("Figure1bc.pdf",useDingbats=FALSE)       





p <- ggplot(permutedDistribution, aes(x = Length, y = NrOfClusters, group = round_any(Length, 5000, floor)))
p + geom_boxplot()+ facet_grid(. ~ DataSet)+
  theme_bw()+
  geom_point(data = trueDistribution, size = 5, aes(Length, NrOfClusters, shape = factor(sigLevel),colour = factor(DataSet2)))
ggsave("Figure1bc.pdf")       

