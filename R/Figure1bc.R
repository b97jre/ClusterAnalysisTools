library(ggplot2)
library(plyr)
library(reshape)



setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/figure1")
setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/figure1/denseDist")


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
PositionMatrix = matrix(nrow  = length(DSfactors), ncol = length(LengthFactors))


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
    PositionMatrix[i,j] = nrOfPermutations- position
    if(floor(-log(pvalue, base = 10)) > 2){
      trueDistribution$sigLevel[trueDistribution$DataSet == DSfactors[i] 
                         & trueDistribution$Length == LengthFactors[j]] = floor(-log(pvalue,base = 10))
    }
  }
}

trueDistribution$pValue = ">= 1e10-3" 
trueDistribution$pValue[trueDistribution$sigLevel>2 & trueDistribution$sigLevel<5] = "< 1e10-3" 
trueDistribution$pValue[trueDistribution$sigLevel>4 ] = "< 1e10-5" 

colnames(PvalueMatrix) = levels(factor(permutedDistribution$Length))
rownames(PvalueMatrix) = levels(factor(permutedDistribution$DataSet))
 

melt(PvalueMatrix)

PM = melt(PvalueMatrix)
colnames(PM ) <- c("Dataset", "Length", "pValue")

points <- ggplot(PM, aes(factor(Length), log(pValue,base = 10)))
points + geom_point(aes(colour = factor(Dataset), shape = factor(Dataset)))  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

lines <- ggplot(PM, aes(x=Length, y=log(pValue,base = 10), group=Dataset))
lines + geom_line(aes(colour = factor(Dataset)))



head(permutedDistribution)
p <- ggplot(permutedDistribution, aes(x = factor(Length), y = NrOfClusters))
p + geom_boxplot(aes(fill = factor(DataSet)))+ 
  theme_bw()+
  geom_point(data = trueDistribution, 
             size = 3, aes(factor(Length), NrOfClusters, 
                           shape = factor(pValue),colour = factor(DataSet),fill = factor(DataSet)))

ggsave("Figure1b.pdf",useDingbats=FALSE)       





p <- ggplot(permutedDistribution, aes(x = Length, y = NrOfClusters, group = round_any(Length, 5000, floor)))
p + geom_boxplot()+ facet_grid(. ~ DataSet)+
  theme_bw()+
  geom_point(data = trueDistribution, size = 5, aes(Length, NrOfClusters, shape = factor(sigLevel),colour = factor(DataSet2)))
ggsave("Figure1bc.pdf")       

