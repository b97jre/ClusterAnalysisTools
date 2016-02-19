
loadGeneInfo(gffFile,arrayInfoFile,OrthoInfoFile)arrayInfoFile
GenesInfo = getGenes(gffFile)
GenesInfo = addSpottedGenesOnArrayInfo(arrayInfoFile,GenesInfo)
GenesInfo = addOrthoInfo(OrthoInfoFile, GenesInfo)

GenesInfo = removeGenesNotOnArray(GenesInfo)

plotDistrubution <- function(GenesInfo, ExpressionFile, maxDistances)



maxDistances= c(2000,5000,10000,20000,40000)

TrueDistributions <- data.frame() 
for(i in 1:length(ExpressionFile)){
  Genes = addExpressionInfo(ExpressionFile[i],GenesInfo)
  for(j in 1:length(maxDistances)){
    trueCluster = getClustersTrue(Genes = Genes,aboveCutoff = aboveCutoff, belowCutoff = belowCutoff, maxDistance = maxDistances[j],nrOfGOIsPerClusterCutoff = nrOfGenesInCluster, nrOfOrthologs = nrOfOrthologs)
    nrOfGOI = getNrOfGOIs(Genes = Genes,aboveCutoff = aboveCutoff, belowCutoff = belowCutoff)
    nrOfClusters =length(unique(trueCluster$Cluster))
    nrOfGOIinCluster = nrow(trueCluster)

    newRow = data.frame(Dataset=Dataset[i],maxDistance = maxDistances[j] , GOI=nrOfGOI,Clusters=nrOfClusters,ClusteredGOI=nrOfGOIinCluster)
    TrueDistributions = rbind(TrueDistributions,newRow)
  }
}

TrueDistributions$pValue = 1


BackgroundDistributions <- data.frame() 
nrOfTimes = 100000

for(i in 1:length(ExpressionFile)){
  Genes = addExpressionInfo(ExpressionFile[i],GenesInfo)
  nrOfGOI = getNrOfGOIs(Genes = Genes,aboveCutoff = aboveCutoff, belowCutoff = belowCutoff)
  for(j in 1:length(maxDistances)){
    ScrambleInfo = list(Genes = Genes,count = nrOfGOI,dist = maxDistances[j],nrOfGOIsPerClusterCutoff = nrOfGenesInCluster,nrOrt =  nrOfOrthologs)
    files <- lapply(1:nrOfTimes ,getClustersScrambledFromList, listOfVariables = ScrambleInfo)
    Distribution = data.frame(matrix(unlist(files), nrow=nrOfTimes, byrow=T))
    colnames(Distribution) <- c("NrOfClusters","ClusteredGOI") 
    Distribution$maxDistance = maxDistances[j] 
    Distribution$Dataset = Dataset[i]
    Distribution$GOI = nrOfGOI
    Distribution =Distribution[,c(4,3,5,1,2)]
    BackgroundDistributions = rbind(BackgroundDistributions,Distribution)
  }
}


for(i in 1:length(ExpressionFile)){
  for(j in 1:length(maxDistances)){
    Distribution = BackgroundDistributions[BackgroundDistributions$maxDistance == maxDistances[j] & BackgroundDistributions$Dataset == Dataset[i],  ]
    TrueDistribution = TrueDistributions[TrueDistributions$maxDistance == maxDistances[j] & TrueDistributions$Dataset == Dataset[i],  ]
    position = length(which(Distribution$NrOfClusters >=  TrueDistribution$Clusters)) +1
    pValue = position / length(Distribution$NrOfClusters)
    TrueDistributions$pValue[TrueDistributions$maxDistance == maxDistances[j] & TrueDistributions$Dataset == Dataset[i]] = pValue
    
  }
}

TrueDistributions$pValue2[TrueDistributions$pValue >= 0.001] = ">= 1e10-3" 
TrueDistributions$pValue2[TrueDistributions$pValue < 0.001 ] = "< 1e10-3" 
TrueDistributions$pValue2[TrueDistributions$pValue < 0.00001] = "< 1e10-5" 


p <- ggplot(BackgroundDistributions, aes(x = factor(maxDistance), y = NrOfClusters))
p + geom_boxplot(aes(fill = factor(Dataset)))+ 
  theme_bw()+
  geom_point(data = TrueDistributions, size = 3, aes(factor(maxDistance), Clusters,shape = factor(pValue2), colour = factor(Dataset),fill = factor(Dataset)))

ggsave("figure1b.pdf")



g <- ggplot(TrueDistributions, aes(factor(maxDistance), Clusters))
g +  geom_point( )
  )data = TrueDistributions,aes(factor(maxDistance), Clusters))

)
## Genes spotted on array based on blast analysis


