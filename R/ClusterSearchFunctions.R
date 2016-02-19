## Gene info
library(ggplot2)
library(plyr)

getGenes <- function(gffFile){
  gff3OutPut = read.table(gffFile, sep = "\t", header =  FALSE, skip =1, stringsAsFactors=FALSE ) 
  gff3Genes = gff3OutPut[which(gff3OutPut[[3]] == "gene"), ]
  
  
  ll <- unlist(strsplit(gff3Genes$V9, ";", fixed=TRUE))
  ll = ll[seq(from = 1, to = length(ll)-1, by = 2)]
  ll = unlist(strsplit(ll, "=", fixed=TRUE))
  ll = ll[seq(from = 2, to = length(ll), by = 2)]
  gff3Genes$ID = ll
  
  Genes = gff3Genes[,c(1,10,4,5,7)]
  colnames(Genes) <- c("chr","ID","start","stop","strand")
  return(Genes)
}

addSpottedGenesOnArrayInfo <- function(ArrayInfoFile, Genes){
  Genes$array = 0
  arraysMatches = read.table( sep = "\t", file = "Operon_set1_4.0_oligos.original_ACC_match.Athaliana_167_TAIR9.trimmed.matchGff3.blast")
  Genes$array[Genes$ID %in% arraysMatches$V16] = 1
  return (Genes)
}

loadGeneInfo <- function(gffFile,arrayInfoFile,OrthoInfoFile){
  GenesInfo = getGenes(gffFile)
  GenesInfo = addSpottedGenesOnArrayInfo(arrayInfoFile,GenesInfo)
  GenesInfo = addOrthoInfo(OrthoInfoFile, GenesInfo)
  GenesInfo = removeGenesNotOnArray(GenesInfo)
  return(GenesInfo)
}



addOrthoInfo <- function(OrthoInfoFile, Genes){
  ## OrthoMCLInfo
  OrthoInfo = read.table(file = OrthoInfoFile,  sep = "\t")
  head(OrthoInfo)
  OrthoInfo = OrthoInfo[grep("\\.1",OrthoInfo$V2) , ]
  ll <- unlist(strsplit(as.character(OrthoInfo$V2), ".", fixed=TRUE))
  ll = ll[seq(from = 1, to = length(ll)-1, by = 2)]
  OrthoInfo$Genes = ll
  OrthoInfo = OrthoInfo[with(OrthoInfo, order(V2)),]
  
  
  
  OrthoInfoGenes = OrthoInfo[!duplicated(OrthoInfo$Genes),c(1,3)]
  colnames(OrthoInfoGenes) <- c("OrthoMCLCluster","Genes") 
  
  Genes = merge(Genes,OrthoInfoGenes, by.x = "ID", by.y = "Genes", all.x = TRUE)
  Genes$OrthoMCLCluster = as.character(Genes$OrthoMCLCluster)
  Genes$OrthoMCLCluster[which(is.na(Genes$OrthoMCLCluster))] = which(is.na(Genes$OrthoMCLCluster))
  return (Genes)
}

addExpressionInfo <- function(ExpressionFile, Genes){
  Expression = read.table(file = ExpressionFile[1],  sep = "\t")
  Expression = Expression[!duplicated(Expression$V1),]
  if(dim(Expression)[2] > 2){
    Expression[, "max"] <- apply(Expression[, 2: dim(Expression)[2]  ], 1, max)
    Expression[, "min"] <- apply(Expression[, 2: dim(Expression)[2]  ], 1, min)
  }
  else{
    Expression$max = Expression$min = Expression[, 2] 
  }
  
  Expression = Expression[,c("V1","max","min")]
  colnames(Expression) <- c("Genes","max","min")  
  Genes = merge(Genes,Expression,by.x = "ID", by.y = "Genes", all.x = TRUE)  
}


getNrOfGOIs<- function(Genes ,aboveCutoff ,belowCutoff){
  GOIgenes = Genes[union(which(Genes$max >=aboveCutoff),which(Genes$min <= belowCutoff)), ]
  return (nrow(GOIgenes))
}

# getCluster <- function(GOIgenes,maxDistance){
#   GOIgenes =  GOIgenes[with(GOIgenes, order(chr,start)),]
#   nrOfGenes = dim(GOIgenes)[1]
#   compare = nrOfGenes -1 
#   GOIgenes$distanceForward = c(GOIgenes$start[2:nrOfGenes] -GOIgenes$stop[1:compare],maxDistance+1)
#   GOIgenes$distanceReverse = c(maxDistance+1,GOIgenes$distanceForward[1:compare] ) 
#   GOIgenes$nextChr = c(GOIgenes$chr[2:nrOfGenes],GOIgenes$chr[nrOfGenes]) 
#   GOIgenes$distanceForward[which(GOIgenes$chr != GOIgenes$nextChr)] = maxDistance+1
#   GOIgenes$distanceReverse[which(GOIgenes$chr != GOIgenes$nextChr)+1] = maxDistance+1
#   GOIgenes$geneBefore = c("not",GOIgenes$ID[1:compare])
#   GOIgenes$geneAfter = c(GOIgenes$ID[2:nrOfGenes],"not")
#   
#   withinDistance = sort(which(GOIgenes$distanceForward <=maxDistance & GOIgenes$distanceReverse <= maxDistance))
#   GOIgenesDistance = GOIgenes[withinDistance,]
#   
#   ClusteredGenes = union(GOIgenesDistance$ID,GOIgenesDistance$geneBefore)
#   ClusteredGenes = union(ClusteredGenes,GOIgenesDistance$geneAfter)
#   
#   
#   finalSet =  GOIgenes[GOIgenes$ID %in% ClusteredGenes, ]
#   return(finalSet)
#   
# }


getCluster <- function(GOIgenes,maxDistance){
  GOIgenes =  GOIgenes[with(GOIgenes, order(chr,start)),]
  nrOfGenes = dim(GOIgenes)[1]
  compare = nrOfGenes -1 
  GOIgenes$distanceForward = c(GOIgenes$start[2:nrOfGenes] -GOIgenes$stop[1:compare],maxDistance+1)
  GOIgenes$nextChr = c(GOIgenes$chr[2:nrOfGenes],GOIgenes$chr[nrOfGenes]) 
  GOIgenes$distanceForward[which(GOIgenes$chr != GOIgenes$nextChr)] = maxDistance+1
  GOIgenes$geneAfter = c(GOIgenes$ID[2:nrOfGenes],"not")
  
  withinDistance = sort(which(GOIgenes$distanceForward <=maxDistance))
  GOIgenesDistance = GOIgenes[withinDistance,]
  
  ClusteredGenes = union(GOIgenesDistance$ID,GOIgenesDistance$geneAfter)
  finalSet =  GOIgenes[GOIgenes$ID %in% ClusteredGenes, ]
  return(finalSet)
  
}

IdentifyClusters <- function(finalSet,nrOfGOIsPerClusterCutoff,maxDistance){
  nrOfGOIs = nrow(finalSet)
  if(nrOfGOIs == 0){
    return (finalSet)
  }
  breaks =  c(which(finalSet$distanceForward[1:(nrOfGOIs-1)] >  maxDistance),nrOfGOIs)
  nrOfClusters = length(breaks)
  
  if(nrOfClusters > 1){
    nrOfGOISperCluster = breaks - c(0,breaks[1:nrOfClusters-1])
    finalSet$Cluster = rep.int(1:nrOfClusters, nrOfGOISperCluster)
  }else{
    nrOfGOISperCluster = nrOfGOIs
    finalSet$Cluster = 1
  }
  
  keep = which(nrOfGOISperCluster >= nrOfGOIsPerClusterCutoff)
  finalSetFiltered = finalSet[finalSet$Cluster %in% keep, ]
  nrOfGOIsFiltered = nrow(finalSetFiltered)
  if(nrOfGOIsFiltered == 0){
    return (finalSetFiltered)
  }
  breaks =  c(which(finalSetFiltered$distanceForward[1:(nrOfGOIsFiltered-1)] >  maxDistance),nrOfGOIsFiltered)
  nrOfClusters = length(breaks)
  if(nrOfClusters > 1){
    nrOfGOISperCluster = breaks - c(0,breaks[1:nrOfClusters-1])
    finalSetFiltered$Cluster = rep.int(1:nrOfClusters, nrOfGOISperCluster)
  }else{
    
    finalSetFiltered$Cluster = 1
  }
  return(finalSetFiltered)
}


getOrthologousCluster <- function(ClusteredGenes,nrOfOrthologs){
  out <- split( ClusteredGenes$OrthoMCLCluster , f = ClusteredGenes$Cluster )
  uniqueOrhtologs <- lapply(out ,unique)
  uniqueOrthologsCount <- which(unlist(lapply(uniqueOrhtologs,length))>=nrOfOrthologs)
  Clustered_nonOrthologous_Genes = ClusteredGenes[ClusteredGenes$Cluster %in% uniqueOrthologsCount, ]
  
  return (Clustered_nonOrthologous_Genes)
}



getClustersTrue <- function(Genes ,aboveCutoff ,belowCutoff, maxDistance, nrOfGOIsPerClusterCutoff ,nrOfOrthologs ){
  ## Identify genes that are of interest for study
  GOIgenes = Genes[union(which(Genes$max >=aboveCutoff),which(Genes$min <= belowCutoff)), ]
  
  ## Identify it the GOI that are at least pairs in clusters
  
  finalSet = getCluster(GOIgenes,maxDistance)
  
  ## Identify the Clusters
  finalSetFiltered = IdentifyClusters(finalSet,nrOfGOIsPerClusterCutoff = nrOfGOIsPerClusterCutoff,maxDistance = maxDistance)
  
  
  finalSet_nonOrthologous = getOrthologousCluster(finalSetFiltered,nrOfOrthologs)
  clustersNR = unique(finalSet_nonOrthologous$Cluster)
  finalSet_nonOrthologous$Cluster = mapvalues(finalSet_nonOrthologous$Cluster,clustersNR,1:length(clustersNR))
  
  return(finalSet_nonOrthologous)
  
  
}

getClustersScrambledInfo <- function(Genes ,nrOfGOI ,maxDistance = 10000,nrOfGOIsPerClusterCutoff, nrOfOrthologs){
  GOIgenes = Genes[sort(sample(nrow(Genes), nrOfGOI)), ]
  finalSet = getCluster(GOIgenes,maxDistance)
  finalSet = IdentifyClusters(finalSet,nrOfGOIsPerClusterCutoff = nrOfGOIsPerClusterCutoff,maxDistance = maxDistance)
  finalSet = getOrthologousCluster(finalSet,nrOfOrthologs )
  
  nrOfGOIs = nrow(finalSet)
  nrOfClusters = length(unique(finalSet$Cluster))
  return(c(nrOfClusters, nrOfGOIs))
}


getClustersScrambledFromList <- function(x,listOfVariables){
  Genes2 = listOfVariables[[1]]  
  count2 = listOfVariables[[2]]
  dist2 = listOfVariables[[3]]
  nrOfGenes = listOfVariables[[4]]
  nrOrt = listOfVariables[[5]]
  
  clusterInfo = getClustersScrambledInfo(Genes = Genes2,nrOfGOI = count2, maxDistance = dist2 ,nrOfGOIsPerClusterCutoff =  nrOfGenes, nrOfOrthologs = nrOrt)
  return (clusterInfo)
}

removeGenesNotOnArray <- function(Genes){
  Genes = Genes[Genes$array ==1, ]
}


getDistribution <- function(GenesInfo, ExpressionFile, belowCutoff,aboveCutoff,maxDistances,nrOfGenesInCluster,nrOfOrthologs, nrOfIterations,outFileName ){
  
  
  
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
  
  for(i in 1:length(ExpressionFile)){
    Genes = addExpressionInfo(ExpressionFile[i],GenesInfo)
    nrOfGOI = getNrOfGOIs(Genes = Genes,aboveCutoff = aboveCutoff, belowCutoff = belowCutoff)
    for(j in 1:length(maxDistances)){
      ScrambleInfo = list(Genes = Genes,count = nrOfGOI,dist = maxDistances[j],nrOfGOIsPerClusterCutoff = nrOfGenesInCluster,nrOrt =  nrOfOrthologs)
      files <- lapply(1:nrOfIterations ,getClustersScrambledFromList, listOfVariables = ScrambleInfo)
      Distribution = data.frame(matrix(unlist(files), nrow=nrOfIterations, byrow=T))
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
  
  plotDistribution(TrueDistributions,BackgroundDistributions)
  
  ggsave(paste(outFileName,"figure.pdf",sep = "_"))
  write.table(x = TrueDistributions, file = paste(outFileName,"true_distributions.tab.txt",sep = "_"),sep = "\t")       
  write.table(x = BackgroundDistributions, file = paste(outFileName,"random_distributions.tab.txt",sep = "_"),sep = "\t")       
  
}

plotDistribution <- function(TrueDistributions, BackgroundDistributions,TitleInfo){
  
  p <- ggplot(BackgroundDistributions, aes(x = factor(maxDistance), y = NrOfClusters))
  p + geom_boxplot(aes(fill = factor(Dataset)))+ 
    theme_bw()+ggtitle(TitleInfo)+
    geom_point(data = TrueDistributions, size = 3, aes(factor(maxDistance), Clusters,shape = factor(pValue2), colour = factor(Dataset),fill = factor(Dataset)))
  
}

calculatePvalues <- function(TrueDistributions,BackgroundDistributions){
  Dataset = unique(TrueDistributions$Dataset)
  maxDistances = unique(TrueDistributions$maxDistance)
  for(i in 1:length(Dataset)){
    for(j in 1:length(maxDistances)){
      Distribution = BackgroundDistributions[BackgroundDistributions$maxDistance == maxDistances[j] & BackgroundDistributions$Dataset == Dataset[i],  ]
      TrueDistribution = TrueDistributions[TrueDistributions$maxDistance == maxDistances[j] & TrueDistributions$Dataset == Dataset[i],  ]
      position = length(which(Distribution$NrOfClusters >=  TrueDistribution$Clusters)) +1
      pValue = position / length(Distribution$NrOfClusters)
      TrueDistributions$pValue[TrueDistributions$maxDistance == maxDistances[j] & TrueDistributions$Dataset == Dataset[i]] = pValue
      
    }
  }
  
  TrueDistributions$pValue2 = as.character(TrueDistributions$pValue2)
  
  TrueDistributions$pValue2[TrueDistributions$pValue >= 0.01] = ">= 1e10-2" 
  TrueDistributions$pValue2[TrueDistributions$pValue < 0.01 ] = "< 1e10-2" 
  TrueDistributions$pValue2[TrueDistributions$pValue < 0.001 ] = "< 1e10-3" 
  TrueDistributions$pValue2[TrueDistributions$pValue < 0.0001] = "< 1e10-4" 
  TrueDistributions$pValue2[TrueDistributions$pValue < 0.00001] = "< 1e10-5" 
  
  return(TrueDistributions)
}

getClusters <- function(GenesInfo, ExpressionFile, belowCutoff,aboveCutoff, maxDistance, nrOfGenesInCluster, nrOfOrthologs){
  
  Genes = addExpressionInfo(ExpressionFile,GenesInfo)
  Cluster = getClustersTrue(Genes = Genes,aboveCutoff = aboveCutoff, belowCutoff = belowCutoff, maxDistance = maxDistance,nrOfGOIsPerClusterCutoff = nrOfGenesInCluster, nrOfOrthologs = nrOfOrthologs)
  return(Cluster) 
}





