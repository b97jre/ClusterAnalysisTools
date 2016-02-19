
require(stats)


setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/")
fileName1 = '/Users/johanreimegard/Vetenskap/Data/clusterProject/clusters/Wellmer2004Clusters.txt'
fileName2 = '/Users/johanreimegard/Vetenskap/Data/clusterProject/clusters/Wellmer2007Clusters.txt'


mergedClusters = getMergedClusters(fileName1,fileName2)
length(which(mergedClusters$GOI == 1))


#########################
##   FUNCTIONS START
########################


trim <- function (x) gsub("^\\s+|\\s+$", "", x)



getClusterInfo <-  function(filename){
  GenesInCluster = read.table(file = filename,header = TRUE)
  GenesInCluster$ACCs_of_all_genes_in_cluster = trim(GenesInCluster$ACCs_of_all_genes_in_cluster)
  NrOfClusters = strsplit(as.character(GenesInCluster$ACCs_of_all_genes_in_cluster),",")
  Genes = trim(unlist(strsplit(as.character(GenesInCluster$ACCs_of_all_genes_in_cluster),",")))
  GOI =trim( unlist(strsplit(as.character(GenesInCluster$SEGorNot),"")))
  
  ClusterInfo = data.frame(Genes,GOI)
  ClusterInfo$Cluster =-1
  
  ## adding genecluster to Cluster info
  pointer = 1
  for(i in 1:length(NrOfClusters)){
    for(j in 1:length(NrOfClusters[[i]])){
      ClusterInfo$Cluster[pointer] = i
      pointer = pointer+1 
    }    
  }
  
  return (ClusterInfo)
  
}



getMergedClusters <- function(fileName1, fileName2){   
  
  Cluster1 = getClusterInfo(fileName1)
  colnames(Cluster1) <- c("Gene","GOI","Cluster")
  
  Cluster2 = getClusterInfo(fileName2)
  colnames(Cluster2) <- c("Gene","GOI","Cluster")
  
  Shared1Clusters =  as.numeric(levels(as.factor(Cluster1$Cluster[Cluster1$Gene %in% intersect(Cluster1$Gene, Cluster2$Gene)])))
  Cluster1Shared = Cluster1[which(Cluster1$Cluster %in% Shared1Clusters), ] 
  
  Shared2Clusters =  as.numeric(levels(as.factor(Cluster2$Cluster[Cluster2$Gene %in% intersect(Cluster1$Gene, Cluster2$Gene)])))
  Cluster2Shared = Cluster2[Cluster2$Cluster %in% Shared2Clusters, ] 
  
  mergedClusters = merge(Cluster1Shared,Cluster2Shared, by = "Gene", all = TRUE)
  mergedClusters$Gene = as.character(mergedClusters$Gene)
  mergedClusters = mergedClusters[with(mergedClusters, order(Gene)),]
  
  mergedClusters$newCX = 0
  mergedClusters$newCY = 0
  
  clusterNr = 1
  for(i in unique(mergedClusters$Cluster.x)[!is.na(unique(mergedClusters$Cluster.x))]){
    mergedClusters$newCX[mergedClusters$Cluster.x == i] = clusterNr 
    clusterNr = clusterNr+1
  }
  clusterNr = 1
  for(i in unique(mergedClusters$Cluster.y)[!is.na(unique(mergedClusters$Cluster.y))]){
    mergedClusters$newCY[mergedClusters$Cluster.y == i] = clusterNr 
    clusterNr = clusterNr+1
  }
  
  mergedClusters$Cluster = mergedClusters$newCX
  mergedClusters$Cluster[mergedClusters$newCX == 0] = mergedClusters$newCY[mergedClusters$newCX == 0]
  
  mergedClusters[is.na(mergedClusters)] = 0
  
  mergedClusters$GOI.x = as.numeric(mergedClusters$GOI.x)-1
  mergedClusters$GOI.y = as.numeric(mergedClusters$GOI.y)-1
  
  mergedClusters$GOI = ceiling((mergedClusters$GOI.x + mergedClusters$GOI.y)/2)
  
  finalMergedClusters = mergedClusters[,c("Gene","Cluster","GOI")]
  
  return (finalMergedClusters)
}
