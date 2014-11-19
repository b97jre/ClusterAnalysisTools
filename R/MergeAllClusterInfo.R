
library(stringr)
require(stats)


setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/clusters")


## Gene info
gff3OutPut = read.table("Athaliana_167_gene.gff3", sep = "\t", header =  FALSE, skip =1, stringsAsFactors=FALSE ) 
gff3Genes = gff3OutPut[which(gff3OutPut[[3]] == "gene"), ]


ll <- unlist(strsplit(gff3Genes$V9, ";", fixed=TRUE))
ll = ll[seq(from = 1, to = length(ll)-1, by = 2)]
ll = unlist(strsplit(ll, "=", fixed=TRUE))
ll = ll[seq(from = 2, to = length(ll), by = 2)]
gff3Genes$ID = ll

## Genes spotted on array based on blast analysis
gff3Genes$array = 0

arraysMatches = read.table( sep = "\t", file = "Operon_set1_4.0_oligos.original_ACC_match.Athaliana_167_TAIR9.trimmed.matchGff3.blast")
gff3Genes$array[gff3Genes$ID %in% arraysMatches$V16] = 1




##ClusterInfo
GenesInCluster = read.table(file = 'ClusterNrAndGenesInvolved.txt',  sep = "\t")

NrOfClusters = strsplit(as.character(GenesInCluster$V3),",")
Genes = trim(unlist(strsplit(as.character(GenesInCluster$V3),",")))
GOI =trim( unlist(strsplit(as.character(GenesInCluster$V2),"")))

ClusterInfo = data.frame(Genes,GOI)
ClusterInfo$Cluster =-1

## adding genecluster to Cluster info
pointer = 1
for(i in 1:length(NrOfClusters)){
  for(j in 1:length(NrOfClusters[[i]])){
    ClusterInfo$Cluster[pointer] = GenesInCluster$V1[i]
    pointer = pointer+1 
  }    
}

welmer2004ClusterInfo = getClusterInfo('/Users/johanreimegard/Vetenskap/Data/clusterProject/clusters/Wellmer2004Clusters.txt')
welmer2007ClusterInfo = getClusterInfo('/Users/johanreimegard/Vetenskap/Data/clusterProject/clusters/Wellmer2007Clusters.txt')

welmer2004ClusterInfo2 = welmer2004ClusterInfo[,c(1,3)]
colnames(welmer2004ClusterInfo2) <- c("Gene","Wellmer2004Cluster")

welmer2007ClusterInfo2 = welmer2007ClusterInfo[,c(1,3)]
colnames(welmer2007ClusterInfo2) <- c("Gene","Alves_Fereria2007Cluster")


## OrthoMCLInfo
OrthoInfo = read.table(file = 'orthoMCL.parsed.AthalianaOnly.atLeast2.onePerLine.tab.txt',  sep = "\t")

ll <- unlist(strsplit(as.character(OrthoInfo$V2), ".", fixed=TRUE))
ll = ll[seq(from = 1, to = length(ll)-1, by = 2)]
OrthoInfo$Genes = ll
OrthoInfo = OrthoInfo[with(OrthoInfo, order(V2)),]


OrthoInfoGenes = OrthoInfo[!duplicated(paste(OrthoInfo$V1,OrthoInfo$Genes)),c(1,3)]
colnames(OrthoInfoGenes) <- c("OrthoMCLCluster","Genes") 
head(OrthoInfoGenes)

## Expression 
Wellmer2004 = read.table(file = 'Wellmer2004_stamen.TAIR10.expression',  sep = "\t")
AlvesFerreria2007 = read.table(file = 'Wellmer2007_MS1.TAIR10.expression',  sep = "\t")

AlvesFerreria2007 = AlvesFerreria2007[AlvesFerreria2007$V9 >=2 , c(1,9) ]

colnames(AlvesFerreria2007) <- c("Genes","AlvesFerreria2007")
colnames(Wellmer2004) <- c("Genes","Wellmer2004")




## Merge Data
GenesWithOrtho = merge(gff3Genes,OrthoInfoGenes, by.x = "ID", by.y = "Genes", all.x = TRUE)
head(GenesWithOrtho)

GenesWithOrthoAndCluster = merge(GenesWithOrtho,ClusterInfo,by.x = "ID", by.y = "Genes", all.x = TRUE)

GenesWithOrthoAndClusterAndExpression= merge(GenesWithOrthoAndCluster,Wellmer2004,by.x = "ID", by.y = "Genes", all.x = TRUE)
GenesWithOrthoAndClusterAndExpression= merge(GenesWithOrthoAndClusterAndExpression,AlvesFerreria2007,by.x = "ID", by.y = "Genes", all.x = TRUE)
GenesWithOrthoAndClusterAndExpression$kind = "Other"
GenesWithOrthoAndClusterAndExpression$kind[GenesWithOrthoAndClusterAndExpression$Cluster>0] = "Cluster"
for(i in seq(from=3,to = dim(GenesWithOrthoAndClusterAndExpression)[1]-2,by = 1)){
  if(GenesWithOrthoAndClusterAndExpression$kind[i] =="Other"){
    if(GenesWithOrthoAndClusterAndExpression$kind[i+2] == "Cluster"){GenesWithOrthoAndClusterAndExpression$kind[i] = "Adjacent"}
    if(GenesWithOrthoAndClusterAndExpression$kind[i+1] == "Cluster"){GenesWithOrthoAndClusterAndExpression$kind[i] = "Adjacent"}
    if(GenesWithOrthoAndClusterAndExpression$kind[i-1] == "Cluster"){GenesWithOrthoAndClusterAndExpression$kind[i] = "Adjacent"}
    if(GenesWithOrthoAndClusterAndExpression$kind[i-2] == "Cluster"){GenesWithOrthoAndClusterAndExpression$kind[i] = "Adjacent"}
  }
}
 
GenesWithOrthoAndClusterAndExpression$GOI = "No"
GenesWithOrthoAndClusterAndExpression$GOI[GenesWithOrthoAndClusterAndExpression$Wellmer2004 >1 & GenesWithOrthoAndClusterAndExpression$AlvesFerreria2007 >1] = "yes"
GenesWithOrthoAndClusterAndExpression$GOIsingle = "No"
GenesWithOrthoAndClusterAndExpression$GOIsingle[GenesWithOrthoAndClusterAndExpression$Wellmer2004 >1 | GenesWithOrthoAndClusterAndExpression$AlvesFerreria2007 >1] = "yes"

GenesWithOrthoAndClusterAndExpression$OrthoMCLCluster = as.character(GenesWithOrthoAndClusterAndExpression$OrthoMCLCluster)
GenesWithOrthoAndClusterAndExpression$OrthoMCLCluster[is.na(GenesWithOrthoAndClusterAndExpression$OrthoMCLCluster)] <- "No_cluster"
GenesWithOrthoAndClusterAndExpression$OrthoMCLCluster = as.factor(GenesWithOrthoAndClusterAndExpression$OrthoMCLCluster)
head(GenesWithOrthoAndClusterAndExpression)





finalOutput = GenesWithOrthoAndClusterAndExpression[ ,c(1,2,5,6,8,11,17,14,13,18,15,16,12)] 
FinalColNames = c("Gene","Chr","left","right","strand", "array","kind","Cluster","GOI","GOI_in_one","Wellmer2004","AlvesFerreria2007","OrthoMCLCluster")
colnames(finalOutput) <- FinalColNames

finalOutput = merge(finalOutput,welmer2004ClusterInfo2,by.x = "Gene", by.y = "Gene", all.x = TRUE)
finalOutput = merge(finalOutput,welmer2007ClusterInfo2,by.x = "Gene", by.y = "Gene", all.x = TRUE)


finalOutput$GOI[finalOutput$array == 0 ]=NA
finalOutput$GOI_in_one[finalOutput$array == 0 ]=NA
finalOutput$Wellmer2004[finalOutput$array == 0 ]=NA
finalOutput$AlvesFerreria2007[finalOutput$array == 0 ]=NA

FinalColNames <-colnames(finalOutput)
FinalColNames[9] = "GOI_in_both"
colnames(finalOutput) <- FinalColNames

write.table(finalOutput,file="Merged_Final_Table.tab.txt", sep="\t", quote =  FALSE,row.names=FALSE)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


##ClusterInfo

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
