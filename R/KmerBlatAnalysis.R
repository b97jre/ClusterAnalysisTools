setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/promoterRegions")

EightMers = read.table(file = 'Athaliana_167_gene.gff3.upstream.1000.fa.fa_At_reliable_REGmotifs.fasta.matches')
SevenMers = read.table(file = 'Athaliana_167_gene.gff3.upstream.1000.fa.fa_ATTED-II-7MERS.fa.matches')



ll <- unlist(strsplit(as.character(EightMers$V1), "\\."))
EightMers$Genes = ll[seq(from = 1, to = length(ll)-1, by = 2)]


ll <- unlist(strsplit(as.character(SevenMers$V1), "\\."))
SevenMers$Genes = ll[seq(from = 1, to = length(ll)-1, by = 2)]

class(SevenMers$Genes)



GenesInCluster = read.table(file = 'ClusterNrAndGenesInvolved.txt',  sep = "\t")
NrOfClusters = strsplit(as.character(GenesInCluster$V3),",")
Genes = trim(unlist(strsplit(as.character(GenesInCluster$V3),",")))
GOI = trim(unlist(strsplit(as.character(GenesInCluster$V2),"")))
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

ClusterInfo$GOI = as.character(ClusterInfo$GOI)
ClusterInfo$GOI[ClusterInfo$GOI == 0] = "no"
ClusterInfo$GOI[ClusterInfo$GOI == 1] = "yes"


MerClusterColnames  = c("Genes","mRNA","MER","Dir","start","stop","GOI","Cluster") 

SevenMersCluster = merge(x =SevenMers, y = ClusterInfo, by = "Genes" ,all.y = TRUE)
colnames(SevenMersCluster) <- MerClusterColnames
EightMersCluster = merge(x =EightMers, y = ClusterInfo, by = "Genes" )
colnames(EightMersCluster) <- MerClusterColnames


write.table(SevenMersCluster,file="ATTED_II_CLUSTER_Table.tab.txt", sep="\t", quote =  FALSE,row.names=FALSE)
write.table(EightMersCluster,file="PPDB_CLUSTER_Table.tab.txt", sep="\t", quote =  FALSE,row.names=FALSE)


trim <- function (x) gsub("^\\s+|\\s+$", "", x)
