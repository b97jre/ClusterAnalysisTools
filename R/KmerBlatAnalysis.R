setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/promoterRegions")

EightMers = read.table(file = 'Athaliana_167_gene.gff3.upstream.1000.fa.fa_At_reliable_REGmotifs.fasta.matches')
SevenMers = read.table(file = 'Athaliana_167_gene.gff3.upstream.1000.fa.fa_ATTED-II-7MERS.fa.matches')


head(SevenMers)


ll <- unlist(strsplit(as.character(EightMers$V1), "\\."))
EightMers$Genes = ll[seq(from = 1, to = length(ll)-1, by = 2)]

EightMers$GOI = 0
EightMers$Cluster = 0


ll <- unlist(strsplit(as.character(SevenMers$V1), "\\."))
SevenMers$Genes = ll[seq(from = 1, to = length(ll)-1, by = 2)]




GenesInCluster = read.table(file = 'ClusterNrAndGenesInvolved.txt',  sep = "\t")

NrOfClusters = strsplit(as.character(GenesInCluster$V3),",")
Genes = unlist(strsplit(as.character(GenesInCluster$V3),","))
GOI = unlist(strsplit(as.character(GenesInCluster$V2),""))

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

MerClusterColnames  = c("Genes","mRNA","MER","Dir","start","stop","GOI","Cluster") 
SevenMersCluster = merge(x =SevenMers, y = ClusterInfo, by = "Genes" )
colnames(SevenMersCluster) <- MerClusterColnames
EightMersCluster = merge(x =EightMers, y = ClusterInfo, by = "Genes" )


