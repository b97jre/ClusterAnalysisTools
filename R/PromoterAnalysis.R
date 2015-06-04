
setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/promoterRegions/17clusterAnalysis")



##Reading in information from file regarding ATTEDs found in cluster promoter regions

AttedInfo = read.table("17clusters.promoter_ATTED_number_of_genes_per_motif_with_genenames.transposed.tab.colnames.txt" ,header = TRUE, sep = "\t")

AttedNames = colnames(AttedInfo)

Genes = read.table("17clusters.promoter_ATTED_number_of_genes_per_motif_with_genenames.transposed.tab.withoutHeader.txt",
                   header = FALSE, sep = "\t", col.names = colnames, fill = TRUE)






#### Reading in info from file regarding ATTED found in all promoterRegions
AttedGenome = read.table("Athaliana_167_gene.gff3.upstream.1000.fa.fa_ATTED-II-7MERS.fa.matches", sep = "\t", header = FALSE)
ll <- unlist(strsplit(as.character(AttedGenome$V1), ".", fixed=TRUE))
ll = ll[seq(from = 1, to = length(ll)-1, by = 2)]
AttedGenome$Genes = ll



AttedGenome$Both  = paste(AttedGenome[, 'Genes'],AttedGenome[, 2],sep = "_")
AttedGenomeUnique = AttedGenome[!duplicated(AttedGenome$Both) ,]
ATTED = AttedGenomeUnique[, 2]
GenomeWide  = table(ATTED)

GenomeWidePercent = GenomeWide/27416


#####Merging the data##############


AllClusters = character() 
for(i in 1:length(AttedNames)){
  AttedGenes = Genes[,i]
  AttedInfo = AttedNames[i]
  AttedProb = GenomeWidePercent[AttedNames[i]]
  for(j in 1:17){
    ClusterGenes = finalOutput[ which(finalOutput$Cluster == j), 'Gene' ]
    ClusterGenesInAtted = ClusterGenes[ClusterGenes %in% AttedGenes]
    if(length(ClusterGenes) == length(ClusterGenesInAtted)){
      prob = as.numeric(AttedProb^length(ClusterGenes))
      sep1 = paste(paste("Cluster",j,sep = "_"), paste("prob", prob, sep = " ") , sep = ",")
      AttedInfo = paste(AttedInfo, paste("(", ")",sep = sep1), sep = " ")
      AllClusters = c(AllClusters, paste("(", ")",sep = sep1), sep = " ")
    }
  }
  print(AttedInfo)
}
 
table(AllClusters)





temp = AllClusters[!duplicated(AllClusters)]
FinalClusters = temp[2:length(temp)] 
FinalClusters[order(FinalClusters)]


}