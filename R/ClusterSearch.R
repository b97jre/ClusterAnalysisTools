

#Load Functions that should be in the same as the ClusterSerch.R location
source("/Users/johanreimegard/git/clusterProject/R/ClusterSearchFunctions.R")


# Set the correct parameters for the program. 

aboveCutoff = 2
belowCutoff = -1000
#Change directory to where all the data is.

setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/newSearch/")

gffFile = "Athaliana_167_gene.gff3"
arrayInfoFile = "Operon_set1_4.0_oligos.original_ACC_match.Athaliana_167_TAIR9.trimmed.matchGff3.blast" 
OrthoInfoFile = "orthoMCL.parsed.AthalianaOnly.atLeast2.onePerLine.tab.txt2"

ExpressionFile = c("Wellmer2007_MS1.TAIR10.expression", "Wellmer2004_stamen.TAIR10.expression")
Dataset = c("Alves_2007","Wellmer2004")
maxDistances = c(2000,5000,10000,20000,40000)

GenesInfo = loadGeneInfo(gffFile = gffFile,arrayInfoFile = arrayInfoFile,OrthoInfoFile = OrthoInfoFile)
  

# Set up to determine best distance 
maxDistances = c(2000,5000,10000,20000,40000)
nrOfIterations = 100100



nrOfGenesInCluster= 2
nrOfOrthologs = 2
outFileName = "Distribution_2_2.1"
getDistribution(GenesInfo, ExpressionFile, belowCutoff,aboveCutoff,maxDistances,nrOfGenesInCluster,nrOfOrthologs, nrOfIterations,outFileName)

nrOfGenesInCluster= 3
nrOfOrthologs = 2
outFileName = "Distribution_3_2.1"
getDistribution(GenesInfo, ExpressionFile, belowCutoff,aboveCutoff,maxDistances,nrOfGenesInCluster,nrOfOrthologs, nrOfIterations,outFileName)

nrOfGenesInCluster= 3
nrOfOrthologs = 3
outFileName = "Distribution_3_3.1"
getDistribution(GenesInfo, ExpressionFile, belowCutoff,aboveCutoff,maxDistances,nrOfGenesInCluster,nrOfOrthologs, nrOfIterations,outFileName)



dist2_2_1 = read.table("Distribution_2_2.1_random_distributions.tab.txt" ,sep ="\t", header = TRUE)
dist2_2 = read.table("Distribution_2_2_random_distributions.tab.txt" ,sep ="\t", header = TRUE)
trueDist2_2 = read.table("Distribution_2_2_true_distributions.tab.txt",sep ="\t", header = TRUE)
dist2_2 = rbind(dist2_2,dist2_2_1)
trueDist2_2 = calculatePvalues(TrueDistributions =trueDist2_2, BackgroundDistributions = dist2_2 )
outFileName = "Distribution_2_2"
plotDistribution(TrueDistributions = trueDist2_2, BackgroundDistributions = dist2_2,TitleInfo = "Distributions \n N=2\nO=2")
ggsave(paste(outFileName,"figure.pdf",sep = "_"), useDingbats=FALSE)


dist3_2_1 = read.table("Distribution_3_2.1_random_distributions.tab.txt" ,sep ="\t", header = TRUE)
dist3_2 = read.table("Distribution_3_2_random_distributions.tab.txt" ,sep ="\t", header = TRUE)
trueDist3_2 = read.table("Distribution_3_2_true_distributions.tab.txt",sep ="\t", header = TRUE)
dist3_2 = rbind(dist3_2,dist3_2_1)
trueDist3_2 = calculatePvalues(TrueDistributions =trueDist3_2, BackgroundDistributions = dist3_2 )
outFileName = "Distribution_3_2"
plotDistribution(TrueDistributions = trueDist3_2, BackgroundDistributions = dist3_2,TitleInfo = "Distributions \n N=3\nO=2")
ggsave(paste(outFileName,"figure.pdf",sep = "_"), useDingbats=FALSE)

dist3_3_1 = read.table("Distribution_3_3.1_random_distributions.tab.txt" ,sep ="\t", header = TRUE)
dist3_3 = read.table("Distribution_3_3_random_distributions.tab.txt" ,sep ="\t", header = TRUE)
trueDist3_3 = read.table("Distribution_3_3_true_distributions.tab.txt",sep ="\t", header = TRUE)
dist3_3 = rbind(dist3_3,dist3_3_1)
trueDist3_3 = calculatePvalues(TrueDistributions =trueDist3_3, BackgroundDistributions = dist3_3)
outFileName = "Distribution_3_3"
plotDistribution(TrueDistributions = trueDist3_3, BackgroundDistributions = dist3_3,TitleInfo = "Distributions \n N=3\nO=3" )
ggsave(paste(outFileName,"figure.pdf",sep = "_"), useDingbats=FALSE)



