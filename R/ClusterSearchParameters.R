#Load Functions that should be in the same as the ClusterSerch.R location
source("ClusterSearchFunctions.R")


# Set the correct parameters for the program. 

aboveCutoff = 2
belowCutoff = -1000
maxDistance = 10000
nrOfGenesInCluster= 3
nrOfOrthologs = 2

#Change directory to where all the data is.

setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/newSearch/")

gffFile = "Athaliana_167_gene.gff3"
arrayInfoFile = "Operon_set1_4.0_oligos.original_ACC_match.Athaliana_167_TAIR9.trimmed.matchGff3.blast" 
OrthoInfoFile = "orthoMCL.parsed.AthalianaOnly.atLeast2.onePerLine.tab.txt2"
ExpressionFile = c("Wellmer2007_MS1.TAIR10.expression", "Wellmer2004_stamen.TAIR10.expression")
Dataset = c("Alves_2007","Wellmer2004")
nrOfTimes = 10000

