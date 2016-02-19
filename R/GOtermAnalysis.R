# Importing data from master data table. 
library(stringr)
require(stats)


setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/DAVID")
masterTable = read.table(file="Merged_Final_Table.tab.txt", sep="\t", header =TRUE)
head(masterTable)
dim(masterTable)


spottedGenes = masterTable[masterTable$array ==1,]


AlvesGOI = spottedGenes[!is.na(spottedGenes$AlvesFerreria2007) & spottedGenes$AlvesFerreria2007>2, ]
WellmerGOI = spottedGenes[!is.na(spottedGenes$Wellmer2004) &spottedGenes$Wellmer2004 >2, ]
BothGOI = AlvesGOI[AlvesGOI$Gene %in% intersect(WellmerGOI$Gene,AlvesGOI$Gene ), ]


WellmerClusterGOI = WellmerGOI[!is.na(WellmerGOI$Wellmer2004Cluster)  , ]
AlvesClusterGOI =  AlvesGOI[!is.na(AlvesGOI$Alves_Fereria2007Cluster), ]
BothClusterGOI = AlvesClusterGOI[AlvesClusterGOI$Gene %in% intersect(WellmerClusterGOI$Gene,AlvesClusterGOI$Gene ), ]


write.table(spottedGeneNames, "1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(WellmerGOI$Gene, "2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(AlvesGOI$Gene, "3.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(BothGOI$Gene, "4.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(WellmerClusterGOI$Gene, "5.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(AlvesClusterGOI$Gene, "6.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(BothClusterGOI$Gene, "7.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(union(AlvesClusterGOI$Gene,WellmerClusterGOI$Gene), "5_6.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(union(WellmerGOI$Gene,AlvesGOI$Gene), "2_3.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

WellmerGOI = Wellmer[]

dim(Wellmer)

dim(Alves)

dim(spottedGenes)


ModelGeneGOrelationship = getModelGeneGOrelationship(ModelGeneGOrelationshipFile = "/Users/johanreimegard/Vetenskap/Data/GOtermAnalysis/gene_association.tair",ModelGeneGOrelationshipFileFormat = "GOC")

ModelGeneGOrelationship = ModelGeneGOrelationship[, c("Locus","Function","GOterm")]


######23 vs 7
BackgroundTable = data.frame(Locus= union(WellmerGOI$Gene,AlvesGOI$Gene))
GOItable = data.frame(Locus= BothClusterGOI$Gene)

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "P", GOtermCutoff = 0.001,AnalysisName = "3∪4_vs_8_P") 
printTopGOresults(topGOdata = TopGOdata, outName ="3∪4_vs_8_P")                   

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "F", GOtermCutoff = 0.001,AnalysisName = "3∪4_vs_8_F") 
printTopGOresults(topGOdata = TopGOdata, outName ="3∪4_vs_8_F")                   

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "C", GOtermCutoff = 0.001,AnalysisName = "3∪4_vs_8_C") 
printTopGOresults(topGOdata = TopGOdata, outName ="3∪4_vs_8_C")                   



######56 vs 7
BackgroundTable = data.frame(Locus= BothGOI$Gene)
GOItable = data.frame(Locus= BothClusterGOI$Gene)

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "P", GOtermCutoff = 0.001,AnalysisName = "5_vs_8_P") 
printTopGOresults(topGOdata = TopGOdata, outName ="5_vs_8_P")                   

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "F", GOtermCutoff = 0.001,AnalysisName = "5_vs_8_F") 
printTopGOresults(topGOdata = TopGOdata, outName ="5_vs_8_F")                   

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "C", GOtermCutoff = 0.001,AnalysisName = "5_vs_8_C") 
printTopGOresults(topGOdata = TopGOdata, outName ="5_vs_8_C")         

######56 vs 7
BackgroundTable = data.frame(Locus= spottedGenes$Gene)
GOItable = data.frame(Locus= BothClusterGOI$Gene)

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "P", GOtermCutoff = 0.001,AnalysisName = "2_vs_8_P") 
printTopGOresults(topGOdata = TopGOdata, outName ="2_vs_8_P")                   

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "F", GOtermCutoff = 0.001,AnalysisName = "2_vs_8_F") 
printTopGOresults(topGOdata = TopGOdata, outName ="2_vs_8_F")                   

TopGOdata = getTopGOdata(ModelGeneGOrelationship = ModelGeneGOrelationship,BackgroundTable = BackgroundTable , GOItable = GOItable, GOtermFunction = "C", GOtermCutoff = 0.001,AnalysisName = "2_vs_8_C") 
printTopGOresults(topGOdata = TopGOdata, outName ="2_vs_8_C")         






