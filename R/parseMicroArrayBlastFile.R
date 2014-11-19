library(stringr)
require(stats)

oligoLength = 70



setwd("/Users/johanreimegard/Vetenskap/Data/clusterProject/blast")


gff3OutPut = read.table("Athaliana_167_gene.gff3", sep = "\t", header =  FALSE, skip =1, stringsAsFactors=FALSE ) 
gff3Genes = gff3OutPut[which(gff3OutPut[[3]] == "gene"), ]


ll <- unlist(strsplit(gff3Genes$V9, ";", fixed=TRUE))
ll = ll[seq(from = 1, to = length(ll)-1, by = 2)]
ll = unlist(strsplit(ll, "=", fixed=TRUE))
ll = ll[seq(from = 2, to = length(ll), by = 2)]
gff3Genes$ID = ll




blastOutput = read.table("Operon_set1_4.0_oligos.original_ACC_match.Athaliana_167_TAIR9.blast.tab", sep = "\t", header =  FALSE,stringsAsFactors=FALSE)
blastOutput[[1]] = toupper(str_match(blastOutput[[1]], "^At.g....."))

# calculate how well the thing mapped assuming 70 nt long oligo
blastOutput$match = ((blastOutput$V8-blastOutput$V7+1)*(blastOutput$V3/100))/oligoLength



# only keep the ones with the same AT.G..... name

BlastOutputOrdered = blastOutput[with(blastOutput, order(V1, -match)), ]

BlastOutputOrderedSingle = BlastOutputOrdered[ !duplicated(BlastOutputOrdered$V1), ]

BlastOutputOrderedSingleGood = BlastOutputOrderedSingle[which(BlastOutputOrderedSingle$match >0.5), ]
BlastOutputOrderedSingleGood$location = (BlastOutputOrderedSingleGood[[9]]+BlastOutputOrderedSingleGood[[10]]) /2
BlastOutputOrderedSingleGood$dir = "+"
BlastOutputOrderedSingleGood$dir[which(BlastOutputOrderedSingleGood[[9]]-BlastOutputOrderedSingleGood[[10]] > 0)] = "-"
BlastOutputOrderedSingleGood$TAIR10Name = "noHit"
head(BlastOutputOrderedSingleGood)


for(i in 1: length(BlastOutputOrderedSingleGood[[1]])){
  hits = gff3Genes[which( gff3Genes[[4]] < BlastOutputOrderedSingleGood$location[[i]]), ] 
  hits =  hits[ which(hits[[5]] > BlastOutputOrderedSingleGood$location[[i]]),]  
  hits =  hits[which(hits[[7]] == BlastOutputOrderedSingleGood$dir[[i]]),]
  hits =  hits[which(as.character(hits[[1]]) == as.character(BlastOutputOrderedSingleGood[[2]][[i]])),]
  if(length(hits[[1]]) == 1) {
    BlastOutputOrderedSingleGood$TAIR10Name[[i]] = hits$ID[[1]]
  }else if(length(hits[[1]]) > 1){
    BlastOutputOrderedSingleGood$TAIR10Name[[i]] = "many"
  }
}




BlastOutputOrderedSingleGoodWithHit = BlastOutputOrderedSingleGood[which(BlastOutputOrderedSingleGood$TAIR10Name != "noHit"), ] 
write.table(BlastOutputOrderedSingleGoodWithHit, sep = "\t", quote = FALSE, file = "Operon_set1_4.0_oligos.original_ACC_match.Athaliana_167_TAIR9.trimmed.matchGff3.blast", row.names = FALSE,col.names = FALSE)







# Load expression File 
Wellmer2004 = read.table("/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2004_stamen.expression", sep = "\t", header =  FALSE)
Wellmer2004Hits = Wellmer2004[substr(Wellmer2004[[1]], 0, 9) %in% BlastOutputOrderedSingleGoodWithHit$V1, ]
Wellmer2004Hits$TAIR10Name = paste(BlastOutputOrderedSingleGoodWithHit[ BlastOutputOrderedSingleGoodWithHit$V1 %in% substr(Wellmer2004Hits[[1]], 0, 9),16])
Wellmer2004_TAIR10 = Wellmer2004Hits[ ,c(3,2)]

write.table(Wellmer2004_TAIR10, sep = "\t", quote = FALSE, file = "/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2004_stamen.TAIR10.expression", row.names = FALSE,col.names = FALSE)



Wellmer2007 = read.table("/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2007_MS1.expression", sep = "\t", header =  FALSE)
Wellmer2007Hits = Wellmer2007[toupper(substr(Wellmer2007[[1]], 0, 9)) %in% BlastOutputOrderedSingleGoodWithHit$V1, ]
Wellmer2007Hits_TAIR10 = Wellmer2007Hits
Wellmer2007Hits_TAIR10[[1]] = paste(BlastOutputOrderedSingleGoodWithHit[ BlastOutputOrderedSingleGoodWithHit$V1 %in% toupper(substr(Wellmer2007Hits[[1]], 0, 9)),16])

results = as.matrix(Wellmer2007Hits_TAIR10[,c(2,3,4,5,6,7,8)])
Wellmer2007Hits_TAIR10$max = as.numeric(apply(results, 1, max))
Wellmer2007Hits_TAIR10$min = as.numeric(apply(results, 1, min))

write.table(Wellmer2007Hits_TAIR10, sep = "\t", quote = FALSE, file = "/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2007_MS1.TAIR10.expression", row.names = FALSE,col.names = FALSE)


Wellmer2007Hits_TAIR10_upregulated = Wellmer2007Hits_TAIR10[Wellmer2007Hits_TAIR10$max>=2, ] 
Wellmer2007Hits_TAIR10_downregulated = Wellmer2007Hits_TAIR10[Wellmer2007Hits_TAIR10$min< -2, ] 

AlvesFereria = Wellmer2007Hits_TAIR10_upregulated[,c(1,9)]


blastOutput = read.table("Operon_set1_4.0_oligos.original_ACC_match.fasta.Athaliana_CDS.blast", sep = "\t", header =  FALSE)
blastOutput = read.table("Operon_set1_4.0_oligos.original_ACC_match.fasta.Athaliana_CDS.blast", sep = "\t", header =  FALSE)

blastOutput[[1]] = toupper(str_match(blastOutput[[1]], "^At.g....."))
blastOutput[[2]] = str_match(blastOutput[[2]], "^AT.G.....")

# calculate how well the thing mapped assuming 70 nt long oligo
blastOutput$match = ((blastOutput$V8-blastOutput$V7+1)*(blastOutput$V3/100))/oligoLength

head(blastOutput, n = 50)
# only keep the ones with the same AT.G..... name


BlastOutputOrdered = blastOutput[with(blastOutput, order(V1, -match)), ]
BlastOutputOrderedSingle = BlastOutputOrdered[ !duplicated(BlastOutputOrdered$V1), ]

BlastOutputOrderedSingleSameName = BlastOutputOrderedSingle[which(BlastOutputOrderedSingle[[1]] == BlastOutputOrderedSingle[[2]]),]
BlastOutputOrderedSingleOtherName = BlastOutputOrderedSingle[which(BlastOutputOrderedSingle[[1]] != BlastOutputOrderedSingle[[2]]),]

BlastOutputOrderedSingleOtherName$compare1 = substr(BlastOutputOrderedSingleOtherName$V1, 1 , 8)
BlastOutputOrderedSingleOtherName$compare2 = substr(BlastOutputOrderedSingleOtherName$V2, 1 , 8)

BlastOutputOrderedSingleOtherNameClose = BlastOutputOrderedSingleOtherName[which(BlastOutputOrderedSingleOtherName$compare1 == BlastOutputOrderedSingleOtherName$compare2),]





SameNameBlastOutput = blastOutput[which(blastOutput[[1]] == blastOutput[[2]]),]


SameNameBlastOutputOrdered = SameNameBlastOutput[with(SameNameBlastOutput, order(V1, -match)), ]
SameNameBlastOutputOrderedSingle = SameNameBlastOutputOrdered[ !duplicated(SameNameBlastOutputOrdered$V1), ]

SameNameBlastOutputOrderedSingleGood = SameNameBlastOutputOrderedSingle[which(SameNameBlastOutputOrderedSingle$match >0.5), ]

SameNameBlastOutputOrderedSingleBad = SameNameBlastOutputOrderedSingle[which(SameNameBlastOutputOrderedSingle$match <=0.5), ]

write.table(SameNameBlastOutputOrderedSingleGood, sep = "\t", quote = FALSE, file = "Operon_set1_4.0_oligos.original_ACC_match.fasta.Athaliana_CDS.trimmed.match.blast", row.names = FALSE,col.names = FALSE)


# Load expression File 
Wellmer2004 = read.table("/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2004_stamen.expression", sep = "\t", header =  FALSE)

Wellmer2004_TAIR10 = Wellmer2004[substr(Wellmer2004[[1]], 0, 9) %in% SameNameBlastOutputOrderedSingleGood$V1,]

write.table(Wellmer2004_TAIR10, sep = "\t", quote = FALSE, file = "/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2004_stamen.TAIR10.expression", row.names = FALSE,col.names = FALSE)

Wellmer2007 = read.table("/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2007_MS1.expression", sep = "\t", header =  FALSE)

Wellmer2007_TAIR10 = Wellmer2007[toupper(substr(Wellmer2007[[1]], 0, 9)) %in% SameNameBlastOutputOrderedSingleGood$V1,]

write.table(Wellmer2007_TAIR10, sep = "\t", quote = FALSE, file = "/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2007_MS1.TAIR10.expression", row.names = FALSE,col.names = FALSE)






