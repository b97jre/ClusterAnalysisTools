#!/usr/bin/env Rscript

#source("http://bioconductor.org/biocLite.R") 
#biocLite("topGO")
#install.packages("memisc") 
#biocLite("Rgraphviz")


library(topGO)
library(memisc)
library(Rgraphviz)







getGOtermFunctions <- function(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat){
  ModelGeneGOrelationship = getModelGeneGOrelationship(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat)
  Functions = unique(ModelGeneGOrelationship$Function)
  return (Functions)
}

getModelGeneGOrelationship <- function(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat){
  if(ModelGeneGOrelationshipFileFormat == "blast2GO"){
    # SeqName Hit-Desc        GO-Group        GO-ID   Term
    # HMEL025033-PA   homeobox protein abdominal-a homolog    F       GO:0003700      sequence-specific DNA binding transcription factor activity
    ModelGeneGOrelationship <- read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Name", "Function", "GOterm", "Description")
    
  }else if(ModelGeneGOrelationshipFileFormat == "argot2"){
    # Sequence        Aspect  GO ID   Name    Total Score     Internal Confidence     Information Content
    # HMEL015695-PA   F       GO:0004222      metalloendopeptidase activity   6.728758473824388       0.5     7.3129653461121995
    
    ModelGeneGOrelationship <- read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Function", "GOterm", "Name", "Total Score", "Internal Confidence","Information Content")
  }else if(ModelGeneGOrelationshipFileFormat == "GOC"){
    # From the annotation that can be foudn here http://geneontology.org/page/download-annotations
    ModelGeneGOrelationship <- read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F,comment.char = "!")
    colnames(ModelGeneGOrelationship) <- c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GOterm","DB_reference","Evidence","WithOrFrom","Function",
                                           "Locus","DB_Object_Synonym","DB_Object_Type","taxon","Date","Assigned_by") 
  }else{
    ##Asumes that there are three colums with 1.GeneName 2.Function  3. Goterm and it containts a header
    
    # Locus        Function  GOterm  
    # HMEL015695-PA   F       GO:0004222
    ModelGeneGOrelationship <-  read.table(ModelGeneGOrelationshipFile, header=T, fill=T, sep="\t", quote="", stringsAsFactors=F)
    colnames(ModelGeneGOrelationship) <- c("Locus", "Function", "GOterm")
    
  }
  return (ModelGeneGOrelationship)
}



getTopGOdata <- function(ModelGeneGOrelationship,GOtermFunction,BackgroundTable, GOItable,GOtermCutoff,AnalysisName = "GotermAnalysis"){


  rownames(BackgroundTable) <- rownames(BackgroundTable$Locus)
  BackgroundTable = addGOIinfo(BackgroundTable,GOItable$Locus)

  #############################################################
  # Change the GOterm table into right format
  # and filter for only those present in List 
  ##############################################################
  
  ModelGeneGOrelationship <- getGOtermFormat(ModelGeneGOrelationship,Function=GOtermFunction,LocusList =BackgroundTable$Locus)
  printGOtableInfo(ModelGeneGOrelationship)
  
 
  #############################################################
  # Move the format into the way TopGO wants it. 
  ##############################################################
  
  geneList = getTopGOGenelistFormat(ModelGeneGOrelationship$Locus,BackgroundTable$Locus[BackgroundTable$GOI==1])
  
  gene2GOlist =getTopGOGene2GOrelationshipFormat(ModelGeneGOrelationship)
  
  #############################################################
  # Create a topGOdata class 
  ##############################################################
  
  
  GOtermCutoff = as.character(GOtermCutoff)
  if(GOtermFunction == "F"){ 
    topGOdata <- new("topGOdata", description=paste(AnalysisName,GOtermFunction,GOtermCutoff,sep="_"), ontology="MF", allGenes=geneList, 
                     annot = annFUN.gene2GO, gene2GO = gene2GOlist, nodeSize=0)
  }else if(GOtermFunction == "P"){
    topGOdata <- new("topGOdata", description=paste(AnalysisName,GOtermFunction,GOtermCutoff,sep="_"), ontology="BP", allGenes=geneList, 
                     annot = annFUN.gene2GO, gene2GO = gene2GOlist, nodeSize=0)
  }else if(GOtermFunction == "C"){
    topGOdata <- new("topGOdata", description=paste(AnalysisName,GOtermFunction,GOtermCutoff,sep="_"), ontology="CC", allGenes=geneList, 
                     annot = annFUN.gene2GO, gene2GO = gene2GOlist, nodeSize=0)
  }
  return (topGOdata) 
  
  
}

#statistic tests
printTopGOresults <- function(topGOdata,outName){
  test.statFis <- new("classicCount", testStatistic=GOFisherTest, name="Fisher test")
  resultFis <- getSigGroups(topGOdata, test.statFis)
  
  
  resultFisadj <- resultFis
  Fisadj <- p.adjust(score(resultFis),method="BH") 
  score(resultFisadj) <- Fisadj
  
  #  test.statKS <- new("classicScore", testStatistic=GOKSTest, name="KS tests")
  #  resultKS <- getSigGroups(topGOdata, test.statKS)
  #  resultKSadj <- resultKS
  #  KSadj <- p.adjust(score(resultKS),method="BH")
  #  score(resultKSadj) <- KSadj
  
  #  test.statElim <- new("elimCount", testStatistic=GOFisherTest, name="Fisher test", cutOff=0.01)
  #  resultElim <- getSigGroups(topGOdata, test.statElim)
  
  #  test.statWei <- new("weightCount", testStatistic=GOFisherTest, name="Fisher test", sigRatio="ratio")
  #  resultWei <- getSigGroups(topGOdata, test.statWei)
  
  #list <- list(classic=score(resultFis), KS=score(resultKS), elim=score(resultElim), weight=score(resultWei))
  #listadj <- list(classicadj=score(resultFisadj), KSadj=score(resultKSadj), elim=score(resultElim), weight=score(resultWei))
  #allRes <- GenTable(topGOdata, classicFisher = resultFis, classicKS = resultKS, elimFisher = resultElim, WeightedFisher = resultWei, orderBy="classicFisher", ranksOf="classicFisher", topNodes=30)  
  #allResadj <- GenTable(topGOdata, classicFisheradj = resultFisadj, classicKSadj = resultKSadj, elimFisher = resultElim, WeightedFisher = resultWei, orderBy="WeightedFisher", ranksOf="classicFisheradj", topNodes=30)
  
  
  
  list <- list(classic=score(resultFis))
  listadj <- list(classicadj=score(resultFisadj))
  
  allRes <- GenTable(topGOdata, classicFisher = resultFis, orderBy="classicFisher", ranksOf="classicFisher", topNodes=30)  
  allResadj <- GenTable(topGOdata, classicFisheradj = resultFisadj,  orderBy="classicFisheradj", ranksOf="classicFisheradj", topNodes=30)
  
  GOIDs = allResadj$GO.ID
  
  
  #output of data
  #saves everything in Results order within the folder one is at
  if (!file.exists("results")){
    dir.create(file.path("results"))
  }
  
  write.table(allResadj, file=paste("results/",outName,"_stat_adj.txt",sep=""), row.names=FALSE, col.names=TRUE)
  
  # pdf(paste("results/GO_graphFis_",outName,"_15nodes_adj.pdf",sep=""))
  # showSigOfNodes(topGOdata, score(resultFisadj), firstSigNodes=15, useInfo="all")
  #  dev.off()
  
  write.table(allRes, file=paste("results/",outName,"_stat.txt",sep=""), row.names=FALSE, col.names=TRUE)
  
  pdf(paste("results/P_value_distribution_",outName,".pdf",sep=""))
  par(mfrow = c(2,2))
  for (nn in names(list)){
    p.val <- list[[nn]]
    hist(p.val[p.val < 1], br = 50, xlab = "p values", 
         main = paste("Histogram for method:", nn, sep = " "))
  }
  dev.off()
  
  pdf(paste("results/GO_graphFis_",outName,"_5nodes.pdf",sep=""))
  showSigOfNodes(topGOdata, score(resultFisadj), sigForAll = FALSE,firstSigNodes=5, useInfo="all")
  dev.off()
  
  #pdf(paste("results/GO_graphFis_TEST_NEW_3_15nodes.pdf",sep=""))
  #showSigOfNodes(topGOdata, score(resultFis), firstSigNodes=15, useInfo="all")
  #dev.off()
  
  #pdf(paste("results/GO_graphWeiFis_",outName,"_15nodes.pdf",sep=""))
  #showSigOfNodes(topGOdata, score(resultWei), firstSigNodes=15, useInfo="all")
  #dev.off()
  
}

getGOtermFormat <- function(ModelGeneGOrelationship, Function="P",LocusList=ModelGeneGOrelationship$Locus,GOtermList=ModelGeneGOrelationship$GOterm){
  
  
  if(checkFormat(ModelGeneGOrelationship,c("Locus","GOterm","Function"))){
    cat("GOterm table format is correct. \n")
    #Keep only columns that are interesting for this purpose
    ModelGeneGOrelationship = ModelGeneGOrelationship[,c("Locus","GOterm","Function")]
    
    #Filter based on function 
    ModelGeneGOrelationship = ModelGeneGOrelationship[which(ModelGeneGOrelationship$Function == Function),]
    
    #Filter based on Locus 
    ModelGeneGOrelationship = ModelGeneGOrelationship[ModelGeneGOrelationship$Locus %in% LocusList,]
    
    #Filter based on GOterm 
    ModelGeneGOrelationship = ModelGeneGOrelationship[ModelGeneGOrelationship$GOterm %in% GOtermList,]
    return (ModelGeneGOrelationship)
  }else{
    stop("Something wrong with the GOterm table column names")
  }
  return (NULL)
  
}







printOrthoInfo <- function(ModelGeneOrthologsRightFormat){
  cat("This is the informtion about the orthologs\n")
  print(head(ModelGeneOrthologsRightFormat))
  cat("Nr of unique Locuses:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Locus)))
  cat("\nNr of unique orthologs:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Ortholog)))
  cat("\n\n")
}



printOrthoInfoAdvanced <- function(ModelGeneOrthologsRightFormat){
  cat("This is the informtion about the orthologs\n")
  print(head(ModelGeneOrthologsRightFormat))
  cat("Nr of unique Locuses:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Locus)))
  cat("\nNr of unique orthologs:\n")
  cat(length(ModelGeneOrthologsRightFormat$GOI==1))
  cat("\n\n")
}


printGOtableInfo <- function(ModelGeneOrthologsRightFormat){
  cat("This is the informtion about the GOtable\n")
  print(head(ModelGeneOrthologsRightFormat))
  cat("\nNr of unique Functions:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Function)))
  cat("\nNr of unique Locuses:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$Locus)))
  cat("\nNr of unique GOterms:\n")
  cat(length(unique(ModelGeneOrthologsRightFormat$GOterm)))
  cat("\n\n")
}



checkFormat <- function(dat,Mandatory){
  found = Mandatory %in% colnames(dat)
  return (all(found))
}


getRightOrthoInfo <- function(ModelGeneOrthologs,mRNA2Gene=FALSE){
  
  if(checkFormat(ModelGeneOrthologs,c("Locus","Ortholog"))){
    cat("Ortholog format is correct. Getting right Ortho info\n")
    
    #Keep only interesting columns
    ModelGeneOrthologs <- ModelGeneOrthologs[ , c("Locus","Ortholog")]
    
    #Remove rows where there is no ortholog
    ModelGeneOrthologs = ModelGeneOrthologs[ModelGeneOrthologs$Ortholog>0 , ]
    
    if(mRNA2Gene){
      ##Going from isoform to gene not always applicable
      ModelGeneOrthologs$Ortholog =getLocusFromGeneName(ModelGeneOrthologs$Ortholog)
    }
    
    return (ModelGeneOrthologs)
  }else{
    stop("Something wrong with the Ortholog column names")
  }
  return (NULL)
  
}

getLocusFromGeneName <-function(GeneName){
  #This assumes geneName to be AT1GXXXXX.1
  Raw = unlist(strsplit(GeneName,"[.]") )
  Locus = Raw[seq(1,length(Raw)-1,by=2)]
  return (Locus)
}


addGOIinfo <- function(BackgroundList, GOIList){
  BackgroundList$GOI  <- as.integer(BackgroundList$Locus %in% GOIList) 
  return(BackgroundList)
}



filterGeneList <- function(GenesWithOrthologs, SubsetList=GenesWithOrthologs$Locus, GOIList){
  #GenesWithOrthologs should be a two row data.frame  with names "Locus" "Ortholog"
  
  #First remove all GeneNames that is not in the SubsetList
  FilteredGenesWithOrthologs <- GenesWithOrthologs[GenesWithOrthologs$Locus %in% SubsetList,]
  #Add column if Gene with ortholog is interesting
  FilteredGenesWithOrthologs$GOI  <- as.integer(FilteredGenesWithOrthologs$Locus %in% GOIList) 
  return(FilteredGenesWithOrthologs)
  
}

getTopGOGenelistFormat <- function(Genes,GenesOfInterest){
  #First step
  Genes = unique(Genes)
  GenesGOI <- factor(as.integer(Genes %in% GenesOfInterest))
  names(GenesGOI) <- Genes
  return (GenesGOI)
}

getTopGOGene2GOrelationshipFormat <- function(ModelGeneGOrelationship){
  gene2GOlist_sb <- split(ModelGeneGOrelationship$GOterm,  ModelGeneGOrelationship$Locus)
  gene2GOlist_sb <- lapply(gene2GOlist_sb, unique) #remove duplicates
  return (gene2GOlist_sb)
}


getTopGO2GenerelationshipFormat <- function(ModelGeneGOrelationship){
  gene2GOlist_sb <- split(ModelGeneGOrelationship$Locus,  ModelGeneGOrelationship$GOterm)
  gene2GOlist_sb <- lapply(gene2GOlist_sb, unique) #remove duplicates
  return (gene2GOlist_sb)
}


Functions = getGOtermFunctions(ModelGeneGOrelationshipFile,ModelGeneGOrelationshipFileFormat)

DiffExprFileName = strsplit(DiffExprFile,split="/")[[1]][length(strsplit(DiffExprFile,split="/")[[1]])]

for(i in 1:length(Functions)){
  outName = paste(DiffExprFileName,ModelGeneGOrelationshipFileFormat,Functions[[i]],GOtermCutoff,sep="_")
  topGOdata = getTopGOdata(ModelGeneGOrelationshipFile = ModelGeneGOrelationshipFile,
                           ModelGeneGOrelationshipFileFormat = ModelGeneGOrelationshipFileFormat,
                           GOtermFunction = Functions[[i]],
                           DiffExprProgram = DiffExprProgram,
                           DiffExprFile = DiffExprFile, 
                           PvalueCutoff = PvalueCutoff,
                           GOtermCutoff = GOtermCutoff,
                           Determinant = determinant)
  printTopGOresults(topGOdata,outName)
}

