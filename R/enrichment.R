
#Return gene id correspondance, GO species code and KEGG species code
getSpeciesData<-function(sample.species="Human",updateSpeciesPackage=FALSE){
	require(gage)
	require(AnnotationDbi)
	data(bods)

	species<-list()
	species.data<-data.frame(bods)
	species.index<-which(species.data$species==sample.species)
	if(len(species.index)!=1) stop("Wrong species name, type \ndata(bods)\nbods[,\"species\"]\nto see available species")
	species$package<-as.character(species.data[species.index,"package"])
	species$kegg<-as.character(species.data[species.index,"kegg.code"])
	species$go<-strsplit(as.character(species$package),split = ".",fixed = TRUE)[[1]][2]
	if(updateSpeciesPackage | !(require(species$package,character.only = TRUE))){
		require("BiocManager")
		print(paste0("Downloading species package: ",species.data$package))
		install(species$package, update=FALSE)
	}
	require(species$package,character.only = TRUE)
	suppressMessages(species$GeneIdTable<-AnnotationDbi::select(get(species$package),keys = AnnotationDbi::keys(get(species$package),"SYMBOL") , "ENTREZID","SYMBOL"))
	species$species<-sample.species
	return(species)
}


#Take a gene list, add transcription factor data and long Gene name
detailOnGenes<-function(x,tfDat,speciesDat){
	require(AnnotationDbi)
	if(is.data.frame(x) | is.matrix(x)){
		geneSym<-rn(x)
		res<-data.frame(x)
	}else{
		if(is.null(names(x))){
			res<-data.frame(row.names=x)
			geneSym<-x
		}else{
			geneSym<-names(x)
			res<-data.frame(val=x,row.names=names(x))
		}
	}
	geneNames<-select(get(speciesDat$package),geneSym,"GENENAME","SYMBOL")
	res$LongName<-ConvertKey(geneSym,tabKey = geneNames,colOldKey = "SYMBOL",colNewKey = "GENENAME")
	genesTF<-intersect(rn(tfDat),geneSym)
	res$TFdegree<-0
	res[genesTF,"TFdegree"]<-tfDat[genesTF,"tf_degree"]
	return(res)
}


# 3rd generation enrichment (fgsea algorithm)
#' @param x : vector or dataframe/matrix of one column. Values are used to classify the genes, example, it can be Log2(Fold-Change0). Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-rnorm(n = 4); names(x)<-c("GATA2","SOX17","KLF4","POU5F1")
#' @param corrIdGenes : dataframe of genes id (used to convert genes), automatically computed if not provided
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param minSize : mininmum number of gene in each term
#' @param maxSize : maximum number of gene in each term
#' @param nperm : number of permutation in the GSEA algorithm
#' @param customAnnot : custom annotation database, as a list af gene symbols, named by annotations name
#' @param returnLeadingEdge : return genes that were the most important for the enrichment of term
#' @param keggDisease : retain kegg disease term ?
#' @param species : species, example: "Rat", "Mouse", "Human"
#' @param db_terms : precomputed list of term database, automatically computed if not provided
#' @param ... : additionnal parameters that are passed to fgsea
#' @param speciesData : result of getSpeciesData2 function, automatically gathered if not provided
enrich.fcs<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
											maxSize=500,minSize=2,customAnnot=NULL,returnGenes=FALSE,
											keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL,...){
	require(fgsea)

	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}

	if(class(x)!="numeric") stop("Values must be numeric")
	if(is.null(db_terms)) db_terms<-getDBterms2(geneSym=names(x), corrIdGenes=corrIdGenes,database=database,
																							customAnnot=customAnnot,keggDisease=keggDisease,species=species)
	if(length(db_terms)==0) stop("Error, no term in any database was found")
	res<-list()
	for(db in names(db_terms)){
		res[[db]]<-fgseaMultilevel (db_terms[[db]], x ,minSize=minSize,maxSize=maxSize,eps = 0,...)
		res[[db]]<-res[[db]][order(res[[db]]$padj),]
		res[[db]]$database<-db
		res[[db]]$leadingEdge<-NULL
		if(returnGenes) res[[db]]$genes <- db_terms[[db]][res[[db]]$pathway]
	}
	res<-do.call("rbind", res)
	res$padj<-p.adjust(res$pval,method = "BH")
	return(res)
}

# 2nd generation enrichment (fisher algorithm)
#' @param x : vector or dataframe/matrix of one column. Values are booleans and say if gene is from the list of interest or not. Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-c(TRUE,TRUE,FALSE,FALSE); names(x)<-c("GATA2","SOX17","KLF4","POU5F1"). In this case, GATA2, SOX17, KLF4, POU5F1 are the universe of gene and GATA2 and SOX17 are the genes of interest
#' @param corrIdGenes : dataframe of genes id (used to convert genes), automatically computed if not provided
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param minSize : mininmum number of gene in each term
#' @param maxSize : maximum number of gene in each term
#' @param customAnnot : custom annotation database, as a list af gene symbols, named by annotations name
#' @param returnGenes : return genes of interest that are in the term
#' @param keggDisease : retain kegg disease term ?
#' @param db_terms : precomputed list of term database, automatically computed if not provided
#' @param species : species, example: "Rat", "Mouse", "Human"
#' @param speciesData : result of getSpeciesData2 function, automatically gathered if not provided
enrich.ora<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
												minSize=2,maxSize=500,returnGenes=FALSE, keggDisease=FALSE,species="Human",
												customAnnot=NULL,db_terms=NULL,speciesData=NULL){
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(sum(database%in%validDBs)==0) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")

	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}

	if(class(x)!="logical") stop("Values must be logical (TRUE or FALSE)")
	if(class(names(x))!="character") stop("Values must be named with genes symbol")

	if(is.null(db_terms)) db_terms<-getDBterms2(geneSym=names(x), corrIdGenes=corrIdGenes,database=database,customAnnot=customAnnot,keggDisease=keggDisease,species=species)

	nInterest<-length(which(x))
	nNotInterest<-length(which(!x))

	results<-list()
	for(db in names(db_terms)){
		len_term<-sapply(db_terms[[db]],length)
		db_terms[[db]]<-db_terms[[db]][len_term>=minSize & len_term<=maxSize]

		nGeneByterm<-sapply(db_terms[[db]],length)
		nGeneOfInterestByterm<-sapply( db_terms[[db]],function(term){
			return(length(which(x[term])))
		})
		results[[db]]<-data.frame(row.names = names(db_terms[[db]]))
		results[[db]]$term <- names(db_terms[[db]])
		results[[db]]$pval<-phyper(q = nGeneOfInterestByterm-0.5, m = nInterest,n = nNotInterest, k = nGeneByterm, lower.tail=FALSE)
		results[[db]]$nGeneOfInterest<-nGeneOfInterestByterm
		results[[db]]$nGene<-nGeneByterm
		results[[db]]$database<-db
		if(returnGenes){
			results[[db]]$genes<- db_terms[[db]]
		}
	}
	results<-do.call("rbind", results)
	results$padj<-p.adjust(results$pval,method = "BH")
	return(results)
}


computeActivationScore<-function(expressionMatrix,corrIdGenes=NULL,scaleScores=FALSE,centerScores=TRUE,
																 database=c("kegg","reactom","goBP","goCC","goMF"),
																 maxSize=500,minSize=2,nperm=1000,customAnnot=NULL,
																 keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL){

	if(!class(expressionMatrix)[1]%in%c("data.frame","matrix")) stop("expressionMatrix should be a matrix or a dataframe")
	if(class(rownames(expressionMatrix))!="character") stop("rows of expression matrix should be named with genes symbol")
	if(is.null(db_terms)){
		db_terms<-getDBterms(geneSym=rownames(expressionMatrix), corrIdGenes=corrIdGenes,database=database,
												 customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	}

	if(length(db_terms)==0) stop("Error, no term in any database was found")

	lapply(db_terms,function(database){
		database<-lapply(database,function(genesOfTerm) intersect(genesOfTerm,rownames(expressionMatrix)))
		nGenePerTerm<-sapply(database,length)
		database<-database[nGenePerTerm>minSize & nGenePerTerm<maxSize]
		resPerPathway<-lapply(database,function(genesOfTerm){
			eigengenes(exprMatrix = expressionMatrix,genes = genesOfTerm,returnContribution = TRUE,scale = scaleScores,center = centerScores)
		})
		list(
			eigen=t(sapply(resPerPathway,function(term) term$eigen)),
			contribution = lapply(resPerPathway,function(term) term$contribution)
		)
	})
}


GSDA<-function(geneSetEigens=NULL,expressionMatrix=NULL,colData,contrast, corrIdGenes=NULL,
							 database=c("kegg","reactom","goBP","goCC","goMF"),
							 maxSize=500,minSize=2,customAnnot=NULL,keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL){

	if(is.null(geneSetEigens) & is.null(geneSetEigens)) stop("At least expressionMatrix or geneSetEigens miiust be given")

	if(is.null(db_terms)){
		db_terms<-getDBterms(geneSym=rownames(expressionMatrix), corrIdGenes=corrIdGenes,database=database,
												 customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	}

	if(is.null(geneSetEigens)){
		geneSetEigens<-computeActivationScore(expressionMatrix=expressionMatrix,db_terms=db_terms)
	}

	res<-list()

	for(db in names(db_terms)){
		if(is.list(geneSetEigens[[db]])){
			eigenPerPathway<-geneSetEigens[[db]]$eigen
		}else{
			eigenPerPathway<-geneSetEigens[[db]]
		}

		db_terms[[db]]<-db_terms[[db]][rownames(eigenPerPathway)]
		res[[db]]<-dfres<-data.frame(term=names(db_terms[[db]]),multiLinearModel(eigenPerPathway,colData,contrast),database=db,
																 size=sapply(db_terms[[db]],length),sd=apply(eigenPerPathway,1,sd),row.names = NULL)
	}
	do.call("rbind", res)

}


getDBterms2<-function(geneSym,geneEntrez=NULL, corrIdGenes=NULL, speciesData=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),customAnnot=NULL,
											keggDisease=FALSE,species="Human",returnGenesSymbol=TRUE){
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(!(combineLogical(database%in%validDBs))) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
	if(is.null(speciesData)){
		speciesData<-getSpeciesData2(species)
	}else{
		species<-speciesData$species
	}
	if(is.null(corrIdGenes)) corrIdGenes<-speciesData$GeneIdTable
	options(warn=-1)
	if(is.null(geneEntrez)){
		geneEntrez<-ConvertKey(geneSym,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
		geneEntrez<-geneEntrez[!is.na(geneEntrez)]
	}
	db_terms<-list()
	if(is.list(customAnnot)){
		db_terms$custom<-lapply(customAnnot,function(x){
			new_x<-ConvertKey(x,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
			new_x[!is.na(new_x)]
		})
	}
	if(!(length(database)<=1 & database[1]=="custom")){
		if("reactom"%in%database){
			require("fgsea")
			require("reactome.db")
			db_terms$reactom<- reactomePathways(geneEntrez)
			db_terms$reactom<-db_terms$reactom[unique(names(db_terms$reactom))]
		}
		if("kegg"%in%database){
			require("gage")
			kg.species <- kegg.gsets(speciesData$kegg, id.type="entrez")
			db_terms$kegg<- if(keggDisease) kg.species$kg.sets else kg.species$kg.sets[kg.species$sigmet.idx]
		}
		if("go"%in%substr(database,1,2)){
			require("gage")
			go.species <- go.gsets(tolower(species), id.type="entrez")
			if("goBP"%in%database) db_terms$goBP<-go.species$go.sets[go.species$go.subs$BP]
			if("goMF"%in%database) db_terms$goMF<-go.species$go.sets[go.species$go.subs$MF]
			if("goCC"%in%database) db_terms$goCC<-go.species$go.sets[go.species$go.subs$CC]
		}
	}
	options(warn=0)

	if(returnGenesSymbol){
		lapply(db_terms,function(db) lapply(db,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL"))
	}else{
		db_terms
	}
}


calConsensusRanking<-function(genes,pvalues,logFoldChanges){
	InvPvalues <- 1-pvalues;
	InvPvalues[logFoldChanges<0]<- -InvPvalues[logFoldChanges<0]
	names(InvPvalues)<-genes
	InvPvalues
}
