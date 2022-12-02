#' Return gene id correspondance, GO species code and KEGG species code
#'
#' @param sample.species Character. Shortname of the species as described in `data("bods")`.
#' @param updateSpeciesPackage Logical. Download or update automatically the annotation org.db package corresponding to the species.
#'
#' @return A list describing specific data for the species (gene IDs, annotation package...).
#' @export
#'
#' @examples
#' library(gage)
#' data(bods)
#' bods
#' getSpeciesData("Human")
#' getSpeciesData("Mouse")
getSpeciesData<-function(sample.species="Human",updateSpeciesPackage=FALSE){
	data("bods",package = "gage")
	species<-list()
	species.data<-data.frame(bods)
	species.index<-which(species.data$species==sample.species)
	if(length(species.index)!=1) stop("Wrong species name, type \ndata(bods)\nbods[,\"species\"]\nto see available species")
	species$package<-as.character(species.data[species.index,"package"])
	species$kegg<-as.character(species.data[species.index,"kegg.code"])
	species$go<-strsplit(as.character(species$package),split = ".",fixed = TRUE)[[1]][2]
	if(updateSpeciesPackage | !(require(species$package,character.only = TRUE))){
		print(paste0("Downloading species package: ",species.data$package))
		BiocManager::install(species$package, update=FALSE)
	}
	require(species$package,character.only = TRUE)
	species$GeneIdTable<-AnnotationDbi::select(get(species$package),keys = AnnotationDbi::keys(get(species$package),"ENTREZID") , columns = c("ENTREZID","SYMBOL","ENSEMBL")) |> suppressMessages()
	species$species<-sample.species
	return(species)
}



#' Functional class scoring enrichment (fgsea algorithm)
#'
#' @param x vector or dataframe/matrix of one column. Values are used to classify the genes, example, it can be Log2(Fold-Change0). Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-rnorm(n = 4); names(x)<-c("GATA2","SOX17","KLF4","POU5F1")
#' @param corrIdGenes Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param maxSize maximum number of gene in each term
#' @param minSize Minimum number of gene in each term
#' @param customAnnot custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param returnGenes  return genes that were the most important for the enrichment of term
#' @param keggDisease Logical. Retain kegg disease term ?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' If this argument is not NULL, no additional database are downloaded.
#' @param speciesData  object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#' @param ... Additionnal parameters that are passed to fgsea
#'
#' @return
#' A dataframe with the following columns:
#' - pathway: name of the pathway/term
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' - log2err: the expected error for the standard deviation of the P-value logarithm
#' - ES: enrichment score, same as in Broad GSEA implementation
#' - NES: enrichment score normalized to mean enrichment of random samples of the same size
#' - size: number of gene in the term after removing genes not present
#' - genes (if `returnGenes`). Vector of genes of the term.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
#' resEnrich<-enrich.fcs(fcsScore,database = "kegg",species = "Human")
#' View(resEnrich)
enrich.fcs<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
											maxSize=500,minSize=2,customAnnot=NULL,returnGenes=FALSE,
											keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL,...){
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}

	if(class(x)!="numeric") stop("Values must be numeric")
	if(is.null(db_terms)) db_terms<-getDBterms(geneSym=names(x), corrIdGenes=corrIdGenes,database=database,
																							customAnnot=customAnnot,keggDisease=keggDisease,species=species)
	if(length(db_terms)==0) stop("Error, no term in any database was found")
	res<-list()
	for(db in names(db_terms)){
		res[[db]]<-suppressWarnings(fgsea::fgseaMultilevel (db_terms[[db]], x ,minSize=minSize,maxSize=maxSize,eps = 0,...))
		res[[db]]<-res[[db]][order(res[[db]]$padj),]
		res[[db]]$database<-db
		res[[db]]$leadingEdge<-NULL
		if(returnGenes) res[[db]]$genes <- db_terms[[db]][res[[db]]$pathway]
	}
	res<-do.call("rbind", res)
	res$padj<-p.adjust(res$pval,method = "BH")
	return(res)
}



#' Over Representation Analysis (enrichment, Fischer tests)
#'
#' @param x vector or dataframe/matrix of one column.
#' Values are booleans and say if gene is from the list of interest or not.
#' Genes are contained in the names/rownames of the vector/dataframe.
#' Example of valid x: x<-c(TRUE,TRUE,FALSE,FALSE); names(x)<-c("GATA2","SOX17","KLF4","POU5F1").
#' In this case, GATA2, SOX17, KLF4, POU5F1 are the universe of gene and GATA2 and SOX17 are the genes of interest
#' @param corrIdGenes  Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param minSize Minimum number of gene in each term.
#' @param maxSize Maximum number of gene in each term.
#' @param returnGenes Return genes that were the most important for the enrichment of term.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param customAnnot Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#'
#' @return
#' A dataframe with the following columns:
#' - term: name of the term
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' - nGeneOfInterest: number of gene of interest in the term.
#' - nGene: number of gene in the term after removing genes not present.
#' - genes (if `returnGenes`). Vector of genes of the term.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' vectorIsDE<-DEgenesPrime_Naive$isDE!="NONE";names(vectorIsDE)<-rownames(DEgenesPrime_Naive)
#' resEnrich<-enrich.ora(vectorIsDE,database = "kegg",species = "Human")
#' View(resEnrich)
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

	if(is.null(db_terms)) db_terms<-getDBterms(geneSym=names(x), corrIdGenes=corrIdGenes,database=database,customAnnot=customAnnot,keggDisease=keggDisease,species=species)

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


#' Compute the activation score of gene sets from an expression matrix.
#'
#' @description Perform a PCA for each gene set, from the matrix of [ genes from gene set × all samples ]. Return the first PCs as activations scores of the gene sets.
#'
#' @param expressionMatrix An expression matrix (normalized log2(x+1) counts). Genes as rows and sample as columns. If `db_terms` is not given, must be named by gene symbols.
#' @param corrIdGenes Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param scaleScores Logical. Divide expression of gene by its standard deviation before doing the PCA.
#' @param centerScores Logical. Subtract mean to gene expression before doing the PCA.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param maxSize Maximum number of gene in each term.
#' @param minSize Minimum number of gene in each term.
#' @param customAnnot Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#'
#' @return A list where each element is a database of gene set given as input.
#' For each database, contain a list of activation score (eigen), with gene sets as rows and samples as columns ;
#' and the list of contribution (or weight) to activation score of each gene per gene set.
#'
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' keggDB<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,db_terms = keggDB)
#' #same as
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,database = "kegg")
computeActivationScore<-function(expressionMatrix,corrIdGenes=NULL,scaleScores=FALSE,centerScores=TRUE,
																 database=c("kegg","reactom","goBP","goCC","goMF"),
																 maxSize=500,minSize=2,customAnnot=NULL,
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


#' Gene Set Differential Activation (GSDA)
#'
#' @param geneSetActivScore A list of database with an element "eigen" containing the matrix of gene set activation score (see what `computeActivationScore` returns).
#' If NULL, this is computed automatically from the `expressionMatrix` and the gene set database given via `db_terms` or requested via `database`.
#' @param expressionMatrix An expression matrix (normalized log2(x+1) counts). Genes as rows and sample as columns. If `db_terms` is not given, must be named by gene symbols.
#' @param colData An annotation dataframe. Each column is a feature, each row a sample. Same number of samples than in `expressionMatrix`.
#' @param contrast A vector of 3 character.
#' 1. Name of the experimental variable that have to be used for differential activation. Must be a column name of `colData`.
#' 2. Condition considered as the reference.
#' 3. Condition considered as the target group.
#' @param corrIdGenes Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param maxSize Maximum number of gene in each term.
#' @param minSize Minimum number of gene in each term.
#' @param customAnnot  Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param db_terms A list or NULL. A named list were each element is a database. Inside each database, a list terms, named by the term and containing gene vectors as gene symbols.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#'
#' @return
#' A dataframe with the following columns:
#' - term: name of the term/gene set
#' - baseMean: mean of activation score in the gene set
#' - sd: standard deviation  in the gene set
#' - log2FoldChange: Log(Log Fold Change) of activation score between the two tested groups.
#' - pval: an enrichment p-value
#' - padj: a BH-adjusted p-value
#' - database: origin of the gene set
#' - size: number of gene in the term after removing genes not present.
#' @export
#'
#' @examples
#' data("bulkLogCounts")
#' data("sampleAnnot")
#'
#' keggDB<-getDBterms(rownames(bulkLogCounts),database = "kegg")
#' geneSetActivScore<-computeActivationScore(bulkLogCounts,db_terms = keggDB)
#' resGSDA<-GSDA(geneSetActivScore = geneSetActivScore,colData = sampleAnnot,contrast = c("culture_media","T2iLGO","KSR+FGF2"),db_terms =  keggDB)
#'
#' #or
#' resGSDA<-GSDA(expressionMatrix = bulkLogCounts,colData = sampleAnnot,contrast = c("culture_media","T2iLGO","KSR+FGF2"),database = "kegg")
GSDA<-function(geneSetActivScore=NULL,expressionMatrix=NULL,colData,contrast, corrIdGenes=NULL,
							 database=c("kegg","reactom","goBP","goCC","goMF"),
							 maxSize=500,minSize=2,customAnnot=NULL,keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL){

	if(is.null(geneSetActivScore) & is.null(expressionMatrix)) stop("At least expressionMatrix or geneSetEigens must be given")

	if(is.null(db_terms)){
		db_terms<-getDBterms(geneSym=rownames(expressionMatrix), corrIdGenes=corrIdGenes,database=database,
												 customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	}

	if(is.null(geneSetActivScore)){
		geneSetActivScore<-computeActivationScore(expressionMatrix=expressionMatrix,db_terms=db_terms)
	}

	res<-list()

	for(db in names(db_terms)){
		if(is.list(geneSetActivScore[[db]])){
			eigenPerPathway<-geneSetActivScore[[db]]$eigen
		}else{
			eigenPerPathway<-geneSetActivScore[[db]]
		}

		db_terms[[db]]<-db_terms[[db]][rownames(eigenPerPathway)]
		res[[db]]<-dfres<-data.frame(term=names(db_terms[[db]]),multiLinearModel(eigenPerPathway,colData,contrast),database=db,
																 size=sapply(db_terms[[db]],length),sd=apply(eigenPerPathway,1,sd),row.names = NULL)
	}
	do.call("rbind", res)

}


#' Download database of term/pathway for enrichment analyses.
#'
#' @param geneSym A vector of gene symbols.
#' @param geneEntrez NULL or a vector of gene Entrez ID.
#' @param corrIdGenes Dataframe of gene ID correspondence where each column is a gene ID type. If not NULL `species` and `speciesData` arguments wont be used.
#' @param speciesData object returned by `getSpeciesData`. If not NULL `species` argument wont be used.
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom.
#' @param customAnnot  Custom annotation database, as a list of terms, each element contain a vector of gene symbols.
#' @param keggDisease Logical. Retain kegg disease term in kegg database?
#' @param species Character. Shortname of the species as described in `data("bods")`.
#' @param returnGenesSymbol Logical, return gene symbol in each term instead of Entrez ID.
#'
#' @return
#' A list where each element is a database that contain a list of term with the associated gene symbols.
#' @export
#' @import AnnotationDbi
#' @examples
#' data("bulkLogCounts")
#' enrichDBs<-getDBterms(rownames(bulkLogCounts),species="Human",database=c("kegg","reactom"))
getDBterms<-function(geneSym,geneEntrez=NULL, corrIdGenes=NULL, speciesData=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),customAnnot=NULL,
											keggDisease=FALSE,species="Human",returnGenesSymbol=TRUE){
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(sum(database%in%validDBs)<length(database)) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
	if(is.null(speciesData)){
		speciesData<-getSpeciesData(species)
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
			db_terms$reactom<- fgsea::reactomePathways(geneEntrez)
			db_terms$reactom<-db_terms$reactom[unique(names(db_terms$reactom))]
		}
		if("kegg"%in%database){
			kg.species <- gage::kegg.gsets(speciesData$kegg, id.type="entrez")
			db_terms$kegg<- if(keggDisease) kg.species$kg.sets else kg.species$kg.sets[kg.species$sigmet.idx]
		}
		if("go"%in%substr(database,1,2)){
			go.species <- gage::go.gsets(tolower(species), id.type="entrez")
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


#' Compute a "interest score" for a set of genes. Useful for GSEA.
#'
#' @param genes A vector of gene names
#' @param pvalues A vector of numeric corresponding to p-values.
#' @param logFoldChanges A vector of numeric corresponding to Log2(Fold-change)
#' @param logPval Logical. Compute the p-value score as `-log10(pval)` instead of `1-pval`.
#'
#' @return A vector of numeric corresponding to interest scores, named by genes.
#' @export
#'
#' @examples
#' data("DEgenesPrime_Naive")
#' fcsScore<-fcsScoreDEgenes(rownames(DEgenesPrime_Naive),DEgenesPrime_Naive$pvalue,DEgenesPrime_Naive$log2FoldChange)
fcsScoreDEgenes<-function(genes,pvalues,logFoldChanges,logPval=FALSE){
	if( sum(length(genes) == c(length(pvalues),length(logFoldChanges))) <2) stop("genes, pvalues and logPval should have the same length")
	if(logPval){
		pvalScore<- -log10(pvalues)
	}else{
		pvalScore <- 1-pvalues
	}
	pvalScore[logFoldChanges<0]<- -pvalScore[logFoldChanges<0]
	names(pvalScore)<-genes
	pvalScore
}
