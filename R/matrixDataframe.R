
rowScale<-function(data,center=TRUE,scaled=FALSE){
	data<-t(data)
	data<-t(scale(data,center=center,scale=scaled))
	return(data)
}

chead<-function(x, n=5){
	print(x[1:min(n,nrow(n)),1:min(n,ncol(n))])
}

ConvertKey<-function(keyList, tabKey,colOldKey=1,colNewKey=2){
	#@param keyList: vector of id to be converted
	#@param tabKey : dataframe of x columns with at least two contain old and new keys correspondance
	#@param colOldKey : index of old keys column in tabKey
	#@param colNewKey : index of new keys column in tabKey
	hashCorr<-tabKey[,colNewKey]
	names(hashCorr)<-tabKey[,colOldKey]
	returned<-hashCorr[keyList]
	names(returned)<-NULL
	return(as.character(returned))
}

ConvertKeyMatrix<-function(tab,tabKey,colOldKey=1,colNewKey=2,first=TRUE,dim=1,fun){
	if(dim==2) tab<-t(tab)
	keyList<-rownames(tab)
	newkey<-ConvertKey(keyList, tabKey,colOldKey,colNewKey)
	names(newkey)<-as.character(1:length(newkey))
	if(first){
		newkey<-newkey[which(!is.na(newkey))]
		newkey<-takefirst(newkey)
		tab<-tab[as.numeric(names(newkey)),]
		rownames(tab)<-newkey
	}else{
		tab<-tab[which(!is.na(newkey)),]
		newkey<-newkey[which(!is.na(newkey))]
		tab<-aggregRows(tab,newkey,fun=fun)
	}
	if(dim==2) tab<-t(tab)
	return(tab)
}

supprNAnames<-function(x ,side=1){
	#Vector case
	if(is.null(dim(x))){
		return(x[which(!is.na(names(x)))]);
	}
	#Col case
	if(side==2){
		return(x[,which(!is.na(colnames(x)))]);
	}
	#Row case
	else{
		return(x[which(!is.na(rownames(x))),]);
	}
}


#create empty dataframe
emptyDF<-function(row.names,col.names,defaultValue=NA){
	tempMat<-matrix(data = defaultValue,nrow=length(row.names),ncol = length(col.names),dimnames = list(row.names,col.names))
	data.frame(tempMat)
}


matrixCoord1D_2D<-function(x,Mat){
	matDim<-dim(Mat)
	c((x-1) %% matDim[1] +1 , (x-1) %/% matDim[1] + 1)
}


matrixFromDimnames<-function(row,col,value=0){
	matrix(value,ncol=length(col),nrow=length(row),dimnames=list(row,col))
}


reScale <- function(nonCorrected, corrected) {
	apply(nonCorrected, 1, function(x) max(x) - min(x))/
		(apply(corrected, 1, function(x) max(x)-min(x)))*
		(corrected-apply(corrected, 1, max))+apply(nonCorrected, 1, max)
}


matDist<-function(M1,M2,method="euclidean"){
	if(!is.matrix(M1)) M1<-as.matrix(M1)
	if(!is.matrix(M2)) M2<-as.matrix(M2)

	res<-matrix(apply(expand.grid(1:nrow(M1),1:nrow(M2)),1,function(row){
		as.vector(dist(rbind(M1[row[1],],M2[row[2],]),method = method))
	}),nrow =nrow(M1))

	rownames(res)<-rownames(M1)
	colnames(res)<-rownames(M2)

	res
}

matDistEuclidean<-function(M1,M2){
	if(!is.matrix(M1)) M1<-as.matrix(M1)
	if(!is.matrix(M2)) M2<-as.matrix(M2)

	res<-sqrt(do.call("+",lapply(1:ncol(M1),function(i){
		outer(M1[,i],M2[,i],"-")^2
	})))

	rownames(res)<-rownames(M1)
	colnames(res)<-rownames(M2)

	res
}


intersectionEnrichment<-function(isInGroupMat){
	universe<-nrow(isInGroupMat)
	expected<-prod(apply(isInGroupMat,2,sum)/universe)*universe
	real<-sum(apply(isInGroupMat,1,function(x) sum(x))==ncol(isInGroupMat))
	real/expected
}

formatAnnotFromMeta<-function(annotDataFrame, metaAnnot){
	require(stringr)
	colorScales=list()
	for(feature in rn(metaAnnot)){
		if(! feature %in% colnames(annotDataFrame))
			annotDataFrame[,feature]<- do.call(paste0("as.",metaAnnot[feature,"Type"]),list(x=annotDataFrame[,feature]))
		if(metaAnnot[feature,"colorScale"] != ""){
			colorScale=strsplit(str_remove_all(metaAnnot[feature,"colorScale"]," "),split = ",")[[1]]
			splitted=sapply(colorScale,strsplit,"=")
			colorScales[[feature]]<-sapply(splitted,function(x) x[2])
			names(colorScales[[feature]])<-sapply(splitted,function(x) x[1])
			if(metaAnnot[feature,"Type"]=="factor"){
				annotDataFrame[,feature]<-factor(annotDataFrame[,feature],levels = names(colorScales[[feature]]))
			}
		}
	}
	attr(annotDataFrame,"colorScales")<-colorScales
	annotDataFrame
}


subSampleColumnPerGroup<-function(mat,groupVector,n=NULL){
	lvlTable<-table(groupVector)
	lvls<-names(lvlTable)
	if(is.null(n)){
		n<-min(lvlTable)
	} else {
		if(n > min(lvlTable)) stop("n superior to minimum level freq")
	}
	subSampledData<-lapply(lvls,function(lvl){
		samples <-which(groupVector==lvl)
		selected<-sample(samples,size = n)
		mat[,selected]
	});names(subSampledData)<-lvls
	do.call("cbind",subSampledData)
}


#' Convert factors in a dataframe to strings
#'
#' @param dt A dataframe containing factor column
#'
#' @return A dataframe where all the factor column have been converted to srings
#'
#' @examples
#' a<-data.frame(x=factor(c("a","a","b")),y=1:3)
#' b<-factorAsStrings(a)
#' b$x
factorAsStrings<-function(dt){
	for(i in 1:ncol(dt)){
		if(is.factor(dt[,i])) dt[,i]<-as.character(dt[,i])
	}
	dt
}

