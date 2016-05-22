getTCGA <- function(cancer, datatype = c("methylation", "expression"), platform = c("450k", "27k"), idat=FALSE, idatDir=NULL, verbose=FALSE, n.samples=NULL){
	cancer <- tolower(cancer)	
	datatype <- match.arg(datatype)
	platform <- match.arg(platform)
	#what <- match.arg(what)
	if (idat & is.null(idatDir)){
		stop("idatDir must be specified if idat=TRUE")
	}
	if (datatype != "methylation"){
		stop("Only methylation data are supported at the moment")
	}
	if (datatype=="methylation"){
		object <- .getTCGA.meth(cancer = cancer, platform = platform, verbose = verbose, n.samples = n.samples, idat=idat, idatDir=idatDir)
	} 
	object
}


	






	




	









	

			


