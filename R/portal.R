


getTCGA <- function(cancer, datatype = c("methylation", "expression"), platform = c("450k", "27k"), verbose=FALSE, n.samples=NULL){
	cancer <- tolower(cancer)	
	datatype <- match.arg(datatype)
	platform <- match.arg(platform)
	#what <- match.arg(what)

	if (datatype != "methylation" & datatype != "expression"){
		stop("Only methylation and expression data are supported for the moment")
	}
	if (datatype=="methylation"){
		object <- getTCGA.meth(cancer = cancer, platform = platform, verbose = verbose, n.samples = n.samples)
	} else if (datatype == "expression"){
		#object <- 
	}
	object
}
	





	




	









	

			


