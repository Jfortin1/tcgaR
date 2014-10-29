library(RCurl)
library(downloader)



getTCGA <- function(cancer, datatype = c("methylation", "expression"), platform = c("27k", "450k"), verbose=FALSE){
	cancer <- tolower(cancer)	
	datatype <- match.arg(datatype)
	platform <- match.arg(platform)
	#what <- match.arg(what)

	if (datatype != "methylation" & datatype != "expression"){
		stop("Only methylation and expression data are supported for the moment")
	}
	if (datatype=="methylation"){
		object <- getTCGA.meth(cancer = cancer, platform = platform, verbose = verbose)
	} else if (datatype == "expression"){
		#object <- 
	}
	object
}
	





	




	









	

			


