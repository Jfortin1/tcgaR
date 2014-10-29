library(RCurl)
library(downloader)



getTCGA <- function(cancer, datatype = c("methylation"), platform = c("27k", "450k"), verbose=FALSE, what = c("both", "normal","tumor")){
	cancer <- tolower(cancer)	
	datatype <- match.arg(datatype)
	platform <- match.arg(platform)
	what <- match.arg(what)

	if (datatype != "methylation"){
		stop("Only methylation data are supported for the moment")
	}
	rgset <- getTCGA.meth(cancer = cancer, platform = platform, verbose = verbose, what = what)
	rgset

}
	





	




	









	

			


