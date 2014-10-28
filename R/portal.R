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
	
	
	getClinicalData <- function(cancer = c("acc", "blca", "brca","coad","cesc" ,"dlbc","esca", "gbm","hnsc","kich",
	 		"kirc","kirp","laml","lgg","lihc","luad","lusc", "meso","ov","paad",
	 		"pcpg","prad","read","sarc","skcm","stad","thca","ucec","ucs","uvm")){

		cancer <- match.arg(cancer)
		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		tail <- "/bcr/biotab/clin/nationwidechildrens.org_clinical_patient_"
		url <- paste0(root,cancer, tail, cancer, ".txt")
		clinical.data <- read.csv( text = getURL(url), sep="\t")
		clinical.data 
	}



	getCancers <- function(){
		d <- getURL("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/")
		d <- strsplit(d, split="\n")
		d <- unlist(d)
		d <- d[grepl("<a href",d) & 
				!grepl("tcgafiles",d) &
				!grepl("README",d) &
				!grepl("lost",d)
			  ]

		start.patt <- "<a href="
		end.patt   <- ">"

		start <- regexpr(start.patt,d) + nchar(start.patt)+1 
		stop  <- regexpr(end.patt, d) -3
		substr(d, start, stop)
	}










	

			


