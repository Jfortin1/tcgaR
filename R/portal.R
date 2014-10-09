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

	# download.tcga <- function(datatype = "methylation", platform = c("27k", "450k"), cancers, 
	# 	level = c(1,2,3), mappings = TRUE, destdir = getwd()){

	# 	require(RCurl)
	# 	require(downloader)

	# 	if (level %in% c(2,3)){
	# 		stop("Levels 2 and 3 are not supported for the moment")
	# 	}

	# 	if (datatype != "methylation"){
	# 		stop("Only methylation data are supported for the moment")
	# 	}

	# 	if (!(platform %in% c("27k", "450k") )){
	# 		stop("Platform must be 27k or 450k")
	# 	}

	# 	availableCancers <- cancersList(datatype = dataype, platform = platform)
	# 	if (sum(!(cancers %in% availableCancers)) > 0 ){
	# 		stop("Methylation data are not available for ones of the cancers. Please use cancersList() to obtain the list of available data")
	# 	}
	# 	cat("[tcga.download] Fetching the available data \n]")
	# 	links <- extract.links(datatype = datatype, platform = platform, cancers = cancers)

	# 	if (mappings){
	# 		cat("[tcga.download] Obtaining the phenotypic data \n")
	# 		mappings <- extract.mappings(datatype = datatype, platform = platform, cancers = cancers)
	# 	} else {
	# 		mappings = NULL
	# 	}
	# 	cat("[tcga.download] Downloading the samples \n")
	# 	download.meth.idat(links = links, destdir = destdir)

	# 	return(mappings)
		
	# 	# Add a function that greps the names of the samples downloaded


	# }



	
	# download.meth.idat <- function(links, destdir){
	# 	m <- length(links)
	# 	for (i in 1:m){
	# 		url <- links[i]
	# 		cancer <- substr(url, regexpr("tumor", url) +6, regexpr("cgcc", url)-2 )
	# 		if (platform == "27k"){
	# 			start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1")
	# 		} else {
	# 			start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1")
	# 		}
			
	# 		stop.patt  <- ".tar.gz"
	# 		start <- regexpr(start.patt,url)[1] + nchar(start.patt) + 1
	# 		stop  <- regexpr(stop.patt, url)[1] -1
	# 		name <- paste0(start.patt, substr(url, start, stop), stop.patt)

	# 		download(url, destfile=file.path(destdir, name))
	# 		print(i)
	# 	}

	# }



	# # In the future the list must come from the website
	# cancersList <- function(datatype = "methylation", platform = c("27k", "450k")){

	# 	if (datatype != "methylation"){
	# 		stop("Only methylation data are supported for the moment")
	# 	}

	# 	if (!(platform %in% c("27k", "450k") )){
	# 		stop("Platform must be 27k or 450k")
	# 	}


	# 	if (platform == "27k"){
	# 		cancers <- c("brca", "coad", "gbm", "kirc", "kirp", "laml", 
	# 			"luad", "lusc", "ov", "read", "stad", "ucec")

	# 	} else {
	# 		cancers <- c("acc", "blca", "brca","coad","cesc", "dlbc","esca", "gbm","hnsc","kich",
	# 			"kirc","kirp","laml","lgg","lihc","luad","lusc", "meso","ov","paad",
	# 			"pcpg","prad","read","sarc","skcm","stad","thca","ucec","ucs","uvm")
	# 	}
	# 	return(cancers)


	# }




	# # This function extracts the links to download the tar.gz files containing the samples:
	# extract.links <- function(datatype = "datatype", platform=c("27k", "450k"), cancers){

	
	# 	m <- length(cancers)

	# 	root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
	# 	if (platform == "27k"){
	# 		tail="/cgcc/jhu-usc.edu/humanmethylation27/methylation/"
	# 	} else {
	# 		tail="/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
	# 	}
		
	# 	links <- vector("list", m)

	# 	for (kk in 1:m){

	# 		extract.version <- function(x, cancer){
	# 			if (platform == "27k"){
	# 				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1")
	# 			} else {
	# 				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1")
	# 			}
	# 			stop.patt  <- ".tar.gz"
	# 			start <- regexpr(start.patt,x)[1] + nchar(start.patt) + 1
	# 			stop  <- regexpr(stop.patt, x)[1] -1
	# 			substr(x, start, stop)
	# 		}

	# 		cancer <- cancers[kk]
	# 		url <- paste0(root,cancer,tail)

	# 		d <- getURL(url)
	# 		d <- strsplit(d, split="\n")
	# 		d <- unlist(d)
	# 		if (platform == "27k"){
	# 			d <- d[grepl("HumanMethylation27.Level_1",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
	# 		} else {
	# 			d <- d[grepl("HumanMethylation450.Level_1",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
	# 		}
			

	# 		versions <- unlist(lapply(as.list(d),function(x) extract.version(x, cancer=cancer)))
	# 		versions <- unlist(strsplit(versions, split="[.]"))
	# 		versions <- matrix(as.numeric(versions), ncol=3, byrow=TRUE)

	# 		o <- order(versions[,1], versions[,2], versions[,3])
	# 		versions <- versions[o, ,drop=FALSE]

	# 		ids <- unique(versions[,1])
	# 		retained.versions <- matrix(NA, length(ids), 3)
	# 		for (i in 1:length(ids)){
	# 			indices <- which(versions[,1]==ids[i])
	# 			last.index <- indices[length(indices)]
	# 			retained.versions[i,] <- versions[last.index,]
	# 		}

	# 		versions <- paste(retained.versions[,1], 
	# 						  retained.versions[,2],
	# 						  retained.versions[,3], sep="."
	# 					)

	# 		if (platform=="27k") {
	# 			files <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1.",versions,".tar.gz")
	# 		} else {
	# 			files <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1.",versions,".tar.gz")
	# 		}
			
	# 		links[[kk]] <- paste0(url,files)
	# 		#print(kk)

	# 	}
	# 	names(links) <- cancers
	# 	unlist(links)

	# }






	# extract.mappings <- function(datatype="methylation", platform=c("27k", "450k"), cancers){

	# 	m <- length(cancers)

	# 	root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
	# 	if (platform == "27k"){
	# 		tail="/cgcc/jhu-usc.edu/humanmethylation27/methylation/"
	# 	} else {
	# 		tail="/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
	# 	}


	# 	m <- length(cancers)
	# 	mappings <- vector("list", m)

	# 	for (kk in 1:m){

	# 		cancer <- cancers[kk]
	# 		url <- paste0(root,cancer,tail)

	# 		d <- getURL(url)
	# 		d <- strsplit(d, split="\n")
	# 		d <- unlist(d)
	# 		d <- d[grepl("aux",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]

	# 		extract.version <- function(x, cancer){
	# 			if (platform=="27k"){
	# 				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.aux.")	
	# 			} else {
	# 				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.aux.")
	# 			}
				
	# 			stop.patt  <- ".tar.gz"
	# 			start <- regexpr(start.patt,x)[1] + nchar(start.patt) 
	# 			stop  <- regexpr(stop.patt, x)[1] -1
	# 			substr(x, start, stop)
	# 		}

	# 		versions <- unlist(lapply(as.list(d),function(x) extract.version(x, cancer=cancer)))
	# 		versions <- unlist(strsplit(versions, split="[.]"))
	# 		versions <- matrix(as.numeric(versions), ncol=3, byrow=TRUE)

	# 		o <- order(versions[,1], versions[,2], versions[,3])
	# 		versions <- versions[o, ,drop=FALSE]

	# 		ids <- unique(versions[,1])
	# 		retained.versions <- matrix(NA, length(ids), 3)
	# 		for (i in 1:length(ids)){
	# 			indices <- which(versions[,1]==ids[i])
	# 			last.index <- indices[length(indices)]
	# 			retained.versions[i,] <- versions[last.index,]
	# 		}

	# 		versions <- paste(retained.versions[,1], 
	# 						  retained.versions[,2],
	# 						  retained.versions[,3], sep="."
	# 					)

	# 		if (platform=="27k"){
	# 			dir  <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.aux.",versions,"/")
	# 		} else {
	# 			dir  <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.aux.",versions,"/")
	# 		}
			
	# 		csv.file  <- paste0(url, dir, toupper(cancer),".mappings.csv")
	# 		mappings[[kk]] <- read.csv(text = getURL(csv.file))
	# 		#print(kk)

	# 	}

	# 	# Merging of the mappings: 
	# 	names <- names(mappings[[1]])
	# 	mappings.all <- mappings[[1]]
	# 	n <- length(mappings)
	# 	if (n>1){
	# 		for (j in 2:n){
	# 			map <- mappings[[j]]
	# 			if (!("plate" %in% names(map))){
	# 				n <- nrow(map)
	# 				map$plate <- rep("-",n)
	# 				map$well  <- rep("-",n)
	# 			}
	# 			map <- map[,match(names, names(map))]
	# 			mappings.all <- rbind(mappings.all, map)
	# 		}
	# 	}
		
	# 	mappings <- mappings.all
	# 	mappings
	# }






	# extract.clinical.data <- function(cancers) {

	# 		# To download the mappings:
	# 		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
	# 		tail <- "/bcr/biotab/clin/nationwidechildrens.org_clinical_patient_"

	# 		m <- length(cancers)
			
	# 		clinical.data <- vector("list", m)

	# 		for (kk in 1:m){

	# 			cancer <- cancers[kk]
	# 			url <- paste0(root,cancer, tail, cancer, ".txt")

	# 			clinical.data[[kk]] <- read.csv( text = getURL(url), sep="\t")
	# 			#print(kk)

	# 		}
	# 		names(clinical.data) <- cancers
	# 		clinical.data
	# }










	

			


