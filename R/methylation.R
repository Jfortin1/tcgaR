#library(RCurl)
#library(downloader)
#library(minfi)
#library(illuminaio)
#library(methylumi)


	
	#mappings.meth <- .getMethMappings("brca", "450k")
	#mappings.rna  <- .getRNAMappings("brca", "genes")	
	#common.samples <- intersect(mappings.meth$sample.id, mappings.rna$sample.id)
	

	#new <- mappings.meth[duplicated(mappings.meth$sample.id),]






	getTCGA.number <- function(cancer, platform = c("27k", "450k")){
		platform <- match.arg(platform)
		cancer <- tolower(cancer)
		filenames <- tcgaR:::.getIdatNames(cancer = cancer, platform = platform)
		n <- length(filenames[[1]])
		cat(paste0("[tcga.meth] ", n," samples have been found \n"))
	}

	





	getTCGA.meth <- function(cancer, platform = c("27k", "450k"), what = c("both", "normal", "cancer"), verbose = FALSE, n.samples=NULL){
		platform <- match.arg(platform)
		what <- match.arg(what)
		cancer <- tolower(cancer)

		# Let's see if the cancer exist:
		doesItExist <- .cancer.exists(cancer = cancer, platform = platform)

		filenames <- .getIdatNames(cancer = cancer, platform = platform)
		if (!is.null(n.samples)){
			filenames <- filenames[1:n.samples,, drop="FALSE"]
		}
		mappings  <- .getMethMappings(cancer = cancer, platform = platform)
		mappings  <- mappings[match(filenames$idat.name, mappings$barcode),]
		if (what == "normal"){
			retained.samples <- mappings$barcode[mappings$tissue == "Matched Normal"]
		} else if (what == "tumor"){
			retained.samples <- mappings$barcode[mappings$tissue != "Matched Normal" & mappings$tissue != "Cell Line Control"]
		} else {
			retained.samples <- mappings$barcode
		}
		mappings <- mappings[match(retained.samples, mappings$barcode),]
		indices <- match(retained.samples, filenames[[2]])
		filenames[[1]] <- filenames[[1]][indices]
		filenames[[2]] <- filenames[[2]][indices]
		n <- length(filenames[[1]])
		cat(paste0("[tcga.meth] ", n," samples have been found \n"))

		if (platform=="450k"){
			cat("[getTCGA.meth] Constructing the RGChannelSet \n")
			object <- read.450k.con(basenames = filenames$idat.name, con = filenames$idat.con, verbose = verbose)
			pData(object) <- mappings
		} else {
			cat("[getTCGA.meth] Constructing the MethyLumi object \n")
			object <- read.27k.con(barcodes = filenames$idat.name, con = filenames$idat.con)
			pData(object) <- mappings
		}
		object
	}

	





	.getMethMappings <- function(cancer , platform=c("27k","450k")) {
		cancer <- tolower(cancer)
		platform <- match.arg(platform)
		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		if (platform == "27k"){
			tail="/cgcc/jhu-usc.edu/humanmethylation27/methylation/"
		} else {
			tail="/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
		}


		url <- paste0(root,cancer,tail)

		d <- getURL(url)
		d <- strsplit(d, split="\n")
		d <- unlist(d)
		d <- d[grepl("aux",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]

		extract.version <- function(x, cancer){
			if (platform=="27k"){
				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.aux.")	
			} else {
				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.aux.")
			}
			
			stop.patt  <- ".tar.gz"
			start <- regexpr(start.patt,x)[1] + nchar(start.patt) 
			stop  <- regexpr(stop.patt, x)[1] -1
			substr(x, start, stop)
		}

		versions <- unlist(lapply(as.list(d),function(x) extract.version(x, cancer=cancer)))
		versions <- unlist(strsplit(versions, split="[.]"))
		versions <- matrix(as.numeric(versions), ncol=3, byrow=TRUE)

		o <- order(versions[,1], versions[,2], versions[,3])
		versions <- versions[o, ,drop=FALSE]

		ids <- unique(versions[,1])
		retained.versions <- matrix(NA, length(ids), 3)
		for (i in 1:length(ids)){
			indices <- which(versions[,1]==ids[i])
			last.index <- indices[length(indices)]
			retained.versions[i,] <- versions[last.index,]
		}

		versions <- paste(retained.versions[,1], 
						  retained.versions[,2],
						  retained.versions[,3], sep="."
					)

		if (platform=="27k"){
			dir  <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.aux.",versions,"/")
		} else {
			dir  <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.aux.",versions,"/")
		}
		
		csv.file  <- paste0(url, dir, toupper(cancer),".mappings.csv")
		mappings <- read.csv(text = getURL(csv.file))
	
		if ("X" %in% colnames(mappings)){
			mappings[,"X"] <- NULL
		}
		#barcodes <- .processBarcodes(mappings$TCGA.ID)
		#mappings <- cbind(mappings, barcodes)
		#rownames(mappings) <- mappings$sample.id
		mappings
	}













	.getIdatNames <- function(cancer, platform=c("27k", "450k")){
		cancer <- tolower(cancer)
		platform <- match.arg(platform)
		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		if (platform == "27k"){
			tail="/cgcc/jhu-usc.edu/humanmethylation27/methylation/"
		} else {
			tail="/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
		}

		extract.version <- function(x, cancer){
				if (platform == "27k"){
					start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1")
				} else {
					start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1")
				}
				stop.patt  <- ".tar.gz"
				start <- regexpr(start.patt,x)[1] + nchar(start.patt) + 1
				stop  <- regexpr(stop.patt, x)[1] -1
				substr(x, start, stop)
		}	

		url <- paste0(root,cancer,tail)
		d <- getURL(url)
		d <- strsplit(d, split="\n")
		d <- unlist(d)

		if (platform == "27k"){
			d <- d[grepl("HumanMethylation27.Level_1",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
		} else {
			d <- d[grepl("HumanMethylation450.Level_1",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
		}

		versions <- unlist(lapply(as.list(d),function(x) extract.version(x, cancer=cancer)))
		versions <- unlist(strsplit(versions, split="[.]"))
		versions <- matrix(as.numeric(versions), ncol=3, byrow=TRUE)

		o <- order(versions[,1], versions[,2], versions[,3])
		versions <- versions[o, ,drop=FALSE]

		ids <- unique(versions[,1])
		retained.versions <- matrix(NA, length(ids), 3)
		for (i in 1:length(ids)){
			indices <- which(versions[,1]==ids[i])
			last.index <- indices[length(indices)]
			retained.versions[i,] <- versions[last.index,]
		}

		versions <- paste(retained.versions[,1], 
						  retained.versions[,2],
						  retained.versions[,3], sep="."
					)

		if (platform=="27k") {
				subdirs <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1.",versions,"/")
		} else {
				subdirs <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1.",versions,"/")
		}


		idat.con <- c()
		idat.name <- c()

		extract.idat.filename <- function(x){
			if (platform == "27k"){
				start.patt <- 'a href=\\\"'
			} else {
				start.patt <- 'a href=\\\"'
			}
			stop.patt  <- ".idat"
			start <- regexpr(start.patt,x)[1] + nchar(start.patt) -1
			stop  <- regexpr(stop.patt, x)[1] -5
			substr(x, start, stop)
		}


		for (kk in 1:length(subdirs)){
			url <- paste0(root,cancer,tail, subdirs[kk])
			d <- getURL(url)
			d <- strsplit(d, split="\n")
			d <- unlist(d)
			d <- d[grepl(".idat",d)]

			d <- unique(extract.idat.filename(d))
			n <- length(d)
			idat.con <- c(idat.con, rep(url, n))
			idat.name <- c(idat.name, d)
			#print(kk)
			Sys.sleep(2) # Otherwise TCGA Portal complains. 
		}
		return(list(idat.con = idat.con, idat.name = idat.name))

	}















	read.450k.con <- function(basenames, con, extended = FALSE, verbose = FALSE) {
		td = tempdir()
	    basenames <- sub("_Grn\\.idat$", "", basenames)
	    basenames <- sub("_Red\\.idat$", "", basenames)
	    G.files <- paste(basenames, "_Grn.idat", sep = "")
	    R.files <- paste(basenames, "_Red.idat", sep = "")
	    G.con   <- file.path(con, G.files)
	    R.con   <- file.path(con, R.files)
	    names(G.files) <- basename(basenames)
	    names(R.files) <- basename(basenames)
	 	

	    stime <- system.time({
	        G.idats <- lapply(G.con, function(xx) {
	            if(verbose) cat("[read.450k] Reading", basename(xx), "\n")
	            tf = tempfile(tmpdir=td, fileext=".idat")
	        	download(xx, destfile=tf, quiet=TRUE)
	            readIDAT(tf)
	        })
	        R.idats <- lapply(R.con, function(xx) {
	            if(verbose) cat("[read.450k] Reading", basename(xx), "\n")
	            tf = tempfile(tmpdir=td, fileext=".idat")
	        	download(xx, destfile=tf, quiet=TRUE)
	            readIDAT(tf)
	        })
	    })[3]
	  
	    ptime1 <- proc.time()
	    GreenMean <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[, "Mean"]))
	    RedMean <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[, "Mean"]))
	    if(extended) {
	        GreenSD <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[, "SD"]))
	        RedSD <- do.call(cbind, lapply(R.idats, function(xx) xx$Quants[, "SD"]))
	        NBeads <- do.call(cbind, lapply(G.idats, function(xx) xx$Quants[, "NBeads"]))
	    }
	    ptime2 <- proc.time()
	    stime <- (ptime2 - ptime1)[3]
	    ptime1 <- proc.time()
	    if(extended) {
	        out <- new("RGChannelSetExtended", Red = RedMean, Green = GreenMean,
	                   RedSD = RedSD, GreenSD = GreenSD, NBeads = NBeads)
	    } else {
	        out <- new("RGChannelSet", Red = RedMean, Green = GreenMean)
	    }
	    featureNames(out) <- rownames(G.idats[[1]]$Quants)
	    annotation(out) <- c(array = "IlluminaHumanMethylation450k", annotation = minfi:::.default.450k.annotation)
	    ptime2 <- proc.time()
	    stime <- (ptime2 - ptime1)[3]
	    out
	}












	read.27k.con <- function(barcodes, con){
		td = tempdir()

		myIDATsToDFs <- function(barcodes, con, fileExts=list(Cy3="Grn.idat", Cy5="Red.idat")) { 
			names(barcodes) = as.character(barcodes)
			#listOfDFs = lapply(barcodes, methylumi:::IDATtoDF, fileExts=fileExts, idatPath=getwd())
			listOfDFs <- list()
			n <- length(barcodes)
			for (kk in 1:n){
				listOfDFs[[kk]] <- myIDATtoDF(barcode = barcodes[kk], con = con[kk])
			}
			names(listOfDFs) = as.character(barcodes)
			return(listOfDFs)
		} 

		myIDATtoDF <- function(barcode,fileExts=list(Cy3="Grn.idat", Cy5="Red.idat"), con) { 
		  processed = lapply(fileExts, function(chan) {
		  	con.file <- file.path(con, paste(barcode, chan, sep="_"))
		  	tf <- tempfile(tmpdir = td, fileext=".idat")
		  	download(con.file, destfile=tf, quiet=TRUE)
		    dat = readMethyLumIDAT(tf)
		    return(list(Quants=as.data.frame(dat$Quants), 
		                RunInfo=dat$RunInfo,
		                ChipType=dat$ChipType))
		  })
		  probe.data = as.data.frame(lapply(processed, function(x) x[['Quants']]))
		  attr(probe.data, 'RunInfo') = processed[[1]][['RunInfo']]
		  attr(probe.data, 'ChipType') = processed[[1]][['ChipType']]
		  return(probe.data)
		} 


		mlumi <- methylumi:::NChannelSetToMethyLumiSet(
			methylumi:::DFsToNChannelSet(
			myIDATsToDFs(barcodes = barcodes, con = con), 
			IDAT=TRUE), n=FALSE, oob=TRUE)

		return(mlumi[ sort(featureNames(mlumi)), ])
	}












	# Old function:
	# This function extracts the links to download the tar.gz files containing the samples:
	extract.links <- function(platform=c("27k", "450k")){

		if (platform == "27k"){
			cancers <- c("brca", "coad", "gbm", "kirc", "kirp", "laml", 
				"luad", "lusc", "ov", "read", "stad", "ucec")

		} else {
			cancers <- c("acc", "blca", "brca","coad","cesc", "dlbc","esca", "gbm","hnsc","kich",
				"kirc","kirp","laml","lgg","lihc","luad","lusc", "meso","ov","paad",
				"pcpg","prad","read","sarc","skcm","stad","thca","ucec","ucs","uvm")
		}

		m <- length(cancers)

		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		if (platform == "27k"){
			tail="/cgcc/jhu-usc.edu/humanmethylation27/methylation/"
		} else {
			tail="/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
		}
		
		links <- vector("list", m)

		for (kk in 1:m){

			extract.version <- function(x, cancer){
				if (platform == "27k"){
					start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1")
				} else {
					start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1")
				}
				stop.patt  <- ".tar.gz"
				start <- regexpr(start.patt,x)[1] + nchar(start.patt) + 1
				stop  <- regexpr(stop.patt, x)[1] -1
				substr(x, start, stop)
			}

			cancer <- cancers[kk]
			url <- paste0(root,cancer,tail)

			d <- getURL(url)
			d <- strsplit(d, split="\n")
			d <- unlist(d)
			if (platform == "27k"){
				d <- d[grepl("HumanMethylation27.Level_1",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
			} else {
				d <- d[grepl("HumanMethylation450.Level_1",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
			}
			

			versions <- unlist(lapply(as.list(d),function(x) extract.version(x, cancer=cancer)))
			versions <- unlist(strsplit(versions, split="[.]"))
			versions <- matrix(as.numeric(versions), ncol=3, byrow=TRUE)

			o <- order(versions[,1], versions[,2], versions[,3])
			versions <- versions[o, ,drop=FALSE]

			ids <- unique(versions[,1])
			retained.versions <- matrix(NA, length(ids), 3)
			for (i in 1:length(ids)){
				indices <- which(versions[,1]==ids[i])
				last.index <- indices[length(indices)]
				retained.versions[i,] <- versions[last.index,]
			}

			versions <- paste(retained.versions[,1], 
							  retained.versions[,2],
							  retained.versions[,3], sep="."
						)

			if (platform=="27k") {
				files <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1.",versions,".tar.gz")
			} else {
				files <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1.",versions,".tar.gz")
			}
			
			links[[kk]] <- paste0(url,files)
			print(kk)

		}
		names(links) <- cancers
		links

	}














	# # Old function
	# extract.mappings <- function(platform=c("27k","450k")) {

	# 	if (platform == "27k"){
	# 		cancers <- c("brca", "coad", "gbm", "kirc", "kirp", "laml", 
	# 			"luad", "lusc", "ov", "read", "stad", "ucec")

	# 	} else {
	# 		cancers <- c("acc", "blca", "brca","coad","cesc", "dlbc","esca", "gbm","hnsc","kich",
	# 			"kirc","kirp","laml","lgg","lihc","luad","lusc", "meso","ov","paad",
	# 			"pcpg","prad","read","sarc","skcm","stad","thca","ucec","ucs","uvm")
	# 	}

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
	# 		print(kk)

	# 	}
	# 	names(mappings) <- cancers
	# 	mappings
	# }





	addClinical <- function(object){
		if (!is(object,"RGChannelSet") & 
			!is(object, "MethylSet") & 
			!is(object, "RatioSet") &
			!is(object, "GenomicRatioSet") &
			!is(object, "GenomicMethylSet")
		) {stop("Object must be an RGChannelSet or a (Genomic)RatioSet or (Genomic)MethylSet")}

		pd <- pData(object)
		if (!("diseaseabr" %in% colnames(pd))){
			stop("diseaseabr must be provided in the phenotype data")
		}
		cancer <- tolower(unique(pd$diseaseabr))
		clinicalData <- tcgaR:::getClinicalData(cancer)
		pd.tcga.id <- substr(pd$TCGA.ID, 1,12)
		pd.clinical <- clinicalData[match(pd.tcga.id, clinicalData$bcr_patient_barcode),]
		new.var <- colnames(pd.clinical)
		if (sum(new.var %in% colnames(pd))>0){
			stop("Automatic addition of clinical data has already been done or cannot be done.")
		}
		if (ncol(pd.clinical) !=0){
			pd <- cbind(pd, pd.clinical)
		}
		n <- ncol(pd.clinical)
		cat(paste0(n, " clinical variables have been added \n"))
		pData(object) <- pd
		object
	}









	# destdir : directory where the .tar.gz files will be saved
	download.tcga.idat <- function(cancer.types, destdir, platform=c("27k","450k")){
		cat("Updating the downloads from the TCGA portal \n")
		if (platform=="27k"){
			links <- extract.links(platform="27k")
		} else {
			links <- extract.links(platform="450k")
		}
		
		indices <- match(cancer.types, names(links))
		links <- unlist(links[indices])
		m <- length(links)
		for (i in 1:m){
			url <- links[i]
			cancer <- substr(url, regexpr("tumor", url) +6, regexpr("cgcc", url)-2 )
			if (platform == "27k"){
				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation27.Level_1")
			} else {
				start.patt <- paste0("jhu-usc.edu_",toupper(cancer),".HumanMethylation450.Level_1")
			}
			
			stop.patt  <- ".tar.gz"
			start <- regexpr(start.patt,url)[1] + nchar(start.patt) + 1
			stop  <- regexpr(stop.patt, url)[1] -1
			name <- paste0(start.patt, substr(url, start, stop), stop.patt)

			download(url, destfile=file.path(destdir, name))
			print(i)
		}
	}









