

	#Tumor types range from 01 - 09, 
	#normal types from 10 - 19 and 
	#control samples from 20 - 29. 
	.processBarcodes <- function(barcodes){
		tss <- substr(barcodes, 6, 7)
		participant <- substr(barcodes, 9, 12)
		sample.type <- substr(barcodes, 14, 15)
		sample.type <- as.numeric(as.character(sample.type))
		center <- substr(barcodes, 27, 28)
		vial <- substr(barcodes, 16, 16)
		portion <- substr(barcodes, 18,19)
		portion <- as.numeric(as.character(portion))
		sampleType <- read.csv(file.path("/Users/Jean-Philippe/tcgaR/data/", "sampleType.csv"))
		sample.code <- sample.type
		sample.type <- as.character(sampleType$Definition[match(as.character(sample.type), sampleType$Code)])
		sample.id   <- paste(tss, participant, sample.code, vial,portion, sep="-")
		status <- sample.code
		status[status >= 1 & status <=9]  <- "Tumor"
		status[status >=10 & status <=19] <- "Normal"
		status[status >=20 & status <=29] <- "Control"

		data.frame(tss=tss, participant = participant, 
			sample.type = sample.type, 
			sample.code = sample.code, 
			status = status,
			center = center, 
			vial = vial,
			portion = portion,
			sample.id = sample.id
		)
	}

	





	.cancer.exists <- function(cancer, platform=c("27k", "450k")){
			cancer <- tolower(cancer)
			platform <- match.arg(platform)
			root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
			if (platform == "27k"){
				tail="/cgcc/jhu-usc.edu/humanmethylation27/methylation/"
			} else {
				tail="/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
			}
			url <- paste0(root,cancer,tail)
			exist <- url.exists(url)
			exist
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


	getClinicalData <- function(cancer = c("acc", "blca", "brca","coad","cesc" ,"dlbc","esca", "gbm","hnsc","kich",
	 		"kirc","kirp","laml","lgg","lihc","luad","lusc", "meso","ov","paad",
	 		"pcpg","prad","read","sarc","skcm","stad","thca","ucec","ucs","uvm")){

		cancer <- match.arg(cancer)
		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		tail <- "/bcr/biotab/clin/nationwidechildrens.org_clinical_patient_"
		url <- paste0(root,cancer, tail, cancer, ".txt")
		clinical.data <- read.csv(text = getURL(url), sep="\t")
		clinical.data <- clinical.data[-c(1,2),]
		clinical.data 
	}



	a <- getClinicalData("gbm")


	# cancers <- getCancers()
	# cancers.27k <- unlist(lapply(cancers, function(x) {
	# 	Sys.sleep(5)
	# 	.cancer.exists(x)
	# }))

	# cancers.450k <- unlist(lapply(cancers, function(x) {
	# 	Sys.sleep(5)
	# 	.cancer.exists(x, "450k")
	# }))

	# data <- cbind(cancers, cancers.27k, cancers.450k)
	# bad <- c("cntl", "fppp", "lnnh", "misc")
	# data <- data[-match(bad, cancers),]



	# sleep() is necessary otherwise TCGA detects a robot.
	# cancers <- getCancers()
	# exist27k <- exist450k <- c()
	# for (i in 1:length(cancers)){
	# 	exist27k[i] <- .cancer.exists(cancers[i], "27k")
	# 	exist450k[i] <- .cancer.exists(cancers[i], "450k")
	# 	Sys.sleep(5)
	# 	print(i)
	# }

	












