library(RCurl)
library(downloader)




data <- getTCGA.expression("kirc", "genes")

link <- "https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=Tissue%20Source%20Site"


.processBarcode <- function(barcodes){
	tss <- substr(barcodes, 6, 7)
	participant <- substr(barcodes, 9, 12)
	sample.type.num <- substr(barcodes, 14, 15)
	sample.type <- 
	center <- substr(barcodes, 27, 28)
	data.frame(tss=tss, participant = participant, sample.type = sample.type, center = center)

}

a <- .processBarcode(barcode.meth)

getTCGA.expression <- function(cancer, platform = c("genes","junctions","isoforms","genes.normalized","isoforms.normalized", "exons"), what = c("both", "normal", "cancer"), verbose = FALSE){
		platform <- match.arg(platform)
		what <- match.arg(what)
		cancer <- tolower(cancer)

		# Let's see if the cancer exist:
		doesItExist <- .cancer.exist(cancer = cancer, platform = platform)

		filenames <- .getRNANames(cancer = cancer, platform = platform)
		sampleNames <- substr(filenames, )
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

		if (platform=="genes"){
			cat("[getTCGA.expression] Constructing of the gene expression matrix \n")
			jj =1 
			temp <- read.table(text = getURL(filenames[jj]), head=TRUE)
			raw <- temp[, "raw_count", drop=FALSE]
			scaled <- temp[,"scaled_estimate", drop=FALSE]
			for (jj in 1:length(filenames)){
				temp <- read.table(text = getURL(filenames[jj]), head=TRUE)
				raw <- cbind(raw, temp$raw_count)
				scaled <- cbind(scaled, temp$scaled_estimate)
				print(jj)
			}
			rownames(raw) <- rownames(scaled) <- temp$gene_id
			
		}
		object
	}


	.getRNANames <- function(cancer, platform = c("genes","junctions","isoforms","genes.normalized","isoforms.normalized", "exons")){
		cancer <- tolower(cancer)
		platform <- match.arg(platform)
		if (platform=="genes" | platform== "isoforms" | platform=="genes.normalized" |
			platform== "isoforms.normalized" ){
				platform <- paste0(platform, ".results")
		} else {
			platform <- paste0(platform, "_quantification")
		}
		root <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		tail <- "/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/"
		

		extract.version <- function(x, cancer){
				
				start.patt <- paste0("unc.edu_",toupper(cancer),".IlluminaHiSeq_RNASeqV2.Level_3")
				
				stop.patt  <- ".tar.gz"
				start <- regexpr(start.patt,x)[1] + nchar(start.patt) + 1
				stop  <- regexpr(stop.patt, x)[1] -1
				substr(x, start, stop)
		}	

		url <- paste0(root,cancer,tail)
		d <- getURL(url)
		d <- strsplit(d, split="\n")
		d <- unlist(d)
		d <- d[grepl("IlluminaHiSeq_RNASeqV2.Level_3",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]
		
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

		
		subdirs <- paste0("unc.edu_",toupper(cancer),".IlluminaHiSeq_RNASeqV2.Level_3.",versions,"/")

		extract.rna.filename <- function(x){
			start.patt <- 'a href=\\\"'
			stop.patt  <- platform
			start <- regexpr(start.patt,x)[1] + nchar(start.patt) -1
			stop  <- regexpr(stop.patt, x)[1] + nchar(stop.patt)-1
			substr(x, start, stop)
		}
		
		base.url <- paste0(root,cancer,tail)
		rna.names <- list()
		for (kk in 1:length(subdirs)){
			url <- paste0(root,cancer,tail, subdirs[kk])
			d <- getURL(url)
			d <- strsplit(d, split="\n")
			d <- unlist(d)
			d <- d[grepl(platform,d)]
			d <- unique(extract.rna.filename(d))
			n <- length(d)
			rna.names[[kk]] <- paste0(base.url,subdirs[kk], d)
			#print(kk)
			Sys.sleep(2) # Otherwise TCGA Portal complains. 
		}
		rna.names <- unlist(rna.names)
		return(rna.names)
	}





# barcode <- mappings$Comment..TCGA.Barcode.



	.getRNAMappings <- function(cancer, platform = c("genes","junctions","isoforms","genes.normalized","isoforms.normalized", "exons")) {
		cancer <- tolower(cancer)
		platform <- match.arg(platform)	
		
		root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		tail <- "/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/"
		url <- paste0(root,cancer,tail)

		d <- getURL(url)
		d <- strsplit(d, split="\n")
		d <- unlist(d)
		d <- d[grepl("mage-tab",d) & grepl("tar.gz",d) &!grepl("tar.gz.md5",d)]

		extract.version <- function(x, cancer){
			start.patt <- paste0("unc.edu_",toupper(cancer),".IlluminaHiSeq_RNASeqV2.mage-tab.")
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


		dir  <- paste0("unc.edu_",toupper(cancer),".IlluminaHiSeq_RNASeqV2.mage-tab.",versions,"/")		
		file <- paste0("unc.edu_",toupper(cancer),".IlluminaHiSeq_RNASeqV2.",versions, ".sdrf.txt")
		file <- paste0(url, dir, file)
		mappings <- read.csv(text = getURL(file), sep="\t")
	
		if ("X" %in% colnames(mappings)){
			mappings[,"X"] <- NULL
		}
		mappings
	}




	getTCGA.number <- function(cancer, platform = c("27k", "450k")){
		platform <- match.arg(platform)
		cancer <- tolower(cancer)
		filenames <- tcgaR:::.getIdatNames(cancer = cancer, platform = platform)
		n <- length(filenames[[1]])
		cat(paste0("[tcga.meth] ", n," samples have been found \n"))
	}


	# Need to be modified for RNA:
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






