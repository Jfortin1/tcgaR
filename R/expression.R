library(RCurl)
library(downloader)





	getTCGA.expression <- function(cancer, 
		platform = c("genes","junctions","isoforms","genes.normalized","isoforms.normalized", "exons"), 
		verbose = FALSE, n=NULL){

		platform <- match.arg(platform)
		what <- match.arg(what)
		cancer <- tolower(cancer)

		# Let's see if the cancer exist:
		#doesItExist <- .cancer.exists(cancer = cancer, platform = platform)

		file.info <- .getRNANames(cancer=cancer, platform=platform)
		filenames <- file.info$rna.files
		filepaths <- file.info$rna.con
		if (is.null(n)){ n <- length(filenames)}
		
		filenames <- filenames[1:n]
		filepaths <- filepaths[1:n]
		
		mappings <- .getRNAMappings(cancer=cancer, platform=platform)
		mappings <- mappings[match(filenames, mappings$Derived.Data.File),]
		sampleNames <- as.character(mappings$Comment..TCGA.Barcode.)

		cat(paste0("[tcga.expression] ", n," samples have been found \n"))

		if (platform=="genes"){
			cat("[getTCGA.expression] Constructing the gene expression matrix \n")
			jj =1 
			temp <- read.table(text = getURL(filepaths[jj]), head=TRUE)
			raw <- temp[, "raw_count", drop=FALSE]
			scaled <- temp[,"scaled_estimate", drop=FALSE]
			for (jj in 2:length(filepaths)){
				temp <- read.table(text = getURL(filepaths[jj]), head=TRUE)
				raw <- cbind(raw, temp$raw_count)
				scaled <- cbind(scaled, temp$scaled_estimate)
				print(jj)
			}
			rownames(raw) <- rownames(scaled) <- temp$gene_id
			
		}
		colnames(raw) <- colnames(scaled) <- sampleNames
		list(raw=raw, scaled=scaled)
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
		rna.files <- list()
		for (kk in 1:length(subdirs)){
			url <- paste0(root,cancer,tail, subdirs[kk])
			d <- getURL(url)
			d <- strsplit(d, split="\n")
			d <- unlist(d)
			d <- d[grepl(platform,d)]
			d <- unique(extract.rna.filename(d))
			n <- length(d)
			rna.names[[kk]] <- paste0(base.url,subdirs[kk], d)
			rna.files[[kk]] <- d
			#print(kk)
			Sys.sleep(2) # Otherwise TCGA Portal complains. 
		}
		rna.names <- unlist(rna.names)
		rna.files <- unlist(rna.files)
		return(list(rna.con = rna.names, rna.files = rna.files))
	}







	.getRNAMappings <- function(cancer, platform = c("genes","junctions","isoforms","genes.normalized","isoforms.normalized", "exons")) {
		cancer <- tolower(cancer)
		platform <- match.arg(platform)
		if (platform=="genes" | platform== "isoforms" | platform=="genes.normalized" |
			platform== "isoforms.normalized" ){
				platform <- paste0(platform, ".results")
		} else {
			platform <- paste0(platform, "_quantification")
		}	
		
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
		mappings <- mappings[, c("Comment..TCGA.Barcode.", "Derived.Data.File")]
		mappings <- mappings[grepl(platform,mappings$Derived.Data.File),]
		barcodes <- .processBarcodes(mappings$Comment..TCGA.Barcode.)
		mappings <- cbind(mappings, barcodes)
		mappings
	}






	getTCGA.number <- function(cancer, platform = c("27k", "450k")){
		platform <- match.arg(platform)
		cancer <- tolower(cancer)
		filenames <- tcgaR:::.getIdatNames(cancer = cancer, platform = platform)
		n <- length(filenames[[1]])
		cat(paste0("[tcga.meth] ", n," samples have been found \n"))
	}







