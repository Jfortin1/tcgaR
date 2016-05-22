cancer.exists.meth <- function(cancer, platform=c("450k", "27k")){
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

cancer.exists <- function(cancer){
		cancers <- getCancers()
		cancer %in% cancers
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
