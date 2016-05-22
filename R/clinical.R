getClinicalData <- function(cancer){
	stopifnot(cancer.exists(cancer))
	root="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
	tail <- "/bcr/biotab/clin/nationwidechildrens.org_clinical_patient_"
	url <- paste0(root,cancer, tail, cancer, ".txt")
	clinical.data <- read.csv(text = getURL(url), sep="\t")
	clinical.data <- clinical.data[-c(1,2),]
	clinical.data 
}