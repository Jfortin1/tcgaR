


#data <- getTCGA.expression("kirc", "genes")

.processBarcodes <- function(barcodes){
	tss <- substr(barcodes, 6, 7)
	participant <- substr(barcodes, 9, 12)
	sample.type <- substr(barcodes, 14, 15)
	sample.type <- as.numeric(as.character(sample.type))
	center <- substr(barcodes, 27, 28)
	sampleType <- read.csv(file.path("/Users/Jean-Philippe/tcgaR/data/", "sampleType.csv"))
	sample.code <- sample.type
	sample.type <- as.character(sampleType$Definition[match(as.character(sample.type), sampleType$Code)])
	data.frame(tss=tss, participant = participant, sample.type = sample.type, sample.code = sample.code, center = center)
}

#Tumor types range from 01 - 09, 
#normal types from 10 - 19 and 
#control samples from 20 - 29. 

#a <- .processBarcode(barcode.meth)

mappings <- .getMethMappings("brca", "27k")
barcodes <- mappings$TCGA.ID 
barcodes <- .processBarcodes(barcodes)







