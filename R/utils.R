


#data <- getTCGA.expression("kirc", "genes")

.processBarcode <- function(barcodes){
	tss <- substr(barcodes, 6, 7)
	participant <- substr(barcodes, 9, 12)
	sample.type.num <- substr(barcodes, 14, 15)
	sample.type <- 
	center <- substr(barcodes, 27, 28)
	data.frame(tss=tss, participant = participant, sample.type = sample.type, center = center)

}

#a <- .processBarcode(barcode.meth)


