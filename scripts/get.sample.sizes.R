








cancers <- getCancers()
cancers.27k <- unlist(lapply(cancers, function(x) {
	Sys.sleep(5)
	.cancer.exists(x)
}))

cancers.450k <- unlist(lapply(cancers, function(x) {
	Sys.sleep(5)
	.cancer.exists(x, "450k")
}))




cancers.27k <- cancers[cancers.27k]
normal <- tumor <- c()
for (i in 1:length(cancers.27k)){
	cancer <- cancers.27k[i]
	a <- .getIdatNames(cancer, "27k")
	b <- .getMethMappings(cancer, "27k")
	hist <- b$histology[match(a$idat.name,b$barcode)]
	n <- length(a$idat.name)
	n.control <- sum(grepl("Control", hist))
	n.normal  <- sum(grepl("Normal", hist))
	n.tumor <- n - n.normal - n.control
	normal[i] <- n.normal
	tumor[i]  <- n.tumor
	print(i)
}

data.27k <- cbind(cancers.27k, normal, tumor)

cancers.27k normal tumor
 [1,] "brca"      "27"   "318"
 [2,] "coad"      "37"   "166"
 [3,] "gbm"       "0"    "296"
 [4,] "kirc"      "199"  "219"
 [5,] "kirp"      "5"    "16" 
 [6,] "laml"      "0"    "194"
 [7,] "luad"      "24"   "127"
 [8,] "lusc"      "27"   "134"
 [9,] "ov"        "14"   "603"
[10,] "read"      "5"    "68" 
[11,] "stad"      "59"   "82" 
[12,] "ucec"      "1"    "117"





cancers.450k <- cancers[cancers.450k]
normal <- tumor <- c()
for (i in 1:length(cancers.450k)){
	cancer <- cancers.450k[i]
	a <- .getIdatNames(cancer, "450k")
	b <- .getMethMappings(cancer, "450k")
	hist <- b$histology[match(a$idat.name,b$barcode)]
	n <- length(a$idat.name)
	n.control <- sum(grepl("Control", hist))
	n.normal  <- sum(grepl("Normal", hist))
	n.tumor <- n - n.normal - n.control
	normal[i] <- n.normal
	tumor[i]  <- n.tumor
	print(i)
}