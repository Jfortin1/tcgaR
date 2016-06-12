\name{getTCGA}
\alias{getTCGA}
\title{
  Import TCGA data into R and create the appropriate R object
}
\description{
  This function is the main user-level function in the tcgaR package. It downloads files from the TCGA portal for methylation and expression data and create the corresponding R objects via the minfi package. 
}

\usage{
getTCGA(cancer, datatype = c("methylation", "expression"), platform = c("450k", "27k"), idat=FALSE, idatDir=NULL, verbose=FALSE, n.samples=NULL, return=TRUE)
}

\arguments{
  \item{cancer}{An object of class \code{RGChannelSet}.}
  \item{datatype}{String indicating what data type should be imported. Should be either methylation or expression.}
  \item{platform}{String indicating which methylation platform should be used. Should be either 450k or 27k.}
  \item{idat}{Should the IDAT files be downloaded and saved on the disk?}
  \item{idatDir}{Directory in which the IDAT files will be saved if \code{idat=TRUE}}
  \item{verbose}{Should the function be verbose?}
  \item{n.samples}{Maximum number of samples downloaded for the cancer type. Mostly used for testing.}
  \item{return}{Should an R object be returned?}
}

\details{
  
This function implements functional normalization preprocessing for
Illumina methylation microarrays. Functional normalization extends the
idea of quantile normalization by adjusting for known covariates
measuring unwanted variation. For the 450k array, the first k principal
components of the internal control probes matrix play the role of the
covariates adjusting for technical variation. The number k of principal
components can be set by the argument \code{nPCs}. By default
\code{nPCs} is set to 2, and have been shown to perform consistently
well across different datasets.  This parameter should only be modified
by expert users.  The normalization procedure is applied to the Meth and
Unmeth intensities separately, and to type I and type II signals
separately. For the probes on the X and Y chromosomes we normalize males
and females separately using the gender information provided in the
\code{sex} argument. For the Y chromosome, standard quantile
normalization is used due to the small number of probes, which results
in instability for functional normalization. If sex is unspecified
(\code{NULL}), a guess is made using by the \code{getSex} function using
copy number information. Note that this algorithm does not rely on any
assumption and therefore can be be applicable for cases where global
changes are expected such as in cancer-normal comparisons or tissue
differences.

}
\value{If \code{return} is TRUE, an object of class \code{RGChannelSet} for 450k array data, and an object of class \code{MethylSet} for 27k array data. If \code{idat} is TRUE, the raw files are saved to the disk.}


\author{
  Jean-Philippe Fortin \email{jfortin@jhsph.edu},
}

\examples{
\dontrun{
  obj <- getTCGA("coad", datatype="methylation", platform="27k", n.samples=10)
}
}
