# tcgaR
### TCGA R Portal for methylation array data.

---------

**Creator**: Jean-Philippe Fortin, jeanphi@mail.med.upenn.edu

##### Software status

| Resource:      | Travis CI     |
| -------------  |  ------------- |
| Platform:      | Linux       |
| R CMD check    | <a href="https://travis-ci.org/Jfortin1/tcgaR"><img src="https://travis-ci.org/Jfortin1/tcgaR.svg?branch=master" alt="Build status"></a> |

The package contains functions to import Illumina Methylation arrays data from The Cancer Genome Atlas (TCGA) into the R softare. For Illumina 450k data, `minfi` objects are created (`RGChannelSet`), while `methylumi` objects are created for 27k data. The use of `methylumi` for 27k data is temporary; we are in the process of adapting `minfi` to support 27k array data. 

### Installation

In R, type the following commands to install `tcgaR`:
```{r}
library(devtools)
install_github("jfortin1/tcgaR")
```

Note that the following Bioconductor packages will need to be installed prior to installing `tcgaR`: `minfi`, `methylumi` and `illuminaio`.

### 450k array data

The different TCGA cancers for which 450k data are available are listed in the section "Available cancers for 27k and 450k platforms" below. To import all 450k methylation data for Colon adenocarcinoma (abbreviation BLCA) from the 450k platform, for both normal and tumor samples, we first load the package and use the `getTCGA` function: 
```{r}
library(tcgaR)
rgset <- getTCGA(cancer="blca", platform="450k")
```
We note that the function is not case-sensitive. This creates an `RGChannelSet` (see `minfi` package), the starting object in `minfi` that contains the raw array data. The option `idat=TRUE` will download the raw .idat files associated with each sample (one for the green channel and one for the red channel), and the option `idatDir` specifies the directory in which the .idat files will be saved. To obtain the metadata associated with the methylation data, we have created the function `getMethMappings`. For instance, to download the metadata for the BLCA dataset, we use
```{r}
metadata <- getMethMappings(cancer="blca", platform="450k")
```
This returns a `data.frame` with the sample names, the histology of the samples, and the TCGA barcode. One can add the metadata to the `RGChannelSet` using 
```{r}
pData(rgset) <- metadata
```
For further processing of the data, we recommend to use the function `preprocessFunnorm` which uses functional normalization to process data. Functional normalization is an extension of quantile normalization for data that show global epigenetic change between two conditions, as in the case of a normal/cancer comparison. To apply functional normalization, we use the following command: 

```{r}
grset <- preprocessFunnorm(rgset)
```



### 27k array data


For post-normalization statistical analyses, please read the `minfi` package vignette.
 
### Available cancers for 27k and 450k platforms:

| Abbreviation      | Full name   | 27k | 450k
| -------------  |  ------------- | ---- | ------ |
  ACC | Adrenocortical carcinoma || :ballot_box_with_check: 
  BLCA | Bladder Urothelial Carcinoma || :ballot_box_with_check:
  BRCA | Breast invasive carcinoma | :ballot_box_with_check:| :ballot_box_with_check:
  CESC | Cervical squamous cell carcinoma || :ballot_box_with_check:
  CHOL | Cholangiocarcinoma || :ballot_box_with_check:
COAD | Colon adenocarcinoma | :ballot_box_with_check:| :ballot_box_with_check:
DLBC | Large B-cell Lymphoma || :ballot_box_with_check:
ESCA | Esophageal carcinoma || :ballot_box_with_check: 
GBM | Glioblastoma multiforme | :ballot_box_with_check:| :ballot_box_with_check:
HNSC | Head and Neck squamous cell carcinoma  || :ballot_box_with_check:
KICH | Kidney Chromophobe|| :ballot_box_with_check:
KIRP| Kidney renal papillary cell carcinoma| :ballot_box_with_check:| :ballot_box_with_check:
KIRC |Kidney renal clear cell carcinoma| :ballot_box_with_check:| :ballot_box_with_check:
LAML | Acute Myeloid Leukemia | :ballot_box_with_check:| :ballot_box_with_check:
LCML | Chronic Myelogenous Leukemia | | 
LGG | Lower grade glioma || :ballot_box_with_check: 
LIHC | Liver hepatocellular carcinoma || :ballot_box_with_check:
LUAD | Lung adenocarcinoma| :ballot_box_with_check:| :ballot_box_with_check:
LUSC | Lung squamous cell carcinoma| :ballot_box_with_check:| :ballot_box_with_check:
MESO | Mesothelioma || :ballot_box_with_check:
OV | Ovarian serous cystadenocarcinoma| :ballot_box_with_check:| :ballot_box_with_check:
PAAD | Pancreatic adenocarcinoma || :ballot_box_with_check:
PCPG | Pheochromocytoma and Paraganglioma || :ballot_box_with_check:
PRAD | Prostate adenocarcinoma|| :ballot_box_with_check: 
READ | Rectum adenocarcinoma| :ballot_box_with_check:| :ballot_box_with_check:
SARC | Sarcoma || :ballot_box_with_check:
SKCM | Skin Cutaneous Melanoma || :ballot_box_with_check:
STAD | Stomach adenocarcinoma| :ballot_box_with_check:| :ballot_box_with_check:
TGCT | Testicular Germ Cell Tumors || :ballot_box_with_check:
THCA | Thyroid carcinoma|| :ballot_box_with_check: 
THYM | Thymoma|| :ballot_box_with_check:  
UCEC | Uterine Corpus Endometrial Carcinoma | :ballot_box_with_check:| :ballot_box_with_check:
UCS | Uterine Carcinosarcoma || :ballot_box_with_check:
UVM | Uveal Melanoma || :ballot_box_with_check:



