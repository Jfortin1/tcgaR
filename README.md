# tcgaR
### TCGA R Portal for methylation array data.

---------

**Creator**: Jean-Philippe Fortin, jeanphi@mail.med.upenn.edu

##### Software status

| Resource:      | Travis CI     |
| -------------  |  ------------- |
| Platform:      | Linux       |
| R CMD check    | <a href="https://travis-ci.org/Jfortin1/tcgaR"><img src="https://travis-ci.org/Jfortin1/tcgaR.svg?branch=master" alt="Build status"></a> |

The package contains functions to import Illumina Methylation arrays data into `minfi` objects in the R software from the The Cancer Genome Atlas (TCGA) data portal. 

### Getting started

```{r}
library(tcgaR)
```

### Available cancers

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



