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
  ACC | Adrenocortical carcinoma || :large_blue_circle: 
  BLCA | Bladder Urothelial Carcinoma || :large_blue_circle:
  BRCA | Breast invasive carcinoma | :large_blue_circle:| :large_blue_circle:
  CESC | Cervical squamous cell carcinoma || :large_blue_circle:
  CHOL | Cholangiocarcinoma || :large_blue_circle:
COAD | Colon adenocarcinoma | :large_blue_circle:| :large_blue_circle:
DLBC | Large B-cell Lymphoma || :large_blue_circle:
ESCA | Esophageal carcinoma || :large_blue_circle: 
GBM | Glioblastoma multiforme | :large_blue_circle:| :large_blue_circle:
HNSC | Head and Neck squamous cell carcinoma  || :large_blue_circle:
KICH | Kidney Chromophobe|| :large_blue_circle:
KIRP| Kidney renal papillary cell carcinoma| :large_blue_circle:| :large_blue_circle:
KIRC |Kidney renal clear cell carcinoma| :large_blue_circle:| :large_blue_circle:
LAML | Acute Myeloid Leukemia | :large_blue_circle:| :large_blue_circle:
LCML | Chronic Myelogenous Leukemia | | 
LGG | Lower grade glioma || :large_blue_circle: 
LIHC | Liver hepatocellular carcinoma || :large_blue_circle:
LUAD | Lung adenocarcinoma| :large_blue_circle:| :large_blue_circle:
LUSC | Lung squamous cell carcinoma| :large_blue_circle:| :large_blue_circle:
MESO | Mesothelioma || :large_blue_circle:
OV | Ovarian serous cystadenocarcinoma| :large_blue_circle:| :large_blue_circle:
PAAD | Pancreatic adenocarcinoma || :large_blue_circle:
PCPG | Pheochromocytoma and Paraganglioma || :large_blue_circle:
PRAD | Prostate adenocarcinoma|| :large_blue_circle: 
READ | Rectum adenocarcinoma| :large_blue_circle:| :large_blue_circle:
SARC | Sarcoma || :large_blue_circle:
SKCM | Skin Cutaneous Melanoma || :large_blue_circle:
STAD | Stomach adenocarcinoma| :large_blue_circle:| :large_blue_circle:
TGCT | Testicular Germ Cell Tumors || :large_blue_circle:
THCA | Thyroid carcinoma|| :large_blue_circle: 
THYM | Thymoma|| :large_blue_circle:  
UCEC | Uterine Corpus Endometrial Carcinoma | :large_blue_circle:| :large_blue_circle:
UCS | Uterine Carcinosarcoma || :large_blue_circle:
UVM | Uveal Melanoma || :large_blue_circle:



