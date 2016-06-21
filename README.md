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
|ACC | Adrenocortical carcinoma | \checkmark | |



