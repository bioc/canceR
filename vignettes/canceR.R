## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
BiocStyle::markdown()

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("RCurl", "XML")
#  install.packages(c("cgdsr","tkrplot","Formula", "RSvgDevice","RCurl" ))

## ----eval=FALSE---------------------------------------------------------------
#  library(biocManager)
#  biocManegr::install("GSEABase", "GSEAlm","geNetClassifier","Biobase", "phenoTest")
#  BiocManager::install("canceR")

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  devtools::install_git("kmezhoud/canceR")

## ----eval=FALSE---------------------------------------------------------------
#  library(canceR)
#  canceR()

