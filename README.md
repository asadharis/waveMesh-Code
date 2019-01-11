# Wavelet regression and additive models for irregularly spaced data
### Asad Haris, Ali Shojaie and Noah Simon

This repository contains accompanying source code for the manuscript ["Wavelet regression and additive models for irregularly spaced data"](http://papers.nips.cc/paper/8112-wavelet-regression-and-additive-models-for-irregularly-spaced-data), NeurIPS 2018, Montreal, QC. 

We briefly describe the relavent files for implementing the main methods/algorithms followed by a description of generating each figure/table of the final manuscript. 

## Main Methods

Function | File (Dependencies) | Description
------------ | ------------- | --------------
`waveMesh` | waveMesh.R (helpers.cpp) | Main function for univariate non-parametric regression using the waveMesh proposal
`cv.waveMesh` | waveMesh.R (helpers.cpp) | K-fold cross validation for univariate waveMesh
`fit.additive` | addWaveMesh.R (waveMesh.R, helpers.cpp) | Main function for fitting (sparse) additive models
`cv.addWaveMesh` |  addWaveMesh.R (waveMesh.R, helpers.cpp) | K-fold cross validation for (sparse) additive waveMesh

## Figures and Tables in Main Manuscript
The table below describes the procedure for generating each of the figures and tables. The numerical experiments were conducted on a linux computing cluster with a Sun Grid Engine (SGE) scheduler (now Oracle Grid Engine).

ID |  Description
------------ | ------------- | 
Table 1  | Run executable file **RUN_SECOND.csh** ---> Run file **resultsOtherMethods.R** ---> *Output*: latex table.
Figure 1 | Run R file **dataUni.R** --->  *Output*: Four plots of Figure 1.
Table 2 | Run executable file **RUN_ADDITIVE.csh** ---> Run file **resultsADD.R** ---> *Output*: 6 vectors for the 6 six columns in Table 2

