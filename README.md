seqPatch
========
R package detecting DNA modifications by kinetic data of SMRT sequencing


System requirement
==================
Linux system. >=32Gb RAM is highly recommanded.


Pre-requisition
===============
seqPatch depends on Rcpp and pbh5. User need to install them first.


Install Rcpp
============
launch R console and type
>install.packages('Rcpp', dependencies=TRUE)


Install pbh5
============
see http://www.smrtcommunity.com/Share/Code?id=a1q70000000GrfbAAC&strRecordTypeName=Code


Install seqPatch
================
R CMD INSTALL downloaded_folder_name


A real example
==============
https://xfengz02.u.hpc.mssm.edu/example_release.tar.gz

Upcoming features
=================
1. Estimating proportion of modifiled DNA molecules at a given genome loci.
2. Empirical Bayesian FDR estimation.


 






