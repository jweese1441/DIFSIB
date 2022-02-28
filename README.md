# DIFSIB

This package performs differential item functioning (DIF) and differential bundle functioning (DBF) using SIBTEST, Crossing-SIBTEST, and POLYSIBTEST procedures.
This package was created by translating the code directly from the open source Fortran code. Some modifications have been made to reduce the number of functions needed
to perform the SIBTEST, Crossing-SIBTEST, and POLYSIBTEST procedures. Further, the modified Crossing-SIBTEST statistics (Chalmers, 2018) has been implemented.

The original code can be found here: https://psychometrics.onlinehelp.measuredprogress.org/tools/dif/ 

The user manual can be found here: [DIFSIB_1.0.2.pdf](https://github.com/jweese1441/DIFSIB/files/8155285/DIFSIB_1.0.2.pdf)

In the below code, using force =TRUE will make the installation occur even if an update has not occured and you want to overwrite your previous install.

# INSTALLATION INSTRUCTIONS

To install this package run the following lines of code. 

install.packages("devtools")

library(devtools)

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)

devtools::install_github("jweese1441/DIFSIB", force =TRUE)
