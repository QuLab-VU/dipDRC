# dipDRC
The R code in the dipDRC repository has now been compiled into a complete R package called `diprate` that includes some basic documentation of the functions. The package contains functions that fall into three classes:
(1) data handling, (2) DIP rate extraction and quality control, and (3) plotting and dose-response curve analysis.

## Downloading/installing
Downloading and installing the `diprate` package requires the `devtools` package, which should be available from any CRAN mirror. You should be able to type `install.packages("devtools")` in the R console to install it.

The `diprate` package is available in this Git repository as a compiled tarball package (code that should run on any machine).

To install, use:
`devtools::install_github("QuLab-VU/dipDRC", subdir="diprate", dependencies=TRUE)`

Questions? Contact d.tyson at vanderbilt.edu
