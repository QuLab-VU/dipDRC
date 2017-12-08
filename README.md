# dipDRC
The dipDRC code has now been compiled into a complete R package including rudimentary documentation of the functions. The package contains functions that fall into three classes:
(1) data handling, (2) DIP rate extraction and quality control, and (3) plotting and dose-response curve analysis.

## Downloading/installing
Downloading and installing the `diprate` package requires the `devtools` package, which should be available from any CRAN mirror. You should be able to type `install.packages("devtools")` in the R console to install it.

The `diprate` package is available in this Git repository as a compiled tarball package (code that should run on any machine).

There are at least three different ways for installing the package, once you have the `devtools` package installed:
1) Download the tarball file, uncompress the file to unpack the `diprate` directory, and type `devtools::load_all(<path_to_diprate>)` in the R console.

2) If you want to install the package into your R framework, perform step 1) but use the R command: 
`install.packages(<path_to_diprate>, repos = NULL, type="source", dependencies=TRUE)`

3) Alternatively, the package can be installed directly from GitHub using R and your user authentication. You can modify the following R code to your system:
`devtools::install_github("QuLab-VU/dipDRC", subdir="diprate", user="<GitHub username>", auth_token="<GitHub auth token>", dependencies=TRUE)`

If you do not install `diprate` into your R framework, you will need to install the `drc` package manually with: `install.packages("drc")`

Follow the links to learn more about the [install_github](https://www.rdocumentation.org/packages/devtools/versions/1.13.3/topics/install_github) function in the `devtools` package and [GitHub authentication tokens](https://github.com/settings/tokens) 

Questions? Contact d.tyson at vanderbilt.edu
