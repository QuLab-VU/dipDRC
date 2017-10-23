# dipDRC
The dipDRC code has now been compiled into a complete R package including rudimentary documentation of the functions. The package contains functions that fall into three classes:
(1) data handling, (2) DIP rate extraction and quality control, and (3) plotting and dose-response curve analysis.

## Downloading/installing
The `diprate` package is available as a compiled tarball package (code that should run on any machine).

Download the file, uncompress the `diprate` directory, and use `devtools::load_all(<path_to_diprate>)`

Alternatively, the package can be installed directly from GitHub using R and your user authentication. This requires the `devtools` package, which should be available from any CRAN mirror. You can modify the following R code to your system:
`devtools::install_github("QuLab-VU/dipDRC", subdir="diprate", user="<GitHub username>", auth_token="<GitHub auth token>")`

Follow the links to learn more about the [install_github](https://www.rdocumentation.org/packages/devtools/versions/1.13.3/topics/install_github) function in the `devtools` package and [GitHub authentication tokens](https://github.com/settings/tokens) 

