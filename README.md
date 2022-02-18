# brainnet
Machine learning on structural, surface-based neuroimaging data for R.

## About

This is an R package for running statistical and machine learning models on surface-based neuroimaging data derived from three-dimensional magnetic resonance images (sMRI) of the human brain. The data required for the methods provided by this package are typically obtained by running [FreeSurfer](https://freesurfer.net/), CAT12 or similar software on raw T1-weighted sMRI scans and consist of per-vertex or region-based measurements on meshes representing the human cortex (e.g., cortical thickness at each vertex or aggregated within Desikan atlas regions).

The [client scripts](./client) that come with the package demonstrate how to run group comparisons and similar tasks using the general linear model (GLM). They can be used as templates that you can adapt for your own analyses.

This package is intended mainly for internal use, and for student projects in our working group. It will not be published on CRAN as it is very specific to what we do.

## State

This is still work in progress, come back another day.

## Installation

You will need a recent R version installed on your system. We recommend to install R for your operating system from the [official R website](https://www.r-project.org/). We also recommend to install the free personal edition of the [rstudio IDE](https://www.rstudio.com/products/rstudio/) after you have installed R itself, but that is optional.

To install this package from an R session (switch to the *console* tab in Rstudio):

```R
install.packages("remotes");
devtools::install_github("dfsp-spirit/brainnet");
```
### System dependencies

Note: *This section only applies if you compile from source (typically only under Linux). You can completely ignore it if you are not doing that. If you are under MacOS or Windows and are unsure, you are almost certainly NOT doing it.*

If you install from source (the default under Linux), make sure to install the system dependencies if the installation fails, then re-run the command above.

The only system dependency of this package is `libglpk-dev`, and it is only required if you want to use the automated sample matching functions with the `glpk` solver (which is part of the `Rglpk` R package).

Under Ubuntu and other Debian-based distros, this should get you started:

```shell
sudo apt install libglpk-dev
```

If you want to use the 3D visualization functions (highly recommended), you should also install the system dependencies for the `fsbrain` package, as explained on the [fsbrain project page](https://github.com/dfsp-spirit/fsbrain). The `fsbrain` package is a dependency of this package and R will try to install it automatically when you install `brainnet`.

If you want to compile manually under MacOS instead of using the available binary packages (advanced users), you probably know that you can get the required libraries using Homebrow or MacPorts. If you want to compile under Windows we cannot provide any assistance, we recommend to use the pre-compiled packages for that OS unless you know what you are doing.

## Usage

The best way to get started is to load the package and have a look at the demo scripts in the [client directory](./client). 

You can get help for all package functions by running `help(package="brainnet")` in R.

Note: Not all examples are complete, and you will have to adapt them substantially to your research question and sample. That is intended. You do NOT have to use these scripts, or the exact methods, functions and parameters in it: they are nothing but examples. Feel free to write your own analysis pipeline from scratch if you prefer that.
