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

## Usage

The best way to get started is to load the package and have a look at the demo scripts in the [client directory](./client).

You can get help for all package functions by running `help(package="brainnet")` in R.
