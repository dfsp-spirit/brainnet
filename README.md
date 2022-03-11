# brainnet
Some helper functions for Machine learning on structural, surface-based neuroimaging data for R. 

Ignore this unless you are a student in our working group.

## About

This is an R package for running statistical and machine learning models on surface-based neuroimaging data derived from three-dimensional magnetic resonance images (sMRI) of the human brain. The data required for the methods provided by this package are typically obtained by running [FreeSurfer](https://freesurfer.net/), CAT12 or similar software on raw T1-weighted sMRI scans and consist of per-vertex or region-based measurements on meshes representing the human cortex (e.g., cortical thickness at each vertex or aggregated within Desikan atlas regions).

The [client scripts](./client) that come with the package demonstrate how to run group comparisons and similar tasks using the general linear model (GLM). They can be used as templates that you can adapt for your own analyses.

The package also provides functions and aggregated data that make it very easy to work with the publicly available [ABIDE](https://fcon_1000.projects.nitrc.org/indi/abide/) and [IXI](https://brain-development.org/ixi-dataset/) data sets. These data sets are also used in the client scripts to demonstrate analysis pipelines.

This package is intended mainly for internal use, and for student projects in our working group. This package will **not** be published on [CRAN](https://cran.r-project.org/) as it is very specific to what we do.


## Installation

**This is still work in progress, come back another day. The following instructions are for developers and *those who dare*:**

You will need R installed on your system. We recommend to install a recent R version for your operating system from the [official R website](https://www.r-project.org/). We also suggest to install the free personal edition of the [Rstudio IDE](https://www.rstudio.com/products/rstudio/) after you have installed R itself, but that is optional.

To install this package from an R session (switch to the *console* tab in Rstudio):

```R
install.packages("remotes");
remotes::install_github("dfsp-spirit/brainnet", dependencies=TRUE);
```

#### System dependencies

Note: *This section only applies if you compile from source (typically only under Linux). You can completely ignore it if you are not doing that. If you are under MacOS or Windows and are unsure, you are almost certainly NOT doing it.*

If you install from source (the default under Linux), make sure to install the system dependencies if the installation fails, then re-run the command above.

Under Ubuntu and other Debian-based distros, this should get you started:

```shell
sudo apt-get install libglpk-dev libnlopt-dev
```

If you want to use the 3D visualization functions (highly recommended), you should also install the system dependencies for the `fsbrain` package, as explained on the [fsbrain project page](https://github.com/dfsp-spirit/fsbrain). The `fsbrain` package is a dependency of this package and R will try to install it automatically when you install `brainnet`. At the time of writing, this will do it for Debian-based systems:

```shell
sudo apt-get install libmagick++-dev libx11-dev libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev libfreetype6-dev libxml2-dev libssh-dev libcurl4-openssl-dev libgfortran4
```

If you want to compile manually under MacOS instead of using the available binary packages (advanced users), you probably know that you can get the required libraries using Homebrow or MacPorts. Don't forget about X11, see the [fsbrain project page](https://github.com/dfsp-spirit/fsbrain) for details. If you want to compile under Windows we cannot provide any assistance, we recommend to use the pre-compiled packages for that OS unless you know what you are doing.

## Usage

Once you have installed the `brainnet` package (see above), I would suggest to download a copy of one of the demo scripts to your computer. Create a new directory for all your analyses, place the file in there, and open it in Rstudio to edit it. 

E.g., download the file [client/demo_ABIDE_regionbased.Rmd](https://raw.githubusercontent.com/dfsp-spirit/brainnet/main/client/demo_ABIDE_regionbased.Rmd) or [browse all scripts in the client directory](./client).

If you downloaded the [demo_ABIDE_regionbased.Rmd](https://raw.githubusercontent.com/dfsp-spirit/brainnet/main/client/demo_ABIDE_regionbased.Rmd) notebook, you can open it in Rstudio and click the small *Knit* button above the *source* panel to run it.


You can get help for all brainnet package functions by running `help(package="brainnet")` in R.

#### Notes:

* Not all example scripts in the client directory are complete, and you will have to adapt them substantially to your research question and sample.
* Of course you do **not** have to use these scripts, or the exact methods, functions and parameters in it when working in our group: they are nothing but examples. Feel free to write your own analysis pipeline from scratch if you prefer that.
* To make our research as reproducible as possible, we prefer a fully or mostly automated analysis pipeline. If some steps cannot be automated (e.g., QC involving manual inspection of the images), make sure to document these steps properly. Also make it clear in the scripts what input data are used, and how they were generated. Referring to external documents is fine, e.g., make it clear that a list of subjects that are removed from the full data set is excluded because they were found to be of bad quality during manual QC. You will also have to explain this in your report.

## References

The functions in these package use data from the [ABIDE I](https://fcon_1000.projects.nitrc.org/indi/abide/) and [IXI](https://brain-development.org/ixi-dataset/) datasets. Make sure to cite the relevant publications for the data sets you use in your report.

## License

* This package contains aggregated data derived from the ABIDE I data set. That data is published under a Creative Commons, Attribution-NonCommercial-Share Alike License as explained in the [ABIDE I Usage Agreement](https://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html).
* This package contains aggregated data derived from the IXI data set. That data is published under the CC BY-SA 3.0 license, as explained on the [IXI website](https://brain-development.org/ixi-dataset/).
* Regarding the source code of the `brainnet` package itself: this is WIP and we have not decided on a license yet: this software has [no license](https://choosealicense.com/no-permission/) at this time. Simply put: do not use it yet.
