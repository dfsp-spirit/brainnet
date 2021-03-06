% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset_abide.R
\name{load_ABIDE_metadata}
\alias{load_ABIDE_metadata}
\title{Load ABIDE demographics and brain stats for all subjects with structural imaging data.}
\usage{
load_ABIDE_metadata(
  impute_data = TRUE,
  exclude_bad_quality = c("UM_1_0050272")
)
}
\arguments{
\item{impute_data}{logical, whether to impute missing data with mice. The filtered data (has_structural_img) contains 69 subjects with missing IQ values.}

\item{exclude_bad_quality}{vector of character strings, the manually excluded subjects, e.g. due to bad quality}
}
\value{
named list with entries 'brainstats', 'demographics' and 'subjects_list'. The first two are data.frames, the last one is a vector of character strings.
}
\description{
Load ABIDE demographics and brain stats for all subjects with structural imaging data. This includes both male and female subjects.
}
\seealso{
Other ABIDE functions: 
\code{\link{aparcstats_files_ABIDE}()},
\code{\link{aparcstats_files_IXI}()},
\code{\link{load_ABIDE_metadata_males}()}
}
\concept{ABIDE functions}
