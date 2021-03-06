% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mass_glm.R
\name{slm_t}
\alias{slm_t}
\title{Compute t map for mass-univariate GLM analysis.}
\usage{
slm_t(X, Y, model.term, output = c("t", "p", "p.adjust"))
}
\arguments{
\item{X}{numerical matrix, the design or model matrix, typically created from the demographics data using \code{\link[stats]{model.matrix}}.}

\item{Y}{numerical matrix, the target value, typically neuroimaging data}

\item{model.term}{character string, the name of the factor in the model for which to compute the stats.}

\item{output}{vector of pre-defined character strings, defined what values to return. Leave alone if in doubt.}
}
\value{
named list with entries according to the \code{output} parameter. By default \code{t}= the t value map, \code{p} = the uncorrected p value map, \code{p.adjust} = the FDR-corrected p value map.
}
\description{
Compute t map for mass-univariate GLM analysis.
}
\author{
C Ecker, documentation by T Schaefer
}
