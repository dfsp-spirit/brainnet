

#' @title Read aparcstats or segstats table file created by FreeSurfer tool 'aparcstats2table' or 'asegstats2table'.
#'
#' @param filepath character string, the path to the aparcstats or asegstats file.
#'
#' @return data.frame with seg/parc stats
#'
#' @export
read.segstats.table <- function(filepath) {
    return(read.table(filepath, sep='\t', header=TRUE, stringsAsFactors=FALSE));
}

