

#' @title Read aparcstats or segstats table file created by FreeSurfer tools.
#'
#' @param filepath character string, the path to the file.
#'
#' @return data.frame with seg/parc stats
#'
#' @export
read.segstats.table <- function(filepath) {
    return(read.table(filepath, sep='\t', header=TRUE, stringsAsFactors=FALSE));
}

