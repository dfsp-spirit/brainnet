
#' @title Get operating system name.
#'
#' @author T Schaefer
#'
#' @return One of 'osx' or 'linux'. Others are not supported.
#'
#' @keywords internal
get_os <- function(){
    sysinf = Sys.info();
    if (!is.null(sysinf)){
        os = sysinf['sysname'];
        if (os == 'Darwin') {
            os = "osx";
        }
    } else { ## mystery machine
        os = .Platform$OS.type;
        if (grepl("^darwin", R.version$os)) {
            os = "osx";
        }
        if (grepl("linux-gnu", R.version$os)) {
            os = "linux";
        }
    }
    return(tolower(unname(os)));
}
