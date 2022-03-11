
#' @title Search for the ABIDE data directory on my machine.
#'
#' @author T Schaefer
#'
#' @note Ignore this internal function. You should simply hard-code the path to the data on your machine in your script. This function only exists to make the demo script short (and still allow it to run on all my machines during development).
#'
#' @return character string, the path to the ABIDE data, pre-processed with FreeSurfer (the recon-all output dir or SUBJECTS_DIR).
#'
#' @keywords internal
get_ABIDE_path_on_tims_machines <- function() {
    if(get_os() == "linux") {
        study_dir = "~/nas/projects/abide";
        if(! dir.exists(study_dir)) {
            study_dir = sprintf("/media/%s/science/data/abide", Sys.getenv("LOGNAME"));
        }
    } else {
        study_dir = "/Volumes/shared/projects/abide";
    }

    if(dir.exists("~/data/abide_min")) {
        study_dir = "~/data/abide_min"; # use local minimal data dir if available, loading from a local SSD is much faster than loading via LAN from the NAS.
    }

    subjects_dir = file.path(study_dir, "structural"); # the FreeSurfer SUBJECTS_DIR containing the neuroimaging data.
    return(subjects_dir);
}
