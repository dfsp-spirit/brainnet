
#' @title Post-process IXI demographiocs data
#'
#' @param demographics the IXI demographics data.frame
#'
#' @return fixed df
#'
#' @export
postproc_IXI_demographics <- function(demographics) {
    demographics$DOB = as.Date(demographics$DOB, format="%Y-%m-%d");
    demographics$STUDY_DATE = as.Date(demographics$STUDY_DATE, format="%Y-%m-%d");
    demographics$DATE_AVAILABLE = NULL; # Delete useless column, the date column would be NA for such dates.

    # Compute scanner site column
    list_of_token_vectors = strsplit(demographics$subject_data_dirname, split="-"); # Each entry is split into 3 parts, which look like this:  "IXI639" "Guys"   "1088". The 2nd one ("Guys" here) is the site we are interested in.
    nr = nrow(demographics);
    sites = rep("?", nr);
    for(row_idx in seq.int(nr)) {
        sites[row_idx] = list_of_token_vectors[[row_idx]][2];
    }
    demographics$site = as.factor(sites);

    # Remap sex columns
    demographics$sex = as.factor(demographics$`SEX_ID (1=m, 2=f)`);
    levels(demographics$sex) = list("male"="1", "female"="2");

    # Remap and clean the ID columns: ETHNIC_ID
    demographics$ethnic = as.factor(demographics$ETHNIC_ID);
    levels(demographics$ethnic) = list("white"="1", "asian_br"="3", "black_br"="4", "chinese"="5", "other"="6");
    demographics$ethnic[which(is.na(demographics$ethnic))] = "other"; # Fix those with ID 2, which is missing from the legend. It is unknown what it is supposed to mean so we set it to "other".

    # Remap and clean the ID columns: MARITAL_ID
    demographics$marital = as.factor(demographics$MARITAL_ID);
    levels(demographics$marital) = list("single"="1", "married"="2", "cohabiting"="3", "divorced_sep"="4", "widowed"="5");

    # Remap and clean the ID columns: OCCUPATION_ID
    demographics$occupation = as.factor(demographics$OCCUPATION_ID);
    levels(demographics$occupation) = list("out_fulltime"="1", "out_parttime"="2", "study"="3", "fulltime_housework"="4", "retired"="5", "unemployed"="6", "work_at_home"=7, "other"=8);

    # Remap and clean the ID columns: QUALIFICATION_ID. This one seems to have an ordering, so we split into an integer column for the value and add a convenience columns that contains the description string.
    demographics$qualification = as.integer(demographics$QUALIFICATION_ID);
    demographics$qualification_desc = as.factor(demographics$QUALIFICATION_ID);
    levels(demographics$qualification_desc) = list("none"="1", "O-levels"="2", "A-levels"="3", "further_edu"="4", "university"="5");

    # Delete the unused fields.
    demographics$`SEX_ID (1=m, 2=f)` = NULL;
    demographics$ETHNIC_ID = NULL;
    demographics$MARITAL_ID = NULL;
    demographics$OCCUPATION_ID = NULL;
    demographics$QUALIFICATION_ID = NULL;
    demographics$subject_id_padded = NULL;

    return(demographics);
}


#' @title Load IXI brain stats and filter by subjects_list.
#'
#' @param subjects_list vector of character strings, the subjects to filter. Pass NULL to get full brain stats.
#'
#' @return data.frame with brain stats
#'
#' @keywords internal
get_IXI_brainstats_for_subjects <- function(subjects_list) {
    brainstats_file = system.file("extdata", "IXI_brainstats_with_pial.txt", package = "brainnet", mustWork = TRUE);
    brainstats = utils::read.table(brainstats_file, header=TRUE, sep = " ");

    if(is.null(subjects_list)) {
        return(brainstats);
    } else {
        filtered_brainstats = brainstats[brainstats$subject %in% subjects_list,];
        cat(sprintf("Received file with brains stats for %d subjects, kept %d which occur in subjects_list of length %d.\n", nrow(brainstats), nrow(filtered_brainstats), length(subjects_list)));
        return(filtered_brainstats);
    }
}


#' @title Load IXI demographics and brain stats
#'
#' @return named list with entries 'brainstats', 'demographics' and 'subjects_list'. The first two are data.frames, the last one is a vector of character strings.
#'
#' @importFrom readxl read_excel
#' @export
load_IXI_metadata <- function() {
    demographics_file = system.file("extdata", "IXI_demographics_filtered_fixed.xlsx", package = "brainnet", mustWork = TRUE);
    demographics = readxl::read_excel(demographics_file);
    demographics = postproc_IXI_demographics(demographics);
    subjects_list = demographics$subject_data_dirname;
    brainstats = get_IXI_brainstats_for_subjects(subjects_list);
    return(list("brainstats"=brainstats, "demographics"=demographics, "subjects_list"=subjects_list));
}


