#' @title Load ABIDE demographics and brain stats.
#'
#' @description Load ABIDE demographics and brain stats and filter so that only subjects from the DiMartino study which have structural imaging data remain.
#'
#' @param impute_data logical, whether to impute missing data with mice. The filtered data (in_DiMartino_study AND has_structural_img) contains 30 subjects with missing IQ values.
#'
#' @param add_merged logical, whether to add return a single merged data.frame to the return value list with key 'merged'.
#'
#' @return named list with entries 'brainstats', 'demographics' and 'subjects_list'. The first two are data.frames, the last one is a vector of character strings.
#'
#' @importFrom mice mice complete
#' @importFrom utils read.table
#'
#' @export
load_ABIDE_metadata <- function(impute_data = TRUE, add_merged=TRUE) {
    abide_metadata = list();


    demographics_file = system.file("extdata", "ABIDE_Phenotypic_V1_0b_preprocessed1.csv", package = "brainnet", mustWork = TRUE);
    md_raw = utils::read.table(demographics_file, header=TRUE, sep = ",", na.strings = c("", "-9999"), comment.char = "");

    # Load the ABIDE brainstats
    brainstats_file = system.file("extdata", "ABIDE_brainstats_with_pial.txt", package = "brainnet", mustWork = TRUE);
    brainstats = utils::read.table(brainstats_file, header=TRUE, sep = " ");
    cat(sprintf("Received brains stats for %d subjects.\n", nrow(brainstats)));

    # filter subjects: use only those from diMartino study that have structural data
    md = subset(md_raw, md_raw$SUB_IN_SMP == 1);
    md = subset(md, md$FILE_ID != "no_filename");

    cat(sprintf("Received demographics for %d subjects, %d of them are in diMartino study and have structural data.\n", nrow(md_raw), nrow(md)));


    demographics = data.frame("subject_id"=md$FILE_ID, "group"=md$DX_GROUP, "site"=md$SITE_ID, "gender"=md$SEX, "age"=md$AGE_AT_SCAN, "iq"=md$FIQ);
    # Rewrite some stuff to strings instead of numbers: group
    demographics$group[md$DX_GROUP == 1] = "asd";
    demographics$group[md$DX_GROUP == 2] = "control";

    ## Rewrite some stuff to strings instead of numbers: gender
    demographics$gender[md$SEX == 1] = "male";
    demographics$gender[md$SEX == 2] = "female";

    if(impute_data) {
        demographics_and_brainstats = base::merge(demographics, brainstats, by.x="subject_id", by.y="subject");
        #cat(sprintf("Merged demographics with %d rows and brainstats with %d, results has %d rows.\n", nrow(demographics), nrow(brainstats), nrow(demographics_and_brainstats)));

        #print(head(demographics_and_brainstats));

        ## Check for NA values in important columns: This shows that 30 subjects have an NA value for "iq".
        #sapply(demographics_and_brainstats, function(x) sum(is.na(x)));
        num_subjects_with_missing_data = sum(as.integer(!complete.cases(demographics_and_brainstats)));
        cat(sprintf("Imputing IQ data for %d subjects with missing data (of %d subjects total).\n", num_subjects_with_missing_data, nrow(demographics_and_brainstats)));

        ## Impute the missing iq data.
        ##  - We only impute data in columns we are interested in, otherwise it takes ages and the results for some columns make no sense anyways.
        ##  - We use the brainstats (computed from FreeSurfer parcellation stats using ExtractBrainMeasuresTS.bash script) to impute the values.

        init = mice::mice(demographics_and_brainstats, maxit=0L);
        imputation_res = mice::mice(demographics_and_brainstats, m=5, maxit=50, method='cart', seed=500);
        demographics_and_brainstats_imputed = mice::complete(imputation_res);
        ## Check again for NAs:
        #sapply(demographics_and_brainstats_imputed, function(x) sum(is.na(x)));
        demographics$iq = demographics_and_brainstats_imputed$iq; # actually replace the NA IQ values with the imputed ones.
    } else {
        message("Not imputing data, some subjects will have missing IQ values.");
    }

    subjects_list = as.character(as.vector(demographics$subject_id));

    abide_metadata$brainstats = brainstats;
    abide_metadata$demographics = demographics;
    abide_metadata$subjects_list = subjects_list;
    if(add_merged) {
        abide_metadata$merged = base::merge(demographics, brainstats, by.x="subject_id", by.y="subject");
    }

    return(abide_metadata);
}

