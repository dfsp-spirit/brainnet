#' @title Post-process ABIDE demographiocs data
#'
#' @param demographics the ABID demographics data.frame
#'
#' @return fixed df
#'
#' @export
load_ABIDE_metadata <- function() {
    abide_metadata = list();


    demographics_file = system.file("extdata", "ABIDE_Phenotypic_V1_0b_preprocessed1.csv", package = "brainnet", mustWork = TRUE);
    md_raw = read.table(demographics_file, header=TRUE, sep = ",", na.strings = c("", "-9999"), comment.char = "");

    # Load the ABIDE brainstats
    brainstats_file = path.expand("~/develop/research/scripts/abide_tools/brainstats_abide_with_pial.txt");
    brainstats = read.table(brainstats_file, header=TRUE, sep = " ")
    cat(sprintf("Received brains stats for %d subjects.\n", nrow(brainstats)));

    # filter subjects: use only those from diMartino study that have structural data
    md = subset(md_raw, md_raw$SUB_IN_SMP == 1);
    md = subset(md, md$FILE_ID != "no_filename");

    cat(sprintf("Received demographics for %d subjects, %d of them are in diMartino study and have structural data.\n", nrow(md_raw), nrow(md)));

    do_quick_test_only = TRUE;
    if(do_quick_test_only) {
        num_test = 100;
        warning(sprintf("!!!!!!!!!!!!  Using the first %d subjects only to speedup script  !!!!!!!!!!!!!!!", num_test));
        md = md[c(1:num_test),];
    }


    demographics = data.frame("subject_id"=md$FILE_ID, "group"=md$DX_GROUP, "site"=md$SITE_ID, "gender"=md$SEX, "age"=md$AGE_AT_SCAN, "iq"=md$FIQ);
    # Rewrite some stuff to strings instead of numbers: group
    demographics$group[md$DX_GROUP == 1] = "asd";
    demographics$group[md$DX_GROUP == 2] = "control";

    ## Rewrite some stuff to strings instead of numbers: gender
    demographics$gender[md$SEX == 1] = "male";
    demographics$gender[md$SEX == 2] = "female";

    impute_data = TRUE;
    if(impute_data) {
        demographics_and_brainstats = merge(demographics, brainstats, by.x="subject_id", by.y="subject");
        cat(sprintf("Merged demographics with %d rows and brainstats with %d, results has %d rows.\n", nrow(demographics), nrow(brainstats), nrow(demographics_and_brainstats)));

        print(head(demographics_and_brainstats));

        ## Check for NA values in important columns: This shows that 30 subjects have a NA value for "iq".
        #sapply(demographics_and_brainstats, function(x) sum(is.na(x)));



        ## Impute the missing iq data.
        ##  - We only impute data in columns we are interested in, otherwise it takes ages and the results for some columns make no sense anyways.
        ##  - We use the brainstats (computed from FreeSurfer parcellation stats using ExtractBrainMeasuresTS.bash script) to impute the values.


        init = mice::mice(demographics_and_brainstats, maxit=0);
        meth = init$method;
        predM = init$predictorMatrix;
        ## Methods for continues ("norm"), binary ("logreg") and ordinal ("polyreg") variables
        #meth[c("age")] = "norm";
        # Start imputation
        #imputation_res = mice(demographics_and_brainstats, method=meth, predictorMatrix=predM, m=5);
        imputation_res = mice::mice(demographics_and_brainstats, m=1, maxit=50, method='cart', seed=500);
        demographics_and_brainstats_imputed = mice::complete(imputation_res);
        ## Check again for NAs:
        #sapply(demographics_and_brainstats_imputed, function(x) sum(is.na(x)));
        demographics = demographics_and_brainstats_imputed;
    }

    subjects_list = as.character(as.vector(demographics$subject_id));

    abide_metadata$brainstats = brainstats;
    abide_metadata$subjects_list = subjects_list;
    abide_metadata$demographics = demographics;
    abide_metadata$demographics_file = demographics_file;
    abide_metadata$subjects_dir = subjects_dir;
    return(abide_metadata);
}

