#!/usr/bin/env Rscript
#
# Trains various machine learning classifiers on FreeSurfer structural neuroimaging data.
#
# Written by Tim Schäfer, 2021-09-10

library("brainnet");
library("fsbrain");
library("kernlab");
library("readxl");


#### Setup data ####

get_IXI_demographics <- function(study_dir, subjects_dir) {

    ##### Load metadata / demographics #####

    ## This is a bit annoying for the IXI dataset, see comments below.
    demographics_file = file.path(study_dir, "IXI_fixed.xls")
    demographics = readxl::read_excel(demographics_file);

    # Various columns contain 0 values for some subjects (e.g., a height of 0 cm), which really mean NA. Fix this.
    demographics$HEIGHT[demographics$HEIGHT == 0] = NA;
    demographics$WEIGHT[demographics$WEIGHT == 0] = NA;
    demographics$ETHNIC_ID[demographics$ETHNIC_ID == 0] = NA;
    demographics$MARITAL_ID[demographics$MARITAL_ID == 0] = NA;
    demographics$OCCUPATION_ID[demographics$OCCUPATION_ID == 0] = NA;
    demographics$QUALIFICATION_ID[demographics$QUALIFICATION_ID == 0] = NA;

    ## The demographics file for the IXI dataset does NOT contain the subject data directory names
    ## and there is no way to construct them from the data in there. We need to search the directories
    ## in the subjects_dir for ones starting with the expected pattern and match on them. ><
    demographics$subject_id_padded = sprintf("%03d", demographics$IXI_ID); # turn stuff like '12' into '012', like the dir names.

    ixi_dirs_and_files = list.files(path = subjects_dir, pattern="IXI[[:digit:]]{3}-");
    ixi_dirs = ixi_dirs_and_files[file.info(file.path(subjects_dir, ixi_dirs_and_files))$isdir];
    dir_info = data.frame("subject_data_dirname"=ixi_dirs, "subject_id_padded"=substring(ixi_dirs, 4,6)); # construct field to match on.


    metadata = base::merge(demographics, dir_info, by="subject_id_padded");     # Merge data.frames to retain only subjects with a FreeSurfer data directory.
    metadata = metadata[complete.cases(metadata),];                             # Remove all subjects with NA values in the metadata. These are subjects for which the age is unknown.
    return(metadata);
}

##### Helper functions #####


#' Compute metrics for classification model evaluation
#'
#' @param actual vector of actual labels
#'
#' @param predicted vector of predicted labels
#'
#' @return Nothing, called for printing side effect.
evaluate_model <- function(actual, predicted) {
    cm = as.matrix(table(Actual = actual, Predicted = predicted)); # confusion matrix

    n = sum(cm) # number of instances
    nc = nrow(cm) # number of classes
    diag = diag(cm) # number of correctly classified instances per class
    rowsums = apply(cm, 1, sum) # number of instances per class
    colsums = apply(cm, 2, sum) # number of predictions per class
    p = rowsums / n # distribution of instances over the actual classes
    q = colsums / n # distribution of instances over the predicted classes

    # per class
    accuracy = sum(diag) / n
    precision = diag / colsums
    recall = diag / rowsums
    f1 = 2 * precision * recall / (precision + recall)
    print(data.frame(precision, recall, f1))

    # macro-averaged
    macroPrecision = mean(precision)
    macroRecall = mean(recall)
    macroF1 = mean(f1)
    print(data.frame(macroPrecision, macroRecall, macroF1))
    return(invisible(NULL));
}


##### Load data #####

measure = "thickness";

if(brainnet:::get_os() == "linux") {
    study_dir = "~/nas/projects/IXI_dataset";
} else {
    study_dir = "/Volumes/shared/projects/IXI_dataset";
}
subjects_dir = file.path(study_dir, "mri/freesurfer"); # the FreeSurfer SUBJECTS_DIR containing the neuroimaging data.
demographics = get_IXI_demographics(study_dir, subjects_dir);
subjects_list = demographics$subject_data_dirname; # The directory names for the subjects, under the SUBJECTS_DIR, that are actually used for the analysis.



# use a subset only
num_subjects_training = 400;
subjects_training = subjects_list[1:num_subjects_training];
fsbrain:::check.subjectslist(subjects_training, subjects_dir = subjects_dir, report_name = "subjects_training");
sex_training = demographics$`SEX_ID (1=m, 2=f)`[1:num_subjects_training];


##### Train and evaluate model #####

do_use_region_data = TRUE;

if(do_use_region_data) {
    data_training = fsbrain::group.agg.atlas.native(subjects_dir, subjects_training, measure, hemi="both", atlas="aparc");
    data_training$subject = NULL;
    data_training$unknown = NULL;
    data_training$corpuscallosum = NULL;
} else {
    data_training = fsbrain::group.morph.standard(subjects_dir, subjects_training, measure, fwhm = "10", df_t = TRUE);
}

data_training$y = sex_training;

##### Train and evaluate model #####

fit_model = kernlab::gausspr(y ~ ., data = data_training, type = "classification");


#### Validate on test data #####

#num_subjects_testing = 30L;
num_subjects_testing = length(subjects_list) - num_subjects_training;
first_idx_testing = num_subjects_training + 1L;
last_idx_testing = first_idx_testing + num_subjects_testing - 1L;
subjects_testing = subjects_list[first_idx_testing:last_idx_testing];
fsbrain:::check.subjectslist(subjects_testing, subjects_dir = subjects_dir, report_name = "subjects_testing");
sex_testing = demographics$`SEX_ID (1=m, 2=f)`[first_idx_testing:last_idx_testing];

if(do_use_region_data) {
    data_testing = fsbrain::group.agg.atlas.native(subjects_dir, subjects_testing, measure, hemi="both", atlas="aparc");
    data_testing$subject = NULL;
    data_testing$unknown = NULL;
    data_testing$corpuscallosum = NULL;
} else {
    data_testing = fsbrain::group.morph.standard(subjects_dir, subjects_testing, measure, fwhm = "10", df_t = TRUE);
}

cat(sprintf("Trained on %d subject, tested on %d.\n", length(subjects_training), length(subjects_testing)));
res = predict(fit_model, data_testing);
evaluate_model(sex_testing, res);

