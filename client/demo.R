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

get_IXI_demographics <- function(study_dir) {
    subjects_dir = file.path(study_dir, "IXI-T1");

    ##### Load metadata / demographics #####

    ## This is a bit annoying for the IXI dataset, see comments below.
    demographics_file = file.path(study_dir, "IXI.xls")
    demographics = readxl::read_excel(demographics_file);

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


##### Load data #####

measure = "thickness";

if(brainnet:::get_os() == "linux") {
    study_dir = "~/nas/projects/IXI_dataset/";
} else {
    study_dir = "/Volumes/shared/projects/IXI_dataset/";
}
subjects_dir = file.path(study_dir, "IXI-T1"); # the FreeSurfer SUBJECTS_DIR containing the neuroimaging data.
demographics = get_IXI_demographics(study_dir);
subjects_list = demographics$subject_data_dirname; # The directory names for the subjects, under the SUBJECTS_DIR, that are actually used for the analysis.

# use a subset only for quick testing
subjects_training = subjects_list[1:500];
sex_training = demographics$`SEX_ID (1=m, 2=f)`[1:500];


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

subjects_testing = subjects_list[501:length(subjects_list)];
sex_testing = demographics$`SEX_ID (1=m, 2=f)`[501:length(subjects_list)];

if(do_use_region_data) {
    data_testing = fsbrain::group.agg.atlas.native(subjects_dir, subjects_testing, measure, hemi="both", atlas="aparc");
    data_testing$subject = NULL;
    data_testing$unknown = NULL;
    data_testing$corpuscallosum = NULL;
} else {
    data_testing = fsbrain::group.morph.standard(subjects_dir, subjects_testing, measure, fwhm = "10", df_t = TRUE);
}

res = predict(fit_model, data_testing);








