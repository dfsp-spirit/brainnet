#!/usr/bin/env Rscript
#
# Trains various machine learning classifiers on FreeSurfer structural neuroimaging data.
#
# Written by Tim Schäfer, 2021-09-10

library("brainnet");
library("fsbrain");
library("kernlab");
library("readxl");
library("fastglm");


##### Load data #####

measure = "thickness";

if(brainnet:::get_os() == "linux") {
    study_dir = "~/nas/projects/IXI_dataset";
    if(! dir.exists(study_dir)) {
        study_dir = sprintf("/media/%s/science/data/IXI_dataset", Sys.getenv("LOGNAME"));
    }
} else {
    study_dir = "/Volumes/shared/projects/IXI_dataset";
}
subjects_dir = file.path(study_dir, "mri/freesurfer"); # the FreeSurfer SUBJECTS_DIR containing the neuroimaging data.
demographics_file = system.file("extdata", "IXI_demographics_filtered_fixed.xlsx", package = "brainnet", mustWork = TRUE);
demographics = readxl::read_excel(demographics_file);


subjects_list = demographics$subject_data_dirname; # The directory names for the subjects, under the SUBJECTS_DIR, that are actually used for the analysis.



# use a subset only
num_subjects_training = 400L;
subjects_training = subjects_list[1:num_subjects_training];
subjects_training_row_indices = which(demographics$subject_data_dirname %in% subjects_training);
fsbrain:::check.subjectslist(subjects_training, subjects_dir = subjects_dir, report_name = "subjects_training");
sex_training = demographics$`SEX_ID (1=m, 2=f)`[subjects_training_row_indices];


##### Train and evaluate model #####

do_use_region_data = TRUE;

if(do_use_region_data) {
    data_full = fsbrain::group.agg.atlas.native(subjects_dir, subjects_list, measure, hemi="both", atlas="aparc");
    data_full$subject = NULL;
    data_full$unknown = NULL;
    data_full$corpuscallosum = NULL;

    # add demographics data.
    data_full$sex = demographics$`SEX_ID (1=m, 2=f)`
    data_full$sex = data_full$sex - 1.0;
    data_full$age = demographics$AGE;
    data_full$height = demographics$HEIGHT;

    data_training = data_full[subjects_training_row_indices, ];
} else {
    data_training = fsbrain::group.morph.standard(subjects_dir, subjects_training, measure, fwhm = "10", df_t = TRUE);
}

# We want to predict sex
data_training$y = data_training$sex;
data_training$sex = NULL;


##### Train and evaluate model #####

fit_model = kernlab::gausspr(y ~ ., data = data_training, type = "classification");


#### Validate on test data #####

#num_subjects_testing = 30L;
num_subjects_testing = length(subjects_list) - num_subjects_training;
first_idx_testing = num_subjects_training + 1L;
last_idx_testing = first_idx_testing + num_subjects_testing - 1L;
subjects_testing = subjects_list[first_idx_testing:last_idx_testing];
subjects_testing_row_indices = which(demographics$subject_data_dirname %in% subjects_testing);
fsbrain:::check.subjectslist(subjects_testing, subjects_dir = subjects_dir, report_name = "subjects_testing");
sex_testing = demographics$`SEX_ID (1=m, 2=f)`[subjects_testing_row_indices];

if(do_use_region_data) {
    data_testing = data_full[subjects_testing_row_indices, ];
} else {
    data_testing = fsbrain::group.morph.standard(subjects_dir, subjects_testing, measure, fwhm = "10", df_t = TRUE);
}

cat(sprintf("Trained on %d subjects, tested on %d.\n", length(subjects_training), length(subjects_testing)));
res = predict(fit_model, data_testing);
brainnet::evaluate_model(actual = sex_testing, predicted =  res);


##### Run GLM for group comparison: predict sex from data #####
data_full_no_y = data_full;
data_full_no_y$sex = NULL;
X = data_full_no_y;
y = data_full$sex;
X = as.matrix(X);
glm_family = "gaussian";
cat(sprintf("Fitting GLM with family '%s'.\n", glm_family));

#fit1 <- fastglm::fastglm(X, y, family = glm_family);
fit1 <- glm(formula = sex ~ age + height, data = data_full);
ggiraphExtra::ggPredict(fit1,colorAsFactor = TRUE,interactive=TRUE)
summary(fit1);

## Re-run without neuroimaging data, demographics only
X2 = as.matrix(data.frame('age'=data_full$age, 'height'=data_full$height));
fit2 <- fastglm::fastglm(X2, y, family = glm_family);
summary(fit2);


