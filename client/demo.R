#!/usr/bin/env Rscript
#
# Trains various machine learning classifiers on FreeSurfer structural neuroimaging data.
#
# Written by Tim Schäfer, 2021-09-10

library("brainnet");
library("fsbrain");
library("kernlab");
library("readxl");
library("effects");


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


data_full = fsbrain::group.agg.atlas.native(subjects_dir, subjects_list, measure, hemi="both", atlas="aparc");
data_full$subject = NULL;  # delete subject ID, not needed.
data_full$unknown = NULL;  # The 'unknown' and 'corpuscallosum' columns do not contain any valid data.
data_full$corpuscallosum = NULL;

considered_atlas_regions = colnames(data_full);

# add demographics data.
data_full$sex = demographics$`SEX_ID (1=m, 2=f)`
data_full$sex = as.factor(data_full$sex - 1.0);
data_full$age = demographics$AGE;
data_full$height = demographics$HEIGHT;

data_training = data_full[subjects_training_row_indices, ];


#data_training = fsbrain::group.morph.standard(subjects_dir, subjects_training, measure, fwhm = "10", df_t = TRUE);



################################################################################
##################### Run Gaussian Process Classification ######################
################################################################################

##### Train and evaluate model #####

fit_model = kernlab::gausspr(sex ~ ., data = data_training, type = "classification");


#### Validate on test data #####

#num_subjects_testing = 30L;
num_subjects_testing = length(subjects_list) - num_subjects_training;
first_idx_testing = num_subjects_training + 1L;
last_idx_testing = first_idx_testing + num_subjects_testing - 1L;
subjects_testing = subjects_list[first_idx_testing:last_idx_testing];
subjects_testing_row_indices = which(demographics$subject_data_dirname %in% subjects_testing);
fsbrain:::check.subjectslist(subjects_testing, subjects_dir = subjects_dir, report_name = "subjects_testing");

data_testing = data_full[subjects_testing_row_indices, ];
sex_testing = data_testing$sex;

cat(sprintf("Trained on %d subjects, tested on %d.\n", length(subjects_training), length(subjects_testing)));
data_testing$sex = NULL; # delete true labels.
res = predict(fit_model, data_testing);
brainnet::evaluate_model(actual = sex_testing, predicted =  res);


################################################################################
##################### Predict Sex using logistic regression ####################
################################################################################


##### Model comparison: with and without the neuroimaging data

##### No neuroimaging data: predict based on age + height only:
fit1 <- glm(formula = sex ~ age + height, data = data_full, family=binomial(link='logit'));
summary(fit1);

##### Add neuroimaging data for all atlas regions:
fit2 <- glm(formula = sex ~ ., data = data_full, family=binomial(link='logit'));
summary(fit1);

## Look at the AIC values in the output (keep in mind that lower AIC sores are better.)


################################################################################
###### Look at effect of sex for predicting NI data for all regions ############
################################################################################

region_idx = 1L;
region_fits = list();
pvalues_sex = list();
for(region_name in considered_atlas_regions) {
    formula_region = sprintf("%s ~ sex + age", region_name);
    region_fits[[region_name]] = glm(formula = formula_region, data = data_full, family=gaussian());
    cat(sprintf("### Handling Region '%s' (%d of %d). ###\n", region_name, region_idx, length(considered_atlas_regions)));
    pvalues_sex[[region_name]] = unname(coef(summary.glm(region_fits[[region_name]]))[2,4]);
    region_idx = region_idx + 1L;
}

# Now investigate region_fits and pvalues_sex.
effects::allEffects(region_fits[[1]]);


