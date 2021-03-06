#!/usr/bin/env Rscript
#
# Trains various machine learning classifiers on FreeSurfer structural neuroimaging data.
#
# Written by Tim Schäfer, 2021-09-10

library("brainnet");
library("fsbrain");   # loading neuroimaging data and visualizing results. you need a version > 0.5.0, to get that run: devtools::install_github("dfsp-spirit/fsbrain")
library("kernlab");   # Gaussian process classificaiton
library("readxl");    # read Excel demographcis file
library("effects");   # GLM effects
library("emmeans");   # GLM effects
library("MatchIt");   # matching of patient/control groups
library("ggplot2");   # general purpose plots
library("rsq");       # to compute R squared of a GLM

#library("sva"); # batch correction using ComBat. To install: install.packages("BiocManager"); BiocManager::install("sva");

################################################################################
########################## Load data and demographics ##########################
################################################################################

measure = "thickness";

if(brainnet:::get_os() == "linux") {
    study_dir = "~/nas/projects/IXI_dataset";
    if(! dir.exists(study_dir)) {
        study_dir = sprintf("/media/%s/science/data/IXI_dataset", Sys.getenv("LOGNAME"));
    }
} else {
    study_dir = "/Volumes/shared/projects/IXI_dataset";
}

if(dir.exists("~/data/IXI_min")) {
    study_dir = "~/data/IXI_min"; # use local minimal data dir if available, loading from a local SSD is much faster than loading via LAN from the NAS.
}

subjects_dir = file.path(study_dir, "mri/freesurfer"); # the FreeSurfer SUBJECTS_DIR containing the neuroimaging data.
demographics_file = system.file("extdata", "IXI_demographics_filtered_fixed.xlsx", package = "brainnet", mustWork = TRUE);
demographics = readxl::read_excel(demographics_file);
demographics = brainnet::postproc_IXI_demographics(demographics);

cat(sprintf("Loading FreeSurfer data from SUBJECTS_DIR '%s'.\n", subjects_dir ));

subjects_list = demographics$subject_data_dirname; # The directory names for the subjects, under the SUBJECTS_DIR, that are actually used for the analysis.



# use a subset only
num_subjects_training = 400L;
subjects_training = subjects_list[1:num_subjects_training];
subjects_training_row_indices = which(demographics$subject_data_dirname %in% subjects_training);
fsbrain:::check.subjectslist(subjects_training, subjects_dir = subjects_dir, report_name = "subjects_training");
sex_training = demographics$sex[subjects_training_row_indices];



data_full = fsbrain::group.agg.atlas.native(subjects_dir, subjects_list, measure, hemi="split", atlas="aparc");
data_full$subject = NULL;  # delete subject ID, not needed.
#data_full$unknown = NULL;  # The 'unknown' and 'corpuscallosum' columns do not contain any valid data.
#data_full$corpuscallosum = NULL;
data_full$lh_unknown = NULL;
data_full$rh_unknown = NULL;
data_full$lh_corpuscallosum = NULL;
data_full$rh_corpuscallosum = NULL;

considered_atlas_regions = colnames(data_full);

# add demographics data.
data_full_dem = data_full;
data_full_dem$sex = demographics$sex;
data_full_dem$age = demographics$AGE;
data_full_dem$height = demographics$HEIGHT;

data_training = data_full_dem[subjects_training_row_indices, ];


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

data_testing = data_full_dem[subjects_testing_row_indices, ];
sex_testing = data_testing$sex;

cat(sprintf("Trained on %d subjects, tested on %d.\n", length(subjects_training), length(subjects_testing)));
data_testing$sex = NULL; # delete true labels.
res = predict(fit_model, data_testing);
brainnet::evaluate_model(actual = sex_testing, predicted =  res);

classfication_gp = data.frame('x'=res, 'y'=seq(length(res)), 'group'=sex_testing);
ggplot2::ggplot(classfication_gp, ggplot2::aes(classfication_glm$x, classfication_glm$y)) + ggplot2::geom_point(ggplot2::aes(colour = classfication_glm$group)) +
    ggplot2::labs(title = "Classification results using Gaussian process classification", x = "Classification result", y = "Subject ID", color = "True group\n");

################################################################################
##################### Predict Sex using logistic regression ####################
################################################################################


##### Model comparison: with and without the neuroimaging data

##### No neuroimaging data: predict based on age + height only:
fit1 <- glm(formula = sex ~ age + height, data = data_full_dem, family=binomial(link='logit'));
classfication_glm = data.frame('x'=fit1$fitted.values, 'y'=seq(length(fit1$fitted.values)), 'group'=data_full_dem$sex, 'age'=data_full_dem$age);
ggplot2::ggplot(classfication_glm, ggplot2::aes(classfication_glm$x, classfication_glm$y)) + ggplot2::geom_point(ggplot2::aes(colour = classfication_glm$group)) +
    ggplot2::labs(title = "Classification results using logistic regression", x = "Classification result", y = "Subject ID", color = "True group\n") +
    ggplot2::geom_vline(xintercept=c(0.5), linetype="longdash");
summary(fit1);

## Compute the R squared value for the model.
cat(sprintf("The R squared value for the GLM predicting sex using a logistic link function without NI data is '%f'.\n", rsq::rsq(fit1)));

## We could also split data into train and test data and predict the test data using the GLM:
##preds <- predict(fit1, newdat=data_training, type="response");



##### Add neuroimaging data for all atlas regions:
fit2 <- glm(formula = sex ~ ., data = data_full_dem, family=binomial(link='logit'));
summary(fit2);
cat(sprintf("The R squared value for the GLM predicting sex using a logistic link function with NI data is '%f'.\n", rsq::rsq(fit2)));

## Look at the AIC values in the output (keep in mind that lower AIC sores are better.)


################################################################################
############################## Match groups ####################################
################################################################################

# check for group differences in age
t.test(data_full_dem$age[data_full_dem$sex == "male"], data_full_dem$age[data_full_dem$sex == "female"]);

# match sample to remove difference
match = MatchIt::matchit(sex ~ age, data = data_full_dem, method = "nearest", distance = "glm");
data_full_matched = MatchIt::match.data(match);

# make sure it worked out:
glm_data = data_full_matched;
t.test(glm_data$age[glm_data$sex == "male"], glm_data$age[glm_data$sex == "female"]);


################################################################################
###### Look at effect of sex for predicting NI data for atlas regions ##########
################################################################################

region_idx = 1L;
region_fits = list();
pvalues_sex = list();
effect_sizes_sex = list();
for(region_name in considered_atlas_regions) {
    cat(sprintf("### Handling Region '%s' (%d of %d). ###\n", region_name, region_idx, length(considered_atlas_regions)));
    formula_region = sprintf("%s ~ sex + age", region_name);
    fit = glm(formula = formula_region, data = glm_data, family=gaussian());
    region_fits[[region_name]] = fit;
    pvalues_sex[[region_name]] = unname(coef(summary.glm(region_fits[[region_name]]))[2,4]);

    raw_sd_male = sd(glm_data[[region_name]][glm_data$sex == "male"]);
    raw_sd_female = sd(glm_data[[region_name]][glm_data$sex == "female"]);
    raw_sd_pooled = sqrt((raw_sd_male * raw_sd_male + raw_sd_female + raw_sd_female) / 2.0);
    effect_sex_male_mean = effects::effect("sex", fit)$fit[1];
    effect_sex_female_mean = effects::effect("sex", fit)$fit[2];
    cohen_d = (effect_sex_male_mean - effect_sex_female_mean) / raw_sd_pooled;
    effect_sizes_sex[[region_name]] = abs(cohen_d); # we are not interested in direction for effect size.

    region_idx = region_idx + 1L;
}

# Now investigate region_fits and pvalues_sex.
fit = region_fits$lh_bankssts;
summary(fit);
plot(effects::allEffects(fit)); # https://www.jstatsoft.org/article/view/v008i15/effect-displays-revised.pdf

#contrast::contrast(fit, list(sex = "male", age=20), list(sex = "female"), age=20);
#emmeans::emmeans(fit, specs = pairwise ~ sex); # https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/

effect_sizes_by_hemi = fsbrain::hemilist.from.prefixed.list(effect_sizes_sex); # split the single list with lh_ and rh_ prefixes into two lh and rh lists.
fsbrain::vis.region.values.on.subject(fsbrain::fsaverage.path(), 'fsaverage', lh_region_value_list = effect_sizes_by_hemi$lh, rh_region_value_list = effect_sizes_by_hemi$rh, atlas = "aparc", draw_colorbar = T);
