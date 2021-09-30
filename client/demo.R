#!/usr/bin/env Rscript
#
# Trains various machine learning classifiers on FreeSurfer structural neuroimaging data.
#
# Written by Tim Schäfer, 2021-09-10

library("brainnet");
library("fsbrain");   # loading neuroimaging data and visualizing results
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

postproc_demographics <- function(demographics) {
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


subjects_dir = file.path(study_dir, "mri/freesurfer"); # the FreeSurfer SUBJECTS_DIR containing the neuroimaging data.
demographics_file = system.file("extdata", "IXI_demographics_filtered_fixed.xlsx", package = "brainnet", mustWork = TRUE);
demographics = readxl::read_excel(demographics_file);
demographics = postproc_demographics(demographics);

cat(sprintf("Loading FreeSurfer data from SUBJECTS_DIR '%s'.\n", subjects_dir ));

subjects_list = demographics$subject_data_dirname; # The directory names for the subjects, under the SUBJECTS_DIR, that are actually used for the analysis.



# use a subset only
num_subjects_training = 400L;
subjects_training = subjects_list[1:num_subjects_training];
subjects_training_row_indices = which(demographics$subject_data_dirname %in% subjects_training);
fsbrain:::check.subjectslist(subjects_training, subjects_dir = subjects_dir, report_name = "subjects_training");
sex_training = demographics$`SEX_ID (1=m, 2=f)`[subjects_training_row_indices];



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
data_full$sex = demographics$sex;
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

classfication_gp = data.frame('x'=res, 'y'=seq(length(res)), 'group'=sex_testing);
ggplot2::ggplot(classfication_gp, ggplot2::aes(x,y)) + ggplot2::geom_point(ggplot2::aes(colour = group)) +
    ggplot2::labs(title = "Classification results using Gaussian process classification", x = "Classification result", y = "Subject ID", color = "True group\n");

################################################################################
##################### Predict Sex using logistic regression ####################
################################################################################


##### Model comparison: with and without the neuroimaging data

##### No neuroimaging data: predict based on age + height only:
fit1 <- glm(formula = sex ~ age + height, data = data_full, family=binomial(link='logit'));
classfication_glm = data.frame('x'=fit1$fitted.values, 'y'=seq(length(fit1$fitted.values)), 'group'=data_full$sex, 'age'=data_full$age);
ggplot2::ggplot(classfication_glm, ggplot2::aes(x,y)) + ggplot2::geom_point(ggplot2::aes(colour = group)) +
    ggplot2::labs(title = "Classification results using logistic regression", x = "Classification result", y = "Subject ID", color = "True group\n") +
    ggplot2::geom_vline(xintercept=c(0.5), linetype="longdash");
summary(fit1);

## Compute the R squared value for the model.
cat(sprintf("The R squared value for the GLM predicting sex using a logistic link function without NI data is '%f'.\n", rsq::rsq(fit1)));

## We could also split data into train and test data and predict the test data using the GLM:
##preds <- predict(fit1, newdat=data_training, type="response");



##### Add neuroimaging data for all atlas regions:
fit2 <- glm(formula = sex ~ ., data = data_full, family=binomial(link='logit'));
summary(fit1);
cat(sprintf("The R squared value for the GLM predicting sex using a logistic link function with NI data is '%f'.\n", rsq::rsq(fit2)));

## Look at the AIC values in the output (keep in mind that lower AIC sores are better.)


################################################################################
###### Look at effect of sex for predicting NI data for all regions ############
################################################################################

# check for group differences in age
t.test(data_full$age[data_full$sex == "male"], data_full$age[data_full$sex == "female"]);

# match sample to remove difference
match = MatchIt::matchit(sex ~ age, data = data_full, method = "nearest", distance = "glm");
data_full_matched = MatchIt::match.data(match);

# make sure it worked out:
glm_data = data_full_matched;
t.test(glm_data$age[glm_data$sex == "male"], glm_data$age[glm_data$sex == "female"]);


# This function splits a named list with keys starting with 'lh_' and 'rh_' into two lists. The prefixes 'lh_' and 'rh_' get stripped, and the entries placed in the respective new list.
# @return: The two separate lists are returned in a single named list, with keys 'lh' and 'rh' for the inner lists.
split_named_list_by_hemi <- function(some_named_list) {
    lh_list = list();
    rh_list = list();
    for(entry in names(some_named_list)) {
        if(startsWith(entry, 'lh_')) {
            lh_list[[substring(entry, 4L)]] = some_named_list[[entry]];
        } else if(startsWith(entry, 'rh_')) {
            rh_list[[substring(entry, 4L)]] = some_named_list[[entry]];
        } else {
            cat(sprintf("Ignoring entry '%s': does not start with 'lh_' or 'rh_'.\n", entry));
        }
    }
    return(list('lh'=lh_list, 'rh'=rh_list));
}

region_idx = 1L;
region_fits = list();
pvalues_sex = list();
effect_sizes_sex = list();
for(region_name in considered_atlas_regions) {
    formula_region = sprintf("%s ~ sex + age", region_name);
    fit = glm(formula = formula_region, data = glm_data, family=gaussian());
    region_fits[[region_name]] = fit;
    cat(sprintf("### Handling Region '%s' (%d of %d). ###\n", region_name, region_idx, length(considered_atlas_regions)));
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
fit = region_fits[[1]];
summary(fit);
plot(effects::allEffects(fit)); # https://www.jstatsoft.org/article/view/v008i15/effect-displays-revised.pdf

#contrast::contrast(fit, list(sex = "male", age=20), list(sex = "female"), age=20);
#emmeans::emmeans(fit, specs = pairwise ~ sex); # https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/

effect_sizes_by_hemi = split_named_list_by_hemi(effect_sizes_sex); # split the single list with lh_ and rh_ prefixes into two lh and rh lists.
fsbrain::vis.region.values.on.subject(fsbrain::fsaverage.path(), 'fsaverage', lh_region_value_list = effect_sizes_by_hemi$lh, rh_region_value_list = effect_sizes_by_hemi$rh, atlas = "aparc", draw_colorbar = T);

