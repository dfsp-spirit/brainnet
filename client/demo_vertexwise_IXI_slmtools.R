#!/usr/bin/env Rscript
#
# Trains various machine learning classifiers on FreeSurfer structural neuroimaging data.
#
# Written by Tim Schäfer, 2021-09-10

library("brainnet");
library("fsbrain");   # loading neuroimaging data and visualizing results. you need a version > 0.5.0, to get that run: devtools::install_github("dfsp-spirit/fsbrain")
library("readxl");    # read Excel demographcis file
library("MatchIt");   # matching of patient/control groups
library("ggplot2");   # general purpose plots
library("slmtools");  # internal R package of Ecker neuroimaging group for mass-univariate GLM analysis.

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


data_full = fsbrain::group.morph.standard(subjects_dir, subjects_list, measure, hemi="both", fwhm="10", df_t=TRUE); # the per-vertex data
considered_vertexcolnames = colnames(data_full);

glm_data = data.matrix(data_full); # Create the unmatched data.matrix from the purely numeric data.frame.

################################################################################
############################## Match groups ####################################
################################################################################

do_matching = FALSE; # Whether to perform cardinality matching for a matched sample.

if(do_matching) {

    # Add demographics data to data.frame, this is required by the MatchIt package.
    data_full_dem = data_full;
    data_full_dem$sex = as.factor(demographics$sex);
    data_full_dem$age = demographics$AGE;
    data_full_dem$site = as.factor(demographics$site);
    data_full_dem$qualification = as.factor(demographics$qualification);
    data_full_dem$subject_id = demographics$subject_data_dirname; # we may need this to identify which subjects were kept after matching.

    # check for group differences in age
    t.test(data_full_dem$age[data_full_dem$sex == "male"], data_full_dem$age[data_full_dem$sex == "female"]);

    # Match sample to remove difference. TODO: The group variable needs to be ordered, and we need to check that factors are used where appropriate.
    solver = "glpk"; # USe "gurobi" if you have it (requires registering for academic license and manual install), or "glpk" if not.
    match = MatchIt::matchit(sex ~ age + site + qualification, data = data_full_dem, method = "cardinality", solver = solver);
    summary(match);
    data_full_matched = MatchIt::match.data(match);

    t.test(glm_data$age[glm_data$sex == "male"], glm_data$age[glm_data$sex == "female"]);

    # Remove covariates from matched data (only the numerical neuroimaging data must remain).
    data_full_matched$sex = NULL;
    data_full_matched$age = NULL;
    data_full_matched$site = NULL;
    data_full_matched$qualification = NULL;

    # Filter the demographics for the retained subjects
    demographics = subset(demographics, demographics$subject_data_dirname %in% data_full_dem$subject_id);
    data_full_dem$subject_id = NULL; # Remove the subject_id column as well.

    # Also remove the matching weights (which are all 1 anyways for cardinality matching)
    data_full_matched$weights = NULL;

    # Create the matched data.matrix from the (now purely numeric) data.frame.
    glm_data = data.matrix(data_full_matched);
}

################################################################################
###### Look at effect of sex for predicting NI data for atlas regions ##########
################################################################################


# create model matrix using factor, needs to be the same model as used for the group comparison in matlab-script
mm <- model.matrix(~ demographics$sex + demographics$AGE + factor(demographics$site) + factor(demographics$qualification));

predictors <- c("Sex", "Age", "Site", "Qualification");
slm.F <- slmtools::slm_effect_sizes(mm, glm_data, predictors, c("cohens.f", "etasq", "power")); # the slmtools package is in my neuroimaging repo as a tar.gz

#visualization of per vertex effect sizes of group on fsaverage6, change expression in slm.F$cohens.F[] to plot the data of other effects, e.g. SA, site etc.
for(pred_idx in seq.int(length(predictors))) {
    predictor = predictors[pred_idx];
    cm = fsbrain::vis.data.on.subject(subjects_dir, "fsaverage", morph_data_both = slm.F$cohens.f[pred_idx,], views=NULL);
    output_img = sprintf("cohenf_%s.png", predictor);
    cat(sprintf("Writing figure to file: %s\n", output_img));
    fsbrain::export(cm, output_img = output_img, colorbar_legend = predictor);
}


