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
library("parallel");  # run stuff in parallel on several CPU cores.
library("bettermc");  # better replacement for parallel::mclapply

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


data_full = fsbrain::group.morph.standard(subjects_dir, subjects_list, measure, hemi="both", fwhm="10", df_t=TRUE);
considered_vertexcolnames = colnames(data_full);

# add demographics data.
data_full_dem = data_full;
data_full_dem$sex = demographics$sex;
data_full_dem$age = demographics$AGE;
data_full_dem$height = demographics$HEIGHT;
data_full_dem$subject_id = demographics$subject_data_dirname;


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
data_full_matched = NULL; # free some RAM
t.test(glm_data$age[glm_data$sex == "male"], glm_data$age[glm_data$sex == "female"]);


################################################################################
###### Look at effect of sex for predicting NI data for atlas regions ##########
################################################################################


num_verts = length(considered_vertexcolnames);
num_metadata_columns_at_end = 7L; # Columns for group, sex, age, ...
num_columns = length(considered_vertexcolnames) + num_metadata_columns_at_end;

num_cores = 44L;
options("mc.cores" = num_cores);

# IMPORTANT:
# The 'chunk_size' setting below is absolutely critical for the performance! Set to high, your system will run out of
# memory and start swapping to disk, which means it will be ORDERS OF MAGNITUDE slower. Setting it too low is a massive
# waste of CPU resources (cores idling). The correct setting depends on the number of cores and RAM of your machine.
# Use proper system monitoring tools (like 'atop' under Linux) that show you when bottlenecks occur, and keep the settings
# just below the value at which they first show up.
chunk_size = num_cores * 500L;  # This works well on my 48 core Linux machine with 128 GB of RAM.


vert_indices = 1L:num_verts;
chunks = split(vert_indices, ceiling(seq_along(vert_indices)/chunk_size));

cat(sprintf("Starting model fitting (%d verts total, %d chunks with chunk size %d) with %d cores at:\n", num_verts, length(chunks), chunk_size, num_cores));
print(Sys.time());


includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait_for_zombies_to_die <- inline::cfunction(body=code, includes=includes, convention='.C')


all_chunks_effect_sizes = c();
all_metadata_column_indices = seq(((num_columns - num_metadata_columns_at_end) + 1L), num_columns);

for(chunk_idx in seq_along(chunks)) {

  vertex_indices_of_this_chunk = chunks[[chunk_idx]];

  # Add the meta data columns (age, sex, ...) from the end. They are required for the GLM.
  vertex_indices_of_this_chunk = c(vertex_indices_of_this_chunk, all_metadata_column_indices); # Column indices of chunk data in full glm_data.

  glm_data_chunk = glm_data[, vertex_indices_of_this_chunk]; # Matrix with dim: num_subjects x chunk_size

  # Fix column names for metadata columns (copy from glm_data) for the names formula to work.
  num_columns_chunk = ncol(glm_data_chunk);
  chunk_first_metadata_column_index = (num_columns_chunk - num_metadata_columns_at_end) + 1L; # These indices are for glm_data.
  chunk_all_metadata_column_indices = seq(chunk_first_metadata_column_index, num_columns_chunk);
  colnames(glm_data_chunk)[chunk_all_metadata_column_indices] = colnames(glm_data)[all_metadata_column_indices];

  fit_model_effect_size <- function(vertex_idx) {
    vertex_colname = considered_vertexcolnames[vertex_idx];
    formula_vertex = sprintf("%s ~ sex + age", vertex_colname);
    fit = glm(formula = formula_vertex, data = glm_data_chunk, family=gaussian());
    raw_sd_male = sd(glm_data_chunk[[vertex_colname]][glm_data_chunk$sex == "male"]);
    raw_sd_female = sd(glm_data_chunk[[vertex_colname]][glm_data_chunk$sex == "female"]);
    raw_sd_pooled = sqrt((raw_sd_male * raw_sd_male + raw_sd_female + raw_sd_female) / 2.0);
    effect_sex_male_mean = effects::effect("sex", fit)$fit[1];
    effect_sex_female_mean = effects::effect("sex", fit)$fit[2];
    cohen_d = (effect_sex_male_mean - effect_sex_female_mean) / raw_sd_pooled;
    return(abs(cohen_d));
  }

  cat(sprintf(" * Starting model fitting for chunk %d of %d (chunk size: %d) with %d cores at:\n", chunk_idx, length(chunks), chunk_size, num_cores));
  print(Sys.time());
  cat(sprintf("    - Data size for this chunk is: %d subjects x %d vertices (verts include %d metadata columns).\n", dim(glm_data_chunk)[1], dim(glm_data_chunk)[2], num_metadata_columns_at_end));


  res_list_effect_sizes_this_chunk = bettermc::mclapply( 1L:num_verts, fit_model_effect_size, mc.cores = num_cores, mc.progress = TRUE );
  wait_for_zombies_to_die();
  effect_sizes_this_chunk = unlist(res_list_effect_sizes_this_chunk);
  all_chunks_effect_sizes = c(all_chunks_effect_sizes, effect_sizes_this_chunk);

  cat(sprintf("    - Model fitting for chunk %d of %d (chunk size: %d) with %d cores done at:\n", chunk_idx, length(chunks), chunk_size, num_cores));
  print(Sys.time());
}


cat(sprintf("Model fitting (%d verts total, %d chunks with chunk size %d) with %d cores done at:\n", num_verts, length(chunks), chunk_size, num_cores));
print(Sys.time());


effect_sizes_file = "~/effects.mgh"
freesurferformats::write.fs.morph(effect_sizes_file, all_chunks_effect_sizes);
cat(sprintf("Wrote results to file '%s'.\n", effect_sizes_file));

#####

do_run_sequential_version = FALSE;
if(do_run_sequential_version) {
    vertex_idx = 1L;
    vertex_fits = list();
    pvalues_sex = rep(NA, num_verts);
    effect_sizes_sex = rep(NA, num_verts);
    for(vertex_colname in considered_vertexcolnames) {
        cat(sprintf("### Handling vertex '%s' (%d of %d). ###\n", vertex_colname, vertex_idx, length(considered_vertexcolnames)));
        formula_vertex = sprintf("%s ~ sex + age", vertex_colname);
        fit = glm(formula = formula_vertex, data = glm_data, family=gaussian());
        vertex_fits[[vertex_colname]] = fit;
        pvalues_sex[vertex_idx] = unname(coef(summary.glm(vertex_fits[[vertex_colname]]))[2,4]);

        raw_sd_male = sd(glm_data[[vertex_colname]][glm_data$sex == "male"]);
        raw_sd_female = sd(glm_data[[vertex_colname]][glm_data$sex == "female"]);
        raw_sd_pooled = sqrt((raw_sd_male * raw_sd_male + raw_sd_female + raw_sd_female) / 2.0);
        effect_sex_male_mean = effects::effect("sex", fit)$fit[1];
        effect_sex_female_mean = effects::effect("sex", fit)$fit[2];
        cohen_d = (effect_sex_male_mean - effect_sex_female_mean) / raw_sd_pooled;
        effect_sizes_sex[vertex_idx] = abs(cohen_d); # we are not interested in direction for effect size.

        vertex_idx = vertex_idx + 1L;
    }
}

#### Visualize results #####


fsbrain::vis.data.on.fsaverage(morph_data_both = effect_sizes_sex, draw_colorbar = T);



##### Use fastglm instead #####

do_try_fastglm_and_qdec = FALSE;

if(do_try_fastglm_and_qdec) {

  ## here we test for a single vertex
  #ffit = fastglm::fastglm(x=as.matrix(data_full$V1), y=as.integer(data_full_dem$sex));#


  ##### Use QDECR instead. It does not support interaction though, so its use seems limited.
  ##### It seems it uses RcppEigen::fastLm internally, and that may (?) support interactions, gotta double-check.
  demographics$age = demographics$AGE; # rename only
  demographics$subject_id = demographics$subject_data_dirname;
  fitqd <- QDECR::qdecr_fastlm(qdecr_thickness ~ age + sex,
                             data = demographics, # I think it only want the demographics here and loads the NI data itself.
                             id = "subject_id",
                             hemi = "lh",
                             project = "test_project",
                             dir_tmp = "/dev/shm/",
                             dir_subj = subjects_dir,
                             dir_fshome = "~/software/freesurfer",
                             n_cores=40,
                             clobber = TRUE);
}
