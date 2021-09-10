#!/usr/bin/env Rscript


library("fsbrain");
library("kernlab");
library("readxl");


#### Setup data ####

study_dir = "/Volumes/shared/projects/IXI_dataset/";
subjects_dir = file.path(study_dir, "IXI-T1/");

measure = "thickness";


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



##### Load data #####

data = fsbrain::group.morph.standard(subjects_dir, subjects_list, measure, fwhm = "10");






