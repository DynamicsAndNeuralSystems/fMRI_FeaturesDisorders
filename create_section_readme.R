github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
setwd(github_dir)

# Step 1

# Step 2

# Step 3
rmarkdown::render(paste0(github_dir,"3_ROI_wise_analysis/STEP3_README.Rmd"))
file.copy(paste0(github_dir,"3_ROI_wise_analysis/STEP3_README.md"),
          paste0(github_dir,"3_ROI_wise_analysis/README.md"),
          overwrite = T)
# Step 4
rmarkdown::render(paste0(github_dir,"4_Feature_wise_analysis/STEP4_README.Rmd"))
file.copy(paste0(github_dir,"4_Feature_wise_analysis/STEP4_README.md"),
          paste0(github_dir,"4_Feature_wise_analysis/README.md"),
          overwrite = T)
# Step 5
rmarkdown::render(paste0(github_dir,"5_ROI_and_feature_wise_analysis/STEP5_README.Rmd"))
file.copy(paste0(github_dir,"5_ROI_and_feature_wise_analysis/STEP5_README.md"),
          paste0(github_dir,"5_ROI_and_feature_wise_analysis/README.md"),
          overwrite = T)
# Step 6

# Step 7