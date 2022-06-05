github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
setwd(github_dir)

# Step 1

# Step 2

# Step 3 -- catch22
rmarkdown::render(paste0(github_dir,"3_ROI_wise_analysis/catch22/Step3_catch22_README.Rmd"))
file.copy(paste0(github_dir,"3_ROI_wise_analysis/catch22/Step3_catch22_README.md"),
          paste0(github_dir,"3_ROI_wise_analysis/catch22/README.md"),
          overwrite = T)

# Step 3 -- catchaMouse16
rmarkdown::render(paste0(github_dir,"3_ROI_wise_analysis/catchaMouse16/Step3_catchaMouse16_README.Rmd"))
file.copy(paste0(github_dir,"3_ROI_wise_analysis/catchaMouse16/Step3_catchaMouse16_README.md"),
          paste0(github_dir,"3_ROI_wise_analysis/catchaMouse16/README.md"),
          overwrite = T)

# Step 4
rmarkdown::render(paste0(github_dir,"4_Feature_wise_analysis/Step4_README.Rmd"))
file.copy(paste0(github_dir,"4_Feature_wise_analysis/Step4_README.md"),
          paste0(github_dir,"4_Feature_wise_analysis/README.md"),
          overwrite = T)
# Step 5
rmarkdown::render(paste0(github_dir,"5_ROI_and_feature_wise_analysis/Step5_README.Rmd"))
file.copy(paste0(github_dir,"5_ROI_and_feature_wise_analysis/Step5_README.md"),
          paste0(github_dir,"5_ROI_and_feature_wise_analysis/README.md"),
          overwrite = T)
# Step 6

# Step 7