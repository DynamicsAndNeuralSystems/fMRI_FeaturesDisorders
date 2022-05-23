# Step 1
setwd("D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders")
# Step 2

# Step 3
rmarkdown::render("3_ROI_wise_analysis/STEP3_README.Rmd", output_file = "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/3_ROI_wise_analysis/README.md")

# Step 4
rmarkdown::render("4_feature_wise_analysis/STEP4_README.Rmd", output_file = "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/4_feature_wise_analysis/README.md")

# Step 5
rmarkdown::render("5_ROI_and_feature_wise_analysis/STEP5_README.Rmd", output_file = "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/5_ROI_and_feature_wise_analysis/README.md")

# Step 6

# Step 7