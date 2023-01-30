python_to_use <- "~/.conda/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(reticulate)
library(theft)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

simData <- theft::simData
simData_catch22 <- calculate_features(simData,
                                     id_var = "id",
                                     time_var = "timepoint",
                                     values_var = "values",
                                     feature_set = "catch22")

################################################################################
# Robust sigmoid in R with theft
################################################################################

# Normalise with robust sigmoid method in theft
simData_catch22_norm_in_theft <- normalise_feature_frame(simData_catch22,
                                                         names_var = "names",
                                                         values_var = "values",
                                                         method = "RobustSigmoid") %>%
  dplyr::rename("Feature" = "names",
                "RobustSigmoid_theft" = "values") %>%
  dplyr::select(-method)

################################################################################
# Robust sigmoid in R with theft
################################################################################

# Cast data from long to wide for robust sigmoid analysis in python
simData_catch22_wide <- simData_catch22 %>%
  pivot_wider(id_cols = id, names_from = names, values_from = values)

# Extract the sample IDs
sample_IDs <- simData_catch22_wide$id

# Extract the catch22 feature names
catch22_data_feature_columns = colnames(simData_catch22_wide)[2:23]

# Convert data to matrix
simData_catch22_matrix <- as.matrix(simData_catch22_wide %>% dplyr::select(-id))

# Import robust sigmoid function from python script
source_python("~/github/fMRI_FeaturesDisorders/classification_analysis/core_classification_functions.py")

# Apply robust sigmoid transform with sklearn
transformer = RobustSigmoidScaler(unit_variance=TRUE)$fit(simData_catch22_matrix)
simData_catch22_norm_in_python_wide = as.data.frame(transformer$transform(simData_catch22_matrix))
colnames(simData_catch22_norm_in_python_wide) <- catch22_data_feature_columns 

simData_catch22_norm_in_python <- simData_catch22_norm_in_python_wide %>%
  mutate(id = sample_IDs) %>%
  pivot_longer(cols = c(-id), 
               names_to = "Feature", 
               values_to = "RobustSigmoid_python")

################################################################################
# Compare results with R vs Python
################################################################################

# Merge the data
merged_norm <- left_join(simData_catch22_norm_in_theft,
                         simData_catch22_norm_in_python)

merged_norm %>%
  ggplot(data=., mapping=aes(x=RobustSigmoid_theft, y=RobustSigmoid_python)) +
  geom_point() +
  ylab("Normalisation in Python") +
  xlab("Normalisation in R") +
  ggtitle("Robust Sigmoid Normalisation for theft::simData") +
  coord_equal() +
  theme(plot.title = element_text(hjust=0.5))

