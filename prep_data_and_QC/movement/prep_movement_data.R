#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

github_dir <- "~/github/"
data_path <- "~/data/"

# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}

library(tidyverse)
library(glue)

# Define task rest BOLD confounds data dir
confounds_dir <- glue("{data_path}/UCLA_CNP/raw_data/confounds_data/")
movement_dir <- glue("{data_path}/UCLA_CNP/movement_data/")
confounds_files <- list.files(confounds_dir)

for (file in confounds_files) {
  subject <- gsub("_task-rest_bold_confounds.tsv", "", file)
  
  subject_motion_data <- read.table(glue(confounds_dir, file),
                                    header=TRUE) %>%
    dplyr::select(X,Y,Z,RotX,RotY,RotZ) %>%
    mutate(Frame = row_number()) %>%
    mutate(X_diff = X - X[Frame==1],
              Y_diff = Y - Y[Frame==1],
              Z_diff = Z - Z[Frame==1],
              RotX_diff = RotX - RotX[Frame==1],
              RotY_diff = RotY - RotY[Frame==1],
              RotZ_diff = RotZ - RotZ[Frame==1])
}
example_movement_data <- read.table(glue("{movement_dir}/sub-10171_movData.txt"), header = F)

write.table(subject_motion_data, 
            file=glue("{data_path}/UCLA_CNP/raw_data/confounds_data/sub-10159_movement_from_confounds.csv"),
            row.names=FALSE,
            sep=",")