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
movement_dir <- glue("{data_path}/UCLA_CNP/movement_data/fmriprep/")
confounds_files <- list.files(confounds_dir)

for (file in confounds_files) {
  subject <- gsub("_task-rest_bold_confounds.tsv", "", file)
  
  subject_motion_data <- read.table(glue(confounds_dir, file),
                                    header=TRUE) %>%
    dplyr::select(X,Y,Z,RotX,RotY,RotZ)
  
  write.table(subject_motion_data, glue("{movement_dir}/{subject}_movData.txt"), 
              row.names=FALSE, col.names=FALSE, sep="\t")
}