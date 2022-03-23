#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 15 March 2022
#--------------------------------------

library(tidyverse)
library(cowplot) 
theme_set(theme_cowplot())

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Get breakdown of control vs schizophrenia subjects in cohort

get_dx_breakdown <- function(subject_csv) {
  dx_data <- read.csv(subject_csv) %>%
    mutate(diagnosis = str_to_sentence(diagnosis)) %>%
    filter(diagnosis %in% c("Control", "Schz"))
  num_total <- nrow(dx_data)
  num_scz <- nrow(filter(dx_data, diagnosis=="Schz"))
  num_ctrl <- nrow(filter(dx_data, diagnosis=="Control"))
  scz2ctrl <- num_scz/num_ctrl
  
  res_df <- data.frame(var1 = num_ctrl,
                       var2 = num_scz,
                       var3 = num_total,
                       var4 = scz2ctrl)
  colnames(res_df) <- c("Control", "Schizophrenia", "Total", "SCZ2Ctrl")
  print(res_df)
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot subject count by diagnosis vs. sex
plot_dx_vs_sex_count <- function(subject_csv) {
  subject_data <- read.csv(subject_csv) %>%
    mutate(diagnosis = str_to_sentence(diagnosis)) %>%
    filter(diagnosis %in% c("Control", "Schz"))
  
  subject_data %>%
    group_by(gender, diagnosis) %>%
    count() %>%
    ggplot(data=., mapping=aes(x=gender, fill=gender, y=n)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=n), vjust=-0.5) +
    facet_grid(diagnosis ~ ., switch="both") +
    scale_y_continuous(expand=c(0,0,0.15,0)) +
    ylab("# Subjects") +
    xlab("Sex") +
    theme(legend.position="none",
          strip.placement = "outside")
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot age by diagnosis vs. sex
plot_dx_vs_sex_age <- function(subject_csv) {
  subject_data <- read.csv(subject_csv) %>%
    mutate(diagnosis = str_to_sentence(diagnosis)) %>%
    filter(diagnosis %in% c("Control", "Schz"))
  
  subject_data %>%
    ggplot(data=., mapping=aes(x=gender, y=age, fill=gender)) +
    geom_boxplot() +
    facet_grid(diagnosis ~ ., switch="both") +
    ylab("Age") +
    xlab("Sex") +
    theme(legend.position="none",
          strip.placement = "outside")
}