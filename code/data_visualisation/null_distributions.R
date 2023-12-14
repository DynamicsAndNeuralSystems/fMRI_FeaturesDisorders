################################################################################
# Load libraries
################################################################################

python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(reticulate)
library(tidyverse)
library(icesTAF)
library(cowplot)
library(theft)
library(glue)
library(ggridges)
library(scales)
library(patchwork)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- glue("{github_dir}/plots/Manuscript_Draft/methods_supplement/null_visualisations/")
TAF::mkdir(plot_path)

data_path <- "~/data/TS_feature_manuscript/"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/processed_data/"
univariate_feature_set <- "catch25"

study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))


# UNIVARIATE 
# Univariate regional nulls
univariate_nulls <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))  
univariate_regional_nulls <- univariate_nulls %>% filter(Analysis_Type == "Univariate_Brain_Region") %>%
  mutate(Null_Balanced_Accuracy = 100*Null_Balanced_Accuracy) %>%
  group_by(Analysis_Type, Comparison_Group, Study, group_var) %>%
  mutate(mu = mean(Null_Balanced_Accuracy),
         sigma = sd(Null_Balanced_Accuracy))


# Plot in one ridgeline plot
UCLA_CNP_summary_curves <- univariate_regional_nulls %>%
  filter(Study=="UCLA_CNP") %>%
  left_join(., study_group_df) %>%
  group_by(group_var, Group_Nickname) %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD"))) %>%
  reframe(x = Null_Balanced_Accuracy,
          y = dnorm(Null_Balanced_Accuracy, 
                    mean = mean(Null_Balanced_Accuracy),
                    sd = sd(Null_Balanced_Accuracy)))

UCLA_CNP_univariate_null_plot <- univariate_regional_nulls %>%
  left_join(., study_group_df) %>%
  filter(Study=="UCLA_CNP") %>%
  mutate(Group_Nickname = factor(Group_Nickname, levels=c("SCZ", "BP", "ADHD"))) %>%
  ggplot(data=., mapping=aes(fill = group_var, x=Null_Balanced_Accuracy)) +
  geom_histogram(aes(y = after_stat(density)), alpha=0.7) +
  facet_grid(group_var~Group_Nickname, scales="free", switch="y") +
  xlab("Null Balanced Accuracy (%)") +
  ylab("Brain Regions") +
  geom_line(data=UCLA_CNP_summary_curves, aes(x=x, y=y), color="black") +
  geom_vline(xintercept = 50, linetype=2) +
  theme(legend.position="none",
        strip.text.y.left = element_text(angle=0, size=9, lineheight = 0.6, hjust=1),
        strip.background = element_blank(),
        strip.placement="outside",
        panel.spacing = unit(-0.25, "lines"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

# Calculate the normal curve density for each group
ABIDE_summary_curves <- univariate_regional_nulls %>%
  filter(Study=="ABIDE_ASD") %>%
  left_join(., study_group_df) %>%
  group_by(group_var) %>%
  reframe(x = Null_Balanced_Accuracy,
          y = dnorm(Null_Balanced_Accuracy, 
                    mean = mean(Null_Balanced_Accuracy),
                    sd = sd(Null_Balanced_Accuracy)))

ABIDE_univariate_null_plot <- univariate_regional_nulls %>%
  filter(Study=="ABIDE_ASD") %>%
  left_join(., study_group_df) %>%
  ggplot(data=., mapping=aes(fill = group_var, x=Null_Balanced_Accuracy)) +
  geom_histogram(aes(y = after_stat(density)), alpha=0.7) +
  facet_grid(group_var~Group_Nickname, scales="free", switch="y",
             labeller = labeller(group_var = label_wrap_gen(45))) +
  xlab("Null Balanced Accuracy (%)") +
  ylab("Brain Regions") +
  geom_line(data=ABIDE_summary_curves, aes(x=x, y=y, fill=group_var), color="black") +
  geom_vline(xintercept = 50, linetype=2) +
  theme(legend.position="none",
        strip.text.y.left = element_text(angle=0, size=9, lineheight = 0.6, hjust=1),
        strip.background = element_blank(),
        strip.placement="outside",
        panel.spacing = unit(-0.25, "lines"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

wrap_plots(list(UCLA_CNP_univariate_null_plot, ABIDE_univariate_null_plot), widths = c(0.75, 0.25))

ggsave(glue("{plot_path}/univariate_null_distributions.svg"), 
       width = 12, height = 10, units="in", dpi=300)
