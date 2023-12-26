library(tidyverse) 
library(broom)
library(circlize)
library(colorspace)
library(ComplexHeatmap)
library(correctR)
library(cowplot)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(ggnewscale)
library(ggpp)
library(ggpubr)
library(ggraph)
library(ggridges)
library(ggseg)
library(ggsegDefaultExtra)
library(ggsegHO)
library(ggsignif)
library(glue)
library(igraph)
library(LaCroixColoR)
library(patchwork)
library(psych)
library(RColorBrewer)
library(rlist)
library(scales)
library(see)
library(splitstackshape)
library(tidyverse)
library(viridis)

python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)
library(reticulate)
# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE"),
                             Disorder = c("SCZ", "BP", "ADHD", "ASD"))
pairwise_TS_feature_info <- read.csv("feature_info/pairwise_feature_info.csv")

UCLA_CNP_metadata <- pyarrow_feather$read_feather("input_data/UCLA_CNP_sample_metadata_filtered.feather") %>%
  mutate(Study="UCLA_CNP")

ABIDE_metadata <- pyarrow_feather$read_feather("input_data/ABIDE_sample_metadata_filtered.feather") %>%
  mutate(Study="ABIDE")

BP_combo_balacc <-  pyarrow_feather$read_feather("classification_results/balanced_accuracy/UCLA_CNP_BP_Univariate_catch25_Pairwise_pyspi14_Linear_SVM_balanced_accuracy.feather")

all_balanced_accuracy_results <- pyarrow_feather$read_feather("classification_results/all_balanced_accuracy_results.feather")
all_p_values <-  pyarrow_feather$read_feather("classification_results/all_p_values.feather")

pairwise_results <- all_balanced_accuracy_results %>%
  filter(Analysis_Type == "Pairwise_SPI") %>%
  left_join(all_p_values %>%
              filter(Analysis_Type == "Pairwise_SPI") %>%
              dplyr::select(Study:group_var, p_value_HolmBonferroni), 
            by = join_by(group_var, Analysis_Type, Disorder, Study)) %>%
  distinct()

combined_univariate_pairwise_results <- all_balanced_accuracy_results %>% 
  filter(Analysis_Type == "Univariate_Pairwise_Combo") %>%
  left_join(all_p_values %>%
              filter(Analysis_Type == "Univariate_Pairwise_Combo") %>%
              dplyr::select(Study:group_var, p_value_HolmBonferroni), 
            by = join_by(group_var, Analysis_Type, Disorder, Study)) %>%
  distinct()

results_df = plyr::rbind.fill(pairwise_results, 
                              combined_univariate_pairwise_results)

corrected_SPI_T_res <- 1:nrow(study_group_df) %>%
  purrr::map_df(~ run_correctR_group(disorder = study_group_df$Disorder[.x],
                                     study = study_group_df$Study[.x],
                                     metadata = plyr::rbind.fill(UCLA_CNP_metadata, ABIDE_metadata),
                                     results_df = results_df)) %>%
  left_join(., pairwise_TS_feature_info, by=c("SPI"="pyspi_name"))

relevant_p_values <- all_p_values %>%
  filter(Analysis_Type %in% c("Pairwise_SPI", "Univariate_Pairwise_Combo")) %>%
  dplyr::rename("SPI" = group_var) %>%
  semi_join(., corrected_SPI_T_res %>% dplyr::select(SPI, Disorder), by = join_by(Disorder, SPI)) %>%
  left_join(., corrected_SPI_T_res, by = join_by(Disorder, SPI)) %>%
  mutate(Analysis_Type = ifelse(Analysis_Type == "Pairwise_SPI", "FC", "FC + Local\nDynamics"),
         individually_significant = ifelse(p_value_HolmBonferroni < 0.05, "pcorr < 0.05", "Not significant"),
         significant_diff_with_univariate = ifelse(p_value_corr_HolmBonferroni < 0.05, "Sig Diff", "No Sig Diff"),
         Analysis_Type = factor(Analysis_Type, levels=c("FC", "FC + Local\nDynamics"))) %>%
  mutate(Disorder = factor(Disorder, levels = c("SCZ", "BP", "ADHD", "ASD"))) 


relevant_p_values%>%
  ggplot(data=., mapping=aes(x=Analysis_Type, y=100*Balanced_Accuracy_Across_Folds,
                             group = SPI)) +
  geom_line(aes(color = Disorder, 
                alpha = significant_diff_with_univariate), show.legend = FALSE) +
  scale_alpha_manual(values=c("Sig Diff" = 1, "No Sig Diff" = 0.2)) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#9d60a8", 
                              "BP" = "#2F77C0", 
                              "ADHD" = "#e45075", 
                              "ASD" = "#E28328")) +
  facet_wrap(Disorder ~ ., ncol=1, scales="fixed", strip.position = "left") +
  new_scale_colour() +  # start a new scale
  geom_point(aes(color = individually_significant)) +
  scale_color_manual(values = c("gray60", "#5BB67B")) +
  # scale_x_discrete(labels = wrap_format(7)) +
  xlab("Analysis Type") +
  ylab("Mean Balanced Accuracy (%)") +
  scale_y_continuous(breaks=c(45, 55, 65)) +
  scale_x_discrete(expand=c(0,0.2,0,0.2)) +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.text.y.left = element_text(angle=0, face="bold"),
        plot.margin = margin(1,30,1,1, unit="pt"),
        legend.title = element_blank())