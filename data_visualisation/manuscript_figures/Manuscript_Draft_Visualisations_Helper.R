# Plot the distribution of each significant region or feature across n=10 repeats
# Show empirical null distribution as a shaded grey band
plot_boxplot_shaded_null <- function(dataset_ID,
                                     grouping_var_name,
                                     main_data_by_repeat,
                                     fill_color,
                                     null_mean_value,
                                     wrap_length = 20,
                                     null_SD_value) {
  ggplot() +
    geom_boxplot(data = main_data_by_repeat %>%
                   mutate(grouping_var = str_replace_all(grouping_var, "ctx-lh-|Left-", "Left ")) %>%
                   mutate(grouping_var = str_replace_all(grouping_var, "ctx-rh-|Right-", "Right ")) %>%
                   mutate(grouping_var = fct_reorder(grouping_var, 
                                                     balanced_accuracy,
                                                     .fun = mean)), 
                 aes(y=grouping_var, x=balanced_accuracy),
                 fill = fill_color, color="black") +
    geom_rect(data = main_data_by_repeat, 
              xmin = null_mean_value - null_SD_value,
              xmax = null_mean_value + null_SD_value,
              ymin=0, ymax=Inf, fill="gray85") +
    geom_vline(xintercept = null_mean_value) +
    scale_x_continuous(limits = c(0.975*(null_mean_value - null_SD_value),
                                  1.025*max(main_data_by_repeat %>% 
                                             filter(bal_acc_p_adj < 0.05) %>% 
                                             pull(balanced_accuracy))),
                       breaks = scales::breaks_pretty(n = 4)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_length)) +
    xlab("Balanced Accuracy") +
    ylab(grouping_var_name)
}

# Plot significant brain regions in a brain map using ggseg
plot_significant_regions_ggseg <- function(dataset_ID,
                                           atlas_name,
                                           atlas_data,
                                           main_data_for_ggseg,
                                           fill_color) {
  
  sig_region_ggseg_data <- main_data_for_ggseg %>%
    left_join(., atlas_data %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    ungroup() %>%
    dplyr::select(-label)
  
  sig_region_ggseg_data %>%
    ggseg(atlas = atlas_name, mapping = aes(fill = fillyes),
          position = "dispersed", colour = "darkgrey") +
    scale_fill_manual(values = c(fill_color), na.value = NA) +
    labs(fill = "Balanced Accuracy") +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.position = "none") 
}

# Plot univariate region vs combo results in violin plot
plot_univar_region_vs_combo_violin <- function(dataset_ID,
                                               ROI_data_full,
                                               ROI_pvals,
                                               combo_data,
                                               region_fill_color = "#F0224B",
                                               combo_fill_color = "#9B51B4") {
  
  # Get average of all brain regions, first by repeat and then across repeats
  ROI_data <- ROI_data_full %>%
    group_by(grouping_var, Noise_Proc, repeat_number) %>%
    summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
    group_by(grouping_var, Noise_Proc) %>%
    summarise(balanced_accuracy = mean(mean_balacc, na.rm=T)) %>%
    left_join(., ROI_pvals)
  
  ROI_data$Analysis <- "Region-wise"
  combo_data$Analysis <- "Combo-wise"
  
  merged_data <- do.call(plyr::rbind.fill, list(ROI_data, combo_data)) %>%
    mutate( Analysis = factor(Analysis, levels = c("Region-wise", 
                                                   "Combo-wise"))) 
  
  merged_data %>%
    ggplot(data=., mapping = aes(x = Analysis, y = balanced_accuracy)) +
    geom_violin(aes(fill = Analysis)) +
    stat_summary(data = subset(merged_data, Analysis == "Combo-wise"),
                 geom = "crossbar", fun = "mean", aes(color=Analysis), size=1,
                 color = combo_fill_color) +
    geom_boxplot(fill=NA, width=0.1, color="black") +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    scale_fill_manual(values = c(region_fill_color)) +
    ylab("Balanced Accuracy") +
    theme(legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.y.left = element_blank(),
          axis.title.x = element_blank())
}

# Plot pairwise results as a spaghetti scatter plot
# Each SPI with vs without univariate combo info available
plot_SPI_with_without_univar <- function(dataset_ID,
                                         SPI_data_main_repeats,
                                         SPI_combo_data_main_repeats,
                                         SPI_only_color = "chartreuse3",
                                         SPI_univar_color = "darkgoldenrod2") {
  # Merge data for paired T-test across repeats
  merged_SPI_combo_data <- SPI_data_main_repeats %>% 
    mutate(Analysis = "Pairwise SPI-wise") %>%
    plyr::rbind.fill(., SPI_combo_data_main_repeats %>% 
                       mutate(Analysis = "Pairwise SPI + Univariate Combo-wise")) %>%
    mutate(Analysis = factor(Analysis, levels = c("Pairwise SPI-wise",
                                                  "Pairwise SPI + Univariate Combo-wise")))
  
  # Average across repeats for plot only
  merged_SPI_combo_data_avg <- merged_SPI_combo_data %>%
    group_by(grouping_var, Analysis) %>% 
    summarise(mean_balacc = mean(balanced_accuracy, na.rm=T)) %>%
    ungroup()
  
  # Find SPIs with a significant difference with vs without univariate combo info
  sig_SPIs_diff <- merged_SPI_combo_data %>%
    nest(data = -grouping_var) %>%
    mutate(test = map(data, ~ t.test(balanced_accuracy ~ Analysis, data = .x, paired=T)),
           tidied = map(test, tidy)) %>%
    unnest(tidied) %>%
    ungroup() %>%
    mutate(p.adj = p.adjust(p.value, method="BH")) %>%
    filter(p.adj < 0.05) %>%
    pull(grouping_var)
  
  # Plot the results
  merged_SPI_combo_data_avg %>%
    ggplot(data=., mapping=aes(x = Analysis, y = mean_balacc)) +
    geom_point(aes(color = Analysis)) +
    scale_color_manual(values = c(SPI_only_color, SPI_univar_color)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15),
                     expand = c(0, 0.5, 0, 3)) +
    geom_line(data = subset(merged_SPI_combo_data_avg, 
                            grouping_var %in% sig_SPIs_diff),
              aes(x=Analysis, y=mean_balacc, group=grouping_var),
              color="black", alpha=0.5) +
    ggrepel::geom_text_repel(data = subset(merged_SPI_combo_data_avg,
                                           grouping_var %in% sig_SPIs_diff & Analysis == "Pairwise SPI + Univariate Combo-wise"),
                             aes(label = grouping_var),
                             hjust=0,
                             
                             seed         = 42,
                             force        = 0.5,
                             force_pull   = 0,
                             direction    = "y",
                             segment.size = 0.2,
                             box.padding  = 0.5,
                             
                             nudge_x = 0.5) +
    ylab("Balanced Accuracy") +
    theme(legend.position = "none",
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.y.left = element_blank(),
          axis.title.x = element_blank())
}