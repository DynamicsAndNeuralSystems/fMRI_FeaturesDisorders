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
                   ungroup() %>%
                   mutate(group_var = str_replace_all(group_var, "ctx-lh-|Left-", "Left ")) %>%
                   mutate(group_var = str_replace_all(group_var, "ctx-rh-|Right-", "Right ")) %>%
                   mutate(group_var = fct_reorder(group_var, 
                                                  Balanced_Accuracy_Across_Folds,
                                                  .fun = mean,)), 
                 aes(y=group_var, x=Balanced_Accuracy_Across_Folds),
                 fill = fill_color, color="black") +
    geom_rect(data = main_data_by_repeat, 
              xmin = null_mean_value - null_SD_value,
              xmax = null_mean_value + null_SD_value,
              ymin=0, ymax=Inf, fill="gray85") +
    geom_vline(xintercept = null_mean_value) +
    scale_x_continuous(limits = c(0.975*(null_mean_value - null_SD_value),
                                  1.025*max(main_data_by_repeat %>% 
                                              ungroup() %>%
                                              filter(p_value_BH < 0.05) %>% 
                                              pull(Balanced_Accuracy_Across_Folds))),
                       breaks = scales::breaks_pretty(n = 4)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_length)) +
    xlab("Balanced Accuracy\nby CV Repeat") +
    ylab(grouping_var_name)
}

# Plot cortical + subcortical data in a divering gradient
plot_data_with_ggseg_diverging <- function(dataset_ID,
                                 atlas_name,
                                 atlas_data,
                                 data_to_plot,
                                 fill_variable,
                                 min_fill=NULL,
                                 max_fill=NULL,
                                 fill_palette="RdBu") {
  
  ggseg_data <- data_to_plot %>%
    left_join(., atlas_data %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    ungroup() %>%
    dplyr::select(-label)
  
  if (atlas_name == "aseg") {
    ggseg_plot <- ggseg_data %>%
      ggplot() +
      geom_brain(atlas = aseg, mapping = aes_string(fill=fill_variable), 
                 side = "coronal", colour = "darkgrey")  +
      scale_fill_continuous_divergingx(palette = fill_palette, 
                                       mid = 0, 
                                       rev = TRUE, 
                                       limits = c(min_fill, max_fill),
                                       na.value="gray90") +
      theme_void() +
      theme(plot.title = element_blank()) 
  } else {
    ggseg_plot <- ggseg_data %>%
      ggseg(atlas = atlas_name, mapping = aes_string(fill = fill_variable),
            position = "stacked", colour = "darkgrey") +
      scale_fill_continuous_divergingx(palette = fill_palette, 
                                       mid = 0, 
                                       rev = TRUE, 
                                       limits = c(min_fill, max_fill),
                                       na.value="gray90") +
      theme_void() +
      theme(plot.title = element_blank()) 
  }
  
  return(ggseg_plot)
}

# Plot cortical + subcortical data in a continuous gradient
plot_data_with_ggseg <- function(dataset_ID,
                                 atlas_name,
                                 atlas_data,
                                 data_to_plot,
                                 fill_variable,
                                 min_fill=NULL,
                                 max_fill=NULL,
                                 line_color="darkgrey",
                                 fill_colors=c("blue", "red")) {
  
  ggseg_data <- data_to_plot %>%
    left_join(., atlas_data %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    ungroup() %>%
    dplyr::select(-label)
  
  if (atlas_name == "aseg") {
    ggseg_plot <- ggseg_data %>%
      ggplot() +
      geom_brain(atlas = aseg, mapping = aes_string(fill=fill_variable), 
                 side = "coronal", colour = line_color)  +
      scale_fill_gradientn(colours = fill_colors, na.value="gray90") +
      theme_void() +
      theme(plot.title = element_blank()) 
  } else {
    ggseg_plot <- ggseg_data %>%
      ggseg(atlas = atlas_name, mapping = aes_string(fill = fill_variable),
            position = "stacked", colour = line_color) +
      scale_fill_gradientn(colours = fill_colors, na.value="gray90") +
      theme_void() +
      theme(plot.title = element_blank()) 
  }
  
  return(ggseg_plot)
}


# Plot significant brain regions in a brain map using ggseg
plot_significant_regions_ggseg <- function(dataset_ID,
                                           atlas_name,
                                           atlas_data,
                                           main_data_for_ggseg,
                                           min_fill=0,
                                           max_fill=1,
                                           fill_color) {
  
  sig_region_ggseg_data <- main_data_for_ggseg %>%
    left_join(., atlas_data %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    ungroup() %>%
    dplyr::select(-label)
  
  if (atlas_name == "aseg") {
    ggseg_plot <- sig_region_ggseg_data %>%
      ggplot() +
      geom_brain(atlas = aseg, mapping = aes(fill=Balanced_Accuracy_Across_Repeats), 
                 side = "coronal", colour = "darkgrey")  +
      scale_fill_gradientn(colors=c(alpha(fill_color, 0.3), fill_color), 
                           limits = c(min_fill, max_fill),
                           na.value=NA) +
      labs(fill = "Mean Balanced Accuracy (%)") +
      theme_void() +
      guides(fill = guide_colorbar(title.position = "top", 
                                   nrow = 1,
                                   barwidth = 10, 
                                   barheight = 0.75,
                                   title.hjust = 0.5,
                                   label.position = "bottom")) +
      theme(plot.title = element_blank(),
            legend.position = "bottom") 
  } else {
    ggseg_plot <- sig_region_ggseg_data %>%
      ggseg(atlas = atlas_name, mapping = aes(fill = Balanced_Accuracy_Across_Repeats),
            position = "stacked", colour = "darkgrey") +
      scale_fill_gradientn(colors=c(alpha(fill_color, 0.3), fill_color), 
                           limits = c(min_fill, max_fill),
                           na.value=NA) +
      labs(fill = "Mean Balanced Accuracy (%)") +
      theme_void() +
      guides(fill = guide_colorbar(title.position = "top", 
                                   nrow = 1,
                                   barwidth = 10, 
                                   barheight = 0.75,
                                   title.hjust = 0.5,
                                   label.position = "bottom")) +
      theme(plot.title = element_blank(),
            legend.position = "bottom") 
  }
  
  
  return(ggseg_plot)
}

