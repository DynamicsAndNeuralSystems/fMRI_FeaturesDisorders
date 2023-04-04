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

# Plot cortical + subcortical data in a continuous gradient
plot_data_with_ggseg_gradient <- function(dataset_ID,
                                 atlas_name,
                                 atlas_data,
                                 data_to_plot,
                                 fill_variable,
                                 min_fill=NULL,
                                 max_fill=NULL,
                                 line_color="darkgrey",
                                 na_color="white",
                                 fill_colors=c("blue", "red")) {
  
  ggseg_data <- data_to_plot %>%
    left_join(., atlas_data %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    ungroup() %>%
    dplyr::select(-label)
  
  if (atlas_name == "aseg") {
    ggseg_plot <- ggseg_data %>%
      filter(type!="cortical") %>%
      ggplot() +
      geom_brain(atlas = aseg, mapping = aes_string(fill=fill_variable), 
                 side = "coronal", colour = line_color)  +
      theme_void() +
      theme(plot.title = element_blank()) 
  } else {
    ggseg_plot <- ggseg_data %>%
      ggseg(atlas = atlas_name, mapping = aes_string(fill = fill_variable),
            position = "stacked", colour = line_color) +
      theme_void() +
      theme(plot.title = element_blank()) 
  }
  
  if (!is.null(min_fill)) {
    ggseg_plot <- ggseg_plot +
      scale_fill_gradientn(colours = fill_colors, na.value=na_color,
                           limits = c(min_fill, max_fill)) 
  } else {
    ggseg_plot <- ggseg_plot +
      scale_fill_gradientn(colours = fill_colors, na.value=na_color) 
  }
  
  return(ggseg_plot)
}

# Plot cortical + subcortical data in a discretized gradient
plot_data_with_ggseg_discrete <- function(dataset_ID,
                                          atlas_name,
                                          atlas_data,
                                          data_to_plot,
                                          fill_variable,
                                          num_bins=5,
                                          line_color="darkgrey",
                                          na_color="white",
                                          bin_seq = NULL,
                                          fill_colors=NULL) {
  
  ggseg_data <- data_to_plot %>%
    left_join(., atlas_data %>% as_tibble()) %>%
    filter(!is.na(region)) %>%
    ungroup() %>%
    dplyr::select(-label)
  
  if (is.null(fill_colors)) {
    fill_colors = viridis::viridis(n=num_bins)
  }
  
  if (atlas_name == "aseg") {
    ggseg_plot <- ggseg_data %>%
      filter(type!="cortical") %>%
      ggplot() +
      geom_brain(atlas = aseg, mapping = aes_string(fill = fill_variable), 
                 side = "coronal", colour = line_color)
  } else {
    ggseg_plot <- ggseg_data %>%
      ggseg(atlas = atlas_name, mapping = aes_string(fill = fill_variable),
            position = "stacked", colour = line_color)
     
  }
  
  if (is.null(bin_seq)) {
    bin_seq <- seq(min(ggseg_data %>% pull(fill_variable)),
                      max(ggseg_data %>% pull(fill_variable)),
                      length.out = num_bins + 1)
  }
  
  ggseg_plot <- ggseg_plot + 
    theme_void() +
    theme(plot.title = element_blank()) +
    binned_scale(aesthetics = "fill",
                 scale_name = "stepsn", 
                 palette = function(x) fill_colors,
                 breaks = bin_seq,
                 limits = c(min(bin_seq), max(bin_seq)),
                 show.limits = TRUE, 
                 guide = "colorsteps"
    )
  
  return(ggseg_plot)
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



