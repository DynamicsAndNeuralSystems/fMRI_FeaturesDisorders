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
    left_join(., atlas_data) %>%
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
    ) +
    guides(col = guide_legend(override.aes = list(color="black")))
  
  return(ggseg_plot)
}

# Function to plot specific data in a given iteration
plot_network_data <- function(edge_color) {
  
  # connect = dataframe of pairwise correlations between cortical ROIs
  connect <- data.frame("from" = sample(rois$to, 500, replace=T),
                        "to" = sample(rois$to,500, replace=T),
                        value = runif(500, 0, 1)) %>%
    arrange(from, to) %>%
    sample_n(100)
  
  # Collect edges where from = selected ROI and to = ROIs connected to the selected ROI
  from <- match(connect$from, vertices$name)
  to <- match(connect$to, vertices$name)
  
  # mygraph = igraph object linking each cortical ROI
  # convert to a circular dendrogram-shaped ggraph object
  p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    theme_void() +  geom_conn_bundle(data = get_con(from = from, to = to, 
                                                    value=connect$value), 
                                     tension=runif(1, 0.55, 0.85), 
                                     width=1,
                                     aes(color=value,
                                         alpha=value))  +
    labs(edge_width="Pearson\nCorrelation") + 
    geom_node_point(aes(filter = leaf, 
                        x = x*1.05, y=y*1.05),   
                    size=3) +
    scale_edge_color_gradientn(colors=c(alpha(edge_color,0.3), edge_color)) +
    theme_void() + 
    theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
          legend.position="none")
  
  return(p)
}

plot_feature_in_brain <- function(fill_color_gradient, region_label="all") {
  if (region_label=="all") {
    p <- dk %>%
      as_tibble() %>%
      mutate(region_values = runif(nrow(.))) %>%
      ggseg(atlas = "dk", mapping = aes(fill = region_values),
            hemisphere="left",
            view = "lateral",
            position = "stacked", colour = "black") +
      scale_fill_gradientn(colors=c(alpha(fill_color_gradient, 0.3), 
                                    fill_color_gradient), 
                           na.value=NA)
  } else {
    p <- dk %>%
      as_tibble() %>%
      mutate(region_values = ifelse(label==region_label, "1", NA_character_)) %>%
      ggseg(atlas = "dk", mapping = aes(fill = region_values),
            hemisphere="left",
            view = "lateral",
            position = "stacked", colour = "gray40") +
      scale_fill_manual(values=c(fill_color_gradient),
                        na.value="white")
  }

  p <- p  +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.position = "none") 
}


# Plot balanced accuracy in the brain
plot_balacc_in_brain <- function(significant_univariate_region_wise_results, 
                                 study_group_df,
                                 ABIDE_brain_region_info,
                                 color_palette=c("#FFEE75", "#FCA769", "#fb6555", "#D32345", "#401057"),
                                 bin_seq_range=seq(50,75,by=5)) {
  
  # Find max fill and min fill values
  min_fill <- floor(min(significant_univariate_region_wise_results$Balanced_Accuracy_Across_Folds))
  max_fill <- ceiling(max(significant_univariate_region_wise_results$Balanced_Accuracy_Across_Folds))
  
  # Initialize list of ggseg plots
  ggseg_plot_list <- list()
  
  # First plot within brain using ggseg
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Disorder[i]
    
    # Define atlas by study
    atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
    
    # If dataset is ABIDE ASD, convert regions to ggseg regions
    if (dataset_ID == "ABIDE") {
      significant_data_for_ggseg <- significant_univariate_region_wise_results %>%
        filter(Study == dataset_ID,
               Disorder == comparison_group) %>%
        dplyr::rename("Brain_Region" = "group_var") %>%
        left_join(., ABIDE_brain_region_info)
      
    } else {
      # Extract sig results to plot
      significant_data_for_ggseg <- significant_univariate_region_wise_results %>%
        filter(Study == dataset_ID,
               Disorder == comparison_group) %>%
        distinct() %>%
        mutate(label = ifelse(str_detect(group_var, "ctx-"),
                              gsub("-", "_", group_var),
                              as.character(group_var))) %>%
        mutate(label = gsub("ctx_", "", label))
    }
    
    
    # Plot balanced accuracy data in cortex
    dataset_ggseg <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                   atlas_name = atlas,
                                                   atlas_data = get(atlas) %>% as_tibble(),
                                                   data_to_plot = significant_data_for_ggseg,
                                                   fill_variable = "Balanced_Accuracy_Across_Folds",
                                                   fill_colors = color_palette,
                                                   bin_seq = bin_seq_range,
                                                   line_color = "gray30",
                                                   na_color = "white") +
      labs(fill = "Mean Balanced Accuracy (%)") +
      theme(plot.title = element_blank())
    
    # Append to list
    ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
    
    # Add subcortical data for UCLA CNP
    if (dataset_ID == "UCLA_CNP") {
      dataset_ggseg_subctx <- plot_data_with_ggseg_discrete(dataset_ID = dataset_ID,
                                                            atlas_name = "aseg",
                                                            atlas_data = aseg %>% as_tibble(),
                                                            data_to_plot = significant_data_for_ggseg,
                                                            fill_variable = "Balanced_Accuracy_Across_Folds",
                                                            fill_colors = color_palette,
                                                            bin_seq = bin_seq_range,
                                                            line_color = "gray30",
                                                            na_color = "white") +
        labs(fill = "Mean Balanced Accuracy (%)") +
        theme(plot.title = element_blank()) 
      # Append to list
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
    }
  }
  return(ggseg_plot_list)
}


repkfold_ttest <- function(data, n1, n2, k, r){
  
  # Arg checks
  
  '%ni%' <- Negate('%in%')
  
  if("model" %ni% colnames(data) || "values" %ni% colnames(data) || "k" %ni% colnames(data) || "r" %ni% colnames(data)){
    stop("data should contain at least four columns called 'model', 'values', 'k', and 'r'.")
  }
  
  if(!is.numeric(data$values) || !is.numeric(data$k) || !is.numeric(data$r)){
    stop("data should be a data.frame with only numerical values in columns 'values', 'k', and 'r'.")
  }
  
  if(!is.numeric(n1) || !is.numeric(n2) || !is.numeric(k) || !is.numeric(r) ||
     length(n1) != 1 || length(n2) != 1 || length(k) != 1 || length(r) != 1){
    stop("n1, n2, k, and r should all be integer scalars.")
  }
  
  if(length(unique(data$model)) != 2){
    stop("Column 'model' in data should only have two unique labels (one for each model to compare).")
  }
  
  # Calculations
  
  d <- c()
  
  for(i in 1:k){
    for(j in 1:r){
      x <- data[data$k == i, ]
      x <- x[x$r == j, ]
      d <- c(d, x[x$model == unique(x$model)[1], c("values")] - x[x$model == unique(x$model)[2], c("values")]) # Differences
    }
  }
  
  # Catch for when there is zero difference(s) between the models
  
  if (sum(unlist(d)) == 0) {
    tmp <- data.frame(statistic = 0, p.value = 1)
  } else{
    
    statistic <- mean(unlist(d), na.rm = TRUE) / sqrt(stats::var(unlist(d), na.rm = TRUE) * ((1/(k * r)) + (n2/n1))) # Calculate t-statistic
    df <- n1 + n2 - 2
    
    if(statistic < 0){
      p.value <- stats::pt(statistic, (k * r) - 1) # p-value for left tail
    } else{
      p.value <- stats::pt(statistic, (k * r) - 1, lower.tail = FALSE) # p-value for right tail
    }
    
    tmp <- data.frame(statistic = statistic, p.value = p.value, df = df)
  }
  
  return(tmp)
}

run_correctR_group <- function(disorder, study, metadata, results_df) {
  # Find number of subjects for the specified comparison group
  num_subjects <- metadata %>%
    filter(Study == study, 
           Diagnosis %in% c("Control", disorder)) %>%
    distinct(Sample_ID) %>%
    nrow()
  
  # Compute the training and test fold sizes for 10-fold CV
  training_size <- ceiling(0.9*num_subjects)
  test_size <- floor(0.1*num_subjects)
  
  # Prep the resulting balanced accuracies with vs without univariate data
  data_for_correctR <- results_df %>%
    filter(Study == study, 
           Disorder == disorder) %>%
    group_by(group_var) %>%
    filter(any(p_value_HolmBonferroni < 0.05)) %>%
    ungroup() %>%
    dplyr::rename("model" = "Analysis_Type",
                  "SPI" = "group_var",
                  "k" = "Fold",
                  "r" = "Repeat_Number",
                  "values" = "Balanced_Accuracy") %>%
    dplyr::select(model, SPI, k, r, values) %>%
    dplyr::mutate(r = r + 1) %>%
    group_by(SPI) %>%
    group_split()

  res <- data_for_correctR %>%
    purrr::map_df(~ as.data.frame(repkfold_ttest(data = .x %>% dplyr::select(-SPI), 
                                                 n1 = training_size,
                                                 n2 = test_size,
                                                 k = 10,
                                                 r = 10)) %>%
                    mutate(SPI = unique(.x$SPI))) %>%
    ungroup() %>%
    dplyr::rename("p_value_corr"="p.value") %>%
    mutate(p_value_corr_HolmBonferroni = p.adjust(p_value_corr, method="bonferroni"),
           Disorder = disorder)
  
  return(res)
  
}


plot_fold_heatmap_for_dataset <- function(fold_assignments_df, group_var_to_use, plot_title) {
  p <- fold_assignments_df %>%
    filter(group_var == group_var_to_use) %>%
    mutate(Repeat = factor(Repeat),
           Fold = factor(Fold)) %>%
    ggplot(data=., mapping=aes(x = Sample_ID, y = Repeat, fill = Fold )) +
    facet_grid(. ~ Diagnosis, scales="free", space="free") +
    geom_tile() +
    xlab("Samples") +
    ylab("CV-SVM Repeat") +
    labs(fill = "k-fold") +
    ggtitle(plot_title) +
    scale_fill_viridis_d() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0.5),
          legend.position = "bottom") +
    guides(fill = guide_legend(title.position = "top", 
                               nrow = 1,
                               title.hjust = 0.5,
                               label.position = "bottom")) 
  
  return(p)
}

plot_group_vs_control_mean_FD <- function(FD_dataset,
                                          study,
                                          dx,
                                          dx_title,
                                          ymin,
                                          ymax,
                                          group_color) {
  p <- FD_dataset %>%
    filter(Study==study,
           Diagnosis %in% c("Control", dx)) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(dx, "Control"))) %>%
    ggplot(data=., mapping=aes(x=Diagnosis, y=as.numeric(Mean_FD_Power))) +
    geom_violinhalf(aes(fill=Diagnosis), scale="width", 
                    position = position_nudge(x=0.2))  +
    geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
                 fill=NA, color="white",
                 position = position_nudge(x=0.27), coef = 0) +
    geom_point(aes(color = Diagnosis), position = position_jitterdodge(dodge.width = 1,
                                                                       jitter.width = 0.5),
               size = 1, alpha=0.7) +
    geom_signif(test = "wilcox.test",
                comparisons = list(c(dx, "Control")), 
                y_position = c(0.55, 0.6),
                map_signif_level=TRUE) +
    scale_fill_manual(values = c(group_color, control_color)) +
    scale_color_manual(values = c(group_color, control_color)) +
    scale_y_continuous(limits = c(ymin, ymax), expand = c(0,0,0.1,0)) +
    scale_x_discrete(labels = c(dx_title, "Control")) +
    ggtitle(dx_title) +
    ylab("Mean FD (Power)") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5, size=14))
  return(p)
}