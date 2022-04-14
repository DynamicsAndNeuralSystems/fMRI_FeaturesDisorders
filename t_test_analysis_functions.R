#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Very simple t-test
t_test_by_region <- function(rdata_path, 
                             noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                             norm_methods = c("z-score", "RobustSigmoid")) {
  
  # Instantiate list for t-test results
  t_test_list <- list()
  
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc) 
    
    # Load catch22 feature matrix
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                         noise_label))) %>%
      mutate(group = factor(group, levels = c("Schz", "Control")))
    
    
    # Define non-normalised data
    non_norm_data <- feature_matrix %>%
      mutate(Norm_Method="non-normalised")
    
    # Calculate t statistics for non-normalised data
    t_stat_res <- non_norm_data %>%
      group_by(Brain_Region, Norm_Method, names) %>%
      nest() %>%
      mutate(
        test = map(data, ~ t.test(.x$values ~ .x$group)), # S3 list-col
        tidied = map(test, tidy)
      ) %>% 
      unnest(tidied) %>%
      dplyr::select(-data, -test) %>%
      mutate(Noise_Proc = noise_proc)
    
    # Append results to list
    t_test_list <- rlist::list.append(t_test_list, t_stat_res)
    
    for (norm_method in norm_methods) {
      
      # Normalise using the given noise-processing method
      normed <- normalise_feature_frame(feature_matrix, 
                                        names_var = "names", 
                                        values_var = "values", 
                                        method = norm_method)
      
      # Calculate t statisics
      t_stat_res <- normed %>%
        group_by(Brain_Region, names) %>%
        nest() %>%
        mutate(
          test = map(data, ~ t.test(.x$values ~ .x$group)), # S3 list-col
          tidied = map(test, tidy)
        ) %>% 
        unnest(tidied) %>%
        dplyr::select(-data, -test) %>%
        mutate(Noise_Proc = noise_proc,
               Norm_Method = norm_method)
      
      t_test_list <- rlist::list.append(t_test_list, t_stat_res)
    }
    
  }
  
  t_test_res <- do.call(plyr::rbind.fill, t_test_list)
  return(t_test_res)
}