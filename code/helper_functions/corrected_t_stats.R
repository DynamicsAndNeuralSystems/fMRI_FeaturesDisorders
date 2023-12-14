# Taken from https://github.com/hendersontrent/feature-set-classification/blob/main/R/corr_t_test.R

#------------------------------------------
# This script sets defines a function for 
# calculating the correlated t-test statistic
# presented in Nadeau & Bengio (2003)
#------------------------------------------

#------------------------------------------
# Author: Trent Henderson, 21 November 2022
#------------------------------------------

#' Function to compute Nadeau & Bengio (2003) correlated t-statistic p-value for train-test splits
#' @param x vector of classification accuracy values for classifier A
#' @param y vector of classification accuracy values for classifier B
#' @param n integer denoting number of repeat samples
#' @param n1 integer denoting train set size
#' @param n2 integer denoting test set size
#' @author Trent Henderson
#' 

corr_t_test <- function(x, y, n, n1, n2){
  
  d <- x - y # Calculate differences
  d_bar <- mean(d, na.rm = TRUE) # Calculate mean of differences
  sigma_2 <- var(d, na.rm = TRUE) # Calculate variance
  sigma_2_mod <- sigma_2 * (1/n + n2/n1) # Calculate modified variance
  t_stat <- d_bar / sqrt(sigma_2_mod) # Calculate t-statistic
  
  if(t_stat < 0){
    p_val <- pt(t_stat, n - 1) # p-value for left tail
  } else{
    p_val <- pt(t_stat, n - 1, lower.tail = FALSE) # p-value for right tail
  }
  
  return(p_val)
}


#' Function to compute Nadeau & Bengio (2003) correlated t-statistic p-value for k-fold CV
#' @param x vector of classification accuracy values for classifier A
#' @param y vector of classification accuracy values for classifier B
#' @param n integer denoting total sample size
#' @param k number of folds
#' @author Trent Henderson
#' 

corr_t_test_folds <- function(x, y, n, k){
  
  d <- x - y # Calculate differences
  d_bar <- mean(d, na.rm = TRUE) # Calculate mean of differences
  sigma_2 <- var(d, na.rm = TRUE) # Calculate variance
  t_stat <- d_bar / sqrt(sigma_2 * ((1/n + (1/k)) / (1 - 1/k))) # Calculate t-statistic
  
  if(t_stat < 0){
    p_val <- pt(t_stat, n - 1) # p-value for left tail
  } else{
    p_val <- pt(t_stat, n - 1, lower.tail = FALSE) # p-value for right tail
  }
  
  return(p_val)
}


corr_t_test_repeated_cv <- function(x, y, n, k, r, n1, n2) {
  
}
