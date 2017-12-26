#' Modified Chumbley Non-Random
#' 
#' The modification is made in the same-shift, where for the 2nd signature 
#' the window is moved a bit towards the left and right to find if the 
#' correlation of the moved windows are higher than the original same shift (OSS) window
#' If its higher the new correlation values are stored
#' 
#' This function computes the Chumbley U-Statistic on systemically chosen pairs of windows rather 
#' than the original method which selects randomly chosen pairs of windows
#' @param data1 The first tool mark as a 1-column matrix
#' @param data2 The second tool mark as a 1-column matrix
#' @param window_opt size of the window to be used in the optimization step
#' @param window_val Size of the window to be used in the validation step
#' @param coarse smoothing parameter for the normalization smooth
#' @importFrom stats pnorm
#' @importFrom ggplot2 ggplot
#' @export
#' @return list with
#' \itemize{
#' \item {same_shift_n} Number of same shift offsets used
#' \item {diff_shift_n} Number of different shift offsets used
#' \item {U} observed U statistic
#' \item {p_value} Corresponding p-value
#' }
chumbley_non_random_modified <- function(data1, data2, window_opt = 500, window_val = 50, coarse = .25){
  
  unity <- function(x) {x / sqrt(sum(x^2))} ## normalize columns of a matrix to make correlation computation faster
  
  ####################################################
  ##Clean the marks and compute the smooth residuals##
  ####################################################
  
  data1 <- matrix(data1[round((0.01*nrow(data1))):round(0.99*nrow(data1)),], ncol = 1)
  data2 <- matrix(data2[round((0.01*nrow(data2))):round(0.99*nrow(data2)),], ncol = 1)
  
  ##Normalize the tool marks
  y1 <- data1 - lowess(y = data1,  x = 1:nrow(data1), f= coarse)$y
  y2 <- data2 - lowess(y = data2,  x = 1:nrow(data2), f= coarse)$y
  
  
  ############################################
  ##Compute the observed maximum correlation##
  ############################################
  
  #####################
  ##Optimization step##
  #####################
  ##Each column in these matrices corresponds to a window in the respective tool mark
  y1_mat_opt <- matrix(NA, ncol = length(1:(length(y1) - (window_opt - 1))), nrow = window_opt)
  for(l in 1:(length(y1) - (window_opt - 1))){
    y1_mat_opt[,l] <- y1[l:(l+(window_opt - 1))]
  }
  y2_mat_opt <- matrix(NA, ncol = length(1:(length(y2) - (window_opt - 1))), nrow = window_opt)
  for(l in 1:(length(y2) - (window_opt - 1))){
    y2_mat_opt[,l] <- y2[l:(l+(window_opt - 1))]
  }
  
  ##Compute the correlation between all pairs of windows for the two marks
  ##Rows in the following matrix are mark 2, columns are mark 1
  y2_mat_opt <- apply(scale(y2_mat_opt), 2, unity)
  y1_mat_opt <- apply(scale(y1_mat_opt), 2, unity)
  corr_mat_opt <- t(y2_mat_opt) %*% y1_mat_opt ##correlation matrix
  max_corr_opt_loc <- which(corr_mat_opt == max(corr_mat_opt), arr.ind = TRUE) ##pair of windows maximizing the correlation
  
  
  ###################
  ##Validation step##
  ###################
  ##Each column in these matrices corresponds to a window in the respective tool mark
  y1_mat_val <- matrix(NA, ncol = length(1:(length(y1) - (window_val - 1))), nrow = window_val)
  for(l in 1:(length(y1) - (window_val - 1))){
    y1_mat_val[,l] <- y1[l:(l+(window_val - 1))]
  }
  y2_mat_val <- matrix(NA, ncol = length(1:(length(y2) - (window_val - 1))), nrow = window_val)
  for(l in 1:(length(y2) - (window_val - 1))){
    y2_mat_val[,l] <- y2[l:(l+(window_val - 1))]
  }
  
  ##Compute the correlation between all pairs of windows for the two marks
  ##Rows in the following matrix are mark 2, columns are mark 1
  y2_mat_val <- apply(scale(y2_mat_val), 2, unity)
  y1_mat_val <- apply(scale(y1_mat_val), 2, unity)
  corr_mat_val <- t(y2_mat_val) %*% y1_mat_val
 
###### SECTION MODIFIED ######################################   
  ##Pull out the correlations that correspond to windows with the same offset as the largest correlation found in the optimization step
  same_shift <- data.frame(row = NA, col = NA, U = NA)
  
  # Note step size and window size in consideration in the step below is windoe_val
  # Because of the initialization of the vector form subtracts window_val
  # In the first step of the following while loop we fall into place
  # i.e the end of the window of optimization becomes the place
  # from which the stepping size of window_val is used to get subsequent 
  # windows of same shift
  rows <- max_corr_opt_loc[1] + (window_opt - window_val) # Initialization
  cols <- max_corr_opt_loc[2] + (window_opt - window_val) # Initialization
  eps<- 1 #round(window_val * 0.05)
  cols_right2<- cols + 2*eps
  # Moving towards the right end of the signature
  while(rows + window_val < nrow(corr_mat_val) & cols + window_val < ncol(corr_mat_val)){
    
    # Original Same shift
    rows <- rows + window_val # step size is window of validation
    cols <- cols + window_val
    
    # Need a window to left of OSS and Need a window to right of OSS
    # The left movement and / or right movement from OSS of the window needs to be
    # Window_val > Left/ right movement
    # Left/ right movement -+ OSS cols = Left/ right window location
    # Pull out the location from the corr_mat_val :: this is the correlation matrix for 
    # validation window taken as window size and step size as
    
    if (cols_right2 + window_val < ncol(corr_mat_val))
    {
      cols_right1<- cols + eps
      cols_right2<- cols + 2*eps
      cols_left1 <- cols - eps
      cols_left2 <- cols - 2*eps
      
    }else{
      cols_right1<- cols
      cols_right2<- cols
      cols_left1 <- cols - eps
      cols_left2 <- cols - 2*eps
    }
    
    max_wiggle_corr<- 
      max(c(corr_mat_val[rows,cols_left2],corr_mat_val[rows,cols_left1],
            corr_mat_val[rows,cols],
            corr_mat_val[rows,cols_right1],corr_mat_val[rows,cols_right2]))
    
    max_wiggle_corr_loc <- which(corr_mat_val == max_wiggle_corr, arr.ind = TRUE)
    
    
    cols_star<- max_wiggle_corr_loc[2]
    
    # Same shift error
    same_shift <- rbind(same_shift, c(rows, cols_star, corr_mat_val[rows,cols_star]))
    
  } 
  
  rows <- max_corr_opt_loc[1]
  cols <- max_corr_opt_loc[2]
  
  cols_left2 <- cols - 2*eps
  # Moving towards the left end of the signature
  while(rows - window_val > 0 & cols - window_val > 0){
    
    rows <- rows - window_val
    cols <- cols - window_val

    if (cols_left2 - window_val > 0)
    {
      cols_right1<- cols + eps
      cols_right2<- cols + 2*eps
      cols_left1 <- cols - eps
      cols_left2 <- cols - 2*eps
      
    }else{
      cols_right1<- cols + eps
      cols_right2<- cols + 2*eps
      cols_left1 <- cols 
      cols_left2 <- cols 
    }
    
    max_wiggle_corr<- 
      max(c(corr_mat_val[rows,cols_left2],corr_mat_val[rows,cols_left1],
            corr_mat_val[rows,cols],
            corr_mat_val[rows,cols_right1],corr_mat_val[rows,cols_right2]))
    
    max_wiggle_corr_loc <- which(corr_mat_val == max_wiggle_corr, arr.ind = TRUE)
    
    
    cols_star<- max_wiggle_corr_loc[2]
    
    # Same shift error
    same_shift <- rbind(same_shift, c(rows, cols_star, corr_mat_val[rows,cols_star]))
    
  }
  
  ###########SECTION MODIFICATION FINISHED##########################
  
  same_shift <- same_shift[-1,]
  
  ##Pull out the correlations that correspond to windows with different offset as the largest correlation found in the optimization step
  ##along a single anti-diagonal
  diff_shift <- data.frame(row = NA, col = NA, U = NA)
  rows <- max_corr_opt_loc[1] + (window_opt - window_val)
  cols <- max_corr_opt_loc[2]
  while(rows + window_val < nrow(corr_mat_val) & cols - window_val > 0){
    
    rows <- rows + window_val
    cols <- cols - window_val
    diff_shift <- rbind(diff_shift, c(rows, cols, corr_mat_val[rows,cols]))
    
  }
  rows <- max_corr_opt_loc[1]
  cols <- max_corr_opt_loc[2] + (window_opt - window_val)
  while(rows - window_val > 0 & cols + window_val < ncol(corr_mat_val)){
    
    rows <- rows - window_val
    cols <- cols + window_val
    diff_shift <- rbind(diff_shift, c(rows, cols, corr_mat_val[rows, cols]))
    
  }
  diff_shift <- diff_shift[-1,]
  
  ######################################
  ##Compute the Ustatistic if possible##
  ######################################
  if(nrow(same_shift) == 0 | nrow(diff_shift) == 0) {
    
    obs_U <- NA
    n <- length(same_shift$U)
    m <- length(diff_shift$U)
    
  }
  
  if(nrow(same_shift) != 0 & nrow(diff_shift) != 0) {
    
    ranks <- rank(c(same_shift$U, diff_shift$U))
    Rx <- ranks[seq_along(same_shift$U)]
    Ry <- ranks[-(1:length(Rx))]
    n <- length(same_shift$U)
    m <- length(diff_shift$U)
    N <- n + m
    
    t <- sum(Rx) ##Test statistic...sum of sample one ranks
    t1 <- (t - n*((N + 1) / 2)) / sqrt( ((n*m)/(N*(N-1)))*sum(c(Rx^2, Ry^2)) - ((n*m*(N+1)^2) / (4*(N-1)))) ##Standardized test statistics
    obs_U <- t1
    
  }
  pval <- 1 - pnorm(obs_U)
  
  list(same_shift_n = n, ##Number of same shift offsets used
       diff_shift_n = m, ##Number of different shift offsets used
       U = obs_U, ##observed U-statistic
       p_value = pval) ##Corresponding p-value
}
