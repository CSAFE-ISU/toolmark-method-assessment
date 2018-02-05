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
{
  chumbley_non_random_modified <- function(data1, data2, window_opt = 500, window_val = 50, coarse = .25, eps = 1){
  
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
  #eps<- 1 #round(window_val * 0.05)
  
  
  #cols_right2<- cols + 2*eps
  cols_right2<- cols + (eps+1)
  
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
      #cols_right2<- cols + 2*eps
      cols_right2<- cols + (eps+1)
      cols_left1 <- cols - eps
      cols_left2 <- cols - (eps+1)
      
    }else{
      cols_right1<- cols
      cols_right2<- cols
      cols_left1 <- cols - eps
      cols_left2 <- cols - (eps+1)
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
  
  # cols_left2 <- cols - 2*eps
  cols_left2 <- cols - (eps+1)
  
  # Moving towards the left end of the signature
  while(rows - window_val > 0 & cols - window_val > 0){
    
    rows <- rows - window_val
    cols <- cols - window_val

    if (cols_left2 - window_val > 0)
    {
      cols_right1<- cols + eps
      #cols_right2<- cols + 2*eps
      cols_right2<- cols + (eps+1)
      cols_left1 <- cols - eps
      cols_left2 <- cols - (eps+1)
      
    }else{
      cols_right1<- cols + eps
      #cols_right2<- cols + 2*eps
      cols_right2<- cols + (eps+1)
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
  
  cols_left2_diffshift <- cols - (eps+1)
  while(rows + window_val < nrow(corr_mat_val) & cols - window_val > 0){
    
    rows <- rows + window_val
    cols <- cols - window_val
    
   # corresponding to the same shift cols - window_val 
    if (cols_left2_diffshift - window_val > 0)
    {
      cols_right1_diffshift<- cols + eps
      #cols_right2<- cols + 2*eps
      cols_right2_diffshift<- cols + (eps+1)
      cols_left1_diffshift <- cols - eps
      cols_left2_diffshift <- cols - (eps+1)
      
    }else{
      cols_right1_diffshift<- cols + eps
      #cols_right2<- cols + 2*eps
      cols_right2_diffshift<- cols + (eps+1)
      cols_left1_diffshift <- cols 
      cols_left2_diffshift <- cols 
    }
    
    max_wiggle_corr_diffshift<- 
      max(c(corr_mat_val[rows,cols_left2_diffshift],corr_mat_val[rows,cols_left1_diffshift],
            corr_mat_val[rows,cols],
            corr_mat_val[rows,cols_right1_diffshift],corr_mat_val[rows,cols_right2_diffshift]))
    
    max_wiggle_corr_loc_diffshift <- which(corr_mat_val == max_wiggle_corr_diffshift, arr.ind = TRUE)
    
    
    cols_star_diffshift<- max_wiggle_corr_loc_diffshift[2]
    
    # diff shift error
    diff_shift <- rbind(diff_shift, c(rows, cols_star_diffshift, corr_mat_val[rows,cols_star_diffshift]))
    #diff_shift <- rbind(diff_shift, c(rows, cols, corr_mat_val[rows,cols]))
    
    
    
  }
  rows <- max_corr_opt_loc[1]
  cols <- max_corr_opt_loc[2] + (window_opt - window_val)
  
  cols_right2_diffshift<- cols + (eps+1)
  while(rows - window_val > 0 & cols + window_val < ncol(corr_mat_val)){
    
    rows <- rows - window_val
    cols <- cols + window_val
    if (cols_right2_diffshift + window_val < ncol(corr_mat_val))
    {
      cols_right1_diffshift<- cols + eps
      #cols_right2<- cols + 2*eps
      cols_right2_diffshift<- cols + (eps+1)
      cols_left1_diffshift <- cols - eps
      cols_left2_diffshift <- cols - (eps+1)
      
    }else{
      cols_right1_diffshift<- cols
      cols_right2_diffshift<- cols
      cols_left1_diffshift <- cols - eps
      cols_left2_diffshift <- cols - (eps+1)
    }
    
    max_wiggle_corr_diffshift<- 
      max(c(corr_mat_val[rows,cols_left2_diffshift],corr_mat_val[rows,cols_left1_diffshift],
            corr_mat_val[rows,cols],
            corr_mat_val[rows,cols_right1_diffshift],corr_mat_val[rows,cols_right2_diffshift]))
    
    max_wiggle_corr_loc_diffshift <- which(corr_mat_val == max_wiggle_corr_diffshift, arr.ind = TRUE)
    
    
    cols_star_diffshift<- max_wiggle_corr_loc_diffshift[2]
    
    
    diff_shift <- rbind(diff_shift, c(rows, cols_star_diffshift, corr_mat_val[rows,cols_star_diffshift]))
    #diff_shift <- rbind(diff_shift, c(rows, cols, corr_mat_val[rows, cols]))
    
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

  
##############################################
  ##########CHECK
  # library(ggplot2)
  # 
  # sigs_graphics<- readRDS("./data/sigs_generate_window.rds")
  # 
  # sigs_graphics.match<- as.data.frame(sigs_graphics[1])
  # sigs_graphics.nonmatch<- as.data.frame(sigs_graphics[2])
  # colnames(sigs_graphics.match)<- gsub("^.*\\.","", colnames(sigs_graphics.match))
  # colnames(sigs_graphics.nonmatch)<- gsub("^.*\\.","", colnames(sigs_graphics.nonmatch))
  # 
  # # Matching signature
  # #gi<-
  # #sigs_graphics.match<- sigs_graphics.match %>% filter(land_id == c(1,9)) # 1 9 17 19 276
  # #sigs_graphics.match<- sigs_graphics.match %>% select(y, l30, land_id)
  # d1<- sigs_graphics.match[which(sigs_graphics.match$land_id==19,arr.ind = TRUE),]
  # d2<- sigs_graphics.match[which(sigs_graphics.match$land_id==276,arr.ind = TRUE),]
  # d3n<- sigs_graphics.nonmatch[which(sigs_graphics.nonmatch$land_id==3,arr.ind = TRUE),]
  # data1<-d1$l30
  # #data1<- data1$y
  # data2<- d2$l30
  # #data2<- data2$y
  # data3n<- d3n$l30
  # window_opt = 120
  # window_val = 30
  # coarse = 1
  # eps = 1
  # 
  # data1<- matrix(unlist(data1))
  # data2<- matrix(unlist(data2))
  # data3n<- matrix(unlist(data3n))
  # 
  # chumbley_non_random_modified(data1, data3n,  window_opt=window_opt, window_val=window_val, coarse=coarse, eps = eps)
  # 
#################################################
library(toolmaRk)
library(RMySQL)
library(dplyr)
library(dbplyr)
#library(dbConnect)
library(ggplot2)
library(purrr)
###
### DB Connection
###

dbname <- "bullets"
user <- "buser"
password <- readLines("../buser_pass.txt")
host <- "127.0.0.1"
con <- src_mysql(dbname, host, user=user, password=password)
conn<- dbConnect(MySQL(), user= user, password= password , dbname=dbname, host=host)

matches <- dbReadTable(conn, "matches")
profiles <- dbReadTable(conn, "profiles")
sigs <- dbReadTable(conn, "signatures")
#ccf <- dbReadTable(conn, "ccf")

data<- tbl(conn, "data")
metadata_derived <- tbl(conn, "metadata_derived")
#ccf<- tbl(conn, "ccf")
metadata<- tbl(conn, "metadata")
runs<- tbl(conn, "runs")

result.metadata_derived<- metadata_derived %>% collect()
result.metadata<- metadata %>% collect()

meta_sub <- metadata %>% filter(study %in% c("Hamby44", "Hamby252")) %>% collect()
profiles_sub <- profiles %>% filter(land_id %in% unique(meta_sub$land_id))
sigs <- sigs %>% left_join(profiles_sub %>% select(land_id, profile_id), by="profile_id")
sigs <- sigs %>% filter(!is.na(land_id))
sigs <- sigs %>% filter(run_id ==3)

sigs_nest <- sigs %>% group_by(land_id, profile_id, run_id) %>% tidyr::nest()

# sigs_land_id_1_match<- sigs %>% filter(land_id %in% c(1,9,17,19,276))
# sigs_land_id_1_nonmatch<- sigs %>% filter(land_id %in% c(1,2,3,4))
# sig_generate<- list(sigs_land_id_1_match, sigs_land_id_1_nonmatch)
# names(sig_generate)<- c("sigs_land_id_1_match", "sigs_land_id_1_nonmatch")
# saveRDS(sig_generate, file="./data/sigs_generate_window.rds")

wo <- 120
wv <- 30
coarse <- 1
eps<- 4

k <- 1
n <- nrow(sigs_nest)
dframe <- data.frame(matrix(rep(NA, n*(n-1)*5), ncol=5))
names(dframe) <- c("land1_id", "land2_id", "profile1_id", "profile2_id", "chumbley")


for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    dframe$chumbley[k] <- I(list(chumbley_non_random_modified(matrix(sigs_nest$data[[i]]$l30), matrix(sigs_nest$data[[j]]$l30),  window_opt=wo, window_val=wv, coarse=coarse, eps = eps)))
    dframe$land1_id[k] <- sigs_nest$land_id[i]
    dframe$profile1_id[k] <- sigs_nest$profile_id[i]
    dframe$land2_id[k] <- sigs_nest$land_id[j]
    dframe$profile2_id[k] <- sigs_nest$profile_id[j]
    k <- k + 1
  }
}
dframe$wo <- wo
dframe$wv <- wv

dframe <- dframe %>% head(min(which(is.na(dframe$land1_id)))-1) %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
matches$match <- TRUE
dframe <- dframe %>% left_join(matches, by = c("land1_id", "land2_id") ) %>% 
  mutate(match = replace(match, is.na(match), FALSE))
dframe %>% ggplot(aes(x = p_value)) + geom_density(aes(fill=match), alpha = 0.5) +
  geom_vline(xintercept=0.05)

dframe$signif <- dframe$p_value <= .05
xtabs(~signif+match, data =dframe)

#saveRDS(object = dframe, file = "chumbley_wo-280_wv-30")

write.table(dframe, file = "algo-mod_v2_sig_chumbley_wo-120_wv-30_coarse1_eps4.csv", sep = ",", append = FALSE, col.names = TRUE)

}
aa<- read.csv2("./data/Algo_modification/Signatures/algo-mod2_v2_sig_chumbley_wo-120_wv-30_coarse1_eps4.csv",header = TRUE, sep = ",")

#  dframe<- read.csv2("./toolmark-method-assessment/data/chumbley_wo-320_wv-50.csv",header = TRUE, sep = ",")
# dframe$U<- as.numeric(as.character(dframe$U))
# dframe$p_value<- as.numeric(as.character(dframe$p_value))
