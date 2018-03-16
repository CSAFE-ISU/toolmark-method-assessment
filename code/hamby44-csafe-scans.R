# This document is meant for the Hamby 44 scan made at Iowa State
# Resolution of the scans 0.0635 microns per pixel

# Loading Signature selection analytics
# Remove unkonwns
# Bullets coming from same barrel should match
# Lands that belong to these bullets must match? Do lands from the same bullet match,

library(readr)
library(toolmaRk)
hamby44<- readRDS("/opt/hamby44/hamby44.rda")
featureplus<- read_csv("/opt/hamby44/featuresplus.csv",col_names = T)

library(tidyverse)
library(purrr)
hamb<- hamby44 %>% select(-c(x3p,ccdata))

#Filtered out the unknowns from the hamb
hamb<- hamb %>% filter(barrel != "Unknowns")

# Groove removal
hamb$groove_left<- purrr::map(hamb$grooves,.f = function(x) x[["groove"]][1])
hamb$groove_right<- purrr::map(hamb$grooves,.f = function(x) x[["groove"]][2])

marking_proc<- hamb %>% select(proc_smooth)
marking_proc<- marking_proc[,2]


for(i in 1:120){
  if(!unlist(purrr::map( marking_proc[[1]][i], .f = function(x) is_empty(x)))){ # check for empty proc_smooth #96 i.e Barrel 9 BUllet 2 Land 3 alldata missing
    marking_proc[[1]][i]<- marking_proc[[1]][i] %>% purrr::map(.f = function(x) as.data.frame(x) %>% filter(y >= hamb$groove_left[i] & y <= hamb$groove_right[i]))  
    } else{
      marking_proc[[1]][i]<-  NA
  }
  

}

# generating id of the form barrel-bullet-land in hamb
names(marking_proc[[1]])<- paste0(str_split(hamb$barrel, pattern = " ",simplify = T)[,2],
                "-",
                str_split(hamb$bullet, pattern = " ",simplify = T)[,2],
                "-",
                str_split(hamb$land, pattern = " ",simplify = T)[,2]
                )

## Next steps Involve making land id's for hamb and featureplus matches
#head(paste0(f_sub$bulletland1,"-",f_sub$land1))
f_sub<- featureplus %>% select(barrel2, barrel1, bullet2, bullet1, bulletland2, bulletland1, land2, land1, sameland, samesource, known)

f_sub$barrel_bullet_land1<- paste0(f_sub$bulletland1,"-",f_sub$land1)
f_sub$barrel_bullet_land2<- paste0(f_sub$bulletland2,"-",f_sub$land2)

f_sub <- f_sub %>% filter(known == TRUE)
f_sub<- f_sub[,c("barrel_bullet_land1","barrel_bullet_land2","sameland")]#,"samesource")]
#f_sub %>% filter(samesource == FALSE, sameland == TRUE)

wo <- 120
wv <- 30
coarse <- 0.25
k <- 1
n <- nrow(marking_proc)
dframe <- data.frame(matrix(rep(NA, n*(n-1)*3), ncol=3))
names(dframe) <- c("barrel_bullet_land1", "barrel_bullet_land2", "chumbley")


for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    marking_I<- marking_proc[[1]][i] %>% purrr::map(.f = function(x) x["l30"]) %>% unlist() %>% as.numeric()
    marking_II<- marking_proc[[1]][j] %>% purrr::map(.f = function(x) x["l30"]) %>% unlist() %>% as.numeric()
    if(!(length(marking_I) == 1) && !(length(marking_II) == 1)){ # check for empty of non-existing signature/ proc_smooth fields typically will have an NA (owing to previous setup)
      dframe$chumbley[k] <- I(list(chumbley_non_random(matrix(marking_I), matrix(marking_II),  window_opt=wo, window_val=wv, coarse=coarse)))  
      dframe$barrel_bullet_land1[k] <- names(marking_proc[[1]][i])
      dframe$barrel_bullet_land2[k] <- names(marking_proc[[1]][j])#sigs_nest$land_id[j]
    }


    k <- k + 1
    rm(marking_I, marking_II)
  }
}
dframe$wo <- wo
dframe$wv <- wv
dframe$coarse <- coarse

dframe_raw<- dframe

dframe <- dframe %>% filter(!is.na(dframe$barrel_bullet_land1)) %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
####
#dframe <- dframe %>% head(min(which(is.na(dframe$barrel_bullet_land1)))-1) %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
####
#f_sub$match <- TRUE
dframe <- dframe %>% left_join(f_sub, by = c("barrel_bullet_land1", "barrel_bullet_land2") ) %>% mutate(sameland = replace(sameland, is.na(sameland), FALSE))
dframe %>% ggplot(aes(x = p_value)) + geom_density(aes(fill=sameland), alpha = 0.5) +
  geom_vline(xintercept=0.05)

dframe$signif <- dframe$p_value <= .05
xtabs(~signif+sameland, data =dframe)

#saveRDS(object = dframe, file = "chumbley_wo-280_wv-30")

write.table(dframe, file = "csafe_sig_chumbley_wo-120_wv-30_coarse-pt25.csv", sep = ",", append = FALSE, col.names = TRUE)
aa<- read.csv2("csafe_sig_chumbley_wo-120_wv-30_coarse-pt25.csv",header = TRUE, sep = ",")

# Check if grooves are already removed
# Left groove
already_removed<- numeric(length = 120)
for(i in 1:95){
  t1<- marking_proc[[1]][i] %>% purrr::map(.f = function(x) x["y"]) %>% unlist() %>% as.numeric()
  #if(is.na(t1[1])){ qtmp<- t1[2]} else{qtmp<- t1[1]}
  if((hamb$groove_left[i] %>% as.numeric()) > t1[1]){
    already_removed[i]<- FALSE
  } else{already_removed[i]<- TRUE}
  
}

for(i in 97:120){
  t1<- marking_proc[[1]][i] %>% purrr::map(.f = function(x) x["y"]) %>% unlist() %>% as.numeric()
  #if(is.na(t1[1])){ qtmp<- t1[2]} else{qtmp<- t1[1]}
  if((hamb$groove_left[i] %>% as.numeric()) > t1[1]){
    already_removed[i]<- FALSE
  } else{already_removed[i]<- TRUE}
  
}

# Right Groove

already_removed<- numeric(length = 120)
for(i in 1:95){
  t1<- marking_proc[[1]][i] %>% purrr::map(.f = function(x) x["y"]) %>% unlist() %>% as.numeric()
  #if(is.na(t1[1])){ qtmp<- t1[2]} else{qtmp<- t1[1]}
  if((hamb$groove_right[i] %>% as.numeric()) < t1[1]){
    already_removed[i]<- FALSE
  } else{already_removed[i]<- TRUE}
  
}

# Number 96 proc_smooth is empty missing data
for(i in 97:120){
  t1<- marking_proc[[1]][i] %>% purrr::map(.f = function(x) x["y"]) %>% unlist() %>% as.numeric()
  #if(is.na(t1[1])){ qtmp<- t1[2]} else{qtmp<- t1[1]}
  if((hamb$groove_right[i] %>% as.numeric()) < t1[1]){
    already_removed[i]<- FALSE
  } else{already_removed[i]<- TRUE}
  
}

# Chop? No data for y values with the grooves in proc_smooth
