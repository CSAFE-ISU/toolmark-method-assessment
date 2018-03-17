# This document is meant for the Hamby 44 scan made at CSAFE
# Resolution of the scans 0.0645 microns per pixel

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

# Chumbley comparison
wo <- 300
wv <- 120
coarse <- 0.25
k <- 1
n <- nrow(marking_proc)
dframe <- data.frame(matrix(rep(NA, n*(n-1)*3), ncol=3))
names(dframe) <- c("barrel_bullet_land1", "barrel_bullet_land2", "chumbley")


for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    # Extracting the two signatures from the nested data
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

#dframe_raw<- dframe

#dframe<- dframe_raw

dframe <- dframe %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
#dframe <- dframe %>% filter(!is.na(dframe$barrel_bullet_land1)) %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
####
#dframe <- dframe %>% head(min(which(is.na(dframe$barrel_bullet_land1)))-1) %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
####

dframe <- dframe %>% left_join(f_sub, by = c("barrel_bullet_land1", "barrel_bullet_land2") ) %>% mutate(sameland = replace(sameland, is.na(sameland), FALSE))
dframe %>% ggplot(aes(x = p_value)) + geom_density(aes(fill=sameland), alpha = 0.5) +
  geom_vline(xintercept=0.05)

dframe$signif <- dframe$p_value <= .05
xtabs(~signif+sameland, data =dframe)

#saveRDS(object = dframe, file = "chumbley_wo-280_wv-30")

write.table(dframe, file = "csafe_sig_chumbley_wo-300_wv-120_coarse-pt25.csv", sep = ",", append = FALSE, col.names = TRUE)
aa<- read.csv2("csafe_sig_chumbley_wo-300_wv-120_coarse-pt25.csv",header = TRUE, sep = ",")

####################################################################################################################
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


###########################################################################################################
# Hamby 44 for NIST scans

# Choosing one file for comparison
nist_aa120_30<- read_csv("./data/signatures/chumbley-out-wo-120-wv-30.csv")#,header = TRUE, sep = ",")
nist_aa120_30<- as_data_frame(nist_aa120_30)
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
ccf<- tbl(conn, "ccf")
metadata<- tbl(conn, "metadata")
runs<- tbl(conn, "runs")

# result.metadata_derived<- metadata_derived %>% collect()
# result.metadata<- metadata %>% collect()

# Hamby 44 identification
meta_sub <- metadata %>% filter(study %in% c("Hamby44")) %>% collect() #%>% as.data.frame()
temp_ms<- meta_sub
#head(meta_sub,10)
levels(factor(meta_sub$barrel))

# Remove Unkowns
meta_sub<- meta_sub[grep("[[:digit:]]", meta_sub$barrel), ]

# Generating ID similar to CSAFE scans(helps in comparison)
meta_sub$bbl<- paste0(meta_sub$barrel,"-",meta_sub$bullet,"-",meta_sub$land)
profiles_sub <- profiles %>% filter(land_id %in% unique(meta_sub$land_id))

#head(profiles_sub)
#head(profiles)
#sigs <- sigs %>% left_join(profiles_sub %>% select(land_id, profile_id), by="profile_id")
#colnames(meta_sub)

# Psuedo land1_id and land2_id for the purpose of joining
meta_sub$land1_id<- meta_sub$land_id
meta_sub$land2_id<- meta_sub$land_id

# Psuedo barrel bullet land id for purpose of uniformity and comparison
meta_sub$bbl1<- meta_sub$bbl
meta_sub$bbl2<- meta_sub$bbl

# No need for this but for clarity in the joining process
ms1<- meta_sub %>% select(land1_id,bbl1)
ms1$land1_id<- as.integer(ms1$land1_id)

# 1st join
nist_sub <- nist_aa120_30 %>% right_join(ms1, by= c("land1_id"))

ms2<-meta_sub %>% select(land2_id,bbl2)
ms2$land2_id<- as.integer(ms2$land2_id)

# 2nd join
nist_sub <- nist_sub %>% right_join(ms2, by= c("land2_id"))

# 2 entries were NA no values in them in the comparison, removing them 
nist_sub<- nist_sub[which(!is.na(nist_sub$bbl1)),]



#anti_join(data.frame(a = nist_sub$bbl1),data.frame(a = dframe$barrel_bullet_land1), by = "a")
#anti_join(data.frame(a = dframe$barrel_bullet_land1), data.frame(a = nist_sub$bbl1), by = "a")

# Checking for what specific lands are different between NIST and CSAFE
# For a more accurate comparison of error rates
tmp1_setdiff<- setdiff(dframe$barrel_bullet_land1, nist_sub$bbl1)
tmp1_rev_setdiff<- setdiff(nist_sub$bbl1, dframe$barrel_bullet_land1)

#setdiff(dframe$barrel_bullet_land2, nist_sub$bbl2)
#setdiff(nist_sub$bbl2, dframe$barrel_bullet_land2)

# Filtering out lands that are not in the iowa state scan
nist_sub<- nist_sub %>% filter(bbl1 != tmp1_setdiff)
#dframe<- dframe %>% filter(barrel_bullet_land1 == tmp1_rev_setdiff)

xtabs(~signif+match, data =nist_sub)
# The False Positive and False Negative rates for Wo 120 and Wv 30 NIST scans
# > 384/(6514+384)
# [1] 0.05566831
# 
# > 18/(18+40)
# [1] 0.3103448
#
# # wo 300 wv 120 csafe scans
# fp =  0.06636784
# fn =  0.3846154

######################################################################################
# 387/(6540+387) # =0.05586834
# 18/(18+40) #= 0.3103448
#nist_sub[which(!is.na(nist_sub$bbl1))]
#nist_sub<- join(nist_aa120_30, ms1, by= c("land1_id","land2_id"))
######################################################################################
# Chop? No data for y values with the grooves in proc_smooth

# xtabs(~signif+sameland, data =aa) 
# 
# 19/(19+21)
# 
# # Read in 200 and 50 and 0.25
# 
# aa200_50<- read.csv2("./data/csafe-scans/signatures/csafe_sig_chumbley_wo-200_wv-50_coarse-pt25.csv",header = TRUE, sep = ",")
# xtabs(~signif+sameland, data =aa200_50) 
#  
# 408/(408+6562) = 0.05853659
# 18/(18+23)= 0.4390244
# 
# aa1750_50<- read.csv2("./data/csafe-scans/signatures/csafe_sig_chumbley_wo-175_wv-50_coarse-pt25.csv",header = TRUE, sep = ",")
# xtabs(~signif+sameland, data =aa1750_50) 
# 
# 422/(422+6553) = 0.06050179
# 18/(18+23)= 0.4390244
# 
# aa250_50<- read.csv2("./data/csafe-scans/signatures/csafe_sig_chumbley_wo-250_wv-50_coarse-pt25.csv",header = TRUE, sep = ",")
# xtabs(~signif+sameland, data =aa250_50) 
# 
# 412/(412+6555) = 0.05913593
# 16/(16+ 24) = 0.4
# 
# 
# nist_aa120_30<- read.csv2("./data/signatures/chumbley-out-wo-120-wv-30.csv",header = TRUE, sep = ",")
# xtabs(~signif+match, data =nist_aa120_30)
# 4523/(79249 + 4523) = 0.05399179
# 386/(386 + 854) = 0.3112903
#
# # wo 300 wv 120 csafe scans
# fp =  0.06636784
# fn =  0.3846154
#
# random
# 15/(6457 + 15)
# 24/(459+24)
# 386/(79249 + 386)
# 4523/(4523 + 854)
#   ############################################