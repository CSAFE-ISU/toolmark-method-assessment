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
password <- readLines("buser_pass.txt")
host <- "127.0.0.1"
con <- src_mysql(dbname, host, user=user, password=password)
conn<- dbConnect(MySQL(), user= user, password= password , dbname=dbname, host=host)

result.metadata_derived<- metadata_derived %>% collect()
result.metadata<- metadata %>% collect()

matches <- dbReadTable(conn, "matches")
profiles <- dbReadTable(conn, "profiles")
sigs <- dbReadTable(conn, "signatures")
ccf <- dbReadTable(conn, "ccf")

data<- tbl(conn, "data")
metadata_derived <- tbl(conn, "metadata_derived")
ccf<- tbl(conn, "ccf")
metadata<- tbl(conn, "metadata")
runs<- tbl(conn, "runs")

meta_sub <- metadata %>% filter(study %in% c("Hamby44", "Hamby252")) %>% collect()
profiles_sub <- profiles %>% filter(land_id %in% unique(meta_sub$land_id))
sigs <- sigs %>% left_join(profiles_sub %>% select(land_id, profile_id), by="profile_id")
sigs <- sigs %>% filter(!is.na(land_id))
sigs <- sigs %>% filter(run_id ==3)

sigs_nest <- sigs %>% group_by(land_id, profile_id, run_id) %>% tidyr::nest()

wo <- 110
wv <- 30

k <- 1
n <- nrow(sigs_nest)
dframe <- data.frame(matrix(rep(NA, n*(n-1)*5), ncol=5))
names(dframe) <- c("land1_id", "land2_id", "profile1_id", "profile2_id", "chumbley")
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    dframe$chumbley[k] <- I(list(chumbley_non_random(matrix(sigs_nest$data[[i]]$l30), matrix(sigs_nest$data[[j]]$l30),  window_opt=wo, window_val=wv, coarse=1)))
    dframe$land1_id[k] <- sigs_nest$land_id[i]
    dframe$profile1_id[k] <- sigs_nest$profile_id[i]
    dframe$land2_id[k] <- sigs_nest$land_id[j]
    dframe$profile2_id[k] <- sigs_nest$profile_id[j]
    k <- k + 1
  }
}
dframe$wo <- wo
dframe$wv <- wv

dframe <- dframe %>% mutate(chumbley = chumbley %>% purrr::map(.f=function(x) data.frame(x))) %>% tidyr::unnest(chumbley)
matches$match <- TRUE
dframe <- dframe %>% left_join(matches, by = c("land1_id", "land2_id") ) %>% 
  mutate(match = replace(match, is.na(match), FALSE))

dframe$signif <- dframe$p_value <= .05
write.csv(dframe, file=sprintf("chumbley-out-wo:%d-wv:%d.csv", wo, wv), row.names=FALSE)
dframe %>% ggplot(aes(x = p_value)) + geom_density(aes(fill=match), alpha = 0.5) +
  geom_vline(xintercept=0.05)
dframe %>% ggplot(aes(x = p_value)) + geom_histogram(aes(fill=match), binwidth=0.01) +
  facet_grid(match~.) +
  geom_vline(xintercept=0.05)

xtabs(~signif+match, data =dframe)
