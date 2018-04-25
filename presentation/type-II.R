# ```{r type2, fig.height = 3, fig.width = 4.5, out.width='0.6\\textwidth', fig.cap="Type II error rates observed across a range of window sizes for optimization $w_o$. For a window size of $w_o = 130$ we see a minimum in type II error rate across all type I rates considered. Smaller validation sizes $w_v$ are typically associated with a smaller type II error."}
errorrate <- function(data, alpha) {
  summ <- data %>% filter(!is.na(p_value)) %>%
    mutate(signif = p_value < alpha) %>%
    group_by(wv, wo, match, coarse, signif) %>% tally()
  summ$error <- with(summ, match != signif)
  summ$alpha <- alpha

  rates <- summ %>% dplyr::group_by(wv, wo, coarse, match, alpha) %>% dplyr::summarize(
    rate = n[error==TRUE]/sum(n)
  )
  totals <- summ %>% dplyr::group_by(wv, wo, coarse, alpha) %>% dplyr::summarize(
    total= sum(n[error==TRUE])/sum(n)
  )
  rates <- rates %>% spread(match,rate) %>% dplyr::rename(
    actual = `FALSE`,
    beta = `TRUE`
  )
  left_join(rates, totals)
}


if (!file.exists("../data/cvcs2.RDS")) {
  files <- dir("../data/crossvalidate-cs2/")
  cvcs2 <- data.frame()
  for (file in files) {
    tmp <- readRDS(file.path("../data/crossvalidate-cs2/",file))
    cvcs2 <- rbind(cvcs2, tmp)
  }
  cvcs2$coarse <- cvcs2$coarseness
  cvcs2$method <- "CS2"
  saveRDS(cvcs2, file = "../data/cvcs2.RDS")
} else cvcs2 <- readRDS("../data/cvcs2.RDS")

if (!file.exists("../data/cvcs1.RDS")) {
  files <- dir("../data/crossvalidate-cs1/")
  cvcs1 <- data.frame()
  for (file in files) {
    tmp <- readRDS(file.path("../data/crossvalidate-cs1/",file))
    cvcs1 <- rbind(cvcs1, tmp)
  }
  cvcs1$coarse <- cvcs1$coarseness
  cvcs1$method <- "CS1"
  saveRDS(cvcs1, file = "../data/cvcs1.RDS")
} else cvcs1 <- readRDS("../data/cvcs1.RDS")

cvcs <- rbind(cvcs1, cvcs2)

library(plotROC)
errorsCV2 <- rbind(errorrate(cvcs2, 0.001), 
                   errorrate(cvcs2, 0.005), 
                   errorrate(cvcs2, 0.01), 
                   errorrate(cvcs2, 0.05))
errorsCV2$method <- "CS2"

errorsCV1 <- rbind(errorrate(cvcs1, 0.001), 
                   errorrate(cvcs1, 0.005), 
                   errorrate(cvcs1, 0.01), 
                   errorrate(cvcs1, 0.05))
errorsCV1$method <- "CS1"

errorsCV <- rbind(errorsCV1, errorsCV2)

type2_c_15perc<- errorsCV %>% 
  ggplot(aes(x = wo, y = beta, colour=factor(alpha))) + 
  geom_point(aes(shape=factor(wv)), size=2.5) +
  geom_smooth(aes(group=alpha), se=FALSE, method="loess", span=1,  size=.7) +
  theme_bw() +
  scale_colour_brewer(expression("Nominal\ntype I error"~alpha),
                      palette="Set2") +
  scale_shape_discrete(expression("Size of\nvalidation window "~w[v])) +
  xlab(expression("Window size for optimization "~w[o])) +
  ylab("Type II error rate") +
  facet_wrap(~method, labeller="label_both")

###############################################################
################### CS1#########################################

library(tidyverse)
if (!file.exists("../data/hamby-profiles.rds")) {
  files <- dir("../data/profiles_crossvalid_windows/", pattern="rds")
  files <- grep("^chumbley-csafe-profiles", files, value = TRUE)
  
  profiles_CS1 <- data.frame()
  for (file in files) {
    tmp <- readRDS(file=file.path(path="./data/", file)) 
    #   if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
    profiles_CS1 <- rbind(profiles_CS1, tmp)
  }
  saveRDS(profiles_CS1, file="../data/hamby-profiles.rds")
} else {
  profiles_CS1 <- readRDS("../data/hamby-profiles.rds")
}

profiles_CS1$coarse <- 0.3




errors_CS1 <- rbind(errorrate(profiles_CS1, 0.001), 
                    errorrate(profiles_CS1, 0.005), 
                    errorrate(profiles_CS1, 0.01), 
                    errorrate(profiles_CS1, 0.05))
errors_CS1$method <- "CS1"

type2_cs1_c_30perc<- errors_CS1 %>% ggplot(aes(x = wo, y = beta, colour=factor(alpha))) + geom_point(aes(shape=factor(wv)), size=2.5) + geom_smooth(aes(group=alpha), se=FALSE, method="loess", span=1,  size=.7) +  theme_bw() +  scale_colour_brewer(expression("Nominal\ntype I error"~alpha), palette="Set2") + scale_shape_discrete(expression("Size of\nvalidation window "~w[v])) + xlab(expression("Window size for optimization "~w[o])) + ylab("Type II error rate") + facet_wrap(~method, labeller="label_both")
# geom_point() + geom_smooth(aes(group=alpha), se=FALSE)


testfailed <- profiles_CS1 %>% group_by(wo, wv, match) %>% summarize(
  total = n(), 
  failed = sum(is.na(p_value)))
# testfailed %>% ggplot(aes(x = wo, y = failed/total, colour=factor(match), shape=factor(wv))) + geom_point() +
#   ylim(0, NA)

type1_cs1_30perc<- errors_CS1 %>% filter(wv %in% c(30,50)) %>%
  ggplot(aes(x = wo, y = actual, colour=factor(alpha))) + 
  geom_hline(aes(yintercept=alpha), colour="grey30", size=.5) +
  geom_point(aes(shape=factor(wv)), size=2.5) +
  #  facet_grid(.~wo, labeller="label_both") + 
  geom_smooth(aes(group=alpha), size=.8, se=FALSE, method="lm", alpha=.9) +
  theme_bw() +
  #  scale_y_log10(breaks=c(0.001,.005, 0.01, 0.05)) +
  scale_shape_discrete(expression("Size of validation window "~w[v])) +
  ylab("Observed type I error rate") +
  xlab(expression("Window size for optimization, "~w[o])) +
  scale_colour_brewer(expression("Nominal type I error "~alpha), palette="Set2") +
  facet_wrap(~alpha, labeller="label_both", scales="free") +
  theme(legend.position="bottom")


tp2<- grid.arrange(type2_c_15perc, type2_cs1_c_30perc)

if(file.exists( "./imgs/type2.rds")){
  file.remove("./imgs/type2.rds")
  saveRDS(tp2, file = "./imgs/type2.rds")
}else{
  saveRDS(tp2, file = "./imgs/type2.rds")
}

#a<- readRDS("./presentation/imgs/type2.rds")

# ```