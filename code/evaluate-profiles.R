library(tidyverse)
if (!file.exists("/opt/hamby44/hamby-profiles.rds")) {
  files <- dir("/opt/hamby44/", pattern="rds")
  files <- grep("^chumbley-csafe-profiles", files, value = TRUE)
  
  profiles <- data.frame()
  for (file in files) {
    tmp <- readRDS(file=file.path(path="/opt/hamby44/", file)) 
    #   if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
    profiles <- rbind(profiles, tmp)
  }
  saveRDS(profiles, file="/opt/hamby44/hamby-profiles.rds")
} else {
  profiles <- readRDS("/opt/hamby44/hamby-profiles.rds")
}

profiles$coarse <- 0.3

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



errors <- rbind(errorrate(profiles, 0.001), 
                errorrate(profiles, 0.005), 
                errorrate(profiles, 0.01), 
                errorrate(profiles, 0.05))
errors %>% ggplot(aes(x = wo, y = beta, colour=factor(alpha), shape=factor(wv))) +
  geom_point() + geom_smooth(aes(group=alpha), se=FALSE)


testfailed <- profiles %>% group_by(wo, wv, match) %>% summarize(
  total = n(), 
  failed = sum(is.na(p_value)))
testfailed %>% ggplot(aes(x = wo, y = failed/total, colour=factor(match), shape=factor(wv))) + geom_point() +
  ylim(0, NA)

errors %>% filter(wv %in% c(30,50)) %>%
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


errors %>% filter(wv %in% c(30,50)) %>%
  ggplot(aes(x = wo, y = (actual-alpha)/alpha, colour=factor(alpha))) + 
  geom_hline(aes(yintercept=0), colour="grey30", size=.5) +
  geom_point(aes(shape=factor(wv)), size=2.5) +
  geom_smooth(aes(group=alpha), size=.8, se=FALSE, method="lm", alpha=.9) +
  geom_smooth(aes(group=1), size=.8, se=FALSE, method="lm", alpha=.9) +
  theme_bw() +
  scale_shape_discrete(expression("Size of validation window "~w[v])) +
  ylab("Relative deviation from\nnominal type I error rate") +
  xlab(expression("Window size for optimization, "~w[o])) +
  scale_colour_brewer(expression("Nominal type I error "~alpha), palette="Set2") +
#  facet_wrap(~alpha, labeller="label_both", scales="free") +
  theme(legend.position="bottom")
