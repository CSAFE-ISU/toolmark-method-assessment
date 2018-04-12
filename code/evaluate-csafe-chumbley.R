if (!file.exists("/opt/hamby44/csafe-sigs.rds")) {
  files <- dir("/opt/hamby44/", pattern="csv")
  files <- grep("^H44-", files, value = TRUE)
  
  csafe <- data.frame()
  for (file in files) {
    tmp <- read_csv(file=file.path(path="/opt/hamby44/", file)) 
 #   if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
    csafe <- rbind(csafe, tmp)
  }
  saveRDS(csafe, file="/opt/hamby44/csafe-sigs.rds")
} else {
  csafe <- readRDS("/opt/hamby44/csafe-sigs.rds")
}

csafe$coarse <- 1

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



errors <- rbind(errorrate(csafe, 0.001), 
                errorrate(csafe, 0.005), 
                errorrate(csafe, 0.01), 
                errorrate(csafe, 0.05))
errors %>% ggplot(aes(x = wo, y = beta, colour=factor(alpha), shape=factor(wv))) +
  geom_point() + geom_smooth(aes(group=alpha), se=FALSE)
