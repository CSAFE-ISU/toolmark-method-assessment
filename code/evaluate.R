library(tidyverse)
files <- dir("../data", pattern="csv")
all <- data.frame()
for (file in files) {
  tmp <- read.csv(file=file.path(path="../data", file)) 
  if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
  all <- rbind(all, tmp)
}

summ <- all %>% filter(!is.na(p_value)) %>% group_by(wv, wo, match, signif) %>% tally()
summ$error <- with(summ, match != signif)
summ %>% ggplot(aes( x  = factor(wv), weight = n, fill=error)) + geom_bar(position="fill") +
  facet_grid(match~wo, labeller="label_both")
