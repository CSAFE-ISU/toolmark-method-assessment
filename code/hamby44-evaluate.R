library(tidyverse)
source("code/roc-function.R") # for EER and AUC
files <- dir("data/hamby44/", pattern="rds", recursive=TRUE)

h44 <- data.frame()
for (file in files) {
  tmp <- readRDS(file.path("data/hamby44", file))
  tmp$source <- file
  tmp <- tmp %>% mutate(
    match = replace(match, is.na(match), FALSE)
  )
  h44 <- rbind(h44, tmp)
}

h44 <- h44 %>% separate(source, into=c("method", "rest"), sep="/") %>%
  select(-rest)

h44.summ <- h44 %>% group_by(wo, wv, method, coarseness) %>%
  mutate(
    fails = sum(is.na(p_value))/n() 
  ) %>% filter(!is.na(p_value)) %>%
  summarize(
    auc = AUC(1-p_value, match),
    eer = EER(1-p_value, match),
    ci = list(AUC_confint(1-p_value, match)),
    fails = fails[1],
    roc = list(roc(match, 1-p_value))
  )
h44.summ <- h44.summ %>% mutate(
  lower = ci %>% purrr::map_dbl(.f = function(x) x[[1]][1]),
  upper = ci %>% purrr::map_dbl(.f = function(x) x[[1]][2])
)


h44.summ %>% gather(type, value, auc, eer) %>% 
  filter(wv %in% c(75, 125), coarseness==0.125,
         between(wo, 210, 390)) %>% 
  ggplot(aes(x = wo, y = value, colour=method, shape=factor(wv)))  +
  geom_point(size=2.5) +
#  geom_point(aes(y = lower), size=1) +
#  geom_point(aes(y = upper), size=1) +
  facet_wrap(~type, scales="free") +
  theme_bw() +
  scale_colour_brewer(palette="Set1") +
  geom_smooth(aes(linetype = factor(wv), colour=method), se=FALSE, span=2) +
  xlab("Size of optimization window")
