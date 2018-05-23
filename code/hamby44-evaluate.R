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
  dplyr::select(-rest)

h44.summ <- h44 %>% group_by(wo, wv, method, coarseness) %>%
  mutate(
    fails = sum(is.na(p_value))/length(p_value)
  ) %>% filter(!is.na(p_value)) %>%
  summarize(
    AUC = AUC(1-p_value, match),
    EER = EER(1-p_value, match),
    ci = list(AUC_confint(1-p_value, match)),
    fails = fails[1],
    delong = list(pROC::ci.auc(match, 1-p_value))
  )
h44.summ <- h44.summ %>% mutate(
  lower = ci %>% purrr::map_dbl(.f = function(x) x[[1]][1]),
  upper = ci %>% purrr::map_dbl(.f = function(x) x[[1]][2]),
  low_delong = delong %>% purrr::map_dbl(.f = function(x) x[1]),
  upp_delong = delong %>% purrr::map_dbl(.f = function(x) x[3]),
  Method = toupper(method)
)

h44.summ %>% gather(type, value, AUC, EER) %>% 
  filter(wv %in% c(75, 125), coarseness==0.125,
         between(wo, 210, 390)) %>% 
  ggplot(aes(x = wo, y = value, colour=Method, shape=factor(wv)))  +
  geom_point(size=2.5) +
#  geom_point(aes(y = lower), size=1) +
#  geom_point(aes(y = upper), size=1) +
  facet_wrap(~type, scales="free") +
  theme_bw() +
  scale_colour_brewer(palette="Set1", guide = guide_legend(ncol = 1)) +
  geom_smooth(aes(linetype = factor(wv), colour=Method), se=FALSE, span=2) +
  xlab(expression("Size of optimization window "~w[o]~" in pixels")) +
  ylab("Value") +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2,"line")) +
  scale_shape(expression("Size of validation window "~w[v]~" in pixels"), guide = guide_legend(ncol=1))+
  scale_linetype(expression("Size of validation window "~w[v]~" in pixels"), guide = guide_legend(ncol=1)) 


h44.summ %>% filter(wv %in% c(75, 125), coarseness==0.125,
                    between(wo, 210, 390)) %>%
  select(AUC, low_delong, upp_delong) %>% 
  mutate(report = sprintf("%.3f (%.3f, %.3f)", AUC, low_delong, upp_delong)) %>%
  select(report) %>% 
  ungroup() %>%
  mutate(method.wv = paste(method, wv, sep= "-")) %>%
  select(-wv, -method) %>%
  spread(method.wv, report)
