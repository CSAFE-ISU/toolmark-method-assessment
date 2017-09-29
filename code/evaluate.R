library(tidyverse)
files <- dir("../data", pattern="csv")
all <- data.frame()
for (file in files) {
  tmp <- read.csv(file=file.path(path="../data", file)) 
  if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
  all <- rbind(all, tmp)
}


errorrate <- function(data, alpha) {
summ <- data %>% filter(!is.na(p_value)) %>%
  mutate(signif = p_value < alpha) %>%
  group_by(wv, wo, match, signif) %>% tally()
summ$error <- with(summ, match != signif)
summ$alpha = alpha
#summ %>% ggplot(aes( x  = factor(wv), weight = n, fill=error)) + geom_bar(position="fill") +
#  facet_grid(match~wo, labeller="label_both")

rates <- summ %>% group_by(wv, wo, match, alpha) %>% summarize(
  rate = n[error==TRUE]/sum(n)
)
rates %>% spread(match,rate) %>% rename(
  actual = `FALSE`,
  beta = `TRUE`
)
}

errors <- rbind(errorrate(all, 0.001), 
                errorrate(all, 0.005), 
                errorrate(all, 0.01), 
                errorrate(all, 0.05))
errors %>% filter(wv %in% c(30,50)) %>% ggplot(aes(x = wv, y = beta, colour=factor(alpha))) + geom_point() +
  facet_grid(.~wo, labeller="label_both") + 
  geom_line(aes(group=alpha))

errors %>% filter(wv %in% c(30,50)) %>% ggplot(aes(x = wv, y = actual, colour=factor(alpha))) + geom_point() +
  facet_grid(.~wo, labeller="label_both") + 
  geom_line(aes(group=alpha))

greens <- RColorBrewer::brewer.pal(name="Greens", n = 6)
errors %>% filter(wv %in% c(30,50)) %>%
  ggplot(aes(x = wo, y = beta, colour=factor(alpha))) + 
  geom_point(aes(shape=factor(wv))) +
#  facet_grid(.~wo, labeller="label_both") + 
  geom_smooth(aes(group=alpha), se=FALSE, method="loess") +
  theme_bw() +
  scale_colour_manual("Nominal type I error", values=greens[-(1:2)]) +
  scale_shape_discrete("Size of \nvalidation window") +
  xlab("Window size for optimization") +
  ylab("Type II error rate")

# Relative error for type I
errors %>% filter(wv %in% c(30,50)) %>%
  ggplot(aes(x = wo, y = actual/alpha, colour=factor(alpha))) +
  geom_hline(yintercept = 1, colour="grey30", size=0.5) +
  geom_point() +
  #  facet_grid(.~wo, labeller="label_both") + 
  geom_smooth(aes(group=alpha), se=FALSE, method="loess") +
  theme_bw() +
#  scale_y_log10(breaks=c(0.001,.005, 0.01, 0.05)) +
  ylab("Ratio observed/nominal type I error rate") +
  xlab("Window size for optimization") +
  scale_colour_brewer("Nominal type I error", palette="Set2")

errors %>% filter(wv %in% c(30,50)) %>%
  ggplot(aes(x = wo, y = actual, colour=factor(alpha))) + 
  geom_point(aes(shape=factor(wv))) +
  #  facet_grid(.~wo, labeller="label_both") + 
  geom_smooth(aes(group=alpha), se=FALSE, method="lm") +
  theme_bw() +
  scale_y_log10(breaks=c(0.001,.005, 0.01, 0.05)) +
  scale_shape_discrete("Size of \nvalidation window") +
  ylab("Observed type I error rate") +
  xlab("Window size for optimization") +
  scale_colour_brewer("Nominal type I error", palette="Set2")

# focus on wo = 120

errors %>% filter(wo == 120) %>% ggplot(aes(x = wv, y = actual, colour=factor(alpha))) + geom_point() +
  facet_grid(.~wo, labeller="label_both") + 
  geom_line(aes(group=alpha))
# type I error rate increases 

errors %>% filter(wo == 120) %>% ggplot(aes(x = wv, y = beta, colour=factor(alpha))) + geom_point() +
  facet_grid(.~wo, labeller="label_both") + 
  geom_line(aes(group=alpha))

