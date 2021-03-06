
falspostive<- function(sig_alpha){
  if(is.null(all)){
    library(tidyverse)
    files <- dir("./data/signatures", pattern="csv")
    all <- data.frame()
    for (file in files) {
      tmp <- read.csv(file=file.path(path="./data/signatures", file)) 
      if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
      all <- rbind(all, tmp)
    }
  }
    summ <- all %>% filter(!is.na(p_value)) %>%
      mutate(signif = p_value < sig_alpha) %>%
      group_by(wv, wo, match, signif) %>% tally()
    summ$error <- with(summ, match != signif)
    rates <- summ %>% group_by(wv, wo, match) %>% summarize(
      rate = n[error==TRUE]/sum(n)
    )
    combo<- list(summ, rates)
    names(combo)<- c("summ", "rates")
    return(combo)
}

a<- falspostive(0.05)

a$summ %>% ggplot(aes( x  = factor(wv), weight = n, fill=error)) + geom_bar(position="fill") +
  facet_grid(match~wo, labeller="label_both")
a$rates %>% ggplot(aes(x = wv, y = rate)) + geom_point() + facet_grid(match~wo, labeller="label_both") +
  geom_line()
    
# summ %>% ggplot(aes( x  = factor(wv), weight = n, fill=error)) + geom_bar(position="fill") +
#       facet_grid(match~wo, labeller="label_both")
# 
#     rates %>% ggplot(aes(x = wv, y = rate)) + geom_point() +
#       facet_grid(match~wo, labeller="label_both") +
#       geom_line()

    # library(tidyverse)
    # files <- dir("./data", pattern="csv")
    # all <- data.frame()
    # for (file in files) {
    #   tmp <- read.csv(file=file.path(path="./data", file)) 
    #   if (length(grep("x",names(tmp)) > 0)) tmp <- tmp %>% select(-x)
    #   all <- rbind(all, tmp)
    # }
    # 
    # summ <- all %>% filter(!is.na(p_value)) %>%
    #   mutate(signif = p_value < 0.01) %>%
    #   group_by(wv, wo, match, signif) %>% tally()
    # summ$error <- with(summ, match != signif)
    # summ %>% ggplot(aes( x  = factor(wv), weight = n, fill=error)) + geom_bar(position="fill") +
    #   facet_grid(match~wo, labeller="label_both")
    # 
    # rates <- summ %>% group_by(wv, wo, match) %>% summarize(
    #   rate = n[error==TRUE]/sum(n)
    # )
    # rates %>% ggplot(aes(x = wv, y = rate)) + geom_point() +
    #   facet_grid(match~wo, labeller="label_both") + 
    #   geom_line()
    
aa1<- read.csv2("./data/signatures/chumbley-out-wo-80-wv-50.csv",header = TRUE, sep = ",")
xtabs(~signif+match, data =aa1)


aa2<- read.csv2("./data/signatures/chumbley_wo-200_wv-50.csv",header = TRUE, sep = ",")
xtabs(~signif+match, data =aa2)

aa3<- read.csv2("./data/signatures/chumbley_wo-320_wv-50.csv",header = TRUE, sep = ",")
xtabs(~signif+match, data =aa3)
    

