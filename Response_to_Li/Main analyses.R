library(tidyverse)
library(ggplot2)
library(metafor)
library(mgcv)
library(scales)


#######################################################################
#                       Input data
#######################################################################

pop2000 <- read.csv("Data\\2000 prefecture population.csv")  # source: Sichuan Population Census Office, Tabulation on the 2000 Population Census of Sichuan Province (China Statistics Press, Beijing, 2002).

TB.RR.agg <- read.csv("Data\\TB.RR.agg.csv") # source: output from the analyses in Cheng Q, Trangucci R, Nelson KN, Fu W, Collender PA, Head JR, Hoover CM, Skaff NK, Li T, Li X, You Y. Prenatal and early-life exposure to the Great Chinese Famine increased the risk of tuberculosis in adulthood across two generations. Proceedings of the National Academy of Sciences. 2020 Nov 3;117(44):27549-55.

surv.rate <- read.csv("Data\\Cohort shrinkage estimates.csv") # estimated with 1964, 1982, and 2000 census data, see Survival_fraction_estimates.R for the process

extra.births <- read.csv("Data\\Delayed births estimates.csv") # estimated with Sichuan fertility rate, see Excessive_births_estimates.R






#######################################################################
#                           Figure 1
#######################################################################
pop2000 <- pop2000 %>%
  left_join(TB.RR.agg[, c("Pref", "Pref.name")]) %>%
  mutate(Pref.name = factor(Pref.name, levels = TB.RR.agg$Pref.name)) %>%
  mutate(birthyear = 2000 - age) %>%
  group_by(Pref.name, birthyear) %>%
  summarise(pop = sum(pop)) %>%
  left_join(surv.rate[, c("BirthYear", "Surv.overall")], by = c("birthyear" = "BirthYear")) %>%
  mutate(pop.adjusted = pop/Surv.overall)

pop2000 %>%
  ggplot(aes(x = birthyear, y = pop.adjusted/1000, col = as.factor(Pref.name), group = as.factor(Pref.name))) +
  # geom_vline(xintercept = c(1953, 1958, 1963, 1968)) +
  geom_rect(aes(xmin = 1953, xmax = 1958, ymin = -Inf, ymax = Inf),  col = FALSE, fill = "gray85") +
  geom_rect(aes(xmin = 1958, xmax = 1963, ymin = -Inf, ymax = Inf),  col = FALSE, fill = "gray70") +
  geom_rect(aes(xmin = 1963, xmax = 1968, ymin = -Inf, ymax = Inf),  col = FALSE, fill = "gray85") +
  geom_line() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5))  +
  xlab("Year of birth") +
  ylab("Birth cohort size (thousands)") +
  labs(col = "") +
  guides(col = guide_legend(ncol = 1)) +
  annotate("text", x= 1955.5, y = 360, label = "Pre-famine", angle = 90, hjust = 1) + 
  annotate("text", x= 1960.5, y = 360, label = "Famine", angle = 90, hjust = 1) + 
  annotate("text", x= 1965.5, y = 360, label = "Post-famine", angle = 90, hjust = 1) 

# ggsave("Figure 1.pdf", width = 8, height = 6)







#######################################################################
#                       Estimate CSSI
#######################################################################
CSSI <- pop2000 %>%
  filter(birthyear >= 1953 & birthyear <= 1967) %>%
  group_by(Pref.name, birthyear) %>%
  summarise(pop = sum(pop)) %>%
  mutate(FM = ifelse(birthyear <= 1962 & birthyear >= 1958, "Famine", ifelse(birthyear <= 1957, "Prefamine", "Postfamine"))) %>%
  ungroup() %>%
  group_by(Pref.name, FM) %>%
  summarize(pop.mean = mean(pop)) %>%
  spread(FM, pop.mean) %>%
  ungroup() %>%
  mutate(postfamine.more.birth = extra.births$excessive) %>%
  mutate(CSSI.pre_fixed = (Prefamine - Famine)/Prefamine, CSSI.pre_post_mean = ((Postfamine + Prefamine)/2 - Famine)/((Postfamine + Prefamine)/2), CSSI.pre_post_mean_adjusted = ((Postfamine + Prefamine - postfamine.more.birth/5)/2 - Famine)/((Postfamine + Prefamine- postfamine.more.birth/5)/2))  # the excessive birth was divided by 5 since Postfamine and Prefamine is average across 5 yrs



## using only data before the famine
# using 10-yr's data
secular.pre_10_lm <- pop2000 %>%
  filter(birthyear <= 1952 & birthyear >= 1943) %>%
  group_by(Pref.name) %>%
  summarise(cohort.predict.mean = mean(predict(lm(pop.adjusted ~ birthyear), newdata = data.frame(birthyear = 1958:1962)))) %>%
  ungroup() %>%
  mutate(CSSI = (cohort.predict.mean - CSSI$Famine)/cohort.predict.mean)

# using 20-yr's data
secular.pre_20_lm <- pop2000 %>%
  filter(birthyear <= 1952 & birthyear >= 1933) %>%
  group_by(Pref.name) %>%
  summarise(cohort.predict.mean = mean(predict(lm(pop.adjusted ~ birthyear), newdata = data.frame(birthyear = 1958:1962)))) %>%
  ungroup() %>%
  mutate(CSSI = (cohort.predict.mean - CSSI$Famine)/cohort.predict.mean)

## using data before and after the famine
# using 5yrs before and 5yrs after the famine
secular.pre_5_post_5_lm <- pop2000 %>%
  filter((birthyear <= 1952 & birthyear >= 1948) | (birthyear <= 1972 & birthyear >= 1968)) %>%
  group_by(Pref.name) %>%
  summarise(cohort.predict.mean = mean(predict(lm(pop.adjusted ~ birthyear), newdata = data.frame(birthyear = 1958:1962)))) %>%
  ungroup() %>%
  mutate(CSSI = (cohort.predict.mean - CSSI$Famine)/cohort.predict.mean)

# all data before and 5 yrs after GAM
secular.pre53_post_5_gam <- pop2000 %>%
  filter((birthyear <= 1952 & birthyear >= 1903) | (birthyear <= 1972 & birthyear >= 1968)) %>%
  group_by(Pref.name) %>%
  summarise(cohort.predict.mean = mean(predict(lm(pop.adjusted ~ birthyear), newdata = data.frame(birthyear = 1958:1962)))) %>%
  ungroup() %>%
  mutate(CSSI = (cohort.predict.mean - CSSI$Famine)/cohort.predict.mean)


# using 10yrs before and 10 yrs after the famine
secular.pre_10_post_10_lm <- pop2000 %>%
  filter((birthyear <= 1952 & birthyear >= 1943) | (birthyear <= 1977 & birthyear >= 1968)) %>%
  group_by(Pref.name) %>%
  summarise(cohort.predict.mean = mean(predict(lm(pop.adjusted ~ birthyear), newdata = data.frame(birthyear = 1958:1962)))) %>%
  ungroup() %>%
  mutate(CSSI = (cohort.predict.mean - CSSI$Famine)/cohort.predict.mean)

# all data before and 10 yrs after GAM
secular.pre_53_post_10_gam <- pop2000 %>%
  filter((birthyear <= 1952 & birthyear >= 1903) | (birthyear <= 1978 & birthyear >= 1968)) %>%
  group_by(Pref.name) %>%
  summarise(cohort.predict.mean = mean(predict(lm(pop.adjusted ~ birthyear), newdata = data.frame(birthyear = 1958:1962)))) %>%
  ungroup() %>%
  mutate(CSSI = (cohort.predict.mean - CSSI$Famine)/cohort.predict.mean)


CSSI <- CSSI %>%
  mutate(CSSI.pre_10_lm = secular.pre_10_lm$CSSI, CSSI.pre_20_lm = secular.pre_20_lm$CSSI, CSSI.pre_5_post_5_lm = secular.pre_5_post_5_lm$CSSI, CSSI.pre_53_post_5_gam = secular.pre53_post_5_gam$CSSI, CSSI.pre_10_post_10_lm = secular.pre_10_post_10_lm$CSSI, CSSI.pre_53_post_10_gam= secular.pre_53_post_10_gam$CSSI)

CSSI %>%
  gather("CSSI_type", "Value", CSSI.pre_fixed:CSSI.pre_53_post_10_gam) %>%
  mutate(CSSI_type = factor(CSSI_type, levels = unique(CSSI_type))) %>%
  ggplot(aes(x = CSSI_type, y = Value)) +
  geom_boxplot(notch = TRUE) +
  coord_flip() +
  theme_bw()

apply(CSSI[, 6:14], 2, function(x) round(quantile(x, c(0.5, 0.25, 0.75)), 2))



#######################################################################
#                 Estimate metaregression coefficient
#######################################################################

# TB result 
CSSI <- CSSI %>%
  mutate(F1.log.mean = TB.RR.agg$F1.log.mean, F1.log.sd = TB.RR.agg$F1.log.sd)

final.result <- NULL
CSSI.option <- colnames(CSSI)[6:14]
for(i in 1:length(CSSI.option))
{
  meta.reg.TB.F1 <- rma(F1.log.mean, sei = F1.log.sd, mods =  CSSI[,CSSI.option[i]], data = CSSI)
  
  current.result <- data.frame(CSSI = CSSI.option[i], coef = meta.reg.TB.F1$beta[2], coef.lb = meta.reg.TB.F1$ci.lb[2], coef.ub = meta.reg.TB.F1$ci.ub[2])
  final.result <- rbind(final.result, current.result)
}


final.result %>%
  mutate(coef = round(coef, 2), coef.lb = round(coef.lb, 2), coef.ub = round(coef.ub, 2))

