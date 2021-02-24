#==== processing population data ======
pop2000 <- read.csv("Data\\2000 prefecture population.csv") 
TB.RR.agg <- read.csv("Data\\TB.RR.agg.csv") 

pop2000 <- pop2000 %>%
  left_join(TB.RR.agg[, c("Pref", "Pref.name")]) %>%
  mutate(Pref.name = factor(Pref.name, levels = TB.RR.agg$Pref.name)) %>%
  mutate(birthyear = 2000 - age) %>%
  group_by(Pref.name, birthyear) %>%
  summarise(pop = sum(pop))




#====== estimate excessive births ======
fert.rate <- read.csv("Data\\SichuanFertility.csv") # source: A. J. Coale, S. L. Chen, Basic data on fertility in the provinces of China, 1940-82.  (1987).

ggplot(fert.rate, aes(x = BIrthYearOfChild, y = Sichuan)) +
  geom_line() +
  geom_point(col = "red") +
  geom_vline(xintercept = c(1952.5, 1957.5, 1962.5, 1967.5), col = "gray") +
  theme_bw() +
  annotate("text", x= 1955.5, y = 8.5, label = "Pre-famine", angle = 90, hjust = 1, vjust = 0) + 
  annotate("text", x= 1960.5, y = 8.5, label = "Famine", angle = 90, hjust = 1, vjust = 0) + 
  annotate("text", x= 1965.5, y = 8.5, label = "Post-famine", angle = 90, hjust = 1, vjust = 0) 

fert.rate.toues <- fert.rate %>%
  filter(!BIrthYearOfChild %in% 1958:1962)

# only estimate the excessive fertility rate for 1963, since according to the above figure, the rate dropped to normal level in 1964
fert.lm <- lm(Sichuan ~ BIrthYearOfChild + I(BIrthYearOfChild == 1963), data = fert.rate.toues) %>%
  summary()

# plot result
ggplot(fert.rate, aes(x = BIrthYearOfChild, y = Sichuan)) +
  geom_line() +
  geom_point(col = "red") +
  geom_vline(xintercept = c(1952.5, 1957.5, 1962.5, 1967.5), col = "gray") +
  theme_bw() +
  annotate("text", x= 1955.5, y = 8.5, label = "Pre-famine", angle = 90, hjust = 1, vjust = 0) + 
  annotate("text", x= 1960.5, y = 8.5, label = "Famine", angle = 90, hjust = 1, vjust = 0) + 
  annotate("text", x= 1965.5, y = 8.5, label = "Post-famine", angle = 90, hjust = 1, vjust = 0) +
  geom_abline(intercept = fert.lm$coefficients[1,1], slope = fert.lm$coefficients[2,1])





#==== estimate excessive birth ===
Excessive.birth.1963 <- pop2000 %>%
  filter(birthyear == 1963) %>%
  mutate(excessive = pop*fert.lm$coefficients[3,1]/fert.rate.toues$Sichuan[fert.rate.toues$BIrthYearOfChild == 1963])

write.csv(Excessive.birth.1963, "Data\\Delayed births estimates.csv", row.names = FALSE)
