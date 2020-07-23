library(readr)
library(ggplot2)
library(cowplot)
library(raster)
library(tidyverse)
library(tmap)
library(classInt)
library(RColorBrewer)
library(scales)
library(APCG1)  # can be downloaded from https://www.math.uh.edu/~fuw/index_Faculty%20page%20Res_new.html
library(mgcv)
library(fANCOVA)
library(foreach)
library(doSNOW)
library(metafor)
library(grid)


theme_set(theme_cowplot())

source("Functions.R")

#######################################################################
#                       Input data
#######################################################################

disease.dat <- read_csv(".\\Data\\Disease and population data.csv")   # prefecture-level disease and population data
famine.intensity <- read_csv(".\\Data\\Famine intensity by pref.csv")  # prefecture-level famine intensity
SC.pop <- read_csv(".\\Data\\Sichuan population projection.csv")   # population projection of the whole province
pref.thres <- read_csv(".\\Data\\F2 time by prefecture population data.csv")



SC.shp <- shapefile(".\\Data\\Sichuan_City.shp")
colnames(SC.shp@data)[2] <- "Pref"
SC.shp@data$Pref.name <- c("Chengdu", "Zigong","Panzhihua","Luzhou","Deyang","Mianyang","Guangyuan","Suining","Neijiang","Leshan","Nanchong","Meishan","Yibin","Guang'an","Dazhou","Yaan","Bazhong","Ziyang","Aba","Ganzi","Liangshan")




#######################################################################
#                       Figures
#######################################################################

#===========================================================
#    Figure 1 - CSSI and annual incidence rate of PTB
#===========================================================
# calculate annual incidence rate of PTB
TB.inc.avg <- disease.dat %>% group_by(Pref, Pref.name) %>% dplyr::summarise(Population = sum(Population), CaseCount = sum(CaseCount.TB)) %>% mutate(incidence = CaseCount/Population*100000) %>% ungroup() %>% mutate(Pref = as.numeric(Pref)) %>% left_join(famine.intensity) 

SC.shp@data <- left_join(SC.shp@data, TB.inc.avg)

# plot the CSSI
f1.1 <- tm_shape(SC.shp) + tm_fill(c("CSSI"), breaks = classIntervals(SC.shp$CSSI, n = 5, style = "jenks")$brks, title = "", legend.format = list(fun = function(x) round(x, 2))) + tm_borders(col = "black") + tm_text("Pref.name", size = 0.7) + tm_legend(position = c("left","bottom")) + tm_layout(title = "(A) Cohort size shrinkage index", inner.margins = c(0.02, 0.02, 0.08, 0.02), title.size = 1.1, frame = FALSE)

# plot the incidence rate
f1.2 <- tm_shape(SC.shp) + tm_fill(c("incidence"), breaks = classIntervals(SC.shp$incidence, n = 5, style = "jenks")$brks, title = "", legend.format = list(fun = function(x) round(x, 0))) + tm_borders(col = "black") + tm_text("Pref.name", size = 0.7) + tm_legend(position = c("left","bottom")) + tm_layout(title = "(B) Annual mean incidence rate (1/100,000)", inner.margins = c(0.02, 0.02, 0.08, 0.02), title.size = 1.1, frame = FALSE, legend.title.size = 0.9)

tmap_arrange(f1.1, f1.2)
# tmap_save(tmap_arrange(f1.1, f1.2), ".\\Results\\Figure 1 CSSI and TB inc map.pdf", width = 5, height = 8.44, units = "in")


#===========================================================
#                   Figure 2 - heatmap
#===========================================================
# calculate province-level incidence rate by 
case.TB.Sichuan <- disease.dat %>% group_by(Age, DiagnoseYear, Sex) %>% dplyr::summarise(CaseCount = sum(CaseCount.TB)) %>% left_join(SC.pop) %>% ungroup() %>% mutate(incidence = CaseCount/Population*100000)

famine.polygons.F1 <- data.frame(ID = rep(1:5, each = 5), x = c(c(2004.5,2007.5, 2007.5, 2004.5, 2004.5) + rep(seq(0,by = 3, length.out = 4), each = 5), 2016.5,2018.5,2018.5,2016.5,2016.5), y = c(15.5, 15.5, 16.5, 16.5, 15.5) + rep(seq(from = 0, by = 1, length.out = 5), each = 5)-3)
famine.polygons.F2 <- data.frame(ID = rep(1:5, each = 5), x = c(c(2004.5,2007.5, 2007.5, 2004.5, 2004.5) + rep(seq(0,by = 3, length.out = 4), each = 5), 2016.5,2018.5,2018.5,2016.5,2016.5), y = c(8.5, 8.5, 9.5, 9.5, 8.5) + rep(seq(from = 0, by = 1, length.out = 5), each = 5)-3)

qn = quantile(case.TB.Sichuan$incidence[case.TB.Sichuan$Sex == 1], seq(0.1,0.9,length.out = 18), na.rm = TRUE)
qn01 <- sort(rescale(c(qn, range(case.TB.Sichuan$incidence[case.TB.Sichuan$Sex == 1]))) )
f2.1 <- ggplot(case.TB.Sichuan[case.TB.Sichuan$Sex == 1,], aes(DiagnoseYear, Age)) + 
  geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
  geom_tile(aes(fill = incidence), color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) + 
  annotate("text", x = 2017.5, y = 18, label = "F1") + 
  annotate("text", x = 2017.5, y = 11, label = "F2") + 
  geom_polygon(data = famine.polygons.F1, aes(x = x, y= y, group = ID), fill = NA, col = "black") + 
  geom_polygon(data = famine.polygons.F2, aes(x = x, y= y, group = ID), fill = NA, col = "black") +
  labs(x = "Year of diagnosis", y = "Age group", fill = "Incidence\nrate (1/100,000)") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position="bottom", legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 10)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
  scale_y_discrete(expand = c(0, 0))

# female
qn = quantile(case.TB.Sichuan$incidence[case.TB.Sichuan$Sex == 0], seq(0.1,0.9,length.out = 18), na.rm = TRUE)
qn01 <- sort(rescale(c(qn, range(case.TB.Sichuan$incidence[case.TB.Sichuan$Sex == 0]))) )
f2.2 <- ggplot(case.TB.Sichuan[case.TB.Sichuan$Sex == 0,], aes(DiagnoseYear, Age)) + 
  geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
  geom_tile(aes(fill = incidence), color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) + 
  annotate("text", x = 2017.5, y = 18, label = "F1") + 
  annotate("text", x = 2017.5, y = 11, label = "F2") + 
  geom_polygon(data = famine.polygons.F1, aes(x = x, y= y, group = ID), fill = NA, col = "black") + 
  geom_polygon(data = famine.polygons.F2, aes(x = x, y= y, group = ID), fill = NA, col = "black") +
  labs(x = "Year of diagnosis", y = "Age group", fill = "Incidence\nrate (1/100,000)") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position="bottom", legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 10)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
  scale_y_discrete(expand = c(0, 0))

Lexis <- plot_grid(f2.1, f2.2, labels = c("(A) Male","(B) Female"), nrow = 1)
Lexis
# save_plot(".\\Results\\Figure 2 heatmap.pdf", Lexis, base_asp = 1.12, base_height = 7.2)


#===========================================================
#              Figure 3 - APC results PTB
#===========================================================
# data processing to be used for APC model
# incidence rate tables
inc.table.male <- case.TB.Sichuan %>% filter(Sex == 1) %>% dplyr::select(Age, DiagnoseYear, incidence) %>% tidyr::spread(DiagnoseYear, incidence) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,79,3)) %>% column_to_rownames(var = "Age")
inc.table.female <- case.TB.Sichuan %>% filter(Sex == 0) %>% dplyr::select(Age, DiagnoseYear, incidence) %>% tidyr::spread(DiagnoseYear, incidence) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,79,3)) %>% column_to_rownames(var = "Age")

# population size tables
male.pop <- SC.pop %>%
  filter(Sex == 1) %>%
  spread(DiagnoseYear, Population) %>%
  select(-Age, -Sex)
female.pop <- SC.pop %>%
  filter(Sex == 0) %>%
  spread(DiagnoseYear, Population) %>%
  select(-Age, -Sex)

# APC model results
male.apc <- apcglmkfit(r = inc.table.male, n.risk = male.pop, header=TRUE, p0 = 3, k = 3, fam = "qlik", amin = 1, pmin = 2005, cmin = 1927, agapyr = 3, pgapyr = 1) 
female.apc <- apcglmkfit(r = inc.table.female, n.risk = female.pop, header=TRUE, p0 = 3, k = 3, fam = "qlik", amin = 1, pmin = 2005, cmin = 1927, agapyr = 3, pgapyr = 1)

male.apc.result <- male.apc %>% reshape.APC.result() %>% mutate(Sex = "Male")
female.apc.result <- female.apc %>% reshape.APC.result() %>% mutate(Sex = "Female")
result.mf.TB <- rbind(male.apc.result, female.apc.result) %>% mutate(X = as.numeric(X))
result.mf.TB$Sex <- factor(result.mf.TB$Sex, levels = c("Male", "Female"))

# TB age effect
f3.1 <- result.mf.TB %>% 
  filter(Type == "Age" & X < 79) %>% 
  ggplot(aes(x = X, y = parameter, col = Sex)) + 
  geom_line() +
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd, fill = Sex), col = NA, alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("(A) Age effect") + 
  theme(plot.title = element_text(hjust = 0)) +
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Age at diagnosis", y = expression(alpha[i]), col = "", fill = "")+
  scale_color_manual(values = rev(gg_color_hue(2)))+
  scale_fill_manual(values = rev(gg_color_hue(2)))

# period effect
f3.2 <- result.mf.TB %>% 
  filter(Type == "Period") %>% 
  ggplot(aes(x = X, y = parameter, col = Sex)) + 
  geom_line() +
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd, fill = Sex), col = NA, alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("(B) Period effect") + 
  theme(plot.title = element_text(hjust = 0)) +
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Year of diagnosis", y = expression(pi[j]), col = "", fill = "")+
  scale_color_manual(values = rev(gg_color_hue(2)))+
  scale_fill_manual(values = rev(gg_color_hue(2)))

# cohort effect
TB.smooth <- province.cohort.smooth(result.mf.TB, F1.start = 1957, F1.end = 1963, F2.start = 1975, F2.end = 1987, method = "GAM")

line_types <- c("Estimated"=1,"Expected"=2)

f3.3 <- result.mf.TB %>% 
  filter(Type == "Cohort" & X >= 1930 & X <= 2000) %>% 
  ggplot(aes(x = X, y = parameter, col = Sex)) + 
  geom_line(aes(linetype="Estimated")) +
  geom_point(data = result.mf.TB[result.mf.TB$Type == "Cohort" & result.mf.TB$X == 1981,], aes(x = X, y = parameter, col = Sex), size = 2) +
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd, fill = Sex), col = NA, alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_line(data = TB.smooth[TB.smooth$Cohort >= 1930 & TB.smooth$Cohort <= 2000,], aes(x = Cohort, y = cohort.predict, linetype="Expected")) +
  geom_vline(xintercept = 1960) + 
  geom_vline(xintercept = 1981) + 
  geom_point(data = result.mf.TB[result.mf.TB$Type == "Cohort" & result.mf.TB$X == 1960,], aes(x = X, y = parameter, col = Sex), size = 2) + 
  ggtitle("(C) Cohort effect") + 
  theme(plot.title = element_text(hjust = 0), legend.position = c(0.8, 0.9), legend.margin = unit(0, "cm")) +
  labs(x = "Year of birth", y = expression(gamma[k]), col = "", fill = "") + 
  annotate("text", 1960, 0.50, label = "F1", size = 6) +
  annotate("text", 1981, 0.50, label = "F2", size = 6) +
  ylim(-0.55, 0.50) +
  scale_color_manual(values = rev(gg_color_hue(2)))+
  scale_fill_manual(values = rev(gg_color_hue(2)))+ 
  scale_linetype_manual(name = "", values=line_types)

left <- plot_grid(f3.1, f3.2, nrow = 2)
plot_grid(left, f3.3, nrow = 1, rel_widths = c(1,2))

# save_plot(".\\Results\\Figure 3 APC results.pdf", plot_grid(left, f3.3, nrow = 1, rel_widths = c(1,2)), base_height  = 6.5, base_aspect_ratio = 1.5)

#===========================================================
#              Figure 4 - famine effect PTB
#===========================================================
#====================== calculate province-level RR ====================
# record the variance-covariance matrix from IE
result.mf.TB.vcov <- list(male.cov = list(a.vcov = male.apc$a.vcov, p.vcov = male.apc$p.vcov, c.vcov = male.apc$c.vcov), female.cov = list(a.vcov = female.apc$a.vcov, p.vcov = female.apc$p.vcov, c.vcov = female.apc$c.vcov))

#  Province IRR
SC.pop.long <- data.frame(Age = rep(rep(seq(10,79,3), each = 14),2), DiagnoseYear = rep(2005:2018, 24*2), Population = c(unlist(t(male.pop)), unlist(t(female.pop))), Sex = rep(c("Male","Female"), each = 14*24))

Sichuan.TB.case.long <- case.TB.Sichuan %>% dplyr::select(Age, DiagnoseYear, CaseCount, Sex) %>% mutate(Age = rep(seq(10,79,3), each = 14*2), Cohort = DiagnoseYear - Age, Sex = fct_recode(as.factor(Sex), "Male" = "1","Female" = "0"))

# since Monte Carlo simulation is used here, the results may not be exactly the same from run to run; this function paralellize the calculation with the foreach package. The number of cores used can be changed by setting the ncore parameter
TB.RR <- province.RR(result.mf.TB, result.mf.TB.vcov, F1.start = 1957, F1.end = 1963, F2.start = 1975, F2.end = 1987, F2 = 1981, n.rep = 1000, method = "GAM", SC.pop.long, Sichuan.TB.case.long, ncore = 7)

# print the result
TB.RR.SC.agg <- TB.RR %>% gather(key = "Generation", value = "RR", F1.RR:F2.avertcase) %>% group_by(Sex, Generation) %>% summarise(mean = median(RR), low = quantile(RR, 0.025), high = quantile(RR, 0.975)) %>% ungroup() %>% mutate(Disease = "PTB")


#=================== calculate prefecture-level RR ===================
case.TB.nosex <- disease.dat %>% group_by(Pref, Age, DiagnoseYear) %>% dplyr::summarise(CaseCount = sum(CaseCount.TB), Population = sum(Population)) %>% ungroup() %>% mutate(incidence = CaseCount/Population*100000)

case.STI.nosex <- disease.dat %>% group_by(Pref, Age, DiagnoseYear) %>% dplyr::summarise(CaseCount = sum(CaseCount.STBBI), Population = sum(Population)) %>% ungroup() %>% mutate(incidence = CaseCount/Population*100000)

# prefecture IRR for both TB and STI
pref.TB <- NULL
pref.STI <- NULL
pref.TB.vcov <- vector("list", 21)
pref.STI.vcov <- vector("list", 21)
for(i in 1:21)
{
  inc.table.TB <- case.TB.nosex %>% filter(Pref == SC.shp@data$Pref[i]) %>% dplyr::select(Age, DiagnoseYear, incidence) %>% tidyr::spread(DiagnoseYear, incidence) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,79,3)) %>% column_to_rownames(var = "Age")
  
  inc.table.STI <- case.STI.nosex  %>% filter(Pref == SC.shp@data$Pref[i]) %>% mutate(DiagnoseYear = rep(c(rep(seq(2006,2015,3), each = 3), 2018,2018),24), Age = fct_recode(Age, "[69,Inf)" = "[69,72)", "[69,Inf)" = "[72,75)","[69,Inf)" = "[75,78)","[69,Inf)" = "[78,Inf)")) %>% group_by(Age, DiagnoseYear) %>% dplyr::summarise(CaseCount = sum(CaseCount), Population = sum(Population)) %>% mutate(incidence = CaseCount/Population*100000) %>% ungroup() %>% dplyr::select(Age, DiagnoseYear, incidence) %>% spread(DiagnoseYear, incidence) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,70,3)) %>% column_to_rownames(var = "Age")
  
  pop.table <- case.TB.nosex %>% filter(Pref == SC.shp@data$Pref[i]) %>% dplyr::select(Age, DiagnoseYear, Population) %>% mutate(Population = round(Population, 0)) %>% spread(DiagnoseYear, Population) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,79,3)) %>% column_to_rownames(var = "Age")
  
  pop.table.STI <- case.STI.nosex  %>% filter(Pref == SC.shp@data$Pref[i])%>% mutate(DiagnoseYear = rep(c(rep(seq(2006,2015,3), each = 3), 2018,2018),24), Age = fct_recode(Age, "[69,Inf)" = "[69,72)", "[69,Inf)" = "[72,75)","[69,Inf)" = "[75,78)","[69,Inf)" = "[78,Inf)")) %>% group_by(Age, DiagnoseYear) %>% dplyr::summarise(Population = round(sum(Population),0)) %>% ungroup() %>% dplyr::select(Age, DiagnoseYear, Population) %>% spread(DiagnoseYear, Population) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,70,3)) %>% column_to_rownames(var = "Age")
  
  current.APC.TB.model <- apcglmkfit(r = inc.table.TB, n.risk = pop.table, header=TRUE, p0 = 3, k = 3, fam = "qlik", amin = 1, pmin = 2005, cmin = 1927, agapyr = 3, pgapyr = 1)
  
  current.APC.TB <- current.APC.TB.model %>% reshape.APC.result() %>% mutate(Pref = SC.shp@data$Pref[i])

  current.APC.STI.model <- apcglmkfit(r = inc.table.STI, header = TRUE, n.risk = pop.table.STI, fam = "qlik", p0 = 1, k = 1)
  current.APC.STI <- current.APC.STI.model %>% reshape.APC.result() %>% mutate(Pref = SC.shp@data$Pref[i])
  
  pref.TB <- rbind(pref.TB, current.APC.TB)
  pref.STI <- rbind(pref.STI, current.APC.STI)
  
  pref.TB.vcov[[i]] <- list(a.vcov = current.APC.TB.model$a.vcov, p.vcov = current.APC.TB.model$p.vcov, c.vcov = current.APC.TB.model$c.vcov)
  pref.STI.vcov[[i]] <- list(a.vcov = current.APC.STI.model$a.vcov, p.vcov = current.APC.STI.model$p.vcov, c.vcov = current.APC.STI.model$c.vcov)
}


Sichuan.TB.case.long <- case.TB.nosex %>% dplyr::select(Pref, Age, DiagnoseYear, CaseCount) %>% mutate(Age = rep(rep(seq(10,79,3), each = 14),21), Cohort = DiagnoseYear - Age)

Sichuan.pop.long <- case.TB.nosex %>% dplyr::select(Pref, Age, DiagnoseYear, Population) %>% mutate(Age = rep(rep(seq(10,79,3), each = 14),21), Cohort = DiagnoseYear - Age)


# prefecture RR
# since Monte Carlo simulation is used here, the results may not be exactly the same from run to run; this function paralellize the calculation with the foreach package. The number of cores used can be changed by setting the ncore parameter
TB.RR.all <- pref.RR(pref.TB, pref.TB.vcov, F1.start = 1957, F1.end = 1963, F2.start = pref.thres$low.cohort, F2.end = pref.thres$high.cohort, F2 = pref.thres$mid.cohort, n.rep = 1000, method = "GAM", Sichuan.pop.long, Sichuan.TB.case.long, collapse = FALSE) 

# summarize the result
TB.RR.agg <- TB.RR.all %>% group_by(Pref) %>% summarise(F1.log.mean = mean(log(F1.RR)), F1.log.sd = sd(log(F1.RR)), F2.log.mean = mean(log(F2.RR)), F2.log.sd = sd(log(F2.RR)), F1.mean = mean(F1.RR), F2.mean = mean(F2.RR), F1.low = quantile(F1.RR, 0.025), F1.high = quantile(F1.RR, 0.975), F2.low = quantile(F2.RR, 0.025), F2.high = quantile(F2.RR, 0.975)) %>% left_join(famine.intensity[,c("Pref","CSSI")]) %>% left_join(SC.shp@data)

TB.RR.agg

# set the values for the prefectures with inaccurate population prediction to be NA
TB.RR.agg$F2.log.mean[TB.RR.agg$Pref %in% c(5101,5108,5116,5120)] <- NA
TB.RR.agg$F2.log.sd[TB.RR.agg$Pref %in% c(5101,5108,5116,5120)] <- NA
TB.RR.agg$F2.mean[TB.RR.agg$Pref %in% c(5101,5108,5116,5120)] <- NA

#========= meta regression =========
meta.reg.TB.F1 <- rma(F1.log.mean, sei = F1.log.sd, mods =  ~ CSSI, data = TB.RR.agg)
# permutest(meta.reg.TB.F1)
preds.TB.F1 <- predict(meta.reg.TB.F1, newmods=seq(0.25,0.6,0.01))
preds.TB.F1 <- data.frame(preds.TB.F1)
preds.TB.F1$x <- seq(0.25,0.6,0.01)

TB.RR.agg$size1 <- (1/TB.RR.agg$F1.log.sd)/(1/max(TB.RR.agg$F1.log.sd))
TB.RR.agg$size2 <- (1/TB.RR.agg$F2.log.sd)/(1/max(TB.RR.agg$F2.log.sd, na.rm = TRUE))


Q1 <- predict(meta.reg.TB.F1, newmods = 43.82262/100)
Q3 <- predict(meta.reg.TB.F1, newmods = 52.07122/100)

slope.sim <- rnorm(100000, summary(meta.reg.TB.F1)$beta[2], summary(meta.reg.TB.F1)$se[2])
intercept.sim <- rnorm(100000, summary(meta.reg.TB.F1)$beta[1], summary(meta.reg.TB.F1)$se[1])

Q1.RR.sim <- exp(intercept.sim + slope.sim*43.82262/100)
Q3.RR.sim <- exp(intercept.sim + slope.sim*52.07122/100)

IQR.F1 <- (Q3.RR.sim - Q1.RR.sim)/Q1.RR.sim
median(IQR.F1)*100
quantile(IQR.F1, 0.025)*100
quantile(IQR.F1, 0.975)*100






#============= plot the results ==================
# panel A RR
TB.RR.SC.agg$Sex <- factor(TB.RR.SC.agg$Sex, levels = c("Male","Female"))
f4.1 <- TB.RR.SC.agg %>% filter(Generation == "F1.RR") %>% ggplot(aes(x = Sex, y = mean, ymin = low, ymax = high, col = Sex)) + 
  geom_errorbar(aes(width = .2)) + 
  geom_point() + 
  scale_color_manual(name = "",values = rev(gg_color_hue(2))) + 
  ylim(1, 1.3) + 
  ylab("IRR") + 
  geom_hline(yintercept = 1, linetype = 3) + 
  xlab("") + 
  theme(legend.position = c(0.5, 0.95), plot.title = element_text(hjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0.4,0.1,0,0.1), "cm")) +
  ggtitle("(A) F1 IRR by sex")

# Panel B IRR map
SC.shp@data <- left_join(SC.shp@data, TB.RR.agg[,c("Pref","F1.mean", "F2.mean")])

tm_shape(SC.shp) + tm_fill(c("F1.mean"), breaks = classIntervals(SC.shp$F1.mean, n = 5, style = "quantile")$brks, title = "IRR", legend.format = list(fun = function(x) round(x, 2))) + tm_borders(col = "black") + tm_text("Pref.name", size = 0.7) + tm_legend(position = c("left","bottom")) + tm_layout(title = expression(bold("(B) F1 IRR by prefecture")), frame = FALSE, inner.margins = c(0, 0, 0.12, 0), title.size = 1.2)
f4.2 <- grid.grab()  # need to adjust the plot size in Rstudio to make it appropriate for saving

# panel c meta-regression
f4.3 <- ggplot() + 
  geom_ribbon(data = data.frame(preds.TB.F1), aes(x = x, ymin = exp(cr.lb), ymax = exp(cr.ub)), fill = "#74add1", alpha = 0.2) + 
  geom_point(data = TB.RR.agg, aes(x = CSSI, y = F1.mean, size = size1), col = "#313695") + 
  geom_line(data = data.frame(preds.TB.F1), aes(x = x, y = exp(pred), col = "PTB"), size = 1.2) +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0), legend.position = c(0.1, 0.9), plot.margin = unit(c(0.4,0.1,0,0.1), "cm")) + 
  geom_hline(yintercept = 1) + 
  xlab("CSSI") + 
  ylab("F1 IRR") + 
  guides(size = FALSE, col = FALSE) +
  ggtitle("(C) F1 meta regression") + 
  scale_color_manual(name = "", values = c("PTB" = "#313695","STIs" = "#a50026"))


f4 <- plot_grid(f4.1, f4.2, f4.3, rel_widths = c(0.15, 0.25, 0.25), nrow = 1)
f4
# save_plot(".\\Results\\Figure 4 F1 effect.pdf", Figure6, base_height = 4.07, base_width = 11)








#######################################################################
#                       Supplementary figures
#######################################################################

#===========================================================
#    Figure S2 - heatmap PTB ratio
#===========================================================
case.TB.Sichuan$incidence.scale <- ave(case.TB.Sichuan$incidence, paste(case.TB.Sichuan$Age, case.TB.Sichuan$Sex), FUN = function(x) x/x[1])


# male
qn = quantile(case.TB.Sichuan$incidence.scale[case.TB.Sichuan$Sex == 1], seq(0.1,0.9,length.out = 18), na.rm = TRUE)
qn01 <- sort(rescale(c(qn, range(case.TB.Sichuan$incidence.scale[case.TB.Sichuan$Sex == 1]))) )
gg.inc.male <- ggplot(case.TB.Sichuan[case.TB.Sichuan$Sex == 1,], aes(DiagnoseYear, Age)) + 
  geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
  geom_tile(aes(fill = incidence.scale), color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) + 
  annotate("text", x = 2017.5, y = 17.7, label = "F1") + 
  annotate("text", x = 2017.5, y = 10.7, label = "F2") + 
  geom_abline(intercept = -662.6667, slope = 1/3, size = 1) + 
  geom_abline(intercept = -655.6667, slope = 1/3, size = 1) + 
  labs(x = "Year of diagnosis", y = "Age group", fill = "IRR to 2005") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position="bottom", legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 10)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
  scale_y_discrete(expand = c(0, 0))


# female
qn = quantile(case.TB.Sichuan$incidence[case.TB.Sichuan$Sex == 0], seq(0.1,0.9,length.out = 18), na.rm = TRUE)
qn01 <- sort(rescale(c(qn, range(case.TB.Sichuan$incidence[case.TB.Sichuan$Sex == 0]))) )
gg.inc.female <- ggplot(case.TB.Sichuan[case.TB.Sichuan$Sex == 0,], aes(DiagnoseYear, Age)) + 
  geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
  geom_tile(aes(fill = incidence.scale), color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) + 
  annotate("text", x = 2017.5, y = 17.7, label = "F1") + 
  annotate("text", x = 2017.5, y = 10.7, label = "F2") + 
  geom_abline(intercept = -662.6667, slope = 1/3, size = 1) + 
  geom_abline(intercept = -655.6667, slope = 1/3, size = 1) + 
  labs(x = "Year of diagnosis", y = "Age group", fill = "IRR to 2005") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position="bottom", legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 10)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
  scale_y_discrete(expand = c(0, 0))

Lexis <- plot_grid(gg.inc.male, gg.inc.female, labels = c("(A) Male","(B) Female"), nrow = 1)
Lexis
save_plot(".\\Results\\Figure S2 PTB heatmap ratio.pdf", Lexis, base_asp = 1.12, base_height = 7.2)



#===========================================================
#    Figure S3 - prefecture level heatmap PTB
#===========================================================
pref.name <- SC.shp@data[, c("Pref", "Pref.name", "CSSI")]

pref.heatmap <- pref.name[order(pref.name$CSSI, decreasing = TRUE),] %>% filter(! Pref %in% c(5101,5108,5116,5120))


pref.hm <- list()
for(i in 1:nrow(pref.heatmap))
{
  current.case <- case.TB.nosex[case.TB.nosex$Pref == pref.heatmap$Pref[i],]
  qn = quantile(current.case$incidence, seq(0.01,0.99,length.out = 18), na.rm = TRUE)
  qn01 <- sort(rescale(c(qn, range(current.case$incidence))) )
  
  #  famine.polygons.F2 <- data.frame(ID = rep(1:5, each = 5), x = c(c(2004.5,2007.5, 2007.5, 2004.5, 2004.5) + rep(seq(0,by = 3, length.out = 4), each = 5), 2016.5,2018.5,2018.5,2016.5,2016.5), y = c(((2006-pref.thres$high.cohort[i])-10)/3+0.5, ((2006-pref.thres$high.cohort[i])-10)/3+0.5, ((2006-pref.thres$low.cohort[i])-10)/3+1.5, ((2006-pref.thres$low.cohort[i])-10)/3+1.5, ((2006-pref.thres$high.cohort[i])-10)/3+0.5) + rep(seq(from = 0, by = 1, length.out = 5), each = 5))
  famine.polygons.F2 <- data.frame(ID = rep(1:5, each = 5), x = c(c(2004.5,2007.5, 2007.5, 2004.5, 2004.5) + rep(seq(0,by = 3, length.out = 4), each = 5), 2016.5,2018.5,2018.5,2016.5,2016.5), y = c(668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3+1, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3+1, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3) + rep(seq(from = 0, by = 1, length.out = 5), each = 5)-3)
  
  pref.hm[[i]] <- ggplot(current.case, aes(DiagnoseYear, Age)) + 
    geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
    geom_tile(aes(fill = incidence), color = NA, alpha = 0.9) + 
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) +
    annotate("text", x = 2017.5, y = 18, label = "F1") + 
    annotate("text", x = 2017.5, y = 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3 + 2.5, label = "F2") + 
    geom_polygon(data = famine.polygons.F1, aes(x = x, y= y, group = ID), fill = NA, col = "black") + 
    geom_polygon(data = famine.polygons.F2, aes(x = x, y= y, group = ID), fill = NA, col = "black") +
    labs(x = "Year of diagnosis", y = "Age group", fill = "Incidence rate\n(1/100,000)") + 
    coord_fixed(ratio = 1) + 
    theme(legend.position="bottom", legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 8, angle = 45), legend.margin = margin(c(0,0,0,0)), legend.title = element_text(size = 10)) + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
    scale_y_discrete(expand = c(0, 0)) +
    ggtitle(pref.heatmap$Pref.name[i])
}

fs3 <- plot_grid(plotlist = pref.hm, ncol = 6)
# save_plot(".\\Results\\Figure S3 prefecture heatmap PTB.pdf", fs3, base_height = 15, base_aspect_ratio = 1.25)


#===========================================================
#    Figure S4 - prefecture level age effect PTB
#===========================================================
pref.heatmap <- SC.shp@data[order(SC.shp$CSSI, decreasing = TRUE),] %>% filter(! Pref %in% c(5101,5108,5116,5120))

fs4 <- pref.TB %>% filter(!Pref %in% c(5101,5108,5116,5120)) %>% left_join(SC.shp@data[,c("Pref","Pref.name")]) %>% mutate(Pref.name = factor(Pref.name, levels = pref.heatmap$Pref.name)) %>% filter(Type == "Age") %>% mutate(X = as.numeric(X)) %>% 
  ggplot(aes(x = X, y = parameter)) + 
  geom_line(col = "#e31a1c") + 
  facet_wrap(~Pref.name, ncol = 6)+
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd), col = NA, alpha = 0.5, fill = "#fb9a99") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Age at diagnosis", y = expression(alpha[i]), col = "", fill = "") +
  xlim(10,76)
fs4
# ggsave(".\\Results\\Figure S4 prefecture age effect PTB.pdf", fs4, width = 8, height = 5)

#===========================================================
#    Figure S5 - prefecture level period effect PTB
#===========================================================
fs5 <- pref.TB %>% filter(!Pref %in% c(5101,5108,5116,5120)) %>% left_join(SC.shp@data[,c("Pref","Pref.name")]) %>% mutate(Pref.name = factor(Pref.name, levels = pref.heatmap$Pref.name)) %>% filter(Type == "Period") %>% mutate(X = as.numeric(X)) %>% 
  ggplot(aes(x = X, y = parameter)) + 
  geom_line(col = "#33a02c") + 
  facet_wrap(~Pref.name, ncol = 6)+
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd), col = NA, alpha = 0.5, fill = "#b2df8a") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Year of diagnosis", y = expression(pi[j]), col = "", fill = "")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fs5

# ggsave(".\\Results\\Figure S5 prefecture period effect PTB.pdf", fs5, width = 8, height = 5)


#===========================================================
#    Figure S6 - prefecture level cohort effect PTB
#===========================================================
TB.pref.smooth <- pref.cohort.smooth(pref.TB, F2.start = pref.thres$low.cohort, F2.end = pref.thres$high.cohort) %>% filter(Cohort >= 1933 & Cohort <= 2002) %>% left_join(SC.shp@data[, c("Pref","Pref.name")])
TB.pref.smooth$Pref.name <- factor(TB.pref.smooth$Pref.name, levels = pref.heatmap$Pref.name)

fs6 <- pref.TB %>% filter(!Pref %in% c(5101,5108,5116,5120)) %>% left_join(SC.shp@data[,c("Pref","Pref.name")]) %>% mutate(Pref.name = factor(Pref.name, levels = pref.heatmap$Pref.name)) %>% filter(Type == "Cohort") %>% filter(X >= 1933 & X <= 2002) %>% mutate(X = as.numeric(X)) %>% 
  ggplot(aes(x = X, y = parameter)) + 
  geom_line(col = "#1f78b4")+
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd), col = NA, alpha = 0.5, fill = "#a6cee3") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  #   geom_vline(xintercept = 1960) + 
  geom_line(data = TB.pref.smooth[!TB.pref.smooth$Pref %in% c(5101,5108,5116,5120),], aes(x = Cohort, y = cohort.predict), linetype = "dashed") + 
  facet_wrap(~Pref.name, ncol = 6) + 
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Year of birth", y = expression(gamma[k]), col = "", fill = "")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fs6
# ggsave(".\\Results\\Figure S6 prefecture cohort effect PTB.pdf", fs6, width = 8, height = 5)


#===========================================================
#    Figure S7 - STBBI incidence rate map
#===========================================================
STI.inc.avg <- disease.dat %>% group_by(Pref, Pref.name) %>% summarise(Population = sum(Population), CaseCount = sum(CaseCount.STBBI)) %>% mutate(incidence = CaseCount/Population*100000) %>% ungroup() %>% mutate(Pref = as.numeric(Pref))

SC.shp@data <- left_join(SC.shp@data, STI.inc.avg)

fs7 <- tm_shape(SC.shp) + tm_fill(c("incidence"), breaks = classIntervals(SC.shp$incidence, n = 5, style = "jenks")$brks, title = "Annual mean STBBI\nincidence rate (1/100,000)") + tm_borders("black") + tm_text("Pref.name", size = 0.7) + tm_legend(position = c("left","bottom"), title = "") + tm_layout(inner.margins = c(0.02, 0.02, 0.02, 0.02), legend.title.size = 0.8)

# tmap_save(fs7, ".\\Results\\Figure S7 STBBI map.pdf", width = 5.5, height = 5, units = "in")

#===========================================================
#    Figure S8 - heatmap STBBI
#===========================================================
case.STI.Sichuan <- disease.dat %>% group_by(Age, DiagnoseYear, Sex) %>% dplyr::summarise(CaseCount = sum(CaseCount.STBBI)) %>% left_join(SC.pop) %>% ungroup() %>% mutate(incidence = CaseCount/Population*100000)

qn = quantile(case.STI.Sichuan$incidence[case.STI.Sichuan$Sex == 1], seq(0.1,0.9,length.out = 18), na.rm = TRUE)
qn01 <- sort(rescale(c(qn, range(case.STI.Sichuan$incidence[case.STI.Sichuan$Sex == 1]))) )
gg.inc.male <- ggplot(case.STI.Sichuan[case.STI.Sichuan$Sex == 1, ], aes(DiagnoseYear, Age)) + 
  geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
  geom_tile(aes(fill = incidence), color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) + 
  annotate("text", x = 2017.5, y = 18, label = "F1") + 
  annotate("text", x = 2017.5, y = 11, label = "F2") +
  geom_polygon(data = famine.polygons.F1, aes(x = x, y= y, group = ID), fill = NA, col = "black") + 
  geom_polygon(data = famine.polygons.F2, aes(x = x, y= y, group = ID), fill = NA, col = "black") +
  labs(x = "Year of diagnosis", y = "Age group", fill = "Incidence\nrate (1/100,000)") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position="bottom", legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 10), legend.margin = margin(c(0,0,0,0))) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
  scale_y_discrete(expand = c(0, 0))

# female
qn = quantile(case.STI.Sichuan$incidence[case.STI.Sichuan$Sex == 0], seq(0.1,0.9,length.out = 18), na.rm = TRUE)
qn01 <- sort(rescale(c(qn, range(case.STI.Sichuan$incidence[case.STI.Sichuan$Sex == 0]))) )
gg.inc.female <- ggplot(case.STI.Sichuan[case.STI.Sichuan$Sex == 0,], aes(DiagnoseYear, Age)) + 
  geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
  geom_tile(aes(fill = incidence), color = NA, alpha = 0.9) + 
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) + 
  annotate("text", x = 2017.5, y = 18, label = "F1") + 
  annotate("text", x = 2017.5, y = 11, label = "F2") +
  geom_polygon(data = famine.polygons.F1, aes(x = x, y= y, group = ID), fill = NA, col = "black") + 
  geom_polygon(data = famine.polygons.F2, aes(x = x, y= y, group = ID), fill = NA, col = "black") +
  labs(x = "Year of diagnosis", y = "Age group", fill = "Incidence\nrate (1/100,000)") + 
  coord_fixed(ratio = 1) + 
  theme(legend.position="bottom", legend.key.width = unit(0.8, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 10), legend.margin = margin(c(0,0,0,0))) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
  scale_y_discrete(expand = c(0, 0))

Lexis <- plot_grid(gg.inc.male, gg.inc.female, labels = c("(A) Male","(B) Female"), nrow = 1)
Lexis
# save_plot(".\\Results\\Figure S8 STBBI heatmap.pdf", Lexis, base_asp = 1.12, base_height = 7.2)

#===========================================================
#    Figure S9 - prefecture level heatmap PTB
#===========================================================
pref.hm <- list()
for(i in 1:nrow(pref.heatmap))
{
  current.case <- case.STI.nosex[case.STI.nosex$Pref == pref.heatmap$Pref[i],]
  qn = quantile(current.case$incidence, seq(0.01,0.99,length.out = 18), na.rm = TRUE)
  qn01 <- sort(rescale(c(qn, range(current.case$incidence))) )
  
  famine.polygons.F2 <- data.frame(ID = rep(1:5, each = 5), x = c(c(2004.5,2007.5, 2007.5, 2004.5, 2004.5) + rep(seq(0,by = 3, length.out = 4), each = 5), 2016.5,2018.5,2018.5,2016.5,2016.5), y = c(668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3+1, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3+1, 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3) + rep(seq(from = 0, by = 1, length.out = 5), each = 5)-3)
  
  pref.hm[[i]] <- ggplot(current.case, aes(DiagnoseYear, Age)) + 
    geom_abline(intercept = seq(-642.6667,-681.6667,-1), slope = 1/3, size = 0.5) + 
    geom_tile(aes(fill = incidence), color = NA, alpha = 0.9) + 
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, "RdYlBu")))(20),values = qn01) +
    annotate("text", x = 2017.5, y = 18, label = "F1") + 
    annotate("text", x = 2017.5, y = 668.8333-pref.thres$mid.cohort[pref.thres$pref == pref.heatmap$Pref[i]]/3 + 2.5, label = "F2") + 
    geom_polygon(data = famine.polygons.F1, aes(x = x, y= y, group = ID), fill = NA, col = "black") + 
    geom_polygon(data = famine.polygons.F2, aes(x = x, y= y, group = ID), fill = NA, col = "black") +
    labs(x = "Year of diagnosis", y = "Age group", fill = "Incidence rate\n(1/100,000)") + 
    coord_fixed(ratio = 1) + 
    theme(legend.position="bottom", legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.3,"cm"), legend.text = element_text(size = 8, angle = 45), legend.margin = margin(c(0,0,0,0)), legend.title = element_text(size = 10)) + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2017,3), limits = c(2004.5,2018.5)) +  
    scale_y_discrete(expand = c(0, 0)) +
    ggtitle(pref.heatmap$Pref.name[i])
}

fs9 <- plot_grid(plotlist = pref.hm, ncol = 6)
fs9

# save_plot(".\\Results\\Figure S9 prefecture heatmap STBBI.pdf", fs9, base_height = 15, base_aspect_ratio = 1.25)

#===========================================================
#    Figure S10 - APC results STBBI
#===========================================================
inc.table.male <- case.STI.Sichuan %>% filter(Sex == 1) %>% dplyr::select(Age, DiagnoseYear, incidence) %>% tidyr::spread(DiagnoseYear, incidence) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,79,3)) %>% column_to_rownames(var = "Age")
inc.table.female <- case.STI.Sichuan %>% filter(Sex == 0) %>% dplyr::select(Age, DiagnoseYear, incidence) %>% tidyr::spread(DiagnoseYear, incidence) %>% dplyr::select(-1) %>% mutate(Age =  seq(10,79,3)) %>% column_to_rownames(var = "Age")

# APC model results
male.apc <- apcglmkfit(r = inc.table.male, n.risk = male.pop, header=TRUE, p0 = 3, k = 3, fam = "qlik", amin = 1, pmin = 2005, cmin = 1927, agapyr = 3, pgapyr = 1) 
female.apc <- apcglmkfit(r = inc.table.female, n.risk = female.pop, header=TRUE, p0 = 3, k = 3, fam = "qlik", amin = 1, pmin = 2005, cmin = 1927, agapyr = 3, pgapyr = 1)

male.apc.result <- male.apc %>% reshape.APC.result() %>% mutate(Sex = "Male")
female.apc.result <- female.apc %>% reshape.APC.result() %>% mutate(Sex = "Female")
result.mf.STI <- rbind(male.apc.result, female.apc.result) %>% mutate(X = as.numeric(X))
result.mf.STI$Sex <- factor(result.mf.STI$Sex, levels = c("Male", "Female"))

result.mf.STI.vcov <- list(male.cov = list(a.vcov = male.apc$a.vcov, p.vcov = male.apc$p.vcov, c.vcov = male.apc$c.vcov), female.cov = list(a.vcov = female.apc$a.vcov, p.vcov = female.apc$p.vcov, c.vcov = female.apc$c.vcov))

# STBBI age effect
fs10.1 <- result.mf.STI %>% 
  filter(Type == "Age") %>% 
  ggplot(aes(x = X, y = parameter, col = Sex)) + 
  geom_line() +
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd, fill = Sex), col = NA, alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("(A) Age effect") + 
  theme(plot.title = element_text(hjust = 0)) +
  guides(col = FALSE, fill = FALSE) + 
  scale_color_manual(values = rev(gg_color_hue(2)))+
  scale_fill_manual(values = rev(gg_color_hue(2)))+ 
  labs(x = "Age at diagnosis", y = expression(alpha[i]), col = "", fill = "")

# period effect
fs10.2 <- result.mf.STI %>% 
  filter(Type == "Period") %>% 
  ggplot(aes(x = X, y = parameter, col = Sex)) + 
  geom_line() +
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd, fill = Sex), col = NA, alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ggtitle("(B) Period effect") + 
  theme(plot.title = element_text(hjust = 0)) +
  guides(col = FALSE, fill = FALSE) +   
  scale_color_manual(values = rev(gg_color_hue(2)))+
  scale_fill_manual(values = rev(gg_color_hue(2)))+ 
  labs(x = "Year of diagnosis", y = expression(pi[j]), col = "", fill = "")

# cohort effect
STI.smooth <- province.cohort.smooth(result.mf.STI, F1.start = 1957, F1.end = 1963, F2.start = 1975, F2.end = 1987, method = "GAM")

line_types <- c("Estimated"=1,"Expected"=2)

fs10.3 <- result.mf.STI %>% 
  filter(Type == "Cohort" & X >= 1933 & X <= 2000) %>% 
  ggplot(aes(x = X, y = parameter, col = Sex)) + 
  geom_line(aes(linetype="Estimated")) +
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd, fill = Sex), col = NA, alpha = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_line(data = STI.smooth[STI.smooth$Cohort >= 1933 & STI.smooth$Cohort <= 2000,], aes(x = Cohort, y = cohort.predict, linetype="Expected")) +
  geom_point(data = result.mf.STI[result.mf.STI$Type == "Cohort" & result.mf.STI$X == 1960,], aes(x = X, y = parameter, col = Sex), size = 2) + 
  geom_point(data = result.mf.STI[result.mf.STI$Type == "Cohort" & result.mf.STI$X == 1981,], aes(x = X, y = parameter, col = Sex), size = 2) + 
  geom_vline(xintercept = 1960) +
  geom_vline(xintercept = 1981) +
  ggtitle("(C) Cohort effect") + 
  theme(plot.title = element_text(hjust = 0), legend.position = c(0.8, 0.9), legend.margin = unit(0, "cm")) +
  labs(x = "Year of birth", y = expression(gamma[k]), col = "", fill = "") + 
  annotate("text", 1960, 1.2, label = "F1") +
  annotate("text", 1981, 1.2, label = "F2") + 
  ylim(-0.9,1.2)+
  scale_color_manual(values = rev(gg_color_hue(2)))+
  scale_fill_manual(values = rev(gg_color_hue(2)))+ 
  scale_linetype_manual(name = "", values=line_types)

left <- plot_grid(fs10.1, fs10.2, nrow = 2)
plot_grid(left, fs10.3, nrow = 1, rel_widths = c(1,2))

# save_plot(".\\Results\\Figure S10 APC results STBBI.pdf", plot_grid(left, f3.3, nrow = 1, rel_widths = c(1,2)), base_height  = 6.5, base_aspect_ratio = 1.5)




#===========================================================
#    Figure S11 - prefecture level age effect STBBI
#===========================================================
fs11 <- pref.STI %>% filter(!Pref %in% c(5101,5108,5116,5120)) %>% left_join(SC.shp@data[,c("Pref","Pref.name")]) %>% mutate(Pref.name = factor(Pref.name, levels = pref.heatmap$Pref.name)) %>% filter(Type == "Age") %>% mutate(X = as.numeric(X)) %>% 
  ggplot(aes(x = X, y = parameter)) + 
  geom_line(col = "#e31a1c") + 
  facet_wrap(~Pref.name, ncol = 6)+
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd), col = NA, alpha = 0.5, fill = "#fb9a99") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Age at diagnosis", y = expression(alpha[i]), col = "", fill = "") +
  xlim(10,76)
fs11
# ggsave(".\\Results\\Figure S11 prefecture age effect STBBI.pdf", fs11, width = 8, height = 5)



#===========================================================
#    Figure S12 - prefecture level period effect STBBI
#===========================================================
fs12 <- pref.STI %>% filter(!Pref %in% c(5101,5108,5116,5120)) %>% left_join(SC.shp@data[,c("Pref","Pref.name")]) %>% mutate(Pref.name = factor(Pref.name, levels = pref.heatmap$Pref.name)) %>% filter(Type == "Period") %>% mutate(X = as.numeric(X)) %>% 
  ggplot(aes(x = X, y = parameter)) + 
  geom_line(col = "#33a02c") + 
  facet_wrap(~Pref.name, ncol = 6)+
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd), col = NA, alpha = 0.5, fill = "#b2df8a") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Year of diagnosis", y = expression(pi[j]), col = "", fill = "")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fs12
# ggsave(".\\Results\\Figure S12 prefecture period effect STBBI.pdf", fs12, width = 8, height = 5)


#===========================================================
#    Figure S13 - prefecture level cohort effect STBBI
#===========================================================
STI.pref.smooth <- pref.cohort.smooth(pref.STI, F2.start = pref.thres$low.cohort, F2.end = pref.thres$high.cohort) %>% filter(Cohort >= 1933 & Cohort <= 2002) %>%  left_join(SC.shp@data[, c("Pref","Pref.name")])
STI.pref.smooth$Pref.name <- factor(STI.pref.smooth$Pref.name, levels = pref.heatmap$Pref.name)

fs13 <- pref.STI %>% filter(!Pref %in% c(5101,5108,5116,5120)) %>% left_join(SC.shp@data[,c("Pref","Pref.name")]) %>% mutate(Pref.name = factor(Pref.name, levels = pref.heatmap$Pref.name))%>% filter(Type == "Cohort") %>% filter(X >= 1933 & X <= 2002) %>% mutate(X = as.numeric(X)) %>% 
  ggplot(aes(x = X, y = parameter)) + 
  geom_line(col = "#1f78b4")+
  geom_ribbon(aes(x = X, ymin = parameter - 1.96*sd, ymax = parameter + 1.96*sd), col = NA, alpha = 0.5, fill = "#a6cee3") + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_line(data = STI.pref.smooth[!STI.pref.smooth$Pref %in% c(5101,5108,5116,5120),], aes(x = Cohort, y = cohort.predict), linetype = "dashed") + 
  facet_wrap(~Pref.name, ncol = 6) + 
  guides(col = FALSE, fill = FALSE) + 
  labs(x = "Year of birth", y = expression(gamma[k]), col = "", fill = "")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fs13
# ggsave(".\\Results\\Figure S13 prefecture cohort effect STBBI.pdf", fs13, width = 8, height = 5)


#===========================================================
#    Figure S14 - famine effect on F1 STBBI
#===========================================================
Sichuan.pop.long <-  data.frame(Age = rep(rep(seq(10,79,3), each = 14),2), DiagnoseYear = rep(2005:2018, 24*2), Population = c(unlist(t(male.pop)), unlist(t(female.pop))), Sex = rep(c("Male","Female"), each = 14*24))

Sichuan.STI.case.long <- case.STI.Sichuan %>% dplyr::select(Age, DiagnoseYear, CaseCount, Sex) %>% mutate(Age = rep(seq(10,79,3), each = 14*2), Cohort = DiagnoseYear - Age, Sex = fct_recode(as.factor(Sex), "Male" = "1","Female" = "0")) 

STI.RR <- province.RR(result.mf.STI, result.mf.STI.vcov, F1.start = 1957, F1.end = 1963, F2.start = 1975, F2.end = 1987, F2 = 1981, n.rep = 1000, method = "GAM", Sichuan.pop.long, Sichuan.STI.case.long) 

STI.Sichuan.agg.long <- STI.RR %>% gather(key = "Generation", value = "RR", F1.RR:F2.RR) %>% mutate(Generation = fct_recode(Generation, "F1" = "F1.RR", "F2" = "F2.RR")) %>% group_by(Sex, Generation) %>% summarise(mean = median(RR), low = quantile(RR, 0.025), high = quantile(RR, 0.975)) %>% ungroup() %>% mutate(Disease = "STIs")

Sichuan.pop.STI.long <- case.STI.nosex %>% filter(Pref != 5123)  %>% mutate(DiagnoseYear = rep(c(rep(seq(2006,2015,3), each = 3), 2018,2018),504), Age = fct_recode(Age, "[69,Inf)" = "[69,72)", "[69,Inf)" = "[72,75)","[69,Inf)" = "[75,78)","[69,Inf)" = "[78,Inf)")) %>% group_by(Pref, Age, DiagnoseYear) %>% dplyr::summarise(Population = round(sum(Population),0)) %>% ungroup() %>% dplyr::select(Pref,Age, DiagnoseYear, Population) %>% mutate(Age = rep(rep(seq(10,70,3), each = 5), 21))

STI.RR.all <- pref.RR(pref.STI, pref.STI.vcov, F1.start = 1957, F1.end = 1963, F2.start = pref.thres$low.cohort, F2.end = pref.thres$high.cohort, F2 = pref.thres$mid.cohort, n.rep = 1000, method = "GAM", Sichuan.pop.STI.long, Sichuan.STI.case.long, collapse = TRUE) 

STI.RR.agg <- STI.RR.all %>% group_by(Pref) %>% summarise(F1.log.mean = mean(log(F1.RR)), F1.log.sd = sd(log(F1.RR)), F2.log.mean = mean(log(F2.RR)), F2.log.sd = sd(log(F2.RR)), F1.mean = mean(F1.RR), F2.mean = mean(F2.RR), F1.low = quantile(F1.RR, 0.025), F1.high = quantile(F1.RR, 0.975), F2.low = quantile(F2.RR, 0.025), F2.high = quantile(F2.RR, 0.975)) %>% left_join(famine.intensity[,c("Pref","CSSI")]) %>% left_join(SC.shp@data[, c("Pref","Pref.name")])

STI.RR.agg$F2.log.mean[STI.RR.agg$Pref %in% c(5101,5108,5116,5120)] <- NA
STI.RR.agg$F2.log.sd[STI.RR.agg$Pref %in% c(5101,5108,5116,5120)] <- NA
STI.RR.agg$F2.mean[STI.RR.agg$Pref %in% c(5101,5108,5116,5120)] <- NA

# panel A RR
STI.Sichuan.agg.long$Sex <- factor(STI.Sichuan.agg.long$Sex, levels = c("Male","Female"))
fs14.1 <- STI.Sichuan.agg.long %>% filter(Generation == "F1") %>% ggplot(aes(x = Sex, y = mean, ymin = low, ymax = high, col = Sex)) + 
  geom_errorbar(aes(width = .2)) + 
  geom_point() +
  ylim(1,1.55) + 
  scale_color_manual(name = "",values = rev(gg_color_hue(2))) + 
  ylab("IRR") + 
  geom_hline(yintercept = 1, linetype = 3) + 
  xlab("") + 
  theme(legend.position = c(0.4, 0.95), plot.title = element_text(hjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0.5,0.5,0,0.5), "cm")) +
  ggtitle("(A) F1 IRR by sex")

# Panel B IRR map
meta.reg.STI.F1 <- rma(F1.log.mean, sei = F1.log.sd, mods =  ~ CSSI, data = STI.RR.agg, method = "SJ")
# permutest(meta.reg.STI.F1)
preds.STI.F1 <- predict(meta.reg.STI.F1, newmods=seq(0.25,0.6,0.01))
preds.STI.F1 <- data.frame(preds.STI.F1)
preds.STI.F1$x <- seq(0.25,0.6,0.01)

STI.RR.agg$size1 <- (1/STI.RR.agg$F1.log.sd)/(1/max(STI.RR.agg$F1.log.sd))
STI.RR.agg$size2 <- (1/STI.RR.agg$F2.log.sd)/(1/max(STI.RR.agg$F2.log.sd, na.rm = TRUE))

SC.shp@data <- SC.shp@data[,!names(SC.shp@data) %in% c("F1.mean","F2.mean")]
SC.shp@data <- left_join(SC.shp@data, STI.RR.agg[,c("Pref","F1.mean", "F2.mean")])

tm_shape(SC.shp) + tm_fill(c("F1.mean"), breaks = classIntervals(SC.shp$F1.mean, n = 5, style = "quantile")$brks, title = "IRR", legend.format = list(fun = function(x) round(x, 2))) + tm_borders(col = "black") + tm_text("Pref.name", size = 0.7) + tm_legend(position = c("left","bottom")) + tm_layout(title = expression(bold("(B) F1 IRR by prefecture")), frame = FALSE, inner.margins = c(0, 0, 0.12, 0), title.size = 1.2)
fs14.2 <- grid.grab()


fs14.3 <- ggplot() + 
  geom_ribbon(data = data.frame(preds.STI.F1), aes(x = x, ymin = exp(cr.lb), ymax = exp(cr.ub)), fill = "#a50026", alpha = 0.2) + 
  geom_point(data = STI.RR.agg, aes(x = CSSI, y = F1.mean, size = size1), col = "#a50026") + 
  geom_line(data = data.frame(preds.STI.F1), aes(x = x, y = exp(pred)), col = "#a50026", size = 1.2) +
  theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0), legend.position = c(0.1, 0.9), plot.margin = unit(c(0.5,0.5,0,0.5), "cm")) + 
  geom_hline(yintercept = 1) + 
  xlab("CSSI (%)") + 
  ylab("F1 IRR") + 
  guides(size = FALSE, col = FALSE) +
  ggtitle("(C) F1 meta regression")

plot_grid(fs14.1, fs14.2, fs14.3, rel_widths = c(0.15, 0.25, 0.25), nrow = 1)

# save_plot(".\\Results\\Figure S14 F1 effect STBBI.pdf", plot_grid(fs14.1, fs14.2, fs14.3, rel_widths = c(0.15, 0.25, 0.25), nrow = 1), base_height = 4.07, base_width = 11)

#===========================================================
#    Figure S15 - famine effect on F2 STBBI
#===========================================================
fs15 <- tm_shape(SC.shp) + tm_fill(c("F2.mean"), breaks = classIntervals(SC.shp$F2.mean, n = 5, style = "quantile")$brks, title = "IRR", legend.format = list(fun = function(x) round(x, 2)), textNA = "Excluded") + tm_borders(col = "black") + tm_text("Pref.name", size = 0.7) + tm_legend(position = c("left","bottom")) + tm_layout(title = expression(bold("F2 IRR by prefecture")), frame = FALSE, inner.margins = c(0, 0, 0.12, 0), title.size = 1.2)
# tmap_save(fs15, ".\\Results\\Figure S15 F2 effect STBBI.pdf", width = 5, height = 4.5, units = "in")
