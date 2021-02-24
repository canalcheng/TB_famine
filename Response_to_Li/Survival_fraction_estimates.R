#This code reads and processes data from the 1964, 1982, and 2000 national censuses of China in order to estimate fractions of people born in each year surviving to the year 2000

library(readr)
library(OrdMonReg)

Chinapop = read_csv('Data\\ChinaPopByAge.csv') # source: 1964 and 1982 census data: Women's Studies Institute of China,the All-china Women's Federation;Women's Studies Institute of Shaanxi Province, Statistics on Chinese Women (China Statistics Press, Beijing, 1991); 2000 census data: Census Office of the State Council; Population and Social Science and Technology Statistics Division of National Bureau of Statistics, Tabulation on the 2000 Population Census of the People's Republic of China (China Statistics Press, Beijing, 2002).

##The first step is reshaping the data such that each cohort is on the same row for the 1964, 1982, 2000 censuses

#Calculate birth years
Chinapop$birthyear1964 = 1964 - Chinapop$Age
Chinapop$birthyear2000 = 2000 - Chinapop$Age
Chinapop$birthyear1982 = 1982 - Chinapop$Age

#split to data specific to each census
pop1964 = Chinapop[,grepl('1964',names(Chinapop))]
pop2000 = Chinapop[,grepl('2000',names(Chinapop))]
pop1982 = Chinapop[,grepl('1982',names(Chinapop))]

#merge on birth year
names(pop1964)[3] = names(pop1982)[3] = names(pop2000)[3] = 'BirthYear'

Chinapop = merge(pop1964,pop1982,by = 'BirthYear',all = T)
Chinapop = merge(Chinapop,pop2000,by = 'BirthYear',all = T)

#subset to years where there is data in at least two censuses
Chinapop = Chinapop[rowSums(is.na(Chinapop))<3,]

#calculate fractions surviving
Chinapop$SurvFrac1964_1982 = (Chinapop$`1982ChinaFemale`+Chinapop$`1982ChinaMale`)/(Chinapop$`1964ChinaFemale`+Chinapop$`1964ChinaMale`)
Chinapop$SurvFrac1982_2000 = (Chinapop$`2000ChinaFemale`+Chinapop$`2000ChinaMale`)/(Chinapop$`1982ChinaFemale`+Chinapop$`1982ChinaMale`)

#plot
plot(Chinapop$BirthYear,Chinapop$SurvFrac1964_1982,lwd = 2,'l',col = rgb(0,0,0,.5),xlab = 'Birth Year',ylab = 'Fraction Surviving', ylim=range(na.omit(Chinapop$SurvFrac1982_2000)))
abline(1,0,lty=2)
lines(Chinapop$BirthYear,Chinapop$SurvFrac1982_2000,lwd = 2,col = rgb(0.5,0,0,.5),xlab = 'Birth Year',ylab = 'Fraction Surviving')
text(1900,1.1,'Survival 1964-1982',pos = 4,col = rgb(0,0,0,.5))
text(1930,0.6,'Survival 1982-2000',pos = 4,col = rgb(0.5,0,0,.5))

##fit bounded isotonic regression for 18 yr survival by age in order to smooth noise in data
isoreg64.82 = BoundedIsoMean(na.omit(Chinapop$SurvFrac1964_1982),
                             w=rep(1,sum(!is.na(Chinapop$SurvFrac1964_1982))),
                             a = rep(0,sum(!is.na(Chinapop$SurvFrac1964_1982))), 
                             b = rep(1,sum(!is.na(Chinapop$SurvFrac1964_1982))))

isoreg82.00 = BoundedIsoMean(na.omit(Chinapop$SurvFrac1982_2000),
                             w=rep(1,sum(!is.na(Chinapop$SurvFrac1982_2000))),
                             a = rep(0,sum(!is.na(Chinapop$SurvFrac1982_2000))), 
                             b = rep(1,sum(!is.na(Chinapop$SurvFrac1982_2000))))

Chinapop$SurvFrac1964_1982.fit = Chinapop$SurvFrac1982_2000.fit = NA

Chinapop$SurvFrac1964_1982.fit[!is.na(Chinapop$SurvFrac1964_1982)] = isoreg64.82
Chinapop$SurvFrac1982_2000.fit[!is.na(Chinapop$SurvFrac1982_2000)] = isoreg82.00

plot(Chinapop$Age1964,Chinapop$SurvFrac1964_1982.fit,lwd = 2,'l',col = rgb(0,0,0,.5),xlab = 'Initial age',ylab = 'Fraction Surviving', ylim=range(na.omit(Chinapop$SurvFrac1982_2000)),
     xlim = c(0,100))
abline(1,0,lty=2)
lines(Chinapop$Age1982,Chinapop$SurvFrac1982_2000.fit,lwd = 2,col = rgb(0.5,0,0,.5),xlab = 'Birth Year',ylab = 'Fraction Surviving')
text(0,0.6,'Survival 1964-1982',pos = 4,col = rgb(0,0,0,.5))
text(0,0.4,'Survival 1982-2000',pos = 4,col = rgb(0.5,0,0,.5))

Chinapop$Surv.overall <- ifelse(is.na(Chinapop$SurvFrac1964_1982.fit), Chinapop$SurvFrac1982_2000.fit, Chinapop$SurvFrac1964_1982.fit*Chinapop$SurvFrac1982_2000.fit)
Chinapop <- Chinapop[, c("BirthYear", "Surv.overall")]

Chinapop <- rbind(Chinapop, data.frame(BirthYear = 1983:2000, Surv.overall = 1))  # append the survival rate from 1983 to 2000

write_csv(Chinapop,'Data\\Cohort shrinkage estimates.csv')
