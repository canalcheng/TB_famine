# Cheng_et_al_2020_TB_famine

This repository contains data and code for reproducing the analyses in *Prenatal and early life exposure to the Great Chinese Famine increased the risk of pulmonary tuberculosis in adulthood across two generations*. 
 
## Data
### 1. Disease and population data.csv 

The population size and the reported case count of PTB and sexually transmitted blood borne infections (STBBI) by prefecture, age group, diagnosis year, and sex. 
 
* **Pref** represents prefecture-level [administrative division codes](https://en.wikipedia.org/wiki/Administrative_division_codes_of_the_People%27s_Republic_of_China). The first two digits "51" represent Sichuan province, China. The third and fourth digits identify the prefecture.
 
* **Pref.name** is the name of the prefecture.
 
* **Age** represents the age group.

* **DiagnoseYear** represents the year of diagnosis.

* **Sex** represents the sex, with 0 representing females and 1 representing males.

* **Population** represents the projected population size in this group.

* **CaseCount.TB** represents the number of reported active PTB cases.

* **CaseCount.STBBI** represents the number of reported STBBI cases.


### 2. Sichuan population projection.csv

The projected population size by age group, year, and sex.

* **Age** represents the age group.

* **DiagnoseYear** represents the year of diagnosis.

* **Sex** represents the sex, with 0 representing females and 1 representing males.

* **Population** represents the projected population size in this group.


### 3. Famine intensity by pref.csv

The cohort size shrinkage index (CSSI) calculated for each prefecture.

* **Pref.name** is the name of the prefecture.

* **CSSI** is the cohort size shrinkage index.

### 4. F2 time by prefecture population data.csv

The timing for F2 for different prefectures.

* **Pref** represents prefecture-level administrative division code.

* **mid.cohort** is the middle year of the defined F2 cohort group (Table S4).

* **low.cohort** and **high.cohort** represent the lower and higher threshold for the cohort groups to remove when estimating the counterfactural effects. Due the overlap in cohort windows, one cohort group before and one cohort group after the F2 cohort gorup were removed.

### 5. Sichuan_City.shp

The base map of Sichuan Province, China.




## Code
* **Figures.R** contains code to regenerate the figures
* **Functions.R** contains functions used in the analysis





