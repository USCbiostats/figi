library(table1)
library(tidyverse)

rm(list = ls())
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData") # newest Epi data subset

numeric_vars <- c("age_ref", "fiber", "calcium_tot", "fruit", "vegetable", "redmeat", "procmeat")
epidat <- Epi %>% 
  filter(drop == 0 & gxe == 1) %>% 
  mutate_at(numeric_vars, as.numeric) %>% 
  mutate(sex = factor(sex, labels = c('Female', 'Male')),
         alcoholc = fct_shift(factor(alcoholc, labels = c("N/A", ">28g/d", "1-28g/d", "nondrinker")), 1),
         smoke = fct_shift(factor(smoke, labels = c("N/A", "Former smoke", "Never smoker", "Smoker")),1),
         famhx1 = fct_shift(factor(famhx1, labels = c("N/A", "No", "Yes")), 1),
         folate_totqc2 = factor(folate_totqc2, labels = c("N/A", "Q1", "Q2", "Q3", "Q4")))

label(epidat$age_ref) <- "Age"
label(epidat$fiber) <- "Fiber Intake"
label(epidat$calcium_tot) <- "Total Calcium Intake"
label(epidat$fruit) <- "Fruit Intake"
label(epidat$vegetable) <- "Vegetable Intake"
label(epidat$redmeat) <- "Redmeat Intake"
label(epidat$procmeat) <- "Processed Meat Intake"
label(epidat$sex) <- "Sex"
label(epidat$alcoholc) <- "Alcohol Intake"
label(epidat$smoke) <- "Smoking"
label(epidat$famhx1) <- "Family History of CRC"

label(epidat$folate_totqc2) <- "Total Folate"



# render 
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), 
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), 
               function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}
table1(~ age_ref + sex + alcoholc + smoke + famhx1 + calcium_tot + fruit + vegetable + redmeat + procmeat | folate_totqc2 , 
       data=epidat,
       render.continuous=my.render.cont, render.categorical=my.render.cat)


#------ Comments ------

# physical activity - LOTS missing
#