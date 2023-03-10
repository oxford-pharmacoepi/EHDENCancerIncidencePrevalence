library(here)
library(dplyr)
library(tidyr)
library(stringr)

# printing numbers with 1 decimal place and commas 
nice.num<-function(x) {
  trimws(format(round(x,1),
                big.mark=",", nsmall = 1, digits=1, scientific=FALSE))}
# printing numbers with 2 decimal place and commas 
nice.num2<-function(x) {
  trimws(format(round(x,2),
                big.mark=",", nsmall = 2, digits=2, scientific=FALSE))}
# printing numbers with 3 decimal place and commas 
nice.num3<-function(x) {
  trimws(format(round(x,3),
                big.mark=",", nsmall = 3, digits=3, scientific=FALSE))}
# for counts- without decimal place
nice.num.count<-function(x) {
  trimws(format(x,
                big.mark=",", nsmall = 0, digits=1, scientific=FALSE))}



#### Load data -----
#data
prevalence_estimates <- readRDS(here("data","prevalence_estimates.rds"))
prevalence_attrition <- readRDS(here("data","prevalence_attrition.rds"))
incidence_estimates <- readRDS(here("data","incidence_estimates.rds"))
incidence_attrition <- readRDS(here("data","incidence_attrition.rds"))
