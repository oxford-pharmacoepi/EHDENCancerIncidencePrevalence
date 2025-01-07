library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(ggplot2)
library(scales)
library("cowplot")
library("gridExtra")
library("ggpubr")
library(ggh4x)
library(patchwork)
library(readr)
library(rio)
library(tidyverse)


# function from dsr package (no longer supported on R version >4.4) to carry out age standardization
dsr <- function(data, event, fu, subgroup, ..., refdata, mp, method="normal", sig=0.95, decimals ) {
  
  subgroup <- enquo(subgroup)
  event <- enquo(event)
  fu <- enquo(fu)
  stdvars <- quos(...)
  
  
  
  #Compute crude and standardized rates and variances
  all_data_st = data %>%
    left_join(refdata) %>%
    group_by(!!subgroup) %>%
    mutate(n=sum(!!event),
           d=sum(!!fu),
           cr_rate=n/d,
           cr_var=n/d^2,
           wts=pop/sum(pop),
           st_rate=sum(wts*(!!event/!!fu)),
           st_var=sum(as.numeric((wts^2)*(!!event/(!!fu )^2)))) %>%
    distinct(!!subgroup, .keep_all=TRUE) %>%
    select(!!subgroup, n, d, cr_rate, cr_var, st_rate, st_var)
  
  
  
  #Compute Confidence Intervals (CI) according to method. The default is 'normal'
  if(method=="gamma") {
    
    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*qgamma((1-sig)/2, shape=cr_rate^2/(cr_var))/(cr_rate/cr_var),
        c_upper=mp*qgamma(1-((1-sig)/2), shape=1+cr_rate^2/(cr_var))/(cr_rate/cr_var),
        s_rate=mp*st_rate,
        s_lower=mp*qgamma((1-sig)/2, shape=st_rate^2/st_var)/(st_rate/st_var),
        s_upper=mp*qgamma(1-((1-sig)/2), shape=1+(st_rate^2/st_var))/(st_rate/st_var)) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)
    
    
  } else if (method=="normal") {
    
    
    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*(cr_rate+qnorm((1-sig)/2)*sqrt(cr_var)),
        c_upper=mp*(cr_rate-qnorm((1-sig)/2)*sqrt(cr_var)),
        s_rate=mp*st_rate,
        s_lower=mp*(st_rate+qnorm((1-sig)/2)*sqrt(st_var)),
        s_upper=mp*(st_rate-qnorm((1-sig)/2)*sqrt(st_var))) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)
    
  } else if (method=="lognormal") {
    
    
    tmp1 =   all_data_st %>%
      mutate(
        c_rate=mp*cr_rate,
        c_lower=mp*exp((log(cr_rate)+qnorm((1-sig)/2)*sqrt(cr_var)/(cr_rate))),
        c_upper=mp*exp((log(cr_rate)-qnorm((1-sig)/2)*sqrt(cr_var)/(cr_rate))),
        s_rate=mp*st_rate,
        s_lower=mp*exp((log(st_rate)+qnorm((1-sig)/2)*sqrt(st_var)/(st_rate))),
        s_upper=mp*exp((log(st_rate)-qnorm((1-sig)/2)*sqrt(st_var)/(st_rate)))) %>%
      select(!!subgroup, n, d, c_rate, c_lower, c_upper, s_rate, s_lower, s_upper)
    
  }
  
  
  #Clean up and output
  
  tmp1$c_rate  <- round(tmp1$c_rate,  digits=decimals)
  tmp1$c_lower <- round(tmp1$c_lower, digits=decimals)
  tmp1$c_upper <- round(tmp1$c_upper, digits=decimals)
  tmp1$s_rate  <- round(tmp1$s_rate, digits=decimals)
  tmp1$s_lower <- round(tmp1$s_lower, digits=decimals)
  tmp1$s_upper <- round(tmp1$s_upper, digits=decimals)
  
  colnames(tmp1) <-  c('Subgroup', 'Numerator','Denominator',
                       paste('Crude Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (Crude)',sep=''),
                       paste(sig*100,'% UCL (Crude)',sep=''),
                       paste('Std Rate (per ',mp,')',sep=''),
                       paste(sig*100,'% LCL (Std)',sep=''),
                       paste(sig*100,'% UCL (Std)',sep=''))
  
  
  tmp1 <- as.data.frame(tmp1)
  
}

pathResults <- "C:/Users/dnewby/Desktop/"

datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/CancerIncidencePrevalanceShiny/shiny/data"

# read in files
incidence_estimates <- readRDS(paste0(datapath ,"/incidence_estimates.rds")) %>% 
  mutate(database_name = ifelse(database_name == "CPRDGoldUpdate2", "CPRD GOLD", database_name))


# get gold and aurum data overall

# any value NA replace with 5
incidence_estimates_overall <- incidence_estimates %>% 
  filter(outcome_cohort_name == "Colorectal",
         analysis_interval == "years",
         denominator_age_group != "All"
  ) %>% 
  mutate(n_events = ifelse(is.na(n_events), 5, n_events)) 


# perform age standardization

# european standard population 2013
ESP13path <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/ESP13.csv"
ESP13 <- read_csv(ESP13path) # european population standard

#collapse ESP13 and remove ages not used for study
ESP13_updated <- ESP13 %>% 
  filter(Agegroup != "0-4",
         Agegroup != "5-9",
         Agegroup != "10-14",
         Agegroup != "15-19" ) %>% 
  add_row(Agegroup = "18 to 29", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '20-24'| Agegroup == '25-29']))) %>% 
  add_row(Agegroup = "30 to 39", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '30-34'| Agegroup == '35-39']))) %>% 
  add_row(Agegroup = "40 to 49", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '40-44'| Agegroup == '45-49']))) %>% 
  add_row(Agegroup = "50 to 59", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '50-54'| Agegroup == '55-59']))) %>% 
  add_row(Agegroup = "60 to 69", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '60-64'| Agegroup == '65-69']))) %>% 
  add_row(Agegroup = "70 to 79", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '70-74'| Agegroup == '75-79']))) %>% 
  add_row(Agegroup = "80 to 89", ESP2013 = with(ESP13, sum(ESP2013[Agegroup == '80-84'| Agegroup == '85-89']))) %>% 
  filter(Agegroup == "18 to 29" | Agegroup == "30 to 39"| Agegroup == "40 to 49"| Agegroup == "50 to 59" | Agegroup == "60 to 69" |
           Agegroup == "70 to 79" |
           Agegroup == "80 to 89" |
           Agegroup == "90+" ) %>% 
  mutate(Agegroup = replace(Agegroup, Agegroup == "90+", "90 +")) 


# UK standard population 2021
mye21path <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/mye21.csv"
mye2021 <- read_csv(mye21path) 

#collapse mye2021 and remove ages not used for study
mye21_updated <- mye2021 %>% 
  filter(Agegroup != "0-4",
         Agegroup != "5-9",
         Agegroup != "10-14",
         Agegroup != "15-19" ) %>% 
  add_row(Agegroup = "18 to 29", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '20-24'| Agegroup == '25-29']))) %>% 
  add_row(Agegroup = "30 to 39", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '30-34'| Agegroup == '35-39']))) %>% 
  add_row(Agegroup = "40 to 49", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '40-44'| Agegroup == '45-49']))) %>% 
  add_row(Agegroup = "50 to 59", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '50-54'| Agegroup == '55-59']))) %>% 
  add_row(Agegroup = "60 to 69", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '60-64'| Agegroup == '65-69']))) %>% 
  add_row(Agegroup = "70 to 79", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '70-74'| Agegroup == '75-79']))) %>% 
  add_row(Agegroup = "80 to 89", mye21_weights = with(mye2021, sum(mye21_weights[Agegroup == '80-84'| Agegroup == '85-89']))) %>% 
  filter(Agegroup == "18 to 29" | Agegroup == "30 to 39"| Agegroup == "40 to 49"| Agegroup == "50 to 59" | Agegroup == "60 to 69" |
           Agegroup == "70 to 79" |
           Agegroup == "80 to 89" |
           Agegroup == "90+" ) %>% 
  mutate(Agegroup = replace(Agegroup, Agegroup == "90+", "90 +")) %>% 
  select(c(Agegroup,
           mye21_weights
           ))



# needs column to be called "pop"

ESP13_updated <- ESP13_updated %>% 
  rename(pop = ESP2013)

mye21_updated <- mye21_updated %>% 
  rename(pop = mye21_weights)

# age standardization


incidence_estimates_Aurum <- incidence_estimates_overall %>% 
  filter(database_name == "CPRD Aurum") %>% 
  rename(Agegroup = denominator_age_group) #rename this to make it the same as agegroup in standard population


incidence_estimates_gold <- incidence_estimates_overall %>% 
  filter(database_name == "CPRD GOLD") %>% 
  rename(Agegroup = denominator_age_group) #rename this to make it the same as agegroup in standard population


# AURUM ESP2013 #

#create a loop for each cancer (all cancers apart from prostate and breast are for single genders)
agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()

for(i in 1:length(table(incidence_estimates_Aurum$outcome_cohort_name))){
  
  
  # for whole population    
  incidence_estimates_i <- incidence_estimates_Aurum  %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
           denominator_sex == "Both")
  
  
  tryCatch({
    
    agestandardizedinc[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedinc[[i]] <- agestandardizedinc[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
             Sex = "Both") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]), "for both sexes done"))
  
  ##########################
  # for female population    
  #########################
  incidence_estimates_i <- incidence_estimates_Aurum %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
           denominator_sex == "Female")
  
  
  tryCatch({
    
    agestandardizedincF[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincF[[i]] <- agestandardizedincF[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
             Sex = "Female") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]), "for females done"))
  
  
  
  
  ##########################
  # for male population    
  #########################
  incidence_estimates_i <- incidence_estimates_Aurum %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
           denominator_sex == "Male")
  
  
  tryCatch({
    
    agestandardizedincM[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincM[[i]] <- agestandardizedincM[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
             Sex = "Male") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]), "for males done"))
  
  
  
  
  
}

agestandardizedinc_final <- bind_rows(agestandardizedinc,
                                      agestandardizedincF,
                                      agestandardizedincM) %>% 
  mutate(Database = "CPRD Aurum",
         Pop_std = "ESP2013")



# Aurum UK std population

#create a loop for each cancer (all cancers apart from prostate and breast are for single genders)
agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()

for(i in 1:length(table(incidence_estimates_Aurum$outcome_cohort_name))){
  
  
  # for whole population    
  incidence_estimates_i <- incidence_estimates_Aurum  %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
           denominator_sex == "Both")
  
  
  tryCatch({
    
    agestandardizedinc[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = mye21_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedinc[[i]] <- agestandardizedinc[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
             Sex = "Both") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]), "for both sexes done"))
  
  ##########################
  # for female population    
  #########################
  incidence_estimates_i <- incidence_estimates_Aurum %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
           denominator_sex == "Female")
  
  
  tryCatch({
    
    agestandardizedincF[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = mye21_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincF[[i]] <- agestandardizedincF[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
             Sex = "Female") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]), "for females done"))
  
  
  
  
  ##########################
  # for male population    
  #########################
  incidence_estimates_i <- incidence_estimates_Aurum %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
           denominator_sex == "Male")
  
  
  tryCatch({
    
    agestandardizedincM[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = mye21_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincM[[i]] <- agestandardizedincM[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]),
             Sex = "Male") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_Aurum$outcome_cohort_name)[i]), "for males done"))
  
  
  
  
  
}

agestandardizedinc_final1 <- bind_rows(agestandardizedinc,
                                      agestandardizedincF,
                                      agestandardizedincM) %>% 
  mutate(Database = "CPRD Aurum",
         Pop_std = "UK2021")



# gold ESP2013

#create a loop for each cancer (all cancers apart from prostate and breast are for single genders)
agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()

for(i in 1:length(table(incidence_estimates_gold$outcome_cohort_name))){
  
  
  # for whole population    
  incidence_estimates_i <- incidence_estimates_gold  %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
           denominator_sex == "Both")
  
  
  tryCatch({
    
    agestandardizedinc[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedinc[[i]] <- agestandardizedinc[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
             Sex = "Both") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_gold$outcome_cohort_name)[i]), "for both sexes done"))
  
  ##########################
  # for female population    
  #########################
  incidence_estimates_i <- incidence_estimates_gold %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
           denominator_sex == "Female")
  
  
  tryCatch({
    
    agestandardizedincF[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincF[[i]] <- agestandardizedincF[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
             Sex = "Female") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_gold$outcome_cohort_name)[i]), "for females done"))
  
  
  
  
  ##########################
  # for male population    
  #########################
  incidence_estimates_i <- incidence_estimates_gold %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
           denominator_sex == "Male")
  
  
  tryCatch({
    
    agestandardizedincM[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincM[[i]] <- agestandardizedincM[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
             Sex = "Male") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_gold$outcome_cohort_name)[i]), "for males done"))
  
  
  
  
  
}

agestandardizedinc_final2 <- bind_rows(agestandardizedinc,
                                      agestandardizedincF,
                                      agestandardizedincM) %>% 
  mutate(Database = "CPRD GOLD",
         Pop_std = "ESP2013")

# gold UK standard population

#create a loop for each cancer (all cancers apart from prostate and breast are for single genders)
agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()

for(i in 1:length(table(incidence_estimates_gold$outcome_cohort_name))){
  
  
  # for whole population    
  incidence_estimates_i <- incidence_estimates_gold  %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
           denominator_sex == "Both")
  
  
  tryCatch({
    
    agestandardizedinc[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = mye21_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedinc[[i]] <- agestandardizedinc[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
             Sex = "Both") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_gold$outcome_cohort_name)[i]), "for both sexes done"))
  
  ##########################
  # for female population    
  #########################
  incidence_estimates_i <- incidence_estimates_gold %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
           denominator_sex == "Female")
  
  
  tryCatch({
    
    agestandardizedincF[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = mye21_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincF[[i]] <- agestandardizedincF[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
             Sex = "Female") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_gold$outcome_cohort_name)[i]), "for females done"))
  
  
  
  
  ##########################
  # for male population    
  #########################
  incidence_estimates_i <- incidence_estimates_gold %>%
    filter(outcome_cohort_name == names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
           denominator_sex == "Male")
  
  
  tryCatch({
    
    agestandardizedincM[[i]] <- dsr(
      data = incidence_estimates_i,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = mye21_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincM[[i]] <- agestandardizedincM[[i]] %>% 
      mutate(Cancer = names(table(incidence_estimates_gold$outcome_cohort_name)[i]),
             Sex = "Male") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates_gold$outcome_cohort_name)[i]), "for males done"))
  
  
  
  
  
}

agestandardizedinc_final3 <- bind_rows(agestandardizedinc,
                                       agestandardizedincF,
                                       agestandardizedincM) %>% 
  mutate(Database = "CPRD GOLD",
         Pop_std = "UK2021")


# merge all the results together
agestd_all <- bind_rows(
  
  agestandardizedinc_final,
  agestandardizedinc_final1,
  agestandardizedinc_final2,
  agestandardizedinc_final3
  
)



























