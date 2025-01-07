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

pathResults <- "C:/Users/dnewby/Desktop/crc_figures"

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






# making plots

# age standardized crc by sex ESP2013

incidenceFigureData <- agestd_all %>%
  filter(Pop_std == "ESP2013" ) %>%
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Database)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Database, fill = Database),size = 3) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.background = element_blank() ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        strip.background = element_rect(colour = "black", size = 0.5),  # Add a black box around facet labels
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years (ESP 2013)",
       col = "Database",
       shape = "Database",
       fill = "Database" ) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"),
               expand = c(0.06,1)) +
  facet_wrap(~Sex)


plotname <- paste0("FIGURE1_Incidence_age_std_ESP2013_crc.pdf")

pdf(paste0(pathResults ,"/", plotname), width = 10, height = 6)

print(incidenceFigureData, newpage = FALSE)
dev.off()


# age standardized crc by sex UK2021

incidenceFigureData <- agestd_all %>%
  filter(Pop_std == "UK2021" ) %>%
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Database)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Database, fill = Database),size = 3) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.background = element_blank() ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        strip.background = element_rect(colour = "black", size = 0.5),  # Add a black box around facet labels
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years (UK 2021)",
       col = "Database",
       shape = "Database",
       fill = "Database" ) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"),
               expand = c(0.06,1)) +
  facet_wrap(~Sex)


plotname <- paste0("FIGURE2_Incidence_age_std_UK2021_crc.pdf")

pdf(paste0(pathResults ,"/", plotname), width = 10, height = 6)

print(incidenceFigureData, newpage = FALSE)
dev.off()

# crude
incidenceFigureData <- agestd_all %>%
  filter(Pop_std == "ESP2013" ) %>%
  ggplot(aes(x = Subgroup,
             y = `Crude Rate (per 1e+05)`,
             group = Database)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Crude)`, 
                  ymax = `95% UCL (Crude)`, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Database, fill = Database),size = 3) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.background = element_blank() ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        strip.background = element_rect(colour = "black", size = 0.5),  # Add a black box around facet labels
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Crude Incidence rate per 100000 person-years",
       col = "Database",
       shape = "Database",
       fill = "Database" ) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"),
               expand = c(0.06,1)) +
  facet_wrap(~Sex)


plotname <- paste0("FIGURE3_Incidence_crude_crc.pdf")

pdf(paste0(pathResults ,"/", plotname), width = 10, height = 6)

print(incidenceFigureData, newpage = FALSE)
dev.off()


# compared to the national cancer registries

# read in the NCRAS crude and age std results
datapath_std <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data/Incidence_data_for_England_2024.csv"
oxford_ncras <- read_csv(datapath_std) 

df_long <- agestd_all %>%
  pivot_longer(
    cols = c("Crude Rate (per 1e+05)", "95% LCL (Crude)", "95% UCL (Crude)", 
             "Std Rate (per 1e+05)", "95% LCL (Std)", "95% UCL (Std)"),
    values_to = "Value",
    names_to = "name"
  ) %>%
  mutate(
    type = case_when(
      str_detect(name, "Crude") ~ "Crude",
      str_detect(name, "Std") ~ "Age Std",
      TRUE ~ NA_character_
    ),
    name = name %>%
      str_remove_all("\\s*\\(?(Crude|Std)\\)?") %>%        # Remove "Crude", "Std", "(Crude)", "(Std)"
      str_remove_all("95%\\s*") %>%                      # Remove "95% CI"
      str_remove_all("\\s*\\(per 1e\\+05\\)\\s*") %>%      # Remove "(per 1e+05)"
      str_trim()                                            # Remove any trailing/leading white space
  )


df_wide <- df_long %>%
  pivot_wider(
    names_from = name,  # Column names come from the 'name' column
    values_from = Value  # Values come from the 'Value' column
  ) %>% 
  mutate(Database = paste(Database,type, Pop_std) ) %>% 
  rename(Year = Subgroup) %>% 
  dplyr::select(-Numerator, -Denominator)


#  [1] "Subgroup"    "Numerator"   "Denominator" "Cancer"      "Sex"         "type"        "Rate"        "LCL"         "UCL"        
# [10] "Database"   


# ncras
# [1] "Year"             "Gender"           "Age_at_Diagnosis" "Geography_code"   "Geography_name"   "ICD10_code"       "Site_description"
# [8] "Count"            "Type_of_rate"     "Rate"             "LCI"              "UCI"              "Flag"    

ncras <- oxford_ncras %>% 
  rename(Sex = Gender,
         LCL = LCI ,
         UCL = UCI,
         type =  Type_of_rate,
         Cancer = Site_description) %>% 
  mutate(
    type = recode(type, 
                  "Non-standardised" = "Crude", 
                  "Age-standardised" = "Age Std"),
    
    Sex = recode(Sex, 
                 "Persons" = "Both")
    
  ) %>% 
  mutate(Database = paste("NCRAS ",type) ) %>% 
  dplyr::select(-c(
    Age_at_Diagnosis,
    Geography_code,   
    Geography_name, 
    ICD10_code,
    Count,
    Flag
    
  )) %>% 
  mutate(
    Year = as.Date(paste(Year, "01", "01", sep = "-"))  # Converts Year to a date format (YYYY-01-01)
  )


final_comb <- bind_rows(df_wide, ncras) %>% 
  mutate(Cancer = str_replace_all(Cancer, "Malignant neoplasm of colon and rectum", "Colorectal"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of bladder", "Bladder"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of breast", "Breast"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of bronchus and lung", "Lung"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of hypopharynx", "Hypopharynx"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of larynx", "Larynx"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of oropharynx", "Oropharynx"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of nasopharynx", "Nasopharynx"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of lip, oral cavity and pharynx", "Head & Neck"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of oesophagus", "Oesophagus"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of stomach", "Stomach"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of prostate", "Prostate"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of pancreas", "Pancreas"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of base of tongue", "Tongue"), # considered oropharynx
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of other and unspecified parts of tongue", "Unspecified Tongue"),
         
         Cancer = str_replace_all(Cancer, "Malignant neoplasm of liver and intrahepatic bile ducts", "Liver")
         
         
         
  ) %>% 
  filter(Cancer == "Colorectal" ) %>% 
  filter(is.na(Pop_std) | Pop_std == "ESP2013" )


lancet_colors <- c("#00468BFF", "#ED0000FF", "#00BFFF", "#0099B4FF", "#925E9FFF", "#FDAF17FF")

plot <- final_comb %>% 
  filter(type == "Age Std") %>% 
  ggplot(aes(x = Year, y = Rate, color = Database, group = Database)) +
  geom_line(size = 1) +  # Thicker lines for clarity
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Database), alpha = 0.2, linetype = 0) +   # Shaded area for confidence intervals
  scale_color_manual(values = lancet_colors) +  # Apply Lancet-style colors to lines
  scale_fill_manual(values = lancet_colors) +   # Apply Lancet-style colors to ribbons
  labs(
    x = "Calendar year",
    y = "Age Standardized Incidence rate per 100000 person-years (ESP2013)",
    color = "Database",
    fill = "Database"
  ) +
  facet_wrap(~ Sex) +  # Facet by Sex
  theme_classic(base_size = 14) +  # Clean, minimal theme
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.2),  # Light grey gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.title = element_text(face = "bold"),  # Bold axis titles
    legend.position = "top",  # Move the legend to the top
    legend.title = element_text(face = "bold"),  # Bold legend title
    strip.background = element_rect(color = "black"),  # White background for facet labels
    strip.text = element_text(face = "bold"),  # Bold facet labels
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add black border around the entire graph
  )


plotname <- paste0("FIGURE3_Incidence_crude_crc.pdf")

pdf(paste0(pathResults ,"/", plotname), width = 10, height = 6)

print(incidenceFigureData, newpage = FALSE)
dev.off()













