library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(quarto)


# Data prep functions -----
# printing numbers with 3 decimal place and commas 
nice.num3<-function(x) {
  trimws(format(round(x,3),
                big.mark=",", nsmall = 3, digits=3, scientific=FALSE))}

nice.num2<-function(x) {
  trimws(format(round(x,2),
                big.mark=",", nsmall = 2, digits=2, scientific=FALSE))}

#preparing the output and renaming numbers for incidence and prevalence
prepare_output<-function(result){
  result <- result %>%
    
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerHypopharynx", "Hypopharynx")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerLarynx", "Larynx")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerNasalCavitySinus", "Nasal Cavity & Sinus")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerNasopharynx", "Nasopharynx")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerOralCavityPrevalent", "Oral Cavity")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerOropharynx", "Oropharynx")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerSalivaryGland", "Salivary Gland")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerTonguePrevalent", "Tongue")) %>% 
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerOralCavityIncidence", "Oral Cavity")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "HeadNeckSubtypeCancerTongueIncidence", "Tongue")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentProstateCancer", "Prostate")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentLungCancer", "Lung")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentBreastCancer", "Breast")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentColorectalCancer", "Colorectal")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentHeadNeckCancer", "Head & Neck")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentLiverCancer", "Liver")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentPancreaticCancer", "Pancreas")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentStomachCancer", "Stomach")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "IncidentEsophagealCancer", "Oesophagus")) %>%    
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentProstateCancer", "Prostate")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentLungCancer", "Lung")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentBreastCancer", "Breast")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentColorectalCancer", "Colorectal")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentHeadNeckCancer", "Head & Neck")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentLiverCancer", "Liver")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentPancreaticCancer", "Pancreas")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentStomachCancer", "Stomach")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PrevalentEsophagealCancer", "Oesophagus")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantProstateCancer", "Prostate")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantLungCancer", "Lung")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantBreastCancer", "Breast")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantColorectalCancer", "Colorectal")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantHeadNeckCancer", "Head & Neck")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantLiverCancer", "Liver")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantPancreaticCancer", "Pancreas")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "MalignantStomachCancer", "Stomach")) %>%
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "ProstateCancerMaleOnly", "Prostate"))
  
  
  result<- result %>% 
    mutate(denominator_age_group= stringr::str_replace(denominator_age_group, ";", " to ")) %>% 
    mutate(denominator_age_group = replace(denominator_age_group, denominator_age_group == "18 to 150", "All")) %>% 
    mutate(denominator_age_group = replace(denominator_age_group, denominator_age_group == "90 to 150", "90 +")) %>% 
    mutate(denominator_age_group = factor(denominator_age_group,
                                          levels = c("All",
                                                     "18 to 29", "30 to 39", "40 to 49",
                                                     "50 to 59", "60 to 69", "70 to 79",
                                                     "80 to 89", "90 +" )))
  
  result <- result %>%
    mutate(database_name = replace(database_name, database_name == "CPRDAurum", "CPRD Aurum")) %>%
    mutate(database_name = replace(database_name, database_name == "CPRDGoldUpdate3", "CPRD GOLD")) 
  
  #filter out the results for both genders for prostate cancer (as cohort only in male)
  result <- result %>%
    filter(!(outcome_cohort_name == "Prostate" & denominator_sex == "Both")) %>%
    filter(!(outcome_cohort_name == "Prostate" & denominator_sex == "Female")) 
  
  return(result)
} # need to update this for the different files

#preparation the output and renaming numbers for survival
prepare_output_survival <- function(result){
  result <- result %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerHypopharynx", "Hypopharynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerLarynx", "Larynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerNasalCavitySinus", "Nasal Cavity & Sinus")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerNasopharynx", "Nasopharynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerOralCavityPrevalent", "Oral Cavity")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerOropharynx", "Oropharynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerSalivaryGland", "Salivary Gland")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerTonguePrevalent", "Tongue")) %>% 
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerOralCavityIncidence", "Oral Cavity")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerTongueIncidence", "Tongue")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentProstateCancer", "Prostate")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentLungCancer", "Lung")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentBreastCancer", "Breast")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentColorectalCancer", "Colorectal")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentHeadNeckCancer", "Head & Neck")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentLiverCancer", "Liver")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentPancreaticCancer", "Pancreas")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentStomachCancer", "Stomach")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentEsophagealCancer", "Oesophagus")) %>%    
    mutate(Cancer = replace(Cancer, Cancer == "PrevalentProstateCancer", "Prostate")) 
  
  
  result<- result %>% 
    mutate(Age = replace(Age, Age == ">=90", "90 +")) %>% 
    mutate(Age = replace(Age, Age == "<30", "18-29")) %>% 
    mutate(Age= stringr::str_replace(Age, "-", " to ")) %>% 
    mutate(CalenderYearGp= stringr::str_replace(CalenderYearGp, "-", " to ")) %>% 
    mutate(Age = factor(Age,
                        levels = c("All",
                                   "18 to 29", "30 to 39", "40 to 49",
                                   "50 to 59", "60 to 69", "70 to 79",
                                   "80 to 89", "90 +" ))) %>%
    mutate(CalenderYearGp = factor(CalenderYearGp,
                                   levels = c("2000 to 2019",
                                              "2000 to 2021",
                                              "2000 to 2004", 
                                              "2005 to 2009", 
                                              "2010 to 2014",
                                              "2015 to 2019",
                                              "2020 to 2021")))

  
  result <- result %>%
    mutate(Database = replace(Database, Database == "CPRDAurum", "CPRD Aurum")) %>%
    mutate(Database = replace(Database, Database == "CPRDGoldUpdate3", "CPRD GOLD")) %>% 
    mutate(Database = replace(Database, Database == "CPRDGoldUpdate2", "CPRD GOLD"))
  
  result <- result %>%
    mutate(Gender=replace(Gender, Cancer=="Prostate", "Male"))
  
  return(result)
} 

##preparation the output for table 1
prepare_output_table1 <- function(result){
  result <- result %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerHypopharynx", "Hypopharynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerLarynx", "Larynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerNasalCavitySinus", "Nasal Cavity & Sinus")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerNasopharynx", "Nasopharynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerOralCavityPrevalent", "Oral Cavity")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerOropharynx", "Oropharynx")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerSalivaryGland", "Salivary Gland")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerTonguePrevalent", "Tongue")) %>% 
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerOralCavityIncidence", "Oral Cavity")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "HeadNeckSubtypeCancerTongueIncidence", "Tongue")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentProstateCancer", "Prostate")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentLungCancer", "Lung")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentBreastCancer", "Breast")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentColorectalCancer", "Colorectal")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentHeadNeckCancer", "Head & Neck")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentLiverCancer", "Liver")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentPancreaticCancer", "Pancreas")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentStomachCancer", "Stomach")) %>%
    mutate(Cancer = replace(Cancer, Cancer == "IncidentEsophagealCancer", "Oesophagus")) %>%    
    mutate(Cancer = replace(Cancer, Cancer == "PrevalentProstateCancer", "Prostate")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentBreastCancer", "Breast")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentColorectalCancer", "Colorectal")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentEsophagealCancer", "Oesophagus")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentHeadNeckCancer", "Head & Neck")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentLiverCancer", "Liver")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentLungCancer", "Lung")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentPancreaticCancer", "Pancreas")) %>%
  mutate(Cancer = replace(Cancer, Cancer == "PrevalentStomachCancer", "Stomach")) 
  
  
  
  
  
  result <- result %>%
    mutate(Database = replace(Database, Database == "CPRDAurum", "CPRD Aurum")) %>%
    mutate(Database = replace(Database, Database == "CPRDGoldUpdate3", "CPRD GOLD")) %>% 
  mutate(Database = replace(Database, Database == "CPRDGoldUpdate2", "CPRD GOLD"))
  
  
  
  return(result)
} 

# Load, prepare, and merge results -----
results <-list.files(here("networkResults"), full.names = TRUE,
                     recursive = TRUE,
                     include.dirs = TRUE,
                     pattern = ".zip")

#unzip data
for (i in (1:length(results))) {
  utils::unzip(zipfile = results[[i]],
               exdir = here("networkResults"))
}

#grab the results from the folders
results <- list.files(
  path = here("networkResults"),
  pattern = ".csv",
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
)



# merge the prevalence estimates
prevalence_estimates_files<-results[stringr::str_detect(results, ".csv")]
prevalence_estimates_files<-results[stringr::str_detect(results, "prevalence_estimates")]


prevalence_estimates <- list()
for(i in seq_along(prevalence_estimates_files)){
  prevalence_estimates[[i]]<-readr::read_csv(prevalence_estimates_files[[i]], 
                                             show_col_types = FALSE)  
}
prevalence_estimates <- dplyr::bind_rows(prevalence_estimates)
prevalence_estimates <- prepare_output(prevalence_estimates)
prevalence_estimates <- prevalence_estimates %>% 
  mutate("Prevalence (95% CI)"= ifelse(!is.na(prevalence),
                                       paste0(paste0(nice.num3(prevalence*100), "%"), " (",
                                              paste0(nice.num3(prevalence_95CI_lower*100), "%")," to ", 
                                              paste0(nice.num3(prevalence_95CI_upper*100), "%"), ")"),
                                       NA
  ))
saveRDS(prevalence_estimates, 
        here("shiny", "data", "prevalence_estimates.rds"))


# prevalence attrition
prevalence_attrition_files<-results[stringr::str_detect(results, ".csv")]
prevalence_attrition_files<-results[stringr::str_detect(results, "prevalence_attrition")]
prevalence_attrition <- list()
for(i in seq_along(prevalence_attrition_files)){
  prevalence_attrition[[i]]<-readr::read_csv(prevalence_attrition_files[[i]], 
                                             show_col_types = FALSE)  
}
prevalence_attrition <- dplyr::bind_rows(prevalence_attrition)
prevalence_attrition <- prepare_output(prevalence_attrition)
saveRDS(prevalence_attrition, 
        here("shiny", "data", "/prevalence_attrition.rds"))


#merge incidence results together
# incidence estimates
incidence_estimates_files<-results[stringr::str_detect(results, ".csv")]
incidence_estimates_files<-results[stringr::str_detect(results, "incidence_estimates")]
incidence_estimates <- list()
for(i in seq_along(incidence_estimates_files)){
  incidence_estimates[[i]]<-readr::read_csv(incidence_estimates_files[[i]], 
                                            show_col_types = FALSE)  
}
incidence_estimates <- dplyr::bind_rows(incidence_estimates)
incidence_estimates <- prepare_output(incidence_estimates)
saveRDS(incidence_estimates, 
        here("shiny", "data", "/incidence_estimates.rds"))

# incidence attrition
incidence_attrition_files<-results[stringr::str_detect(results, ".csv")]
incidence_attrition_files<-results[stringr::str_detect(results, "incidence_attrition")]
incidence_attrition <- list()
for(i in seq_along(incidence_attrition_files)){
  incidence_attrition[[i]]<-readr::read_csv(incidence_attrition_files[[i]], 
                                            show_col_types = FALSE)  
}
incidence_attrition <- dplyr::bind_rows(incidence_attrition)
incidence_attrition <- prepare_output(incidence_attrition)
saveRDS(incidence_attrition, 
        here("shiny", "data", "/incidence_attrition.rds"))


# merge the survival results together
survival_estimates_files<-results[stringr::str_detect(results, ".csv")]
survival_estimates_files<-results[stringr::str_detect(results, "survival_estimates")]
survival_estimates <- list()
for(i in seq_along(survival_estimates_files)){
  survival_estimates[[i]]<-readr::read_csv(survival_estimates_files[[i]],
                                           show_col_types = FALSE)
}
survival_estimates <- dplyr::bind_rows(survival_estimates)
survival_estimates <- prepare_output_survival(survival_estimates)
saveRDS(survival_estimates,
        here("shiny", "data", "/survival_estimates.rds"))


# calculate follow up time from survival estimates

# - only include whole study period - loop over each database each cancer, whole population and sex and age
#filter out calender year results and age*sex stratification
# #add together total number of people censored (event or left the database) and calculate the time so time * number of people with that time
test <- survival_estimates %>%  
  filter(CalenderYearGp == "2000 to 2021" | CalenderYearGp == "2000 to 2019" 
           ) %>% 
  filter(Stratification != "Age*Gender" ) %>% 
  mutate(totalnp = n.event + n.censor) %>% 
  mutate(folup_agg = time * totalnp )


# create a new dataframe which expands each value of time by number censored = so if there were 4 people censored at time 0.5 there will be 4 rows with 0.5
test1 <- test %>% 
  uncount(totalnp)

survival_median_mean_follow_up <-  test1 %>% 
  group_by(Cancer, Database, Age, Gender, Stratification) %>%
  summarise(across(time, c(median, 
                           Q1=~quantile(., probs = 0.25),
                           Q3=~quantile(., probs = 0.75),
                           mean,
                           sd))) %>% 
  rename(median_followup = time_1) %>% 
  rename(mean_followup = time_4) %>% 
  rename(sd_followup = time_5) %>% 
  rename(lower_IQR = time_Q1) %>% 
  rename(upper_IQR = time_Q3)

saveRDS(survival_median_mean_follow_up,
        here("shiny", "data", "/survival_median_mean_follow_up.rds"))

# merge the risk table results together (whole dataset)
survival_risk_table_files<-results[stringr::str_detect(results, ".csv")]
survival_risk_table_files<-results[stringr::str_detect(results, "risk_table_results")]
survival_risk_table_files <- survival_risk_table_files[!stringr::str_detect(survival_risk_table_files, "risk_table_results_cy")]

survival_risk_table <- list()
for(i in seq_along(survival_risk_table_files)){
  survival_risk_table[[i]]<-readr::read_csv(survival_risk_table_files[[i]],
                                            show_col_types = FALSE) %>%
    mutate_if(is.double, as.character)
  
}

survival_risk_table <- dplyr::bind_rows(survival_risk_table)
survival_risk_table <- prepare_output_survival(survival_risk_table)
saveRDS(survival_risk_table,
        here("shiny", "data", "/survival_risk_table.rds"))



# merge the risk table results together (calender year results)
survival_risk_table_cy_files<-results[stringr::str_detect(results, ".csv")]
survival_risk_table_cy_files<-results[stringr::str_detect(results, "risk_table_results_cy")]
survival_risk_cy_table <- list()
for(i in seq_along(survival_risk_table_cy_files)){
  survival_risk_cy_table[[i]]<-readr::read_csv(survival_risk_table_cy_files[[i]],
                                           show_col_types = FALSE)  %>%
    mutate_if(is.double, as.character)
}
survival_risk_cy_table <- dplyr::bind_rows(survival_risk_cy_table)
survival_risk_cy_table <- prepare_output_survival(survival_risk_cy_table)
saveRDS(survival_risk_cy_table,
        here("shiny", "data", "/survival_risk_table_cy.rds"))


# merge the median results together
survival_median_files<-results[stringr::str_detect(results, ".csv")]
survival_median_files<-results[stringr::str_detect(results, "median_survival_results")]
survival_median_table <- list()
for(i in seq_along(survival_median_files)){
  survival_median_table[[i]]<-readr::read_csv(survival_median_files[[i]],
                                               show_col_types = FALSE)  
}
survival_median_table <- dplyr::bind_rows(survival_median_table)
survival_median_table <- prepare_output_survival(survival_median_table)
survival_median_table <- survival_median_table 

# round the values before turning them into characters
survival_median_table <-  survival_median_table %>% 
  mutate(rmean = nice.num3(rmean)) %>% 
  mutate(`se(rmean)` = nice.num3(`se(rmean)`)) %>%
  mutate(median = nice.num3(median)) %>%
  mutate(`0.95LCL` = nice.num3(`0.95LCL`)) %>%
  mutate(`0.95UCL` = nice.num3(`0.95UCL`)) %>% 
  mutate(median = ifelse(median == "NA", NA, median)) %>% 
  mutate(`0.95LCL` = ifelse(`0.95LCL` == "NA", NA, `0.95LCL`)) %>% 
  mutate(`0.95UCL` = ifelse(`0.95UCL` == "NA", NA, `0.95UCL`)) 
  
  


#if events less than 5 turn the result into NA
survival_median_table <-  
  survival_median_table %>% 
mutate(events = ifelse(events <= 5, "<5", events)) %>% 
mutate(median = ifelse(events == "<5", " ", median)) %>% 
mutate(records = ifelse(events == "<5", " ", records)) %>% 
mutate(n.max = ifelse(events == "<5", " ", n.max)) %>% 
mutate(n.start = ifelse(events == "<5", " ", n.start)) %>% 
mutate(rmean = ifelse(events == "<5", " ", rmean)) %>% 
mutate(`se(rmean)` = ifelse(events == "<5", " ", `se(rmean)`)) %>% 
mutate(`0.95LCL` = ifelse(events == "<5", " ", `0.95LCL`)) %>% 
mutate(`0.95UCL` = ifelse(events == "<5", " ",`0.95UCL`)) 

# put reason for obscuring i.e. median not achieved

survival_median_table <-  
  survival_median_table %>% 
  mutate(median = ifelse(is.na(`0.95UCL`) == TRUE, "Not achieved", median)) %>% 
  mutate(`0.95LCL` = ifelse(is.na(`0.95LCL`) == TRUE, "Not calculated", `0.95LCL`)) %>% 
  mutate(`0.95UCL` = ifelse(is.na(`0.95UCL`) == TRUE, "Not calculated",`0.95UCL`)) %>% 
  mutate(`0.95LCL` = ifelse(median == "Not achieved", "Not calculated",`0.95LCL`)) %>% 
  mutate(`0.95LCL` = ifelse(events == "<5", "Result obscured",`0.95LCL`))%>% 
  mutate(`0.95UCL` = ifelse(events == "<5", "Result obscured",`0.95UCL`)) %>% 
  mutate(median = ifelse(events == "<5", "Result obscured",median))%>% 
  mutate(rmean = ifelse(events == "<5", "Result obscured", rmean)) %>% 
  mutate(`se(rmean)` = ifelse(events == "<5", "Result obscured", `se(rmean)`)) 

saveRDS(survival_median_table,
        here("shiny", "data", "/survival_median_table.rds"))



# merge the 1, 5 and 10 survival rates result together
survival_rates_files<-results[stringr::str_detect(results, ".csv")]
survival_rates_files<-results[stringr::str_detect(results, "one_five_ten_survival_rates")]
survival_rates_table <- list()
for(i in seq_along(survival_rates_files)){
  survival_rates_table[[i]]<-readr::read_csv(survival_rates_files[[i]],
                                              show_col_types = FALSE)  
}
survival_rates_table <- dplyr::bind_rows(survival_rates_table)
survival_rates_table <- prepare_output_survival(survival_rates_table)

# obscure those with n.event and n.censor both less than 5
survival_rates_table  <-  
  survival_rates_table  %>% 
  mutate(n.event = ifelse(n.event <= 5, "<5", n.event)) %>% 
  mutate(n.risk = ifelse(n.risk <= 5, "<5", n.risk)) %>% 
  mutate(n.censor = ifelse(n.censor <= 5, "<5", n.censor)) %>% 
  mutate(n.risk = ifelse(n.risk == "<5" & time == 5 & CalenderYearGp != "2000 to 2019", "-", n.risk)) %>% 
  mutate(n.censor = ifelse(n.risk == "<5" & n.event == "<5" & time == 10 & CalenderYearGp != "2000 to 2019"
                       , "<5", n.censor)) %>% 
  mutate(surv = ifelse(n.risk == "<5" & n.event == "<5" & n.censor == "<5"
                       , NA, surv)) %>% 
  mutate(lower = ifelse(n.risk == "<5" & n.event == "<5" & n.censor == "<5"
                       , NA, lower)) %>%   
  mutate(upper = ifelse(n.risk == "<5" & n.event == "<5" & n.censor == "<5"
                       , NA, upper)) %>% 
  mutate(std.err = ifelse(n.risk == "<5" & n.event == "<5" & n.censor == "<5"
                        , NA, std.err)) %>%  
  mutate(cumhaz = ifelse(n.risk == "<5" & n.event == "<5" & n.censor == "<5"
                        , NA, cumhaz)) %>% 
  mutate(std.chaz = ifelse(n.risk == "<5" & n.event == "<5" & n.censor == "<5"
                        , NA, std.chaz)) %>% 
  mutate(surv = ifelse(n.risk == "-" & n.event == "<5" & n.censor == "<5"
                       , NA, surv)) %>% 
  mutate(lower = ifelse(n.risk == "-" & n.event == "<5" & n.censor == "<5"
                        , NA, lower)) %>%   
  mutate(upper = ifelse(n.risk == "-" & n.event == "<5" & n.censor == "<5"
                        , NA, upper)) %>% 
  mutate(std.err = ifelse(n.risk == "-" & n.event == "<5" & n.censor == "<5"
                          , NA, std.err)) %>%  
  mutate(cumhaz = ifelse(n.risk == "-" & n.event == "<5" & n.censor == "<5"
                         , NA, cumhaz)) %>% 
  mutate(std.chaz = ifelse(n.risk == "-" & n.event == "<5" & n.censor == "<5"
                           , NA, std.chaz)) %>% 
  mutate(surv = ifelse(is.na(std.err), NA, surv)) %>% 
  mutate(lower = ifelse(is.na(std.err), NA, lower)) %>% 
  mutate(upper = ifelse(is.na(std.err), NA, upper)) %>% 
  mutate(cumhaz = ifelse(is.na(std.err), NA, cumhaz)) %>%  
  mutate(std.chaz = ifelse(is.na(std.err), NA, std.chaz)) 
  
saveRDS(survival_rates_table,
        here("shiny", "data", "/survival_rates_table.rds"))


# table 1 BOTH genders
table1_files<-results[stringr::str_detect(results, ".csv")]
table1_files<-results[stringr::str_detect(results, "Table1[C|H]")]

table1_results <- list()

for(i in seq_along(table1_files)){
  table1_results[[i]]<-readr::read_csv(table1_files[[i]], 
                                             show_col_types = FALSE)  
}

table1_results <- dplyr::bind_rows(table1_results)

table1_results <- prepare_output_table1(table1_results)

table1_results <- table1_results %>% 
  distinct()

table1_results <- table1_results %>% 
  mutate("Variable"= ifelse(!is.na(percent),
                                       paste0(n, " (",
                                              paste0(percent, ")")),
                                       NA
  )) %>%
  mutate("Variable1"= ifelse(!is.na(mean),
                            paste0(mean, " (SD ",
                                   paste0(standard_deviation, ")")),
                            NA
  )) %>%
  mutate("Variable2"= ifelse(!is.na(median),
                             paste0(median, " (",
                                    paste0(interquartile_range, ")")),
                             NA
  )) %>%
  mutate(Variable = case_when(n == "<5" ~ "<5", TRUE ~ Variable)) %>%
  mutate(Variable = coalesce(Variable, Variable2)) %>%
  mutate(Variable = coalesce(Variable, n))

# remove rows we dont need
table1_results <- table1_results %>% 
  filter(!grepl("Sex: Female",var)) %>%
  filter(!grepl("Death: Alive",var)) %>%
  filter(!grepl("Prior_history_days_study_start",var))%>%
  #filter(!grepl("Prior_history_years_start",var)) %>%
  filter(!grepl("Death: Dead",var)) %>%
  filter(!grepl("time_days",var)) %>%
  filter(!grepl("time_years",var)) %>%
  select(c(var, Variable, Cancer, Database, analysis ))

# make into wider format for results
table1_results <- table1_results %>% 
  tidyr::pivot_wider(names_from = Database,
                     values_from = Variable, 
                     values_fill = NA
  )

# remove drugs for CPRD
table1_results <- table1_results %>% 
  filter(!grepl("AgentsReninAngiotensinSystem",var)) %>%
  filter(!grepl("Antidepressants",var)) %>%
  filter(!grepl("Antiepileptics",var))%>%
  filter(!grepl("AntiinflammatoryAntirheumatic",var)) %>%
  filter(!grepl("Antineoplastics",var)) %>%
  filter(!grepl("Antipsoriatics",var)) %>%
  filter(!grepl("Antipsychotics",var)) %>%
  filter(!grepl("Antithrombotics",var)) %>%
  filter(!grepl("Anxiolytics",var)) %>%
  filter(!grepl("BetaBlockers",var)) %>%
  filter(!grepl("CalciumChannelBlockers",var)) %>%
  filter(!grepl("Diuretics",var))%>%
  filter(!grepl("DrugsAcidRelatedDisorders",var)) %>%
  filter(!grepl("DrugsDiabetes",var)) %>%
  filter(!grepl("DrugsObstructiveAirwayDiseases",var)) %>%
  filter(!grepl("HypnoticsSedatives",var)) %>%
  filter(!grepl("Immunosuppressants",var)) %>%
  filter(!grepl("LipidModifyingAgents",var)) %>%
  filter(!grepl("Opioids",var)) %>%
  filter(!grepl("Psychostimulants",var)) %>% 
  filter(!grepl("HpyloriGIInfection",var)) %>% 
  mutate(Sex = "Both") %>% 
  mutate(across(everything(), as.character))
  

# saveRDS(table1_results, 
#         here("shiny", "data", "table1_results.rds"))
# 



# females
table1_filesF<-results[stringr::str_detect(results, ".csv")]
table1_filesF<-results[stringr::str_detect(results, "Female")]

table1_resultsF <- list()

for(i in seq_along(table1_filesF)){
  table1_resultsF[[i]]<-readr::read_csv(table1_filesF[[i]], 
                                       show_col_types = FALSE)  
}

table1_resultsF <- dplyr::bind_rows(table1_resultsF)

table1_resultsF <- prepare_output_table1(table1_resultsF)

table1_resultsF <- table1_resultsF %>% 
  distinct()

table1_resultsF <- table1_resultsF %>% 
  mutate("Variable"= ifelse(!is.na(percent),
                            paste0(n, " (",
                                   paste0(percent, ")")),
                            NA
  )) %>%
  mutate("Variable1"= ifelse(!is.na(mean),
                             paste0(mean, " (SD ",
                                    paste0(standard_deviation, ")")),
                             NA
  )) %>%
  mutate("Variable2"= ifelse(!is.na(median),
                             paste0(median, " (",
                                    paste0(interquartile_range, ")")),
                             NA
  )) %>%
  mutate(Variable = case_when(n == "<5" ~ "<5", TRUE ~ Variable)) %>%
  mutate(Variable = coalesce(Variable, Variable2)) %>%
  mutate(Variable = coalesce(Variable, n))

# remove rows we dont need
table1_resultsF <- table1_resultsF %>% 
  #filter(!grepl("Sex: Female",var)) %>%
  filter(!grepl("Death: Alive",var)) %>%
  filter(!grepl("Prior_history_days_study_start",var))%>%
  #filter(!grepl("Prior_history_years_start",var)) %>%
  filter(!grepl("Death: Dead",var)) %>%
  filter(!grepl("time_days",var)) %>%
  filter(!grepl("time_years",var)) %>%
  select(c(var, Variable, Cancer, Database, analysis, Gender ))

# make into wider format for results
table1_resultsF <- table1_resultsF %>% 
  tidyr::pivot_wider(names_from = Database,
                     values_from = Variable, 
                     values_fill = NA
  )

# remove drugs for CPRD
table1_resultsF <- table1_resultsF %>% 
  filter(!grepl("AgentsReninAngiotensinSystem",var)) %>%
  filter(!grepl("Antidepressants",var)) %>%
  filter(!grepl("Antiepileptics",var))%>%
  filter(!grepl("AntiinflammatoryAntirheumatic",var)) %>%
  filter(!grepl("Antineoplastics",var)) %>%
  filter(!grepl("Antipsoriatics",var)) %>%
  filter(!grepl("Antipsychotics",var)) %>%
  filter(!grepl("Antithrombotics",var)) %>%
  filter(!grepl("Anxiolytics",var)) %>%
  filter(!grepl("BetaBlockers",var)) %>%
  filter(!grepl("CalciumChannelBlockers",var)) %>%
  filter(!grepl("Diuretics",var))%>%
  filter(!grepl("DrugsAcidRelatedDisorders",var)) %>%
  filter(!grepl("DrugsDiabetes",var)) %>%
  filter(!grepl("DrugsObstructiveAirwayDiseases",var)) %>%
  filter(!grepl("HypnoticsSedatives",var)) %>%
  filter(!grepl("Immunosuppressants",var)) %>%
  filter(!grepl("LipidModifyingAgents",var)) %>%
  filter(!grepl("Opioids",var)) %>%
  filter(!grepl("Psychostimulants",var)) %>% 
  filter(!grepl("HpyloriGIInfection",var)) %>% 
  rename(Sex = Gender) %>% 
  mutate(across(everything(), as.character))


# Males
table1_filesM<-results[stringr::str_detect(results, ".csv")]
table1_filesM<-results[stringr::str_detect(results, "Male")]

table1_resultsM <- list()

for(i in seq_along(table1_filesM)){
  table1_resultsM[[i]]<-readr::read_csv(table1_filesM[[i]], 
                                        show_col_types = FALSE)  
}

table1_resultsM <- dplyr::bind_rows(table1_resultsM)

table1_resultsM <- prepare_output_table1(table1_resultsM)

table1_resultsM <- table1_resultsM %>% 
  distinct()

table1_resultsM <- table1_resultsM %>% 
  mutate("Variable"= ifelse(!is.na(percent),
                            paste0(n, " (",
                                   paste0(percent, ")")),
                            NA
  )) %>%
  mutate("Variable1"= ifelse(!is.na(mean),
                             paste0(mean, " (SD ",
                                    paste0(standard_deviation, ")")),
                             NA
  )) %>%
  mutate("Variable2"= ifelse(!is.na(median),
                             paste0(median, " (",
                                    paste0(interquartile_range, ")")),
                             NA
  )) %>%
  mutate(Variable = case_when(n == "<5" ~ "<5", TRUE ~ Variable)) %>%
  mutate(Variable = coalesce(Variable, Variable2)) %>%
  mutate(Variable = coalesce(Variable, n))

# remove rows we dont need
table1_resultsM <- table1_resultsM %>% 
  #filter(!grepl("Sex: Female",var)) %>%
  filter(!grepl("Death: Alive",var)) %>%
  filter(!grepl("Prior_history_days_study_start",var))%>%
  #filter(!grepl("Prior_history_years_start",var)) %>%
  filter(!grepl("Death: Dead",var)) %>%
  filter(!grepl("time_days",var)) %>%
  filter(!grepl("time_years",var)) %>%
  select(c(var, Variable, Cancer, Database, analysis, Gender ))

# make into wider format for results
table1_resultsM <- table1_resultsM %>% 
  tidyr::pivot_wider(names_from = Database,
                     values_from = Variable, 
                     values_fill = NA
  )

# remove drugs for CPRD
table1_resultsM <- table1_resultsM %>% 
  filter(!grepl("AgentsReninAngiotensinSystem",var)) %>%
  filter(!grepl("Antidepressants",var)) %>%
  filter(!grepl("Antiepileptics",var))%>%
  filter(!grepl("AntiinflammatoryAntirheumatic",var)) %>%
  filter(!grepl("Antineoplastics",var)) %>%
  filter(!grepl("Antipsoriatics",var)) %>%
  filter(!grepl("Antipsychotics",var)) %>%
  filter(!grepl("Antithrombotics",var)) %>%
  filter(!grepl("Anxiolytics",var)) %>%
  filter(!grepl("BetaBlockers",var)) %>%
  filter(!grepl("CalciumChannelBlockers",var)) %>%
  filter(!grepl("Diuretics",var))%>%
  filter(!grepl("DrugsAcidRelatedDisorders",var)) %>%
  filter(!grepl("DrugsDiabetes",var)) %>%
  filter(!grepl("DrugsObstructiveAirwayDiseases",var)) %>%
  filter(!grepl("HypnoticsSedatives",var)) %>%
  filter(!grepl("Immunosuppressants",var)) %>%
  filter(!grepl("LipidModifyingAgents",var)) %>%
  filter(!grepl("Opioids",var)) %>%
  filter(!grepl("Psychostimulants",var)) %>% 
  filter(!grepl("HpyloriGIInfection",var)) %>% 
  rename(Sex = Gender) %>% 
  mutate(across(everything(), as.character))


table1_results_comb <- dplyr::bind_rows(table1_results, table1_resultsF, table1_resultsM)


saveRDS(table1_results_comb,
        here("shiny", "data", "table1_results.rds"))


