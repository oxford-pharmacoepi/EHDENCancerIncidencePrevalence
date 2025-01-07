# liver cancer extra sensitivity analysis figures

library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(quarto)
library(ggplot2)
library(scales)
library(ggh4x)
library(readr)
library(rio)
library(tidyverse)
library(dsr)
library(frailtypack)

#folder of data
datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data/liver"
#path to european population standard 2013
ESP13path <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/ESP13.csv"

#printing numbers with 3 decimal place and commas
nice.num3<-function(x) {
  trimws(format(round(x,3),
                big.mark=",", nsmall = 3, digits=3, scientific=FALSE))}

nice.num2<-function(x) {
  trimws(format(round(x,2),
                big.mark=",", nsmall = 2, digits=2, scientific=FALSE))}


# read in incidence - prevalence data and process
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
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "ProstateCancerMaleOnly", "Prostate")) %>% 
    mutate(outcome_cohort_name = replace(outcome_cohort_name, outcome_cohort_name == "PLC_with_cholangiocarcinoma", "Liver + Cholangiocarcinoma ")) 

  
  
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
    mutate(database_name = replace(database_name, database_name == "GoldUpdate4_liver", "CPRD GOLD"))
  
  #filter out the results for both genders for prostate cancer (as cohort only in male)
  result <- result %>%
    filter(!(outcome_cohort_name == "Prostate" & denominator_sex == "Both")) %>%
    filter(!(outcome_cohort_name == "Prostate" & denominator_sex == "Female"))
  
  return(result)
} # need to update this for the different files

# Load, prepare, and merge results -----
results <-list.files(datapath, full.names = TRUE,
                     recursive = TRUE,
                     include.dirs = TRUE,
                     pattern = ".zip")

#unzip data
for (i in (1:length(results))) {
  utils::unzip(zipfile = results[[i]],
               exdir = datapath)
}


datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data/liver/GoldUpdate4_liverIPResultsAgeStandardization_extra_liver"

#grab the results from the folders
results <- list.files(
  path = datapath,
  pattern = ".csv",
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
)



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
        paste0(datapath,  "/incidence_estimates.rds"))

#read in the incidence and prevalence files
incidence_estimates <- readRDS(paste0(datapath ,"/incidence_estimates.rds"))
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




ESP13_updated <- ESP13_updated %>% 
  rename(pop = ESP2013)

# INCIDENCE
#get data ready - all cancers

incidence_estimates1 <- incidence_estimates %>% 
  #filter(denominator_sex == "Both") %>% 
  #filter(outcome_cohort_name != "Breast" ) %>% 
  filter(analysis_interval != "overall") %>% 
  filter(database_name == "CPRD GOLD") %>% 
  filter(denominator_age_group != "All") %>% 
  rename(Agegroup = denominator_age_group) #rename this to make it the same as agegroup in standard population


agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()

for(i in 1:length(table(incidence_estimates1$outcome_cohort_name))){
  
  
  # for whole population    
  incidence_estimates_i <- incidence_estimates1 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates1$outcome_cohort_name)[i]),
           denominator_sex == "Both")
  
  
  tryCatch({
    
    agestandardizedinc[[i]] <- dsr::dsr(
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
      mutate(Cancer = names(table(incidence_estimates1$outcome_cohort_name)[i]),
             Sex = "Both") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates1$outcome_cohort_name)[i]), "for both sexes done"))
  
  ##########################
  # for female population    
  #########################
  incidence_estimates_i <- incidence_estimates1 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates1$outcome_cohort_name)[i]),
           denominator_sex == "Female")
  
  
  tryCatch({
    
    agestandardizedincF[[i]] <- dsr::dsr(
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
      mutate(Cancer = names(table(incidence_estimates1$outcome_cohort_name)[i]),
             Sex = "Female") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates1$outcome_cohort_name)[i]), "for females done"))
  
  
  
  
  ##########################
  # for male population    
  #########################
  incidence_estimates_i <- incidence_estimates1 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates1$outcome_cohort_name)[i]),
           denominator_sex == "Male")
  
  
  tryCatch({
    
    agestandardizedincM[[i]] <- dsr::dsr(
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
      mutate(Cancer = names(table(incidence_estimates1$outcome_cohort_name)[i]),
             Sex = "Male") 
    
  }, error = function(e) {
    # If an error occurs, print the error message
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  
  print(paste0("age standardization for ", names(table(incidence_estimates1$outcome_cohort_name)[i]), "for males done"))
  
  
  
  
  
}

agestandardizedinc_final <- bind_rows(agestandardizedinc,
                                      agestandardizedincF,
                                      agestandardizedincM)

#save the results
saveRDS(agestandardizedinc_final, paste0(datapath ,"/incidence_estimates_age_sd_liver.rds"))




# age standardized
incidenceFigureData <- agestandardizedinc_final %>%
  filter(Sex == "Both") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Cancer )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Cancer), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Cancer, fill = Cancer),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
 # geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Cancer",
       shape = "Cancer",
       fill = "Cancer") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) 


plotname <- paste0("/FIGURE1_IRs_Liver_sensitivity_ageadjusted.pdf")


pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()



# crude
incidenceFigureDataa <- agestandardizedinc_final %>%
  filter(Sex == "Both") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Crude Rate (per 1e+05)`,
             group = Cancer )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Crude)`, 
                  ymax = `95% UCL (Crude)`, 
                  fill = Cancer), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Cancer, fill = Cancer),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  # geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Crude Incidence rate per 100000 person-years",
       col = "Cancer",
       shape = "Cancer",
       fill = "Cancer") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) 


plotname <- paste0("/FIGURE2_IRs_Liver_sensitivity_crude.pdf")


pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureDataa, newpage = FALSE)

dev.off()


# sex strat crude
incidenceFigureData1 <- agestandardizedinc_final %>%
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Crude Rate (per 1e+05)`,
             group = Cancer )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Crude)`, 
                  ymax = `95% UCL (Crude)`, 
                  fill = Cancer), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Cancer, fill = Cancer),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  labs(x = "Calendar year",
       y = "Crude Incidence rate per 100000 person-years",
       col = "Cancer",
       shape = "Cancer",
       fill = "Cancer") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Sex, scales = "free_y", ncol = 3)

plotname <- paste0("/FIGURE_IRs_Liver_sensitivity_crude_overall_sex.pdf")
pdf(paste0(datapath , plotname),
    width = 12, height = 5)
print(incidenceFigureData1, newpage = FALSE)

dev.off()



# sex strat crude
incidenceFigureData1 <- agestandardizedinc_final %>%
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Cancer )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Cancer), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Cancer, fill = Cancer),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Cancer",
       shape = "Cancer",
       fill = "Cancer") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Sex, scales = "free_y", ncol = 3)

plotname <- paste0("/FIGURE_IRs_Liver_sensitivity_age_std_overall_sex.pdf")
pdf(paste0(datapath , plotname),
    width = 12, height = 5)
print(incidenceFigureData1, newpage = FALSE)

dev.off()


# survival plots

datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data/liver/GoldUpdate4_liverWholeSurvivalResults"

#grab the results from the folders
results <- list.files(
  path = datapath,
  pattern = ".csv",
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
)


survival_estimates_files<-results[stringr::str_detect(results, ".csv")]
survival_estimates_files<-results[stringr::str_detect(results, "survival_estimates")]
survival_estimates <- list()
for(i in seq_along(survival_estimates_files)){
  survival_estimates[[i]]<-readr::read_csv(survival_estimates_files[[i]],
                                            show_col_types = FALSE)
}
survival_estimates <- dplyr::bind_rows(survival_estimates)
#survival_estimates <- prepare_output(survival_estimates)

saveRDS(survival_estimates,
        paste0(datapath,  "/survival_estimates.rds"))


survivalFigureData <- survival_estimates %>%
  filter(Stratification == "None") %>%
  filter(Age == "All") %>%
  ggplot(aes(x = time,
             y = est,
             group = Cancer,
             col = Cancer )) +
  scale_y_continuous( labels = scales::percent, limits = c(0, NA)) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl, 
                  fill = Cancer), alpha = .15, color = NA, show.legend = FALSE) +
  geom_line(aes(linetype = Cancer),size = 0.5) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash","dotted")) +
  labs(x = "Time (Years)",
       y = "Survival Probability",
       col = "Cancer",
       linetype = "Cancer") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  scale_x_continuous(breaks=seq(0, 22, 2))

plotname <- paste0("/FIGURE_IRs_Liver_sensitivity_KM.pdf")
pdf(paste0(datapath , plotname),
    width = 12, height = 5)
print(survivalFigureData, newpage = FALSE)

dev.off()
