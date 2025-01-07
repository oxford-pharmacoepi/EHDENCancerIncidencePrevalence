# esophagus cancer extra sensitivity analysis figures

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
datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data/esophagus"
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
incidence_estimates_files <- incidence_estimates_files[c(5:8)]
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
  filter(analysis_interval != "overall") %>% 
  filter(denominator_age_group != "All") %>% 
  rename(Agegroup = denominator_age_group) #rename this to make it the same as agegroup in standard population


agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()
agestandardizedinc_region <- list()
agestandardizedincF_region <- list()
agestandardizedincM_region <- list()

# overall
for(i in 1:length(table(incidence_estimates1$outcome_cohort_name))){
  
  incidence_estimates_i <- incidence_estimates1 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates1$outcome_cohort_name)[i]),
           denominator_sex == "Both")
  

  
  for(j in 1:length(table(incidence_estimates1$region))){
  # for whole population    
  incidence_estimates_j <- incidence_estimates_i %>%
    filter(region == names(table(incidence_estimates_i$region)[j]))
  
  
#  tryCatch({
    
    agestandardizedinc_region[[j]] <- dsr::dsr(
      data = incidence_estimates_j,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedinc_region[[j]] <- agestandardizedinc_region[[j]] %>% 
      mutate(Cancer = names(table(incidence_estimates1$outcome_cohort_name)[i]),
             region = names(table(incidence_estimates1$region)[j]),
             Sex = "Both") 
    
  }
  
  agestandardizedinc[[i]] <- agestandardizedinc_region
  
}

# females
for(i in 1:length(table(incidence_estimates1$outcome_cohort_name))){
  
  incidence_estimates_i <- incidence_estimates1 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates1$outcome_cohort_name)[i]),
           denominator_sex == "Female")
  
  
  
  for(j in 1:length(table(incidence_estimates1$region))){
    # for whole population    
    incidence_estimates_j <- incidence_estimates_i %>%
      filter(region == names(table(incidence_estimates_i$region)[j]))
    
    
    #  tryCatch({
    
    agestandardizedincF_region[[j]] <- dsr::dsr(
      data = incidence_estimates_j,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincF_region[[j]] <- agestandardizedincF_region[[j]] %>% 
      mutate(Cancer = names(table(incidence_estimates1$outcome_cohort_name)[i]),
             region = names(table(incidence_estimates1$region)[j]),
             Sex = "Female") 
    
  }
  
  agestandardizedincF[[i]] <- agestandardizedincF_region
  
}

# males
for(i in 1:length(table(incidence_estimates1$outcome_cohort_name))){
  
  incidence_estimates_i <- incidence_estimates1 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates1$outcome_cohort_name)[i]),
           denominator_sex == "Male")
  
  
  
  for(j in 1:length(table(incidence_estimates1$region))){
    # for whole population    
    incidence_estimates_j <- incidence_estimates_i %>%
      filter(region == names(table(incidence_estimates_i$region)[j]))
    
    
    #  tryCatch({
    
    agestandardizedincM_region[[j]] <- dsr::dsr(
      data = incidence_estimates_j,  # specify object containing number of deaths per stratum
      event = n_events,       # column containing number of deaths per stratum 
      fu = person_years , # column containing number of population per stratum person years
      subgroup = incidence_start_date,   
      refdata = ESP13_updated, # reference population data frame, with column called pop
      method = "gamma",      # method to calculate 95% CI
      sig = 0.95,            # significance level
      mp = 100000,           # we want rates per 100.000 population
      decimals = 2) 
    
    agestandardizedincM_region[[j]] <- agestandardizedincM_region[[j]] %>% 
      mutate(Cancer = names(table(incidence_estimates1$outcome_cohort_name)[i]),
             region = names(table(incidence_estimates1$region)[j]),
             Sex = "Male") 
    
  }
  
  agestandardizedincM[[i]] <- agestandardizedincM_region
  
}


agestandardizedinc_final <- bind_rows(
  agestandardizedinc
)

agestandardizedinc_finalF <- bind_rows(
  agestandardizedincF
)

agestandardizedinc_finalM <- bind_rows(
  agestandardizedincM
)

agestandardizedinc_final1 <- bind_rows(
  agestandardizedinc_final,
  agestandardizedinc_finalF,
  agestandardizedinc_finalM
  
) %>% 
mutate(Subgroup = as.Date(Subgroup, format = "%d/%m/%Y"))

#save the results
saveRDS(agestandardizedinc_final, paste0(datapath ,"/incidence_estimates_age_sd_oc.rds"))




# age standardized
incidenceFigureData <- agestandardizedinc_final1 %>%
  filter(Sex == "Both") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = region )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = region), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(colour = region),size = 2) +
  #scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "region",
       shape = "region",
       fill = "region") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~region)

incidenceFigureData

plotname <- paste0("/FIGURE1_IRs_esophageal_sensitivity_ageadjusted.pdf")


pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()



#males

# age standardized
incidenceFigureData <- agestandardizedinc_final1 %>%
  filter(Sex == "Male") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = region )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = region), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(colour = region),size = 2) +
  #scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "region",
       shape = "region",
       fill = "region") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~region)


plotname <- paste0("/FIGURE1_IRs_esophageal_sensitivity_ageadjusted_males.pdf")


pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()

# females
# age standardized
incidenceFigureData <- agestandardizedinc_final1 %>%
  filter(Sex == "Female") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = region )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = region), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(colour = region),size = 2) +
  #scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "region",
       shape = "region",
       fill = "region") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~region)

incidenceFigureData

plotname <- paste0("/FIGURE1_IRs_esophageal_sensitivity_ageadjusted_females.pdf")


pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()






# compare with NCRAS
datapath_std <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data/NCRAS_Oxford_oesophagus.csv"
oxford_ncras <- read_csv(datapath_std) %>% 
  filter(`data type` != "GOLD") %>% 
  rename(Sex = sex,
         Cancer = cancer,
         Rate = estimate,
         LCL = lcl,
           UCL = ucl,
         Year = year,
        region = `data type`) %>% 
  mutate(type = "Age Std") %>% 
  mutate(Database = "NCRAS")## %>% 
  #mutate(Database = paste(Database,type, region) )


df_long <- agestandardizedinc_final %>%
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
  mutate(Database = "CPRD GOLD") %>% 
 # mutate(Database = paste(Database,type, region) ) %>% 
  rename(Year = Subgroup) %>% 
  dplyr::select(-Numerator, -Denominator) %>% 
  mutate(Year = as.numeric(substr(Year, 7, 10)))


df_longf <- agestandardizedinc_finalF %>%
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

df_widef <- df_longf %>%
  pivot_wider(
    names_from = name,  # Column names come from the 'name' column
    values_from = Value  # Values come from the 'Value' column
  ) %>% 
  mutate(Database = "CPRD GOLD") %>% 
  # mutate(Database = paste(Database,type, region) ) %>% 
  rename(Year = Subgroup) %>% 
  dplyr::select(-Numerator, -Denominator) %>% 
  mutate(Year = as.numeric(substr(Year, 7, 10)))


df_longm <- agestandardizedinc_finalM %>%
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

df_widem <- df_longm %>%
  pivot_wider(
    names_from = name,  # Column names come from the 'name' column
    values_from = Value  # Values come from the 'Value' column
  ) %>% 
  mutate(Database = "CPRD GOLD") %>% 
  # mutate(Database = paste(Database,type, region) ) %>% 
  rename(Year = Subgroup) %>% 
  dplyr::select(-Numerator, -Denominator) %>% 
  mutate(Year = as.numeric(substr(Year, 7, 10)))


final_oesophageal_gold_ncras <- bind_rows(
  oxford_ncras,
  df_wide ,
  df_widef,
  df_widem
)

# plots


lancet_colors <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF", "#FDAF17FF")
lancet_colors <- c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")


incidenceFigureData <- final_oesophageal_gold_ncras %>%
  filter(type == "Age Std") %>% 
  filter(Sex == "Both") %>% 
  ggplot(aes(x = Year,
             y = Rate,
             group = Database, colour = Database )) +
  geom_line(size = 1) +  # Thicker lines for clarity
  scale_color_manual(values = lancet_colors) +  # Apply Lancet-style colors to lines
  scale_fill_manual(values = lancet_colors) +   # Apply Lancet-style colors to ribbons
  geom_ribbon(aes(ymin = LCL, 
                  ymax = UCL, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  #geom_point(aes(shape = Database, fill = Database),size = 2) +
  #scale_shape_manual(values = c(24,21)) +
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
       col = "Database",
       fill = "Database") +
  facet_wrap(~ region, scales = "free_y")

incidenceFigureData

plotname <- paste0("/FIGURE_IRs_oesophageal_region_ageadjusted_overall.pdf")

pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()



# males

incidenceFigureData <- final_oesophageal_gold_ncras %>%
  filter(type == "Age Std") %>% 
  filter(Sex == "Male") %>% 
  ggplot(aes(x = Year,
             y = Rate,
             group = Database, colour = Database )) +
  geom_line(size = 1) +  # Thicker lines for clarity
  scale_color_manual(values = lancet_colors) +  # Apply Lancet-style colors to lines
  scale_fill_manual(values = lancet_colors) +   # Apply Lancet-style colors to ribbons
  geom_ribbon(aes(ymin = LCL, 
                  ymax = UCL, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  #geom_point(aes(shape = Database, fill = Database),size = 2) +
  #scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  # geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(title = "Male",
       x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Database",
       fill = "Database") +
  facet_wrap(~ region, scales = "free_y")

incidenceFigureData

plotname <- paste0("/FIGURE_IRs_oesophageal_region_ageadjusted_males.pdf")

pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()






# females

incidenceFigureData <- final_oesophageal_gold_ncras %>%
  filter(type == "Age Std") %>% 
  filter(Sex == "Female") %>% 
  ggplot(aes(x = Year,
             y = Rate,
             group = Database, colour = Database )) +
  geom_line(size = 1) +  # Thicker lines for clarity
  scale_color_manual(values = lancet_colors) +  # Apply Lancet-style colors to lines
  scale_fill_manual(values = lancet_colors) +   # Apply Lancet-style colors to ribbons
  geom_ribbon(aes(ymin = LCL, 
                  ymax = UCL, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  #geom_point(aes(shape = Database, fill = Database),size = 2) +
  #scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  # geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(title = "Female",
    x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Database",
       fill = "Database") +
  facet_wrap(~ region, scales = "free_y")

incidenceFigureData

plotname <- paste0("/FIGURE_IRs_oesophageal_region_ageadjusted_female.pdf")

pdf(paste0(datapath , plotname),
    width = 8, height = 7)
print(incidenceFigureData, newpage = FALSE)

dev.off()
