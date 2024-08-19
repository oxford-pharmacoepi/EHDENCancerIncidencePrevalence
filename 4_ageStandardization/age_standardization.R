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

#dsr package needs a specific version to work 
# install.packages("frailtypack")
# devtools::install_github("socale/frailtypack", ref = "master")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/dsr/dsr_0.2.2.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")


#folder of data
datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/4_ageStandardization/data"
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
    mutate(database_name = replace(database_name, database_name == "CPRDGoldUpdate2", "CPRD GOLD"))

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
        paste0(datapath, "/prevalence_estimates.rds"))

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
prevalence_estimates <- readRDS(paste0(datapath ,"/prevalence_estimates.rds"))
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



# try with one cancer : colorectal by hand without packages
# incidence_estimates_CRC <- incidence_estimates %>% 
#   filter(outcome_cohort_name == "Colorectal") %>% 
#   filter(analysis_interval != "overall") %>% 
#   filter(denominator_sex == "Both") %>% 
#   filter(database_name == "CPRD GOLD") %>% 
#   filter(denominator_age_group != "All") %>% 
#   filter(incidence_start_date == "2006-01-01") %>% 
#   rename(Agegroup = denominator_age_group)
#   
# 
# # merge the two files
# esp_CRC <- incidence_estimates_CRC %>% 
#   left_join(ESP13_updated, by = "Agegroup") %>% 
#   select(n_persons,
#          n_events,
#          incidence_start_date,
#          incidence_100000_pys,
#          incidence_100000_pys_95CI_lower,
#          incidence_100000_pys_95CI_upper,
#          ESP2013) %>% 
#   mutate(espbyir = incidence_100000_pys* ESP2013)
# 
# ESP2013_sum <- sum(esp_CRC$ESP2013)
# irage_summ <- sum(esp_CRC$espbyir)
#   
# #age standardize = IR sum/esr2013 sum
# age_adjust_result <- irage_summ/ESP2013_sum

########################################
# Age standardize using dsr package

#rename ESP column to pop (needs to be pop otherwise will not work)
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

#grab breast cancer results
# incidence_estimates2 <- incidence_estimates %>% 
#   filter(denominator_sex == "Female") %>% 
#   filter(outcome_cohort_name == "Breast" ) %>% 
#   filter(analysis_interval != "overall") %>% 
#   filter(database_name == "CPRD GOLD") %>% 
#   filter(denominator_age_group != "All") %>% 
#   rename(Agegroup = denominator_age_group) #rename this to make it the same as agegroup in standard population

#grab prostate cancer results
incidence_estimates3 <- incidence_estimates %>% 
  filter(outcome_cohort_name == "Prostate" ) %>% 
  filter(analysis_interval != "overall") %>% 
  filter(database_name == "CPRD GOLD") %>% 
  filter(denominator_age_group != "All") %>% 
  rename(Agegroup = denominator_age_group) #rename this to make it the same as agegroup in standard population

#merge the results together
incidence_estimates4 <- bind_rows(incidence_estimates1 #,
                                  #incidence_estimates2,
                                 # incidence_estimates3
                                  )




#create a loop for each cancer (all cancers apart from prostate and breast are for single genders)
agestandardizedinc <- list()
agestandardizedincF <- list()
agestandardizedincM <- list()

for(i in 1:length(table(incidence_estimates4$outcome_cohort_name))){
  

# for whole population    
incidence_estimates_i <- incidence_estimates4 %>%
    filter(outcome_cohort_name == names(table(incidence_estimates4$outcome_cohort_name)[i]),
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
  mutate(Cancer = names(table(incidence_estimates4$outcome_cohort_name)[i]),
         Sex = "Both") 

}, error = function(e) {
  # If an error occurs, print the error message
  cat("An error occurred:", conditionMessage(e), "\n")
})


print(paste0("age standardization for ", names(table(incidence_estimates4$outcome_cohort_name)[i]), "for both sexes done"))

##########################
# for female population    
#########################
incidence_estimates_i <- incidence_estimates4 %>%
  filter(outcome_cohort_name == names(table(incidence_estimates4$outcome_cohort_name)[i]),
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
    mutate(Cancer = names(table(incidence_estimates4$outcome_cohort_name)[i]),
           Sex = "Female") 
  
}, error = function(e) {
  # If an error occurs, print the error message
  cat("An error occurred:", conditionMessage(e), "\n")
})


print(paste0("age standardization for ", names(table(incidence_estimates4$outcome_cohort_name)[i]), "for females done"))




##########################
# for male population    
#########################
incidence_estimates_i <- incidence_estimates4 %>%
  filter(outcome_cohort_name == names(table(incidence_estimates4$outcome_cohort_name)[i]),
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
    mutate(Cancer = names(table(incidence_estimates4$outcome_cohort_name)[i]),
           Sex = "Male") 
  
}, error = function(e) {
  # If an error occurs, print the error message
  cat("An error occurred:", conditionMessage(e), "\n")
})


print(paste0("age standardization for ", names(table(incidence_estimates4$outcome_cohort_name)[i]), "for males done"))





}

agestandardizedinc_final <- bind_rows(agestandardizedinc,
                                      agestandardizedincF,
                                      agestandardizedincM)

#save the results
saveRDS(agestandardizedinc_final, paste0(datapath ,"/incidence_estimates_age_sd.rds"))


agestandardizedinc_final <- readRDS(paste0(datapath ,"/incidence_estimates_age_sd.rds"))

###################################################
#plot the results of new age adjusted results INCIDENCE
# all


incidenceFigureData <- agestandardizedinc_final %>%
  filter(Sex == "Both") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = database_name )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.position='none') +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)


plotname <- paste0("/FIGURE1_IRsWholePop_multipleCancers_ageadjusted.pdf")


pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData, newpage = FALSE)

dev.off()

# females 
incidenceFigureData1 <- agestandardizedinc_final %>%
  filter(Sex == "Female") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = database_name )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.position='none') +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(title = "Females" ,
    x = "Calendar year",
    y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)

plotname <- paste0("/FIGURE1_IRsFEMALESPop_multipleCancers_ageadjusted.pdf")
pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData1, newpage = FALSE)

dev.off()


# males 
incidenceFigureData2 <- agestandardizedinc_final %>%
  filter(Sex == "Male") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = database_name )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.position='none') +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       title = "Males" ,
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)

plotname <- paste0("/FIGURE1_IRsMALESPop_multipleCancers_ageadjusted.pdf")
pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData2, newpage = FALSE)

dev.off()


# png(paste0(datapath , plotname),
#     width = 8, height = 7.5, units = "in", res = 1200)
# print(incidenceFigureData, newpage = FALSE)

pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData, newpage = FALSE)

dev.off()



# for lung cancer males and females
incidenceFigureData4 <- agestandardizedinc_final %>%
  filter(Sex != "Both" ) %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Sex )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Sex), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Sex, fill = Sex),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Sex",
       shape = "Sex",
       fill = "Sex") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)

plotname <- paste0("/FIGURE1_IRsSEXPop_multipleCancers_ageadjusted.pdf")
pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData4, newpage = FALSE)
dev.off()


# lung cancer only
incidenceFigureData5 <- agestandardizedinc_final %>%
  filter(Sex != "Both" ,
         Cancer == "Lung") %>% 
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Sex )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Sex), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Sex, fill = Sex),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Sex",
       shape = "Sex",
       fill = "Sex") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)

plotname <- paste0("/FIGURE1_IRsSEXPop_lungCancers_ageadjusted.pdf")
pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData5, newpage = FALSE)

plotname <- paste0("/FIGURE1_IRsSEXPop_lungCancers_ageadjusted.tiff")
tiff(paste0(datapath , plotname),
     width = 6, height = 5 , units = "in", res = 1200)
print(incidenceFigureData5, newpage = FALSE)
dev.off()




# do a loop for each cancer to do a plot for age standardize

for(i in 1:length(table(incidence_estimates4$outcome_cohort_name))){

  agestandardizedinc_final_i <- agestandardizedinc_final %>%
    filter(Sex != "Both" ,
           Cancer == names(table(incidence_estimates4$outcome_cohort_name)[i]))
    

incidenceFigureData5 <- agestandardizedinc_final_i %>%         
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Sex )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Sex), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Sex, fill = Sex),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Sex",
       shape = "Sex",
       fill = "Sex") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)


plotname <- paste0("/FIGUREX_", names(table(incidence_estimates4$outcome_cohort_name)[i]),"_IRsSEX_ageadjusted.tiff")
tiff(paste0(datapath , plotname),
     width = 6, height = 5 , units = "in", res = 1200)
print(incidenceFigureData5, newpage = FALSE)
dev.off()

}


# for all head and nceck cancers combined

agestandardizedinc_final_han <- agestandardizedinc_final %>%
  filter(Sex != "Both" ,
         Cancer == "Head & Neck" |
           Cancer == "Hypopharynx" |
           Cancer == "Larynx" |
           Cancer == "Nasal Cavity & Sinus" |
           Cancer == "Nasopharynx" |
           Cancer == "Oral Cavity" |
           Cancer == "Oropharynx" |
           Cancer == "Salivary Gland" |
           Cancer == "Tongue"  )


incidenceFigureData5 <- agestandardizedinc_final_han %>%         
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = Sex )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#ED0000FF", "#00468BFF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = Sex), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = Sex, fill = Sex),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed")) +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Sex",
       shape = "Sex",
       fill = "Sex") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)


plotname <- paste0("/FIGUREX_HanSubtypes_IRsSEX_ageadjusted.tiff")
tiff(paste0(datapath , plotname),
     width = 9, height = 7 , units = "in", res = 1200)
print(incidenceFigureData5, newpage = FALSE)
dev.off()

plotname <- paste0("/FIGUREX_HanSubtypes_IRsSEX_ageadjusted.png")
png(paste0(datapath , plotname),
     width = 9, height = 7 , units = "in", res = 1200)
print(incidenceFigureData5, newpage = FALSE)
dev.off()



# males and females sep for breast cancer

incidenceFigureData1 <- agestandardizedinc_final %>%
  filter(Sex == "Female",
         Cancer == "Breast") %>%
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = database_name )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.position='none') +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(title = "Females" ,
       x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)

plotname <- paste0("/Figure_age_std_IRs_female_breast.pdf")
pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData1, newpage = FALSE)

dev.off()


incidenceFigureData2 <- agestandardizedinc_final %>%
  filter(Sex == "Male",
         Cancer == "Breast") %>%
  mutate(database_name = "CPRD GOLD") %>% 
  ggplot(aes(x = Subgroup,
             y = `Std Rate (per 1e+05)`,
             group = database_name )) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = `95% LCL (Std)`, 
                  ymax = `95% UCL (Std)`, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 2) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        #panel.spacing.x = unit(0.1,"line"),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.position='none') +
  geom_vline(xintercept = as.numeric(as.Date("2004-04-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
  labs(title = "Males" ,
       x = "Calendar year",
       y = "Age Standardized Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ Cancer, scales = "free_y", ncol = 3)

plotname <- paste0("/Figure_age_std_IRs_male_breast.pdf")
pdf(paste0(datapath , plotname),
    width = 10, height = 9.5)
print(incidenceFigureData2, newpage = FALSE)

dev.off()



################################## Prevalence #######
#plot for all cancer IRs on one facet
# prevalenceFigureData <- agestandardizedprev_final %>%
#   mutate(database_name = "CPRD GOLD") %>% 
#   ggplot(aes(x = Subgroup,
#              y = `Std Rate (per 1)`,
#              group = database_name )) +
#   geom_line(color = "black", size = 0.25) +
#   scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
#   scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
#   geom_ribbon(aes(ymin = `95% LCL (Std)`,
#                   ymax = `95% UCL (Std)`,
#                   fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
#   geom_point(aes(shape = database_name, fill = database_name),size = 2) +
#   scale_shape_manual(values = c(24,21)) +
#   scale_y_continuous( labels = scales::percent, limits = c(0, NA)) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1),
#         panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
#         strip.background = element_rect(color = "black", size = 0.6) ,
#         panel.background = element_blank() ,
#         axis.line = element_line(colour = "black", size = 0.6) ,
#         #panel.spacing.x = unit(0.1,"line"),
#         panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
#         legend.position='none') +
#   geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), linetype="dotted", colour = "#ED0000FF", size = 0.8) +
#   labs(x = "Calendar year",
#        y = "Prevalence",
#        col = "Database name",
#        shape = "Database name",
#        fill = "Database name") +
#   scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
#                expand = c(0.06,1)) +
#   facet_wrap(~ Cancer, scales = "free_y", ncol = 3)
# 
# 
# plotname <- paste0("/PPsWholePop_multipleCancers_age_adjusted.pdf")
# 
# pdf(paste0(datapath , plotname),
#     width = 10, height = 9.5)
# print(prevalenceFigureData, newpage = FALSE)
# 
# dev.off()









