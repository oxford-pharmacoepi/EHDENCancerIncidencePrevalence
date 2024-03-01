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

# whole pop
survival_estimates_whole <- readRDS(here("data","survival_estimates.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp) %>%
  filter(CalendarYearGp == "2000 to 2019" | CalendarYearGp == "2000 to 2021" ) %>%
  droplevels() %>% 
  rename(Sex = Gender)

# survival_risk_table <- readRDS(here("data","survival_risk_table.rds")) %>%
#   rename(CalendarYearGp = CalenderYearGp)

survival_estimates_whole_GOLD <- survival_estimates_whole %>% 
  filter(Database == "CPRD GOLD")

survival_estimates_whole_Aurum <- survival_estimates_whole %>% 
  filter(Database == "CPRD Aurum")


process_data <- function(data) {
  
  # Aggregate data
  df_events <- aggregate(n.event ~ time_1, data = data, FUN = sum)
  df_censor <- aggregate(n.censor ~ time_1, data = data, FUN = sum)
  df_risk <- data %>% 
    group_by(time_1) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(time_1, n.risk, Age, Sex, Database)
  
  # Merge aggregated data back to original dataframe
  overall <- df_events %>% 
    left_join(df_censor, by = "time_1") %>%
    left_join(df_risk, by = "time_1") %>% 
    relocate(n.risk, .before = n.event)
  
  return(overall)
}

breaks <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, Inf)

# Create labels for the breaks
labels <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)

################################################
# GOLD
# Initialize a list to store overall data for each cancer
database_temp_overall <- list()

# Loop over each cancer
for (i in 1:length(table(survival_estimates_whole$Cancer))) { 
  
  # Filter data for the current cancer
  database_temp <- survival_estimates_whole_GOLD %>% 
    filter(Cancer == unique(survival_estimates_whole_GOLD$Cancer)[i]) %>% 
    mutate(time_1 = cut(time, breaks = breaks, labels = labels, include.lowest = TRUE))
  
  # Filter data for each sex
  database_temp_sex <- database_temp %>% 
    filter(Stratification == "Gender" )
  
  #filter age
  database_temp_age <- database_temp %>% 
    filter(Stratification == "Age" & Age != "All") %>%
    droplevels()
  
  #filter age*sex
  database_temp_age_sex <- database_temp %>% 
    filter(Stratification == "Age*Gender" ) 
  
  
  # Filter for overall
  overall_data <- database_temp %>%
    filter(Stratification == "None") %>%
    process_data() %>%
    mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>%
    dplyr::select(c(n.risk, n.event, n.censor)) %>%
    t() %>%
    as_tibble() %>%
    `colnames<-`(labels) %>%
    dplyr::mutate(Cancer = unique(survival_estimates_whole_GOLD$Cancer)[i],
                  Sex = "Both" ,
                  Age = "All",
                  details = c("n.risk", "n.event", "n.censor")) %>%
    dplyr::relocate(details)
  
  # Store the overall data for the current cancer in the list
  database_temp_overall[[i]] <- list(overall = overall_data)
  
  #################
  # sex
  #################
  # Initialize a list to store data for each sex
  database_temp_sex_list <- list()
  
  # for each sex
  for (j in 1:length(table(database_temp_sex$Sex))) { 
    tryCatch({
      # Filter data for the current sex
      sex_data <- database_temp_sex %>% 
        filter(Sex == names(table(database_temp_sex$Sex)[j])) %>% 
        process_data() %>% 
        mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>% 
        dplyr::select(c(n.risk, n.event, n.censor)) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(labels) %>% 
        dplyr::mutate(Cancer = unique(survival_estimates_whole_GOLD$Cancer)[i],
                      Sex = names(table(database_temp_sex$Sex)[j]) ,
                      Age = "All",
                      details = c("n.risk", "n.event", "n.censor")) %>% 
        dplyr::relocate(details)
      
      # Store the data for the current sex in the list
      database_temp_sex_list[[j]] <- sex_data
    }, error = function(e) {
      # Print error message if an error occurs
      print(paste("Error occurred for sex:", names(table(database_temp_sex$Sex)[j])))
      print(e)
    })
  }
  
  # Store the list of data for each sex in the overall list
  database_temp_overall[[i]]$sex_data <- database_temp_sex_list
  
  ####################
  # age 
  ###################
  # Initialize a list to store data for each age
  database_temp_age_list <- list()
  
  # for each sex
  for (a in 1:length(table(database_temp_age$Age))) { 
    
    tryCatch({
      # Filter data for the current sex
      age_data <- database_temp_age %>% 
        filter(Age == names(table(database_temp_age$Age)[a])) %>% 
        process_data() %>% 
        mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>% 
        dplyr::select(c(n.risk, n.event, n.censor)) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(labels) %>% 
        dplyr::mutate(Cancer = unique(survival_estimates_whole_GOLD$Cancer)[i],
                      Sex = "Both" ,
                      Age = names(table(database_temp_age$Age))[a],
                      details = c("n.risk", "n.event", "n.censor")) %>% 
        dplyr::relocate(details)
      
      # Store the data for the current sex in the list
      database_temp_age_list[[a]] <- age_data
    }, error = function(e) {
      # Print error message if an error occurs
      print(paste("Error occurred for age:", names(table(database_temp_age$Age))[a]))
      print(e)
    })
  }
  
  # Store the list of data for each sex in the overall list
  database_temp_overall[[i]]$age_data <- database_temp_age_list
  
  #######################
  # age * sex 
  #######################
  # Initialize a list to store data for each age*sex
  database_temp_age_sex_list <- list()
  
  # for each sex*age
  for (asex in 1:length(table(database_temp_age_sex$GenderAge))) { 
    
    tryCatch({
      # Filter data for the current sex
      age_sex_data <- database_temp_age_sex %>% 
        filter(GenderAge == names(table(database_temp_age_sex$GenderAge)[asex])) %>% 
        process_data() %>% 
        mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>% 
        dplyr::select(c(n.risk, n.event, n.censor)) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(labels) %>% 
        dplyr::mutate(Cancer = unique(survival_estimates_whole_GOLD$Cancer)[i],
                      GenderAge = names(table(database_temp_age_sex$GenderAge))[asex],
                      details = c("n.risk", "n.event", "n.censor")) %>% 
        dplyr::relocate(details) %>% 
        tidyr::separate(col = "GenderAge",
                        into = c("Sex", "Age"),
                        sep = "_",
                        remove = T)
      
      # Store the data for the current sex in the list
      database_temp_age_sex_list[[asex]] <- age_sex_data
    }, error = function(e) {
      # Print error message if an error occurs
      print(paste("Error occurred for age sex:", names(table(database_temp_age_sex$GenderAge))[asex]))
      print(e)
    })
  }
  
  # Store the list of data for each sex in the overall list
  database_temp_overall[[i]]$age_sex_data <- database_temp_age_sex_list
  
}

final_df <- map_dfr(database_temp_overall, bind_rows) %>% 
  mutate(Database = "CPRD GOLD",
         CalendarYearGp = "2000 to 2021")


######################################################
# AURUM
# Initialize a list to store overall data for each cancer
database_temp_overall_aurum <- list()

# Loop over each cancer
for (i in 1:length(table(survival_estimates_whole$Cancer))) { 
  
  # Filter data for the current cancer
  database_temp <- survival_estimates_whole_Aurum %>% 
    filter(Cancer == unique(survival_estimates_whole_Aurum$Cancer)[i]) %>% 
    mutate(time_1 = cut(time, breaks = breaks, labels = labels, include.lowest = TRUE))
  
  # Filter data for each sex
  database_temp_sex <- database_temp %>% 
    filter(Stratification == "Gender" )
  
  #filter age
  database_temp_age <- database_temp %>% 
    filter(Stratification == "Age" & Age != "All") %>%
    droplevels()
  
  #filter age*sex
  database_temp_age_sex <- database_temp %>% 
    filter(Stratification == "Age*Gender" ) 
  
  
  # Filter for overall
  overall_data <- database_temp %>%
    filter(Stratification == "None") %>%
    process_data() %>%
    mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>%
    dplyr::select(c(n.risk, n.event, n.censor)) %>%
    t() %>%
    as_tibble() %>%
    `colnames<-`(labels) %>%
    dplyr::mutate(Cancer = unique(survival_estimates_whole_Aurum$Cancer)[i],
                  Sex = "Both" ,
                  Age = "All",
                  details = c("n.risk", "n.event", "n.censor")) %>%
    dplyr::relocate(details)
  
  # Store the overall data for the current cancer in the list
  database_temp_overall_aurum[[i]] <- list(overall = overall_data)
  
  #################
  # sex
  #################
  # Initialize a list to store data for each sex
  database_temp_sex_list <- list()
  
  # for each sex
  for (j in 1:length(table(database_temp_sex$Sex))) { 
    tryCatch({
      # Filter data for the current sex
      sex_data <- database_temp_sex %>% 
        filter(Sex == names(table(database_temp_sex$Sex)[j])) %>% 
        process_data() %>% 
        mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>% 
        dplyr::select(c(n.risk, n.event, n.censor)) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(labels) %>% 
        dplyr::mutate(Cancer = unique(survival_estimates_whole_Aurum$Cancer)[i],
                      Sex = names(table(database_temp_sex$Sex)[j]) ,
                      Age = "All",
                      details = c("n.risk", "n.event", "n.censor")) %>% 
        dplyr::relocate(details)
      
      # Store the data for the current sex in the list
      database_temp_sex_list[[j]] <- sex_data
    }, error = function(e) {
      # Print error message if an error occurs
      print(paste("Error occurred for sex:", names(table(database_temp_sex$Sex)[j])))
      print(e)
    })
  }
  
  # Store the list of data for each sex in the overall list
  database_temp_overall_aurum[[i]]$sex_data <- database_temp_sex_list
  
  ####################
  # age 
  ###################
  # Initialize a list to store data for each age
  database_temp_age_list <- list()
  
  # for each sex
  for (a in 1:length(table(database_temp_age$Age))) { 
    
    tryCatch({
      # Filter data for the current sex
      age_data <- database_temp_age %>% 
        filter(Age == names(table(database_temp_age$Age)[a])) %>% 
        process_data() %>% 
        mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>% 
        dplyr::select(c(n.risk, n.event, n.censor)) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(labels) %>% 
        dplyr::mutate(Cancer = unique(survival_estimates_whole_Aurum$Cancer)[i],
                      Sex = "Both" ,
                      Age = names(table(database_temp_age$Age))[a],
                      details = c("n.risk", "n.event", "n.censor")) %>% 
        dplyr::relocate(details)
      
      # Store the data for the current sex in the list
      database_temp_age_list[[a]] <- age_data
    }, error = function(e) {
      # Print error message if an error occurs
      print(paste("Error occurred for age:", names(table(database_temp_age$Age))[a]))
      print(e)
    })
  }
  
  # Store the list of data for each sex in the overall list
  database_temp_overall_aurum[[i]]$age_data <- database_temp_age_list
  
  #######################
  # age * sex 
  #######################
  # Initialize a list to store data for each age*sex
  database_temp_age_sex_list <- list()
  
  # for each sex*age
  for (asex in 1:length(table(database_temp_age_sex$GenderAge))) { 
    
    tryCatch({
      # Filter data for the current sex
      age_sex_data <- database_temp_age_sex %>% 
        filter(GenderAge == names(table(database_temp_age_sex$GenderAge)[asex])) %>% 
        process_data() %>% 
        mutate_at(vars(starts_with("n.")), ~ ifelse(. < 5, "<5", as.character(.))) %>% 
        dplyr::select(c(n.risk, n.event, n.censor)) %>% 
        t() %>% 
        as_tibble() %>% 
        `colnames<-`(labels) %>% 
        dplyr::mutate(Cancer = unique(survival_estimates_whole_Aurum$Cancer)[i],
                      GenderAge = names(table(database_temp_age_sex$GenderAge))[asex],
                      details = c("n.risk", "n.event", "n.censor")) %>% 
        dplyr::relocate(details) %>% 
        tidyr::separate(col = "GenderAge",
                        into = c("Sex", "Age"),
                        sep = "_",
                        remove = T)
      
      # Store the data for the current sex in the list
      database_temp_age_sex_list[[asex]] <- age_sex_data
    }, error = function(e) {
      # Print error message if an error occurs
      print(paste("Error occurred for age sex:", names(table(database_temp_age_sex$GenderAge))[asex]))
      print(e)
    })
  }
  
  # Store the list of data for each sex in the overall list
  database_temp_overall_aurum[[i]]$age_sex_data <- database_temp_age_sex_list
  
}

final_df1 <- map_dfr(database_temp_overall_aurum, bind_rows) %>% 
  mutate(Database = "CPRD Aurum",
         CalendarYearGp = "2000 to 2019")

survival_risk_table <- bind_rows(final_df, final_df1) %>% 
  filter(details != "n.censor")

survival_median_table <- readRDS(here("data","survival_median_table.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp) %>%
  filter(CalendarYearGp == "2000 to 2019" | CalendarYearGp == "2000 to 2021" ) %>%
  droplevels()

survival_rates_table <- readRDS(here("data","survival_rates_table.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp) %>%
  filter(CalendarYearGp == "2000 to 2019" | CalendarYearGp == "2000 to 2021" ) %>%
  droplevels() 

#calendar pop
survival_estimates_calendar <- readRDS(here("data","survival_estimates.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp) %>%
  filter(CalendarYearGp != "2000 to 2019") %>% 
  filter(CalendarYearGp != "2000 to 2021" ) %>%
  droplevels() %>% 
  rename(Sex = Gender)

survival_risk_table_cy <- readRDS(here("data","survival_risk_table_cy.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp)

survival_median_table_cy <- readRDS(here("data","survival_median_table.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp) %>%
  filter(CalendarYearGp != "2000 to 2019") %>% 
  filter(CalendarYearGp != "2000 to 2021" ) %>%
  droplevels()

survival_rates_table_cy <- readRDS(here("data","survival_rates_table.rds")) %>%
  rename(CalendarYearGp = CalenderYearGp) %>%
  filter(CalendarYearGp != "2000 to 2019") %>% 
  filter(CalendarYearGp != "2000 to 2021" ) %>%
  droplevels() %>%
  filter(time != 10) %>%
  filter( !(time == 5 & CalendarYearGp == "2020 to 2021" )) %>% 
  mutate(time = as.character(time)) %>% 
  rename(Sex = Gender)

table_one_results <- readRDS(here("data","table1_results.rds")) %>%
  filter(analysis == "Incidence") %>% 
  filter(!grepl("Prior_history_years",var)) 
  
survival_followup_table <- readRDS(here("data","survival_median_mean_follow_up.rds")) 



