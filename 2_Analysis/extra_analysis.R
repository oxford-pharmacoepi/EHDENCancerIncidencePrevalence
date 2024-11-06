# reviewer revisions

# read in lung cancer cohort
lung_cancer_outcome_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "OutcomeCohorts",
  "extra_analysis"
))

lung_cancer_prevalent_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "PrevalentCohorts",
  "extra_analysis"
  ))

smoking_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "OutcomeCohorts",
  "extra_analysis" ,
  "smoking"
))


#tablenames
lc_outcome_table_name <-paste0(outcome_table_stem,"_olc") # for incidence
lc_prevalent_table_name <-paste0(outcome_table_stem,"_plc") # for prevalence
lc_smoking_table_name <- paste0(outcome_table_stem,"_smoke")

instantiatedCohorts <- FALSE
if (instantiatedCohorts == TRUE) {
  
  cdm <- CDMConnector::cdm_from_con(con = db,
                                    cdm_schema = cdm_database_schema,
                                    write_schema = results_database_schema,
                                    cohort_tables = c(lc_outcome_table_name,
                                                      lc_prevalent_table_name))
  
  
} else {
#instantiate the lung cancer cohorts
# incidence
cdm <- CDMConnector::generateCohortSet(cdm, 
                                       cohortSet = lung_cancer_outcome_cohorts,
                                       name = lc_outcome_table_name,
                                       overwrite = TRUE
)

# prevalence
cdm <- CDMConnector::generateCohortSet(cdm = cdm, 
                                       cohortSet = lung_cancer_prevalent_cohorts,
                                       name = lc_prevalent_table_name,
                                       overwrite = TRUE
                                       
)


# smoking
cdm <- CDMConnector::generateCohortSet(cdm = cdm, 
                                       cohortSet = smoking_cohorts,
                                       name = lc_smoking_table_name,
                                       overwrite = TRUE
                                       
)


}


# take denominator (should already be run)
cdm$denominator <- generateDenominatorCohortSet(
  cdm = cdm,
  startDate = as.Date(studyStartDate),
  endDate = as.Date(studyEndDate),
  ageGroup =list(
    c(18, 150),
    c(18, 29),
    c(30, 39),
    c(40, 49),
    c(50, 59),
    c(60, 69),
    c(70, 79),
    c(80, 89),
    c(90, 150)
  ),
  sex = c("Male", "Female", "Both"),
  daysPriorHistory = 365,
  verbose = TRUE
)


#############################################################################
# stratification by location (wales, scotland, england, NI)

# create a denominator with region added in
cdm$denominator <- cdm$denominator %>%
  left_join(cdm$person %>%
              rename("subject_id"= "person_id") %>%
              select("subject_id", "care_site_id"),
            by = "subject_id") %>%
  # add location_id id
  left_join(cdm$care_site %>%
              select("care_site_id", "location_id"),
            by = "care_site_id") %>%
  # add region
  left_join(cdm$location %>%
              rename("region" = "location_source_value") %>%
              select(c("location_id", "region")),
            by = "location_id") %>% 
  mutate(region_collapsed = case_when(region %in% c("Yorkshire  & The Humber", "East Midlands", 
                                                  "West Midlands", "North East", "North West", "East of England", "London", 
                                                  "South East", "South West") ~ "England",
                                    region == "Northern Ireland" ~ "Northern Ireland",
                                    region == "Scotland" ~ "Scotland",
                                    region == "Wales" ~ "Wales"))


# create a denominator for each region
cdm$denominator_england <- cdm$denominator %>% 
  filter(region_collapsed == "England") 

cdm$denominator_scotland <- cdm$denominator %>% 
  filter(region_collapsed == "Scotland") 

cdm$denominator_wales <- cdm$denominator %>% 
  filter(region_collapsed == "Wales") 

cdm$denominator_ni <- cdm$denominator %>% 
  filter(region_collapsed == "Northern Ireland") 


# run incidence and prevalence per region
# ENGLAND
inc_eng <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_england",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

inc_eng_overall <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_england",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "overall", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

#prev
prev_period_eng <- estimatePeriodPrevalence(
  cdm = cdm,
  denominatorTable = "denominator_england",
  outcomeCohortId = lung_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = lc_prevalent_table_name,
  interval = "years" ,
  completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
  fullContribution = FALSE , # individuals only required to be present for one day in interval
  minCellCount = 5,
  verbose = TRUE
)

#get results and export
study_results_eng <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                  resultList=list(inc_eng,prev_period_eng ),
                                                  databaseName = db.name)

study_results_eng1 <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                  resultList=list(inc_eng_overall),
                                                  databaseName = db.name)

exportIncidencePrevalenceResults(result=study_results_eng,
                                 zipName= paste0(db.name, "IPResults_england"),
                                 outputFolder=here::here("Results", db.name))

exportIncidencePrevalenceResults(result=study_results_eng1,
                                 zipName= paste0(db.name, "IPResults_england_overall"),
                                 outputFolder=here::here("Results", db.name))


# SCOTLAND
inc_scot <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_scotland",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

inc_scot_overall <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_scotland",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "overall", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

#prev
prev_period_scot <- estimatePeriodPrevalence(
  cdm = cdm,
  denominatorTable = "denominator_scotland",
  outcomeCohortId = lung_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = lc_prevalent_table_name,
  interval = "years" ,
  completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
  fullContribution = FALSE , # individuals only required to be present for one day in interval
  minCellCount = 5,
  verbose = TRUE
)


#get results and export
study_results_scot <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                      resultList=list(inc_scot,prev_period_scot ),
                                                      databaseName = db.name)

study_results_scot1 <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                       resultList=list(inc_scot_overall),
                                                       databaseName = db.name)

exportIncidencePrevalenceResults(result=study_results_scot,
                                 zipName= paste0(db.name, "IPResults_scotland"),
                                 outputFolder=here::here("Results", db.name))

exportIncidencePrevalenceResults(result=study_results_scot1,
                                 zipName= paste0(db.name, "IPResults_scotland_overall"),
                                 outputFolder=here::here("Results", db.name))


# WALES
inc_wales <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_wales",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

inc_wales_overall <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_wales",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "overall", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

#prev
prev_period_wales <- estimatePeriodPrevalence(
  cdm = cdm,
  denominatorTable = "denominator_wales",
  outcomeCohortId = lung_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = lc_prevalent_table_name,
  interval = "years" ,
  completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
  fullContribution = FALSE , # individuals only required to be present for one day in interval
  minCellCount = 5,
  verbose = TRUE
)


#get results and export
study_results_wales <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                       resultList=list(inc_wales,prev_period_wales ),
                                                       databaseName = db.name)

study_results_wales1 <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                        resultList=list(inc_wales_overall),
                                                        databaseName = db.name)

exportIncidencePrevalenceResults(result=study_results_wales,
                                 zipName= paste0(db.name, "IPResults_wales"),
                                 outputFolder=here::here("Results", db.name))

exportIncidencePrevalenceResults(result=study_results_wales1,
                                 zipName= paste0(db.name, "IPResults_wales_overall"),
                                 outputFolder=here::here("Results", db.name))

# NI
inc_ni <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_ni",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

inc_ni_overall <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_ni",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "overall", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

#prev
prev_period_ni <- estimatePeriodPrevalence(
  cdm = cdm,
  denominatorTable = "denominator_ni",
  outcomeCohortId = lung_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = lc_prevalent_table_name,
  interval = "years" ,
  completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
  fullContribution = FALSE , # individuals only required to be present for one day in interval
  minCellCount = 5,
  verbose = TRUE
)

#get results and export
study_results_ni <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                      resultList=list(inc_ni,prev_period_ni ),
                                                      databaseName = db.name)

study_results_ni1 <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                       resultList=list(inc_ni_overall),
                                                       databaseName = db.name)

exportIncidencePrevalenceResults(result=study_results_ni,
                                 zipName= paste0(db.name, "IPResults_ni"),
                                 outputFolder=here::here("Results", db.name))

exportIncidencePrevalenceResults(result=study_results_ni1,
                                 zipName= paste0(db.name, "IPResults_ni_overall"),
                                 outputFolder=here::here("Results", db.name))


###########################################################################

# smoking
# make distinct cases
# filter out those whose smoking status first event is pass the end of the study period (i.e. first smoking record > 31/12/21)
cdm[[lc_smoking_table_name]] <- cdm[[lc_smoking_table_name]] %>% 
  distinct() %>% 
  filter(cohort_start_date < as.Date(studyEndDate))

# # add on column to the denominator with those with a code for smoking
# dont put distinct as different people contribute to multiple denominators
# only includes people who have ever smoking code before or at the cohort end date for the denominator - 
cdm$denominator_smokers <- cdm$denominator %>%
  inner_join(cdm[[lc_smoking_table_name]] %>%
               rename("smoking_date"= "cohort_start_date") %>%
               select(subject_id, smoking_date), 
               by=c("subject_id"), copy = TRUE) %>% 
  filter(smoking_date <= cohort_end_date )
  

# # stratification by smoking status (any time)
# # inc
inc_smoking <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_smokers",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "years",
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

inc_smoking_overall <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_smokers",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = "overall",
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 0,
  returnParticipants = FALSE,
  verbose = TRUE
)

#prev
prev_period_smoking <- estimatePeriodPrevalence(
  cdm = cdm,
  denominatorTable = "denominator_smokers",
  outcomeCohortId = lung_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0,
  outcomeTable = lc_prevalent_table_name,
  interval = "years" ,
  completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
  fullContribution = FALSE , # individuals only required to be present for one day in interval
  minCellCount = 5,
  verbose = TRUE
)

# gather results
study_results_smoking <- gatherIncidencePrevalenceResults(cdm =cdm,
                                                  resultList=list(inc_smoking,prev_period_smoking ),
                                                  databaseName = db.name)

study_results_smoking_overall<- gatherIncidencePrevalenceResults(cdm =cdm,
                                                         resultList=list(inc_smoking_overall),
                                                         databaseName = db.name)

# export results
exportIncidencePrevalenceResults(result=study_results_smoking,
                                 zipName= paste0(db.name, "IPResults_smoking"),
                                 outputFolder=here::here("Results", db.name))

exportIncidencePrevalenceResults(result=study_results_smoking_overall,
                                 zipName= paste0(db.name, "IPResults_overall_smoking"),
                                 outputFolder=here::here("Results", db.name))


###########################################################################
# table ones stratified by region and by smoking

#instanstiate conditions
info(logger, "- getting feature for diseases definitions")

disease_cohorts <- readCohortSet(here(
  "1_InstantiateCohorts",
  "ConditionFeatureCohorts"
))

info(logger, "- getting features: diseases")
  
  cdm <- generateCohortSet(cdm, 
                           cohortSet = disease_cohorts,
                           name = feature_disease_table_name,
                           overwrite = TRUE
  )
  
info(logger, "- got features for diseases")

#instantiate feature cohorts (disease)
info(logger, "- getting feature for diseases definitions")

disease_cohorts_liver <- readCohortSet(here(
  "1_InstantiateCohorts",
  "ConditionFeatureCohortsLiver"
))

info(logger, "- getting features: diseases liver related")
  
  cdm <- generateCohortSet(cdm, 
                           cohortSet = disease_cohorts_liver,
                           name = feature_disease_liver_table_name,
                           overwrite = TRUE
  )
  


info(logger, "- got features for diseases liver related")


# get denominator to get participants
cdm$denominatordemo <- generateDenominatorCohortSet(
  cdm = cdm,
  startDate = as.Date(studyStartDate),
  endDate = as.Date(studyEndDate),
  ageGroup =list(
    c(18, 150)),
  sex = c("Both"),
  daysPriorHistory = 365,
  verbose = TRUE
)

# get overall incidence 
inc_overall_participants <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominatordemo",
  outcomeTable = lc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
  interval = c("overall"),
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 5,
  returnParticipants = TRUE,
  tablePrefix = outcome_table_stem,
  verbose = TRUE
)

# #extract settings for survival from incidence results
settings_surv <- settings(inc_overall_participants) 

pops <- list()

for (i in 1:length(settings_surv$analysis_id)){
  #extract the participants for each cancer
  
  pops[[i]] <-cdm$person %>%
    inner_join(participants(inc_overall_participants, analysisId = as.numeric(settings_surv$analysis_id[i])) %>% filter(!is.na(outcome_start_date)),
               by = c("person_id" = "subject_id" ), copy = TRUE) %>%
    select(person_id,gender_concept_id,
           year_of_birth, month_of_birth, day_of_birth,
           cohort_start_date,
           cohort_end_date,
           outcome_start_date) %>%
    left_join(cdm$observation_period %>%
                select("person_id",  "observation_period_start_date", "observation_period_end_date") %>%
                distinct(),
              by = "person_id") %>%
    left_join(cdm$death %>%
                select("person_id",  "death_date") %>%
                distinct(),
              by = "person_id") %>%
    collect()
  
  pops[[i]] <- pops[[i]]  %>%
    mutate(outcome_cohort_name = settings_surv$outcome_cohort_name[i]) %>%
    mutate(outcome_cohort_id = settings_surv$outcome_cohort_id[i])
  
}

Pop <- dplyr::bind_rows(pops)

# format data -----

# add region

Pop <-
  Pop %>%
  left_join(
    Pop %>%
   select("person_id") %>%
  left_join(cdm$person %>%
                  select("person_id", "care_site_id"),
                by = "person_id", `copy` = TRUE) %>%
  left_join(cdm$care_site %>%
              select("care_site_id", "location_id"),
            by = "care_site_id", `copy` = TRUE) %>%
  # add region
  left_join(cdm$location %>%
              rename("region" = "location_source_value") %>%
              select(c("location_id", "region")),
            by = "location_id", `copy` = TRUE ) %>% 
  mutate(region_collapsed = case_when(region %in% c("Yorkshire  & The Humber", "East Midlands", 
                                                    "West Midlands", "North East", "North West", "East of England", "London", 
                                                    "South East", "South West") ~ "England",
                                      region == "Northern Ireland" ~ "Northern Ireland",
                                      region == "Scotland" ~ "Scotland",
                                      region == "Wales" ~ "Wales")) 
  ) %>% 
  select(-c(care_site_id, location_id))

#add age -----
Pop$age<- NA
if(sum(is.na(Pop$day_of_birth))==0 & sum(is.na(Pop$month_of_birth))==0){
  # if we have day and month
  Pop <-Pop %>%
    mutate(age=floor(as.numeric((ymd(outcome_start_date)-
                                   ymd(paste(year_of_birth,
                                             month_of_birth,
                                             day_of_birth, sep="-"))))/365.25))
} else {
  Pop <- Pop %>%
    mutate(age= lubridate::year(outcome_start_date)-year_of_birth)
}

# # age age groups ----
Pop <- Pop %>%
  mutate(age_gr=ifelse(age<30,  "18-29",
                       ifelse(age>=30 &  age<=39,  "30-39",
                              ifelse(age>=40 & age<=49,  "40-49",
                                     ifelse(age>=50 & age<=59,  "50-59",
                                            ifelse(age>=60 & age<=69, "60-69",
                                                   ifelse(age>=70 & age<=79, "70-79",
                                                          ifelse(age>=80 & age<=89, "80-89",
                                                                 ifelse(age>=90, ">=90",
                                                                        NA))))))))) %>%
  mutate(age_gr= factor(age_gr,
                        levels = c("18-29","30-39","40-49", "50-59",
                                   "60-69", "70-79","80-89",">=90")))
table(Pop$age_gr, useNA = "always")

# # reformat gender
# # add gender -----
# #8507 male
# #8532 female
Pop <-Pop %>%
  mutate(gender= ifelse(gender_concept_id==8507, "Male",
                        ifelse(gender_concept_id==8532, "Female", NA ))) %>%
  mutate(gender= factor(gender,
                        levels = c("Male", "Female")))
table(Pop$gender, useNA = "always")

# # if missing (or unreasonable) age or gender, drop ----
Pop <-Pop %>%
  filter(!is.na(age)) %>%
  filter(age>=18) %>%
  filter(age<=110) %>%
  filter(!is.na(gender))
#
# # create sex:agegp categorical variables
Pop <- Pop %>%
  unite('genderAgegp', c(gender,age_gr), remove = FALSE) %>%
  mutate(genderAgegp= factor(genderAgegp,
                             levels = c("Female_18-29","Female_30-39","Female_40-49", "Female_50-59",
                                        "Female_60-69", "Female_70-79","Female_80-89","Female_>=90",
                                        "Male_18-29","Male_30-39","Male_40-49", "Male_50-59",
                                        "Male_60-69", "Male_70-79","Male_80-89","Male_>=90")))

# # drop if missing observation period end date ----
Pop <-Pop %>%
  filter(!is.na(observation_period_end_date))

# Prior history from date of diagnosis from start of observation period
Pop  <- Pop %>% 
  mutate(Prior_history_days = as.numeric(outcome_start_date - observation_period_start_date ))

# Prior history from date of diagnosis
Pop  <- Pop %>% 
  mutate(Prior_history_days_study_start = as.numeric(as.Date(cohort_start_date) - observation_period_start_date ))

# calculate Follow up - calculate end of observation period
Pop <- Pop %>%
  mutate(endOfObservation = ifelse(observation_period_end_date >= studyEndDate, studyEndDate, NA)) %>%
  mutate(endOfObservation = as.Date(endOfObservation) ) %>%
  mutate(endOfObservation = coalesce(endOfObservation, observation_period_end_date))

# calculate follow up in years and days
Pop <-Pop %>%
  mutate(time_days=as.numeric(difftime(endOfObservation,
                                       outcome_start_date,
                                       units="days"))) %>%
  mutate(time_years=time_days/365.25)

# binary death outcome (for survival) ---
# need to take into account follow up
# if death date is > database end data set death to 0
Pop <-Pop %>%
  mutate(status= ifelse(!is.na(death_date), 2, 1 )) %>%
  mutate(status= ifelse(death_date > endOfObservation , 1, status )) %>%
  mutate(status= ifelse(is.na(status), 1, status )) %>%
  mutate(Death = recode(status, 
                        "1" = "Alive", 
                        "2" = "Dead"))

# add in smoking status # former, non or smoker in anytime prior history
# If smoking status was not recorded in year, the methods of last observation carried forward (LOCF)
# Records of non-smoking after those for current or previous smoking were amended to former-smoking. 
# As it was not possible to ascertain smoking status prior to CPRD records, we refer to non- rather than never-smoking.

# non smoking 5 years
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter( observation_concept_id == 4222303 |
                             observation_concept_id ==  4144272 
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(1826) ) %>% # removes anyone with observation more than 5 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("non_smoker_date5y"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()

# non smoking 10 years
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter( observation_concept_id == 4222303 |
                             observation_concept_id ==  4144272 
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(3650) ) %>% # removes anyone with observation more than 10 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("non_smoker_date10y"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()


# current non smoking 5 yr
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter(observation_concept_id == 4052464
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(1826) ) %>% # removes anyone with observation more than 5 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("current_non_smoker_date5yr"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()

# current non smoking 10 yr
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter(observation_concept_id == 4052464
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(3650) ) %>% # removes anyone with observation more than 10 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("current_non_smoker_date10yr"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()

# Smoker any time 5 yr
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter(
                     observation_concept_id == 4005823 |
                       observation_concept_id == 4298794 | 
                       observation_concept_id == 762499 |
                       observation_concept_id == 762498 |
                       observation_concept_id == 37395605 |
                       observation_concept_id == 4246415 |
                       observation_concept_id ==  4218917 |
                       observation_concept_id == 4209006 |
                       observation_concept_id ==  4276526 |
                       observation_concept_id == 4209585 |
                       observation_concept_id == 764104 |
                       observation_concept_id == 764103 |
                       observation_concept_id == 4058138 |
                       observation_concept_id == 4042037 |
                       observation_concept_id == 4044777 |
                       observation_concept_id == 4044776 |
                       observation_concept_id == 4044775 |
                       observation_concept_id == 4041511 |
                       observation_concept_id == 4052029 |
                       observation_concept_id == 4058136 |
                       observation_concept_id == 4044778 |
                       observation_concept_id == 4052030 |
                       observation_concept_id == 4144273 |
                       observation_concept_id == 4052947 |
                       observation_concept_id == 4209006 |
                       observation_concept_id == 4204653 |
                       observation_concept_id == 44789712 |
                       observation_concept_id == 40486518 |
                       observation_concept_id ==4046886 |
                       observation_concept_id == 4058137 |
                       observation_concept_id == 4216174 |
                       observation_concept_id == 4215409 |
                       observation_concept_id == 4190573 |
                       observation_concept_id == 4052948 
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(1826) ) %>% # removes anyone with observation more than 5 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("Smoker_date_5yr"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()


# Smokers current or code within 10 years
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter(
                     observation_concept_id == 4005823 |
                       observation_concept_id == 4298794 | 
                       observation_concept_id == 762499 |
                       observation_concept_id == 762498 |
                       observation_concept_id == 37395605 |
                       observation_concept_id == 4246415 |
                       observation_concept_id ==  4218917 |
                       observation_concept_id == 4209006 |
                       observation_concept_id ==  4276526 |
                       observation_concept_id == 4209585 |
                       observation_concept_id == 764104 |
                       observation_concept_id == 764103 |
                       observation_concept_id == 4058138 |
                       observation_concept_id == 4042037 |
                       observation_concept_id == 4044777 |
                       observation_concept_id == 4044776 |
                       observation_concept_id == 4044775 |
                       observation_concept_id == 4041511 |
                       observation_concept_id == 4052029 |
                       observation_concept_id == 4058136 |
                       observation_concept_id == 4044778 |
                       observation_concept_id == 4052030 |
                       observation_concept_id == 4144273 |
                       observation_concept_id == 4052947 |
                       observation_concept_id == 4209006 |
                       observation_concept_id == 4204653 |
                       observation_concept_id == 44789712 |
                       observation_concept_id == 40486518 |
                       observation_concept_id ==4046886 |
                       observation_concept_id == 4058137 |
                       observation_concept_id == 4216174 |
                       observation_concept_id == 4215409 |
                       observation_concept_id == 4190573 |
                       observation_concept_id == 4052948 
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(3650) ) %>% # removes anyone with observation more than 10 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("Smoker_date_10yr"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()

# former smokers 5 year
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter(
                     observation_concept_id == 4052032 | # stopped smoking
                       observation_concept_id == 44802805 # recently stopped smoking
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(1826) ) %>% # removes anyone with observation more than 5 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("PrevSmoker_date5y"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()


# former smokers 10 year
Pop <-
  Pop %>%
  left_join(
    Pop %>%
      select("person_id", "outcome_start_date") %>%
      inner_join(cdm$observation %>%
                   filter(
                     observation_concept_id == 4052032 | # stopped smoking
                       observation_concept_id == 44802805 # recently stopped smoking
                   ) ,
                 by=c("person_id"), copy = TRUE)  %>%
      filter(observation_date < outcome_start_date) %>% # removes anyone with a observation after outcome
      filter(observation_date > outcome_start_date - lubridate::days(3650) ) %>% # removes anyone with observation more than 10 years before the outcome
      select(c(person_id, observation_date, outcome_start_date)) %>%
      distinct() %>%
      group_by(person_id, outcome_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("PrevSmoker_date10y"="observation_date"),
    by= c("person_id", "outcome_start_date")) %>%
  compute()


# rules for updating smoking status
# if only smoking codes == SMOKER, non smoker == NON SMOKER
# Never-smokers were re-assigned as former-smokers if there were contradictory preceding smoking codes
# previous smokers were assigned smokers if date of smoking was after previous smoking code 

# Smoking status using 5 years worth of previous data
# create a new variable for previous and if those who say non smoking convert to former/previous smoker
Pop <- Pop %>% mutate(PrevSmoker_date5y1 = ifelse(!is.na(PrevSmoker_date5y), non_smoker_date5y, PrevSmoker_date5y))
Pop <- Pop %>% mutate(PrevSmoker_date5y1 = replace(PrevSmoker_date5y1, !is.na(PrevSmoker_date5y1), "2Former"))
Pop <- Pop %>% mutate(PrevSmoker_date5y1 = replace(PrevSmoker_date5y1, !is.na(PrevSmoker_date5y), "2Former"))

# remove previous entries if the have a smoking date after previous
Pop <- Pop %>% mutate(Smokernow = Smoker_date_5yr > PrevSmoker_date5y )
Pop <- Pop %>% mutate(PrevSmoker_date5y1 = replace(PrevSmoker_date5y1, Smokernow == TRUE, NA))

# convert smoking into a character variable
Pop <- Pop %>% mutate(Smoker_date_5yr1 = as.character(Smoker_date_5yr))
Pop <- Pop %>% mutate(Smoker_date_5yr1 = replace(Smoker_date_5yr1, !is.na(Smoker_date_5yr1), "3Smoker"))

# convert non smokers into non smoker
Pop <- Pop %>% mutate(non_smoker_date5y1 = as.character(non_smoker_date5y))
Pop <- Pop %>% mutate(non_smoker_date5y1 = replace(non_smoker_date5y1, !is.na(non_smoker_date5y1), "1Non Smoker"))

Pop <- Pop %>% mutate(current_non_smoker_date5yr1 = as.character(current_non_smoker_date5yr))
Pop <- Pop %>% mutate(current_non_smoker_date5yr1 = replace(current_non_smoker_date5yr1, !is.na(current_non_smoker_date5yr1), "1Non Smoker"))

# merge all to get smoking status in past 5 years
Pop <- Pop %>% mutate(smoking_status5yr = coalesce(PrevSmoker_date5y1, Smoker_date_5yr1))
Pop <- Pop %>% mutate(smoking_status5yr = coalesce(smoking_status5yr, non_smoker_date5y1))
Pop <- Pop %>% mutate(smoking_status5yr = coalesce(smoking_status5yr, current_non_smoker_date5yr1))
# replace those with missing observations with "Missing"
Pop <- Pop %>% mutate(smoking_status5yr = replace(smoking_status5yr, is.na(smoking_status5yr), "4Missing"))



# Smoking status 10 years previous
Pop <- Pop %>% mutate(PrevSmoker_date10y1 = ifelse(!is.na(PrevSmoker_date10y), non_smoker_date10y, PrevSmoker_date10y))
Pop <- Pop %>% mutate(PrevSmoker_date10y1 = replace(PrevSmoker_date10y1, !is.na(PrevSmoker_date10y1), "2Former"))
Pop <- Pop %>% mutate(PrevSmoker_date10y1 = replace(PrevSmoker_date10y1, !is.na(PrevSmoker_date10y), "2Former"))

# remove previous entries if the have a smoking date after previous
Pop <- Pop %>% mutate(Smokernow = Smoker_date_10yr > PrevSmoker_date10y )
Pop <- Pop %>% mutate(PrevSmoker_date10y1 = replace(PrevSmoker_date10y1, Smokernow == TRUE, NA))

# convert smoking into a character variable
Pop <- Pop %>% mutate(Smoker_date_10yr1 = as.character(Smoker_date_10yr))
Pop <- Pop %>% mutate(Smoker_date_10yr1 = replace(Smoker_date_10yr1, !is.na(Smoker_date_10yr1), "3Smoker"))

# convert non smokers into non smoker
Pop <- Pop %>% mutate(non_smoker_date10y1 = as.character(non_smoker_date10y))
Pop <- Pop %>% mutate(non_smoker_date10y1 = replace(non_smoker_date10y1, !is.na(non_smoker_date10y1), "1Non Smoker"))

Pop <- Pop %>% mutate(current_non_smoker_date10yr1 = as.character(current_non_smoker_date10yr))
Pop <- Pop %>% mutate(current_non_smoker_date10yr1 = replace(current_non_smoker_date10yr1, !is.na(current_non_smoker_date10yr1), "1Non Smoker"))

# merge all to get smoking status in past 10 years
Pop <- Pop %>% mutate(smoking_status10yr = coalesce(PrevSmoker_date10y1, Smoker_date_10yr1))
Pop <- Pop %>% mutate(smoking_status10yr = coalesce(smoking_status10yr, non_smoker_date10y1))
Pop <- Pop %>% mutate(smoking_status10yr = coalesce(smoking_status10yr, current_non_smoker_date10yr1))
# replace those with missing observations with "Missing"
Pop <- Pop %>% mutate(smoking_status10yr = replace(smoking_status10yr, is.na(smoking_status10yr), "4Missing"))

# conditions (any time in history)
for(i in seq_along(disease_cohorts$cohort_definition_id)){
  
  working_name <- glue::glue("{disease_cohorts$cohort_name[[i]]}")
  working_id <- disease_cohorts$cohort_definition_id[[i]]
  Pop <-
    Pop %>%
    left_join(
      Pop %>%
        select("person_id", "outcome_start_date") %>%
        inner_join(cdm[[feature_disease_table_name]] %>%
                     rename("feature_start_date"="cohort_start_date") %>%
                     filter(cohort_definition_id== working_id ) %>%
                     select(!c(cohort_definition_id,
                               cohort_end_date)),
                   by=c("person_id" = "subject_id"), copy = TRUE) %>%
        filter(feature_start_date < outcome_start_date) %>%
        select(person_id) %>%
        distinct() %>%
        mutate(!!working_name:=1),
      by="person_id")  %>%
    compute()
  
}

# conditions (liver related)
for(i in seq_along(disease_cohorts_liver$cohort_definition_id)){
  
  working_name <- glue::glue("{disease_cohorts_liver$cohort_name[[i]]}")
  working_id <- disease_cohorts_liver$cohort_definition_id[[i]]
  Pop <- 
    Pop %>% 
    left_join(
      Pop %>% 
        select("person_id", "outcome_start_date") %>% 
        inner_join(cdm[[feature_disease_liver_table_name]] %>% 
                     rename("feature_start_date"="cohort_start_date") %>% 
                     filter(cohort_definition_id== working_id ) %>% 
                     select(!c(cohort_definition_id,
                               cohort_end_date)),
                   by=c("person_id" = "subject_id"), copy = TRUE) %>% 
        filter(feature_start_date < outcome_start_date) %>% 
        select(person_id) %>% 
        distinct() %>% 
        mutate(!!working_name:=1),
      by="person_id")  %>% 
    compute()
  
}

# gender splits for characterisation
PopFemale <- Pop %>% 
  filter(gender == "Female")
PopMale <- Pop %>% 
  filter(gender == "Male")

# tidy up results for table 1
get_summary_characteristics <-function(data){
  
  summary_characteristics<- bind_rows(
    data %>% 
      count() %>% 
      mutate(var="N"),
    
    data %>% 
      summarise(mean=nice.num(mean(age)),
                standard_deviation = nice.num(sd(age)),
                median = nice.num(median(age)),
                interquartile_range=paste0(nice.num.count(quantile(age,probs=0.25)),  " to ",
                                           nice.num.count(quantile(age,probs=0.75)))) %>% 
      mutate(var="age"),
    
    data %>% 
      group_by(age_gr) %>% 
      summarise(n=n(),
                percent=paste0(nice.num((n/nrow(data))*100),  "%")) %>%   
      rename("var"="age_gr") %>% 
      mutate(var=paste0("Age group: ", var)),
    
    data %>% 
      mutate(gender=factor(gender, levels=c("Male", "Female"))) %>% 
      group_by(gender) %>% 
      summarise(n=n(),
                percent=paste0(nice.num((n/nrow(data))*100),  "%")) %>%   
      rename("var"="gender") %>% 
      mutate(var=paste0("Sex: ", var)),
    
    data %>% 
      mutate(Death=factor(Death, levels=c("Alive", "Dead"))) %>% 
      group_by(Death) %>% 
      summarise(n=n(),
                percent=paste0(nice.num((n/nrow(data))*100),  "%")) %>%   
      rename("var"="Death") %>% 
      mutate(var=paste0("Death: ", var)),
    
    data %>%
      mutate(smoking_status5yr=factor(smoking_status5yr, levels=c("1Non Smoker", "2Former", "3Smoker", "4Missing"))) %>%
      group_by(smoking_status5yr) %>%
      summarise(n=n(),
                percent=paste0(nice.num((n/nrow(data))*100),  "%")) %>%
      rename("var"="smoking_status5yr") %>%
      mutate(var=paste0("smoking_status5yr: ", var)),
    
    data %>%
      mutate(smoking_status10yr=factor(smoking_status10yr, levels=c("1Non Smoker", "2Former", "3Smoker", "4Missing"))) %>%
      group_by(smoking_status10yr) %>%
      summarise(n=n(),
                percent=paste0(nice.num((n/nrow(data))*100),  "%")) %>%
      rename("var"="smoking_status10yr") %>%
      mutate(var=paste0("smoking_status10yr: ", var)),
    
    data %>% 
      summarise(mean=nice.num.count(mean(Prior_history_days)),
                standard_deviation = nice.num(sd(Prior_history_days)),
                median = nice.num(median(Prior_history_days)),
                interquartile_range=paste0(nice.num.count(quantile(Prior_history_days,probs=0.25)),  " to ",
                                           nice.num.count(quantile(Prior_history_days,probs=0.75)))) %>% 
      mutate(var="Prior_history_days"),
    
    data %>% 
      summarise(mean=nice.num.count(mean((Prior_history_days/365.25))),
                standard_deviation = nice.num(sd((Prior_history_days/365.25))),
                median = nice.num(median((Prior_history_days/365.25))),
                interquartile_range=paste0(nice.num.count(quantile((Prior_history_days/365.25),probs=0.25)),  " to ",
                                           nice.num.count(quantile((Prior_history_days/365.25),probs=0.75)))) %>% 
      mutate(var="Prior_history_years"),
    
    data %>% 
      summarise(mean=nice.num.count(mean(Prior_history_days_study_start)),
                standard_deviation = nice.num(sd(Prior_history_days_study_start)),
                median = nice.num(median(Prior_history_days_study_start)),
                interquartile_range=paste0(nice.num.count(quantile(Prior_history_days_study_start,probs=0.25)),  " to ",
                                           nice.num.count(quantile(Prior_history_days_study_start,probs=0.75)))) %>% 
      mutate(var="Prior_history_days_study_start"),
    
    data %>% 
      summarise(mean=nice.num.count(mean((Prior_history_days_study_start/365.25))),
                standard_deviation = nice.num(sd((Prior_history_days_study_start/365.25))),
                median = nice.num(median((Prior_history_days_study_start/365.25))),
                interquartile_range=paste0(nice.num.count(quantile((Prior_history_days_study_start/365.25),probs=0.25)),  " to ",
                                           nice.num.count(quantile((Prior_history_days_study_start/365.25),probs=0.75)))) %>% 
      mutate(var="Prior_history_years_start"),
    
    
    data %>% 
      summarise(mean=nice.num.count(mean(time_days)),
                standard_deviation = nice.num(sd(time_days)),
                median = nice.num(median(time_days)),
                interquartile_range=paste0(nice.num.count(quantile(time_days,probs=0.25)),  " to ",
                                           nice.num.count(quantile(time_days,probs=0.75)))) %>% 
      mutate(var="time_days"),
    
    data %>% 
      summarise(mean=nice.num.count(mean((time_years))),
                standard_deviation = nice.num(sd((time_years))),
                median = nice.num(median((time_years))),
                interquartile_range=paste0(nice.num.count(quantile((time_years),probs=0.25)),  " to ",
                                           nice.num.count(quantile((time_years),probs=0.75)))) %>% 
      mutate(var="time_years"))
  
  #disease features
  for(i in seq_along(disease_cohorts$cohort_name)){
    working_id_name <- glue::glue("{disease_cohorts$cohort_name[[i]]}")
    summary_characteristics <- bind_rows(summary_characteristics,
                                         data %>%
                                           summarise(n=sum(!is.na(!!rlang::sym(working_id_name))),
                                                     percent=paste0(nice.num((n/nrow(data))*100),  "%"))%>%
                                           mutate(var=working_id_name)
    )
  }
  # #other cancers
  # for(i in seq_along(outcome_cohorts$cohort_name)){
  #   working_name <- glue::glue("{outcome_cohorts$cohort_name[[i]]}")
  #   summary_characteristics <- bind_rows(summary_characteristics,
  #                                        data %>%
  #                                          summarise(n=sum(!is.na(!!rlang::sym(working_name))),
  #                                                    percent=paste0(nice.num((n/nrow(data))*100),  "%"))%>%
  #                                          mutate(var=working_name)
  #   )
  # }
  
  #liver cancer features
  for(i in seq_along(disease_cohorts_liver$cohort_name)){
    working_id_name <- glue::glue("{disease_cohorts_liver$cohort_name[[i]]}")
    summary_characteristics <- bind_rows(summary_characteristics,
                                         data %>% 
                                           summarise(n=sum(!is.na(!!rlang::sym(working_id_name))),
                                                     percent=paste0(nice.num((n/nrow(data))*100),  "%"))%>% 
                                           mutate(var=working_id_name)
    )
  }
  
  
  # filter any less than 5
  summary_characteristics <- summary_characteristics %>% 
    mutate(mean=ifelse(!is.na(n) & n<5, NA, mean)) %>% 
    mutate(percent=ifelse(!is.na(n) & n<5, NA, percent)) %>% 
    mutate(interquartile_range=ifelse(!is.na(n) & n<5, NA, interquartile_range)) %>% 
    mutate(standard_deviation=ifelse(!is.na(n) & n<5, NA, standard_deviation)) %>% 
    mutate(n=ifelse(!is.na(n) & n<5, "<5", n))
  
  return(summary_characteristics %>% 
           relocate(any_of(c("var", "n", "percent",
                             "mean", "standard_deviation",
                             "median", "interquartile_range"))))
  
}

# get a list to put results into
table1Characteristics <- list()

#create a loop that puts table 1 for each outcome
for(j in seq_along(table(Pop$region_collapsed))){
  
  table1Characteristics[[j]] <- get_summary_characteristics(Pop %>% 
                                                              filter(region_collapsed== names(table(Pop$region_collapsed))[[j]])) %>%
    mutate(Cancer = lung_cancer_outcome_cohorts$cohort_name,
           Region = names(table(Pop$region_collapsed))[[j]])
  
  table1Characteristics[[j]]$n <- as.character(table1Characteristics[[j]]$n)
  
  print(paste0("Getting table 1 for ", lung_cancer_outcome_cohorts$cohort_name , " (" , names(table(Pop$region_collapsed))[[j]], ")")) 
  
}

#bind all results together
table1Characteristics_region <- bind_rows(table1Characteristics) %>%
  mutate(Database = db.name, analysis = "Incidence_Region")

#save the results
write_csv(table1Characteristics_region, here::here(paste0("Results/",db.name,"/Table1_by_region",db.name,".csv")))



# get table one for smoking status

# get a list to put results into
table1Characteristics1 <- list()

#create a loop that puts table 1 for each smoking outcome
for(j in seq_along(table(Pop$smoking_status5yr))){
  
  table1Characteristics1[[j]] <- get_summary_characteristics(Pop %>% 
                                                              filter(smoking_status5yr == names(table(Pop$smoking_status5yr))[[j]])) %>%
    mutate(Cancer = lung_cancer_outcome_cohorts$cohort_name,
           Smoking_status = names(table(Pop$smoking_status5yr))[[j]])
  
  table1Characteristics1[[j]]$n <- as.character(table1Characteristics1[[j]]$n)
  
  print(paste0("Getting table 1 for ", lung_cancer_outcome_cohorts$cohort_name , " (" , names(table(Pop$smoking_status5yr))[[j]], ")")) 
  
}

#bind all results together
table1Characteristics_smoking <- bind_rows(table1Characteristics1) %>%
  mutate(Database = db.name, analysis = "Incidence_smoking")

#save the results
write_csv(table1Characteristics_smoking, here::here(paste0("Results/",db.name,"/Table1_by_5yrsmokingstatus",db.name,".csv")))



# survival analysis
Pop <-Pop %>%
  filter(time_days != 0)

# split them up here

#OUTPUT data for whole dataset and strata based on calender year
PopAll <- DataExtraction(dataset = Pop)

### KAPLAIN MEIER CODE ####
# INPUT
# dataset is a dataframe containing the data as processed above
# outcome cohort is a dataframe containing information about cancers
# OUTPUT
# list containing 4 data frames 1) survival estimates 2) risktables 3) median survival 4) survival probabilities

SurAnalysis <- function(dataset, outcomeCohort) {
  
  # whole population
  observedkm <- list()
  observedmedianKM <- list()
  observedsurprobsKM <- list()
  observedrisktableKM <- list()
  
  # loop to carry out for each cancer
  for(j in 1:nrow(outcomeCohort)) { 
    
    #subset the data by cancer type
    data <- dataset %>%
      filter(outcome_cohort_id == j)
    
    tryCatch({
    #carry out km estimate
    observedkm[[j]] <- survfit (Surv(time_years, status) ~ 1, data=data) %>%
      tidy() %>%
      mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Age = "All", Gender = "Both") 
    

    
    print(paste0("KM for observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
    
    # get the risk table ---
    if(ceiling(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")) <=5) {   
      grid <- seq(0,floor(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")),by=1) 
    } else {
      grid <- seq(0,floor(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")),by=2) 
    }
    
    observedrisktableKM[[j]] <- RiskSetCount(grid,data$time_years) %>%
      rbind(grid) %>% as.data.frame() %>%
      `colnames<-`(grid) %>%
      mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Age = "All", Gender = "Both" ) %>%
      slice(1)
    
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
    
    print(paste0("Extract risk table ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
    
    tryCatch({
    # KM median survival---
    modelKM <- survfit(Surv(time_years, status) ~ 1, data=data) %>%
      summary()
    
    observedmedianKM[[j]] <- modelKM$table %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%  
      pivot_longer(-rowname) %>% 
      pivot_wider(names_from=rowname, values_from=value) %>%
      mutate(Method = "Kaplan-Meier", 
             Cancer = outcomeCohort$cohort_name[j],
             Gender = "Both" ,
             Age = "All" ) %>%
      select(-name)
    
    
    print(paste0("Median survival from KM from observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
    
    #grab survival probabilities 1,5,10 years
    sprob <- survfit(Surv(time_years, status) ~ 1, data=data) %>% 
      summary(times = c(1,5,10), extend = TRUE)
    
    cols <- lapply(c(2:15) , function(x) sprob[x])
    observedsurprobsKM[[j]] <- do.call(data.frame, cols) %>%
      mutate(Method = "Kaplan-Meier", 
             Cancer = outcomeCohort$cohort_name[j],
             Gender = "Both" ,
             Age = "All" )
    
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
    
    print(paste0("survival probabilites for 1,5,10 years from KM from observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
    
    
  }
  
  
tryCatch({
  
  # take the results from a list (one element for each cancer) and put into dataframe for KM survival
  observedkmcombined <- dplyr::bind_rows(observedkm) %>%
    rename(est = estimate , ucl = conf.high, lcl = conf.low ) %>%
    mutate(Stratification = "None")
  
}, error = function(e) {
  # Print a custom message
  cat("An error occurred. Skipping over.\n")
  # Print the error message
  print(e)
})
  
  tryCatch({
  medkmcombined <- dplyr::bind_rows(observedmedianKM) %>%
    mutate(Stratification = "None")
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
  
  
  tryCatch({
  # generate the risk table and remove entries < 5 patients
  risktableskm <- dplyr::bind_rows(observedrisktableKM) %>%
    mutate(across(everything(), ~replace(., . ==  0 , NA))) %>%
    mutate(across(everything(), ~replace(., .  <=  5 , "<5"))) %>%
    replace(is.na(.), 0) %>%
    relocate(Cancer) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(Stratification = "None")
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
  
  
  tryCatch({
  #generate probabilities
  sprobkmcombined <- dplyr::bind_rows(observedsurprobsKM) %>%
    mutate(Stratification = "None")
  
  
}, error = function(e) {
  # Print a custom message
  cat("An error occurred. Skipping over.\n")
  # Print the error message
  print(e)
})

  
  info(logger, 'KM analysis for whole population COMPLETE')
  
  # GENDER STRATIFICATION-----
  
  observedkm_gender <- list()
  observedmedianKM_gender <- list()
  observedsurprobsKM_gender <- list()
  observedrisktableKM_gender <- list()
  
  # loop to carry out for each cancer
  for(j in 1:nrow(outcomeCohort)) { 
    
    #subset the data by cancer type
    data <- dataset %>%
      filter(outcome_cohort_id == j) 
    
    tryCatch({
    
    # get the risk table ---
    if(ceiling(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")) <=5) {   
      grid <- seq(0,floor(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")),by=1) 
    } else {
      grid <- seq(0,floor(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")),by=2) 
    }
    
    #creates a test that determines if both genders in the data
    genderlevels <- data %>%
      group_by(gender) %>% summarise(count = n()) %>% tally()
    
    # analysis wont run if only 1 gender present
    if(genderlevels == 2){
      
      # get the risk table ---
      
      observedrisktableKM_gender[[j]] <- RiskSetCount(grid,data$time_years[data$gender == "Male"])%>%
        rbind(grid) %>% as.data.frame() %>%
        `colnames<-`(grid) %>%
        mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Age = "All") %>%
        slice(1) %>%
        rbind(RiskSetCount(grid,data$time_years[data$gender == "Female"]))%>%
        mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Age = "All", Gender = c("Male", "Female"))
      
      print(paste0("Extract risk table ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
      
      tryCatch({
      #carry out km estimate
      observedkm_gender[[j]] <- survfit(Surv(time_years, status) ~ gender, data=data) %>%
        tidy() %>%
        rename(Gender = strata) %>%
        mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Age = "All", Gender = str_replace(Gender, "gender=Male", "Male"), Gender = str_replace(Gender,"gender=Female", "Female")) %>%
        filter(n.risk >= 5) #remove entries with less than 5 patients
      
      print(paste0("KM for observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
      
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
      
      tryCatch({
      # KM median survival ---
      modelKM <- survfit(Surv(time_years, status) ~ gender, data=data) %>%
        summary()
      
      # median survival ---
      observedmedianKM_gender[[j]] <- modelKM$table %>%
        as.data.frame() %>%
        mutate(Method = "Kaplan-Meier", 
               Cancer = outcomeCohort$cohort_name[j], 
               Age = "All" ,
               Gender = c("Male", "Female"))
      
      print(paste0("Median survival from KM from observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
      
      #grab survival probabilities 1,5,10 years
      sprob <- survfit(Surv(time_years, status) ~ gender, data=data) %>%
        summary(times = c(1,5,10), extend = TRUE)
      
      cols <- lapply(c(2:16) , function(x) sprob[x])
      observedsurprobsKM_gender[[j]] <- do.call(data.frame, cols) %>%
        rename(Gender = strata) %>%
        mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Age = "All", Gender = str_replace(Gender, "gender=Male", "Male"), Gender = str_replace(Gender,"gender=Female", "Female"))
      
      print(paste0("survival probabilites for 1,5,10 years from KM from observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
      
      }, error = function(e) {
        # Print a custom message
        cat("An error occurred. Skipping over.\n")
        # Print the error message
        print(e)
      })
      
      
      
    } else {
      
      print(paste0("Gender stratification KM analysis not carried out for ", outcomeCohort$cohort_name[j], " due to only 1 gender present " , Sys.time()))
      
    }
    
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
    
  } # this closes the loop on the analysis containing both genders
  
  tryCatch({
  # take the results from a list (one element for each cancer) and put into dataframe for KM survival
  observedkmcombined_gender <- dplyr::bind_rows(observedkm_gender) %>%
    rename(est = estimate ,ucl = conf.high, lcl = conf.low ) %>%
    mutate(Stratification = "Gender")
  
  medkmcombined_gender <- dplyr::bind_rows(observedmedianKM_gender) %>%
    mutate(Stratification = "Gender")
  
  #generate the risk table and remove entries < 5 patients
  risktableskm_gender <- dplyr::bind_rows(observedrisktableKM_gender) 
  risktableskm_gender <- risktableskm_gender %>%
    mutate_at(.vars = c(1:(ncol(risktableskm_gender)-4)), funs(ifelse(.== 0, NA, .))) %>%  
    mutate_at(.vars = c(1:(ncol(risktableskm_gender)-4)), funs(ifelse(.<= 5, "<5", .))) %>%
    replace(is.na(.), 0) %>%
    relocate(Cancer) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(Stratification = "Gender")
  
  #generate 1,5,10 probabilities
  sprobkmcombined_gender <- dplyr::bind_rows(observedsurprobsKM_gender) %>%
    mutate(Stratification = "Gender")
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
  
  info(logger, 'KM analysis for gender stratification COMPLETE')
  
  ###########################
  # AGE STRATIFICATION
  ##########################
  
  observedkm_age <- list()
  observedmedianKM_age <- list()
  observedsurprobsKM_age <- list()
  observedrisktableKM_age <- list()
  
  # loop to carry out for each cancer
  for(j in 1:nrow(outcomeCohort)) { 
    
    #subset the data by cancer type
    data <- dataset %>%
      filter(outcome_cohort_id == j) 
    
    tryCatch({
    # get the risk table ---
    if(ceiling(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")) <=5) {   
      grid <- seq(0,floor(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")),by=1) 
    } else {
      grid <- seq(0,floor(lubridate::time_length(difftime(max(data$outcome_start_date), min(data$outcome_start_date)), "years")),by=2) 
    }
      
    observedrisktableKM_age[[j]] <- RiskSetCount(grid,data$time_years[data$age_gr == "18-29"]) %>%
      rbind(grid) %>% as.data.frame() %>%
      `colnames<-`(grid) %>%
      slice(1) %>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == "30-39"]))%>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == "40-49"]))%>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == "50-59"]))%>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == "60-69"]))%>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == "70-79"]))%>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == "80-89"]))%>%
      rbind(RiskSetCount(grid,data$time_years[data$age_gr == ">=90"]))%>%
      mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], Gender = "Both", Age = c("18-29" ,"30-39", "40-49" ,"50-59" ,"60-69", "70-79", "80-89" ,">=90")) 
    
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
    
    
    tryCatch({
    #carry out km estimate
    observedkm_age[[j]] <- survfit (Surv(time_years, status) ~ age_gr, data=data) %>%
      tidy() %>%
      rename(Age = strata) %>%
      mutate(Method = "Kaplan-Meier", Cancer = outcomeCohort$cohort_name[j], 
             Age = str_replace(Age, "age_gr=18-29", "18-29"),
             Age = str_replace(Age, "age_gr=30-39", "30-39"),
             Age = str_replace(Age, "age_gr=40-49", "40-49"),
             Age = str_replace(Age, "age_gr=50-59", "50-59"),
             Age = str_replace(Age, "age_gr=60-69", "60-69"),
             Age = str_replace(Age, "age_gr=70-79", "70-79"),
             Age = str_replace(Age, "age_gr=80-89", "80-89"),
             Age = str_replace(Age, "age_gr=>=90", ">=90"),
             Gender = "Both")
    
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
    
    print(paste0("KM for observed data age strat ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))

    tryCatch({
    # KM median survival---
    modelKM <- survfit(Surv(time_years, status) ~ age_gr, data=data) %>%
      summary()
    
    # to create an age label 
    agelabel <- table(data$age_gr) %>% as.data.frame() %>% filter(Freq != 0)
    
    observedmedianKM_age[[j]] <- modelKM$table %>%
      as.data.frame() %>%
      mutate(Method = "Kaplan-Meier", 
             Cancer = outcomeCohort$cohort_name[j], 
             Gender = "Both" ,
             Age = agelabel$Var1)
    
    print(paste0("Median survival from KM from observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
    
    
    #grab survival probabilities 1,5,10 years
    sprob <- survfit(Surv(time_years, status) ~ age_gr, data=data) %>%
      summary(times = c(1,5,10), extend = TRUE)
    
    cols <- lapply(c(2:16) , function(x) sprob[x])
    observedsurprobsKM_age[[j]] <- do.call(data.frame, cols) %>%
      rename(Age = strata) %>%
      mutate(Method = "Kaplan-Meier", 
             Cancer = outcomeCohort$cohort_name[j], 
             Age = str_replace(Age, "age_gr=18-29", "18-29"),
             Age = str_replace(Age, "age_gr=30-39", "30-39"),
             Age = str_replace(Age, "age_gr=40-49", "40-49"),
             Age = str_replace(Age, "age_gr=50-59", "50-59"),
             Age = str_replace(Age, "age_gr=60-69", "60-69"),
             Age = str_replace(Age, "age_gr=70-79", "70-79"),
             Age = str_replace(Age, "age_gr=80-89", "80-89"),
             Age = str_replace(Age, "age_gr=>=90", ">=90"),
             Gender = "Both")
    
    
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
    
    print(paste0("survival probabilites for 1,5,10 years from KM from observed data ", Sys.time()," for ",outcomeCohort$cohort_name[j], " completed"))
    
  }
  
  
  tryCatch({
  # take the results from a list (one element for each cancer) and put into dataframe ----
  observedkmcombined_age <- dplyr::bind_rows(observedkm_age) %>%
    rename(est = estimate ,ucl = conf.high, lcl = conf.low ) %>%
    mutate(Stratification = "Age")
  
  medkmcombined_age <- dplyr::bind_rows(observedmedianKM_age) %>%
    mutate(Stratification = "Age")
  
  #generate the risk table and obscure entries < 5 patients
  risktableskm_age <- dplyr::bind_rows(observedrisktableKM_age)
  risktableskm_age <- risktableskm_age %>%
    mutate_at(.vars = c(1:(ncol(risktableskm_age)-4)), funs(ifelse(.== 0, NA, .))) %>%  
    mutate_at(.vars = c(1:(ncol(risktableskm_age)-4)), funs(ifelse(.<= 5, "<5", .))) %>%
    replace(is.na(.), 0) %>%
    relocate(Cancer) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(Stratification = "Age")
  
  #generate 1,5,10 probabilities
  sprobkmcombined_age <- dplyr::bind_rows(observedsurprobsKM_age) %>%
    mutate(Stratification = "Age")
  
  info(logger, 'KM analysis for AGE stratification COMPLETE')
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
  
  ##################################################################
  # AGE*GENDER STRATIFICATION
  ##########################
  # not carried out
 

  tryCatch({
  # combine all the survival results -----
  survivalResults <- bind_rows(
    observedkmcombined , # all 
    observedkmcombined_gender , # gender strat 
    observedkmcombined_age , # age strat
  ) %>%
    mutate(Database = db.name, CalenderYearGp = paste0(min(lubridate::year(lubridate::ymd(data$cohort_end_date))),"-",
                                                       max(lubridate::year(lubridate::ymd(data$cohort_end_date)))))
  
}, error = function(e) {
  # Print a custom message
  cat("An error occurred. Skipping over.\n")
  # Print the error message
  print(e)
})
  
  tryCatch({
  #risk table # error with characters and double formats
  riskTableResults <- bind_rows(
    risktableskm , # all
    risktableskm_gender , # gender strat
    risktableskm_age , # age strat
  ) %>%
    mutate(Database = db.name, CalenderYearGp = paste0(min(lubridate::year(lubridate::ymd(data$cohort_end_date))),"-",
                                                       max(lubridate::year(lubridate::ymd(data$cohort_end_date)))))
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })

  
  tryCatch({
  #median results
  medianKMResults <- bind_rows( 
    medkmcombined , # all
    medkmcombined_gender , # gender
    medkmcombined_age , # age strat
  ) %>%
    mutate(Database = db.name, CalenderYearGp = paste0(min(lubridate::year(lubridate::ymd(data$cohort_end_date))),"-",
                                                       max(lubridate::year(lubridate::ymd(data$cohort_end_date)))))
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })

  
  tryCatch({
  #1,5,10 survival probabilites results
  SurvProb1510KMResults <- bind_rows( 
    sprobkmcombined , # all
    sprobkmcombined_gender , # gender
    sprobkmcombined_age , # age strat
  ) %>%
    mutate(Database = db.name, CalenderYearGp = paste0(min(lubridate::year(lubridate::ymd(data$cohort_end_date))),"-",
                                                       max(lubridate::year(lubridate::ymd(data$cohort_end_date)))))
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
  
  
  tryCatch({
  # put results all together in a list
  survival_study_results <- list(survivalResults ,
                                 riskTableResults,
                                 medianKMResults,
                                 SurvProb1510KMResults)
  
  names(survival_study_results) <- c(paste0("survival_estimates"),
                                     paste0("risk_table_results"),
                                     paste0("median_survival_results"),
                                     paste0("one_five_ten_survival_rates")
  )
  
  
  print(paste0("Survival Analysis completed"))
  
  return(survival_study_results)
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })
  
}

# create a loop which carries the analysis out on the number of calender year groups (all data plus the calender time splits)
# for each region
SurResults <- list()
SurResults_region <- list()

for(region in seq_along(table(Pop$region_collapsed))) {
  
  Pop_temp <- purrr::map(PopAll, ~ .x[.x$region_collapsed == names(table(Pop$region_collapsed))[[region]], ])

for(l in 1:length(Pop_temp)) {
  
  SurResults[[l]] <- SurAnalysis(dataset = Pop_temp[[l]],
                                 outcomeCohort = lung_cancer_outcome_cohorts)
  
  
}
  
  SurResults_region[[region]] <- SurResults

}

# get the results per region

# extract results for the whole population eng
whole_pop_results_eng <- list(
  SurResults_region[[1]][[1]]$survival_estimates %>% 
    mutate(Region = "England"),
  SurResults_region[[1]][[1]]$risk_table_results %>% 
    mutate(Region = "England"),
  SurResults_region[[1]][[1]]$median_survival_results%>% 
    mutate(Region = "England"),
  SurResults_region[[1]][[1]]$one_five_ten_survival_rates %>% 
    mutate(Region = "England")
)

names(whole_pop_results_eng) <- c(paste0("survival_estimates", db.name),
                              paste0("risk_table_results", db.name),
                              paste0("median_survival_results", db.name),
                              paste0("one_five_ten_survival_rates", db.name))


# extract results for the whole population ni
whole_pop_results_ni <- list(
  SurResults_region[[2]][[1]]$survival_estimates %>% 
    mutate(Region = "Northern Ireland"),
  SurResults_region[[2]][[1]]$risk_table_results %>% 
    mutate(Region = "Northern Ireland"),
  SurResults_region[[2]][[1]]$median_survival_results%>% 
    mutate(Region = "Northern Ireland"),
  SurResults_region[[2]][[1]]$one_five_ten_survival_rates %>% 
    mutate(Region = "Northern Ireland")
)

names(whole_pop_results_ni) <- c(paste0("survival_estimates", db.name),
                                  paste0("risk_table_results", db.name),
                                  paste0("median_survival_results", db.name),
                                  paste0("one_five_ten_survival_rates", db.name))




# extract results for the whole population scotland
whole_pop_results_scot <- list(
  SurResults_region[[3]][[1]]$survival_estimates %>% 
    mutate(Region = "Scotland"),
  SurResults_region[[3]][[1]]$risk_table_results %>% 
    mutate(Region = "Scotland"),
  SurResults_region[[3]][[1]]$median_survival_results%>% 
    mutate(Region = "Scotland"),
  SurResults_region[[3]][[1]]$one_five_ten_survival_rates %>% 
    mutate(Region = "Scotland")
)

names(whole_pop_results_scot) <- c(paste0("survival_estimates", db.name),
                                  paste0("risk_table_results", db.name),
                                  paste0("median_survival_results", db.name),
                                  paste0("one_five_ten_survival_rates", db.name))

# extract results for the whole population wales
whole_pop_results_wales <- list(
  SurResults_region[[4]][[1]]$survival_estimates %>% 
    mutate(Region = "Wales"),
  SurResults_region[[4]][[1]]$risk_table_results %>% 
    mutate(Region = "Wales"),
  SurResults_region[[4]][[1]]$median_survival_results%>% 
    mutate(Region = "Wales"),
  SurResults_region[[4]][[1]]$one_five_ten_survival_rates %>% 
    mutate(Region = "Wales")
)

names(whole_pop_results_wales) <- c(paste0("survival_estimates", db.name),
                                  paste0("risk_table_results", db.name),
                                  paste0("median_survival_results", db.name),
                                  paste0("one_five_ten_survival_rates", db.name))



#whole database eng
exportSurvivalResults(result=whole_pop_results_eng,
                      zipName= paste0(db.name, "WholeSurvivalResults_eng"),
                      outputFolder=here::here("Results", db.name))

#whole database ni
exportSurvivalResults(result=whole_pop_results_ni,
                      zipName= paste0(db.name, "WholeSurvivalResults_ni"),
                      outputFolder=here::here("Results", db.name))


#whole database scot
exportSurvivalResults(result=whole_pop_results_scot,
                      zipName= paste0(db.name, "WholeSurvivalResults_scot"),
                      outputFolder=here::here("Results", db.name))


#whole database wales
exportSurvivalResults(result=whole_pop_results_wales,
                      zipName= paste0(db.name, "WholeSurvivalResults_wales"),
                      outputFolder=here::here("Results", db.name))


# extract calender year results
surres <- list()
rtres <- list()
msres <- list()
oftsrres <- list()

calenderyr_results <- list()

for(w in 1:length(SurResults_region)) {
  
  SurResults_region_temp <- SurResults_region[[w]] 
  
# extract information for calender year (element 1 is whole population so start from 2:n)
for(q in 2:length(PopAll)) {
  
  tryCatch({
    
  surres[[q]] <- SurResults_region_temp[[q]]$survival_estimates %>% 
    mutate(Region = names(table(Pop$region_collapsed))[w] )
  
  rtres[[q]] <- SurResults_region_temp[[q]]$risk_table_results %>% 
    mutate(Region = names(table(Pop$region_collapsed))[w] )
  
  msres[[q]] <- SurResults_region_temp[[q]]$median_survival_results %>% 
    mutate(Region = names(table(Pop$region_collapsed))[w] )
  
  oftsrres[[q]] <- SurResults_region_temp[[q]]$one_five_ten_survival_rates %>% 
    mutate(Region = names(table(Pop$region_collapsed))[w] )
  
  }, error = function(e) {
    # Print a custom message
    cat("An error occurred. Skipping over.\n")
    # Print the error message
    print(e)
  })

}
  

# bind the results for calender years
survival_results_cy <- bind_rows(surres)
risk_table_cy <- bind_rows(rtres)
med_surv_results_cy <- bind_rows(msres)
survival_prob_cy <- bind_rows(oftsrres)

calenderyr_results <- list(
  survival_results_cy,
  risk_table_cy,
  med_surv_results_cy,
  survival_prob_cy)

names(calenderyr_results) <- c(paste0("survival_estimates_cy", db.name),
                               paste0("risk_table_results_cy", db.name),
                               paste0("median_survival_results_cy", db.name),
                               paste0("one_five_ten_survival_rates_cy", db.name))


#calender year stratification
exportSurvivalResults(result=calenderyr_results,
                      zipName= paste0(db.name, "CalenderYrSurvivalResults", names(table(Pop$region_collapsed))[w]),
                      outputFolder=here::here("Results", db.name))

}



#################################################################
# smoking status

# create a loop which carries the analysis out on the number of calender year groups (all data plus the calender time splits)
# for each smoking group
SurResults <- list()
SurResults_smoke <- list()

for(smoking_status in seq_along(table(Pop$smoking_status5yr))) {
  
  Pop_temp <- purrr::map(PopAll, ~ .x[.x$smoking_status5yr == names(table(Pop$smoking_status5yr))[[smoking_status]], ])
  
  for(l in 1:length(Pop_temp)) {
    
    SurResults[[l]] <- SurAnalysis(dataset = Pop_temp[[l]],
                                   outcomeCohort = lung_cancer_outcome_cohorts)
    
    
  }
  
  SurResults_smoke[[smoking_status]] <- SurResults
  
}

# get the results per smoking status

# extract results for the whole population non smoker
whole_pop_results_nonsmoker <- list(
  SurResults_smoke[[1]][[1]]$survival_estimates %>% 
    mutate(smoking_status = "1Non Smoker"),
  SurResults_smoke[[1]][[1]]$risk_table_results %>% 
    mutate(smoking_status = "1Non Smoker"),
  SurResults_smoke[[1]][[1]]$median_survival_results %>% 
    mutate(smoking_status = "1Non Smoker"),
  SurResults_smoke[[1]][[1]]$one_five_ten_survival_rates %>% 
    mutate(smoking_status = "1Non Smoker")
)

names(whole_pop_results_nonsmoker) <- c(paste0("survival_estimates", db.name),
                                  paste0("risk_table_results", db.name),
                                  paste0("median_survival_results", db.name),
                                  paste0("one_five_ten_survival_rates", db.name))



# extract results for the whole population former smoker
whole_pop_results_fsmoker <- list(
  SurResults_smoke[[2]][[1]]$survival_estimates %>% 
    mutate(smoking_status = "2Former"),
  SurResults_smoke[[2]][[1]]$risk_table_results %>% 
    mutate(smoking_status = "2Former"),
  SurResults_smoke[[2]][[1]]$median_survival_results %>% 
    mutate(smoking_status = "2Former"),
  SurResults_smoke[[2]][[1]]$one_five_ten_survival_rates %>% 
    mutate(smoking_status = "2Former")
)

names(whole_pop_results_fsmoker) <- c(paste0("survival_estimates", db.name),
                                        paste0("risk_table_results", db.name),
                                        paste0("median_survival_results", db.name),
                                        paste0("one_five_ten_survival_rates", db.name))



whole_pop_results_smoker <- list(
  SurResults_smoke[[3]][[1]]$survival_estimates %>% 
    mutate(smoking_status = "3Smoker"),
  SurResults_smoke[[3]][[1]]$risk_table_results %>% 
    mutate(smoking_status = "3Smoker"),
  SurResults_smoke[[3]][[1]]$median_survival_results %>% 
    mutate(smoking_status = "3Smoker"),
  SurResults_smoke[[3]][[1]]$one_five_ten_survival_rates %>% 
    mutate(smoking_status = "3Smoker")
)

names(whole_pop_results_smoker) <- c(paste0("survival_estimates", db.name),
                                      paste0("risk_table_results", db.name),
                                      paste0("median_survival_results", db.name),
                                      paste0("one_five_ten_survival_rates", db.name))



whole_pop_results_smoker_miss <- list(
  SurResults_smoke[[4]][[1]]$survival_estimates %>% 
    mutate(smoking_status = "4Missing"),
  SurResults_smoke[[4]][[1]]$risk_table_results %>% 
    mutate(smoking_status = "4Missing"),
  SurResults_smoke[[4]][[1]]$median_survival_results %>% 
    mutate(smoking_status = "4Missing"),
  SurResults_smoke[[4]][[1]]$one_five_ten_survival_rates %>% 
    mutate(smoking_status = "4Missing")
)

names(whole_pop_results_smoker_miss) <- c(paste0("survival_estimates", db.name),
                                     paste0("risk_table_results", db.name),
                                     paste0("median_survival_results", db.name),
                                     paste0("one_five_ten_survival_rates", db.name))




#whole database nonsmoker
exportSurvivalResults(result=whole_pop_results_nonsmoker,
                      zipName= paste0(db.name, "WholeSurvivalResults_nonsmoker"),
                      outputFolder=here::here("Results", db.name))

#whole database former smoker
exportSurvivalResults(result=whole_pop_results_fsmoker,
                      zipName= paste0(db.name, "WholeSurvivalResults_fsmoker"),
                      outputFolder=here::here("Results", db.name))


#whole database smoker
exportSurvivalResults(result=whole_pop_results_smoker,
                      zipName= paste0(db.name, "WholeSurvivalResults_smoker"),
                      outputFolder=here::here("Results", db.name))


#whole database smoking missing
exportSurvivalResults(result=whole_pop_results_smoker_miss,
                      zipName= paste0(db.name, "WholeSurvivalResults_smoke_missing"),
                      outputFolder=here::here("Results", db.name))


# extract calender year results
surres <- list()
rtres <- list()
msres <- list()
oftsrres <- list()

calenderyr_results <- list()

for(w in 1:length(SurResults_smoke)) {
  
  SurResults_smoke_temp <- SurResults_smoke[[w]] 
  
  # extract information for calender year (element 1 is whole population so start from 2:n)
  for(q in 2:length(PopAll)) {
    
    tryCatch({
      
      surres[[q]] <- SurResults_smoke_temp[[q]]$survival_estimates %>% 
        mutate(Region = names(table(Pop$smoking_status5yr))[w] )
      
      rtres[[q]] <- SurResults_smoke_temp[[q]]$risk_table_results %>% 
        mutate(Region = names(table(Pop$smoking_status5yr))[w] )
      
      msres[[q]] <- SurResults_smoke_temp[[q]]$median_survival_results %>% 
        mutate(Region = names(table(Pop$smoking_status5yr))[w] )
      
      oftsrres[[q]] <- SurResults_smoke_temp[[q]]$one_five_ten_survival_rates %>% 
        mutate(Region = names(table(Pop$smoking_status5yr))[w] )
      
    }, error = function(e) {
      # Print a custom message
      cat("An error occurred. Skipping over.\n")
      # Print the error message
      print(e)
    })
    
  }
  
  
  # bind the results for calender years
  survival_results_cy <- bind_rows(surres)
  risk_table_cy <- bind_rows(rtres)
  med_surv_results_cy <- bind_rows(msres)
  survival_prob_cy <- bind_rows(oftsrres)
  
  calenderyr_results <- list(
    survival_results_cy,
    risk_table_cy,
    med_surv_results_cy,
    survival_prob_cy)
  
  names(calenderyr_results) <- c(paste0("survival_estimates_cy", db.name),
                                 paste0("risk_table_results_cy", db.name),
                                 paste0("median_survival_results_cy", db.name),
                                 paste0("one_five_ten_survival_rates_cy", db.name))
  
  
  #calender year stratification
  exportSurvivalResults(result=calenderyr_results,
                        zipName= paste0(db.name, "CalenderYrSurvivalResults", names(table(Pop$smoking_status5yr))[w]),
                        outputFolder=here::here("Results", db.name))
  
}



