# reviewer revisions

# read in oesophageal cancer cohort
oesophageal_cancer_outcome_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "OutcomeCohorts",
  "extra_analysis_oesophageal"
))

oesophageal_cancer_prevalent_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "PrevalentCohorts",
  "extra_analysis_oesophageal"
  ))


#tablenames
oc_outcome_table_name <-paste0(outcome_table_stem,"_olc") # for incidence
oc_prevalent_table_name <-paste0(outcome_table_stem,"_plc") # for prevalence

instantiatedCohorts <- FALSE
if (instantiatedCohorts == TRUE) {
  
  cdm <- CDMConnector::cdm_from_con(con = db,
                                    cdm_schema = cdm_database_schema,
                                    write_schema = results_database_schema,
                                    cohort_tables = c(oc_outcome_table_name,
                                                      oc_prevalent_table_name))
  
  
} else {
#instantiate the oesophageal cancer cohorts
# incidence
cdm <- CDMConnector::generateCohortSet(cdm, 
                                       cohortSet = oesophageal_cancer_outcome_cohorts,
                                       name = oc_outcome_table_name,
                                       overwrite = TRUE
)

# prevalence
cdm <- CDMConnector::generateCohortSet(cdm = cdm, 
                                       cohortSet = oesophageal_cancer_prevalent_cohorts,
                                       name = oc_prevalent_table_name,
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
  outcomeTable = oc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = oesophageal_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
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
  outcomeCohortId = oesophageal_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = oc_prevalent_table_name,
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

exportIncidencePrevalenceResults(result=study_results_eng,
                                 zipName= paste0(db.name, "IPResults_england_oc"),
                                 outputFolder=here::here("Results", db.name))



# SCOTLAND
inc_scot <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_scotland",
  outcomeTable = oc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = oesophageal_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
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
  outcomeCohortId = oesophageal_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = oc_prevalent_table_name,
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


exportIncidencePrevalenceResults(result=study_results_scot,
                                 zipName= paste0(db.name, "IPResults_scotland_oc"),
                                 outputFolder=here::here("Results", db.name))




# WALES
inc_wales <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_wales",
  outcomeTable = oc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = oesophageal_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
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
  outcomeCohortId = oesophageal_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = oc_prevalent_table_name,
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


exportIncidencePrevalenceResults(result=study_results_wales,
                                 zipName= paste0(db.name, "IPResults_wales_oc"),
                                 outputFolder=here::here("Results", db.name))

# NI
inc_ni <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator_ni",
  outcomeTable = oc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = oesophageal_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
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
  outcomeCohortId = oesophageal_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = oesophageal_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = oc_prevalent_table_name,
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


exportIncidencePrevalenceResults(result=study_results_ni,
                                 zipName= paste0(db.name, "IPResults_ni_oc"),
                                 outputFolder=here::here("Results", db.name))


