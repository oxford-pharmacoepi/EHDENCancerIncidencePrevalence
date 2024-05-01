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


#tablenames
lc_outcome_table_name <-paste0(outcome_table_stem,"_olc") # for incidence
lc_prevalent_table_name <-paste0(outcome_table_stem,"_plc") # for prevalence

instantiatedCohorts <- TRUE
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
#%>%
#compute()

class(cdm$denominator) <- c("IncidencePrevalenceDenominator", "cohort_reference" , class(cdm$denominator))


# create a denominator for each region
cdm$denominator_england <- cdm$denominator %>% 
  filter(region_collapsed == "England") 
#%>% 
 # compute()
#update class
#class(cdm$denominator_england) <- c("IncidencePrevalenceDenominator", "cohort_reference" , class(cdm$denominator_england))

cdm$denominator_scotland <- cdm$denominator %>% 
  filter(region_collapsed == "Scotland") 
#%>% 
 # compute()
#update class
#class(cdm$denominator_scotland) <- c("IncidencePrevalenceDenominator","cohort_reference" ,  class(cdm$denominator_scotland))

cdm$denominator_wales <- cdm$denominator %>% 
  filter(region_collapsed == "Wales") 
#%>% 
 # compute()
#update class
#class(cdm$denominator_wales) <- c("IncidencePrevalenceDenominator", "cohort_reference" ,  class(cdm$denominator_wales))

cdm$denominator_ni <- cdm$denominator %>% 
  filter(region_collapsed == "Northern Ireland") 
#%>% 
#  compute()
#update class
#class(cdm$denominator_ni) <- c("IncidencePrevalenceDenominator", "cohort_reference" , class(cdm$denominator_ni))

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
# # add on column to the denominator with those with a code for smoking and filter
cdm$denominator_smokers <- cdm$denominator %>%
  left_join(
    cdm$denominator %>%
      select("subject_id", "cohort_start_date") %>%
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
                 by=c("subject_id"="person_id" ), copy = TRUE)  %>%
      distinct() %>%
      group_by(person_id, cohort_start_date) %>%
      filter(observation_date == max(observation_date)) %>%
      rename("Smoker_date"="observation_date"),
    by= c("person_id" = "subject_id", "cohort_start_date")) 

# filter denominator to only include those with a date of smoking
cdm$denominator_smokers <- cdm$denominator_smokers %>%
  filter(!is.na(Smoker_date))
# 
# # stratification by smoking status (any time)
# # inc
# inc_smoking <- estimateIncidence(
#   cdm = cdm,
#   denominatorTable = "denominator_smokers",
#   outcomeTable = lc_outcome_table_name,
#   denominatorCohortId = NULL,
#   outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
#   outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
#   interval = "years", 
#   outcomeWashout = NULL,
#   repeatedEvents = FALSE,
#   completeDatabaseIntervals = TRUE,
#   minCellCount = 0,
#   returnParticipants = FALSE,
#   verbose = TRUE
# )
# 
# inc_smoking_overall <- estimateIncidence(
#   cdm = cdm,
#   denominatorTable = "denominator_smokers",
#   outcomeTable = lc_outcome_table_name,
#   denominatorCohortId = NULL,
#   outcomeCohortId = lung_cancer_outcome_cohorts$cohort_definition_id,
#   outcomeCohortName = lung_cancer_outcome_cohorts$cohort_name,
#   interval = "overall", 
#   outcomeWashout = NULL,
#   repeatedEvents = FALSE,
#   completeDatabaseIntervals = TRUE,
#   minCellCount = 0,
#   returnParticipants = FALSE,
#   verbose = TRUE
# )
# 
# #prev
# prev_period_smoking <- estimatePeriodPrevalence(
#   cdm = cdm,
#   denominatorTable = "denominator_smokers",
#   outcomeCohortId = lung_cancer_prevalent_cohorts$cohort_definition_id,
#   outcomeCohortName = lung_cancer_prevalent_cohorts$cohort_name,
#   outcomeLookbackDays = 0, 
#   outcomeTable = prevalent_table_name,
#   interval = "years" ,
#   completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
#   fullContribution = FALSE , # individuals only required to be present for one day in interval
#   minCellCount = 5,
#   verbose = TRUE
# )
# 
# # gather results
# study_results_smoking <- gatherIncidencePrevalenceResults(cdm =cdm, 
#                                                   resultList=list(inc_smoking,prev_period_smoking ),
#                                                   databaseName = db.name)
# 
# study_results_smoking_overall<- gatherIncidencePrevalenceResults(cdm =cdm,
#                                                          resultList=list(inc_overall_smoking),
#                                                          databaseName = db.name)
# 
# # export results
# exportIncidencePrevalenceResults(result=study_results_smoking,
#                                  zipName= paste0(db.name, "IPResults_smoking"),
#                                  outputFolder=here::here("Results", db.name))
# 
# exportIncidencePrevalenceResults(result=study_results_smoking_overall,
#                                  zipName= paste0(db.name, "IPResults_overall_smoking"),
#                                  outputFolder=here::here("Results", db.name))