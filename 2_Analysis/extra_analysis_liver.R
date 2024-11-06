# reviewer revisions

# read in lung cancer cohort
liver_cancer_outcome_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "OutcomeCohorts",
  "extra_analysis_liver"
))

liver_cancer_prevalent_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "PrevalentCohorts",
  "extra_analysis_liver"
  ))


#tablenames
livc_outcome_table_name <-paste0(outcome_table_stem,"_olivc") # for incidence
livc_prevalent_table_name <-paste0(outcome_table_stem,"_plivc") # for prevalence

instantiatedCohorts <- FALSE
if (instantiatedCohorts == TRUE) {
  
  cdm <- CDMConnector::cdm_from_con(con = db,
                                    cdm_schema = cdm_database_schema,
                                    write_schema = results_database_schema,
                                    cohort_tables = c(livc_outcome_table_name,
                                                      livc_prevalent_table_name))
  
  
} else {
#instantiate the lung cancer cohorts
# incidence
cdm <- CDMConnector::generateCohortSet(cdm, 
                                       cohortSet = liver_cancer_outcome_cohorts,
                                       name = livc_outcome_table_name,
                                       overwrite = TRUE
)

# prevalence
cdm <- CDMConnector::generateCohortSet(cdm = cdm, 
                                       cohortSet = liver_cancer_prevalent_cohorts,
                                       name = livc_prevalent_table_name,
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

inc <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator",
  outcomeTable = livc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = liver_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = liver_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 5,
  returnParticipants = FALSE,
  verbose = TRUE
)


# Estimate period prevalence ---------
prev_period <- estimatePeriodPrevalence(
  cdm = cdm,
  denominatorTable = "denominator",
  outcomeCohortId = liver_cancer_prevalent_cohorts$cohort_definition_id,
  outcomeCohortName = liver_cancer_prevalent_cohorts$cohort_name,
  outcomeLookbackDays = 0, 
  outcomeTable = livc_prevalent_table_name,
  interval = "years" ,
  completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
  fullContribution = FALSE , # individuals only required to be present for one day in interval
  minCellCount = 5,
  verbose = TRUE
)


# Get the results ----------------


study_results <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                  resultList=list(inc,prev_period ),
                                                  databaseName = db.name)



# Export the results -----


exportIncidencePrevalenceResults(result=study_results,
                                 zipName= paste0(db.name, "IPResults_extra_liver"),
                                 outputFolder=here::here("Results", db.name))


# age standardization

# for age standardization we need results without obscuring therefore if
# age standardization required it will run this code and save the results

if (agestandardization == TRUE) {
  
  inc <- estimateIncidence(
    cdm = cdm,
    denominatorTable = "denominator",
    outcomeTable = livc_outcome_table_name,
    denominatorCohortId = NULL,
    outcomeCohortId = liver_cancer_outcome_cohorts$cohort_definition_id,
    outcomeCohortName = liver_cancer_outcome_cohorts$cohort_name,
    interval = "years", 
    outcomeWashout = NULL,
    repeatedEvents = FALSE,
    completeDatabaseIntervals = TRUE,
    minCellCount = 0,
    returnParticipants = FALSE,
    verbose = TRUE
  )
  
  # Estimate period prevalence ---------
  prev_period <- estimatePeriodPrevalence(
    cdm = cdm,
    denominatorTable = "denominator",
    outcomeCohortId = liver_cancer_prevalent_cohorts$cohort_definition_id,
    outcomeCohortName = liver_cancer_prevalent_cohorts$cohort_name,
    outcomeLookbackDays = 0, 
    outcomeTable = livc_prevalent_table_name,
    interval = "years" ,
    completeDatabaseIntervals = TRUE, # prev only estimate for intervals where db captures all of the interval
    fullContribution = FALSE , # individuals only required to be present for one day in interval
    minCellCount = 0,
    verbose = TRUE
  )
  
  # Get the results ----------------
  study_results1 <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                     resultList=list(inc,prev_period ),
                                                     databaseName = db.name)
  
  exportIncidencePrevalenceResults(result=study_results1,
                                   zipName= paste0(db.name, "IPResultsAgeStandardization_extra_liver"),
                                   outputFolder=here::here("Results", db.name))
  

  
  
  
}