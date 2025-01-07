# reviewer revisions


# 1 # analysis by different regions on crc ie rectal and colon

# read in crc cancer cohort colon and rectal
crc_cancer_outcome_cohorts <- CDMConnector::readCohortSet(here::here(
  "1_InstantiateCohorts",
  "OutcomeCohorts",
  "extra_analysis_crc"
))

#tablenames
crc_outcome_table_name <-paste0(outcome_table_stem,"_ocrc") # for incidence


#instantiate the crc cancer cohorts
# incidence
cdm <- CDMConnector::generateCohortSet(cdm, 
                                       cohortSet = crc_cancer_outcome_cohorts,
                                       name = crc_outcome_table_name,
                                       overwrite = TRUE
)


# 2 # analysis with smaller age groups (add into denominator)

# take denominator 
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
    c(90, 150),
    
    c(18, 24),
    c(25, 29),
    c(30, 34),
    c(35, 39),
    c(40, 44),
    c(45, 49),
    c(50, 54),
    c(55, 59),
    c(60, 64),
    c(65, 69),
    c(70, 74),
    c(75, 79),
    c(80, 84),
    c(85, 89)
    
    
    
  ),
  sex = c("Male", "Female", "Both"),
  daysPriorHistory = 365,
  verbose = TRUE
)


# run incidence

inc <- estimateIncidence(
  cdm = cdm,
  denominatorTable = "denominator",
  outcomeTable = crc_outcome_table_name,
  denominatorCohortId = NULL,
  outcomeCohortId = crc_cancer_outcome_cohorts$cohort_definition_id,
  outcomeCohortName = crc_cancer_outcome_cohorts$cohort_name,
  interval = "years", 
  outcomeWashout = NULL,
  repeatedEvents = FALSE,
  completeDatabaseIntervals = TRUE,
  minCellCount = 5,
  returnParticipants = FALSE,
  verbose = TRUE
)



# Get the results ----------------


study_results <- gatherIncidencePrevalenceResults(cdm =cdm, 
                                                  resultList=list(inc),
                                                  databaseName = db.name)



# Export the results -----


exportIncidencePrevalenceResults(result=study_results,
                                 zipName= paste0(db.name, "IPResults_extra_crc"),
                                 outputFolder=here::here("Results", db.name))


# survival analysis
