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


#instantiate the lung cancer cohorts
# incidence
cdm <- CDMConnector::generateCohortSet(cdm, 
                                       cohortSet = outcome_cohorts,
                                       name = outcome_table_name,
                                       overwrite = TRUE
)

# prevalence
cdm <- CDMConnector::generateCohortSet(cdm = cdm, 
                                       cohortSet = prevalent_cohorts,
                                       name = prevalent_table_name,
                                       overwrite = TRUE
                                       
)


# stratification by smoking status (5 years)






# stratification by location (wales, scotland, england, NI)