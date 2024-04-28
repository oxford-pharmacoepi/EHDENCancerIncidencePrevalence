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


# stratification by smoking status (5 years)






# stratification by location (wales, scotland, england, NI)