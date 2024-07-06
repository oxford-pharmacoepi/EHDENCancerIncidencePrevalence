library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(ggplot2)
library(scales)
library("cowplot")
library("gridExtra")
library("ggpubr")
library(ggh4x)

pathResults <- "C:/Users/dnewby/Desktop/Results"
datapath <- "C:/Users/dnewby/OneDrive - Nexus365/Documents/GitHub/EHDENCancerIncidencePrevalence/CancerIncidencePrevalanceShiny/shiny/data"

# read in files
prevalence_estimates <- readRDS(paste0(datapath ,"/prevalence_estimates.rds")) %>% 
  mutate(database_name = ifelse(database_name == "CPRDGoldUpdate2", "CPRD GOLD", database_name))

prevalence_attrition <- readRDS(paste0(datapath ,"/prevalence_attrition.rds")) %>% 
  mutate(database_name = ifelse(database_name == "CPRDGoldUpdate2", "CPRD GOLD", database_name))

incidence_estimates <- readRDS(paste0(datapath ,"/incidence_estimates.rds")) %>% 
  mutate(database_name = ifelse(database_name == "CPRDGoldUpdate2", "CPRD GOLD", database_name))

incidence_attrition <- readRDS(paste0(datapath ,"/incidence_attrition.rds")) %>% 
  mutate(database_name = ifelse(database_name == "CPRDGoldUpdate2", "CPRD GOLD", database_name))

survival_estimates <- readRDS(paste0(datapath ,"/survival_estimates.rds"))%>% 
  rename(CalendarYearGp = CalenderYearGp ) %>% 
  mutate(Database = ifelse(Database == "CPRDGoldUpdate2", "CPRD GOLD", Database))

# KM for prostate cancer
survivalFigureData <- survival_estimates %>%
  filter(Stratification == "None") %>%
  filter(CalendarYearGp == "2000 to 2019" |
           CalendarYearGp == "2000 to 2021"  ) %>%
  filter(Cancer == "Prostate") %>%
  filter(Age == "All") %>%
  ggplot(aes(x = time,
             y = est,
             group = Database,
             col = Database )) +
  scale_y_continuous( labels = scales::percent, limits = c(0, NA)) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl, 
                  fill = Database), alpha = .15, color = NA, show.legend = FALSE) +
  geom_line(aes(linetype = Database),size = 0.5) +
  scale_linetype_manual(values = c("solid", "dashed", "twodash","dotted")) +
  labs(x = "Time (Years)",
       y = "Survival Probability",
       col = "Database name",
       linetype = "Database name") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  scale_x_continuous(breaks=seq(0, 22, 2))


survival_risk_table <- read.csv(here("3_Plots", "/survival_risk_table_prostate.csv"))


survival_risk_table_prostate <- survival_risk_table %>% 
  filter(details == "n.risk") %>% 
  select(-c(details, Cancer, Sex, Age, Calendar.Year, Database)) 

colnames(survival_risk_table_prostate) <- NULL
#rownames(survival_risk_table_prostate) <- c("CPRD GOLD", "CPRD Aurum")

rownames(survival_risk_table_prostate) <- c("", " ")

table_theme <- ttheme_minimal(
  core = list(fg_params = list(cex = 0.8, fontface = "plain")),  # Change core text size and fontface
  colhead = list(fg_params = list(cex = 0.8, fontface = "plain")),  # Change column header text size and fontface
  rowhead = list(fg_params = list(cex = 0.8, fontface = "plain"))#  ,
  #,  # Change row header text size and fontface
 # padding = unit(c(2, 3), "mm")
  
  
)

# Define column widths based on plot x-axis ranges
# Adjust these ranges according to your actual data
table_grob <- gridExtra::tableGrob(survival_risk_table_prostate, theme = table_theme)

#table_grob$widths <- unit_widths

table_grob$widths <- rep(max(table_grob$widths), (length(table_grob$widths)))

# 
# # Define relative column widths based on plot width

# Combine plot and table using ggarrange
combined_plot <- ggarrange(
  plotlist = list(survivalFigureData, table_grob),
  ncol = 1, nrow = 2,
  heights = c(2, 0.2) 
)


combined_plot <- combined_plot +
  annotation_custom(grob = grid::textGrob("Number at Risk", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.28, xmax = 0, ymin = 0.18, ymax = 0) +
  annotation_custom(grob = grid::textGrob("GOLD", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.125, ymax = 0) +
  annotation_custom(grob = grid::textGrob("Aurum", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.06, ymax = 0)

# Print the combined plot
#print(combined_plot)



#756 x 517

plotname <- paste0("FIGURE5_KMSurvival_prostate.pdf")

pdf(paste0(datapath ,"/", plotname), width = 8, height = 8)

print(combined_plot, newpage = FALSE)
dev.off()








################################################

################################################

# incidence whole population
incidenceFigureData <- incidence_estimates %>%
  filter(denominator_sex == "Male",
         denominator_age_group == "All",
         analysis_interval == "years" ,
         outcome_cohort_name == "Prostate") %>%
  ggplot(aes(x = incidence_start_date,
             y = incidence_100000_pys,
             group = database_name)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = incidence_100000_pys_95CI_lower, 
                  ymax = incidence_100000_pys_95CI_upper, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 3.5) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.background = element_blank() ,
        #axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name" ) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1))


plotname <- paste0("FIGURE1_Incidence_Males_Prostate.pdf")

pdf(paste0(datapath ,"/", plotname), width = 6, height = 5)

print(incidenceFigureData, newpage = FALSE)
dev.off()

# prevalence whole population
prevalenceFigureData <- prevalence_estimates %>%
  filter(denominator_sex == "Male",
         denominator_age_group == "All",
         outcome_cohort_name == "Prostate") %>%
  ggplot(aes(x = prevalence_start_date,
             y = prevalence,
             group = database_name)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = prevalence_95CI_lower, 
                  ymax = prevalence_95CI_upper, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 3.5) +
  scale_shape_manual(values = c(24,21)) +
  scale_y_continuous( labels = scales::percent, limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.background = element_blank() ,
        #axis.line = element_line(colour = "black", size = 0.6) ,
        panel.border = element_rect(colour = "black", fill=NA, size=0.6),
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Prevalence",
       col = "Database name",
       shape = "Database name",
       fill = "Database name" ) +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("2 years"),
               expand = c(0.06,1))


plotname <- paste0("FIGURE3_Prevalence_Males_Prostate.pdf")

pdf(paste0(datapath ,"/", plotname), width = 6, height = 5)

print(prevalenceFigureData, newpage = FALSE)
dev.off()

#################################################
# prostate cancer IR for different age groups
incidenceFigureData <- incidence_estimates %>%
  filter(denominator_sex == "Male",
         analysis_interval == "years" ,
         denominator_age_group != "All", 
         denominator_age_group != "18 to 29", 
         denominator_age_group != "30 to 39", 
         outcome_cohort_name == "Prostate") %>%
  ggplot(aes(x = incidence_start_date,
             y = incidence_100000_pys,
             group = database_name)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = incidence_100000_pys_95CI_lower, 
                  ymax = incidence_100000_pys_95CI_upper, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name),size = 1.5) +
  scale_shape_manual(values = c(24,21)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Incidence rate per 100000 person-years",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ denominator_age_group, scales = "free", ncol = 2)

plotname <- paste0("FIGURE2_Incidence_AgeStrat_Males_Prostate.pdf")
pdf(paste0(datapath ,"/", plotname), width = 6, height = 7)

print(incidenceFigureData , newpage = FALSE)
dev.off()

###################################################
######################################################

prevalenceFigureData <- prevalence_estimates %>%
  filter(denominator_sex == "Male",
         denominator_age_group != "All", 
         denominator_age_group != "18 to 29", 
         denominator_age_group != "30 to 39", 
         outcome_cohort_name == "Prostate") %>%
  ggplot(aes(x = prevalence_start_date,
             y = prevalence,
             group = database_name)) +
  geom_line(color = "black", size = 0.25) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark read, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_ribbon(aes(ymin = prevalence_95CI_lower, 
                  ymax = prevalence_95CI_upper, 
                  fill = database_name), alpha = .15, color = NA, show.legend = FALSE) +
  geom_point(aes(shape = database_name, fill = database_name ),size = 1.5) +
  scale_shape_manual(values = c(24,21)) +
  scale_y_continuous( labels = scales::percent, limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  labs(x = "Calendar year",
       y = "Prevalence",
       col = "Database name",
       shape = "Database name",
       fill = "Database name") +
  scale_x_date(labels = date_format("%Y"), breaks = date_breaks("4 years"),
               expand = c(0.06,1)) +
  facet_wrap(~ denominator_age_group, scales = "free", ncol = 2)

plotname <- paste0("FIGURE4_Prevalence_AgeStrat_Males_Prostate.pdf")
pdf(paste0(datapath ,"/", plotname), width = 6, height = 7)
print(prevalenceFigureData , newpage = FALSE)
dev.off()

########################################################################

survivalFigureData <- survival_estimates %>%
  filter(Age == "All") %>%
  filter(Cancer == "Prostate") %>% 
  filter(CalendarYearGp != "2000 to 2019") %>%
  filter(CalendarYearGp != "2000 to 2021") %>%
  
  ggplot(aes(x = time,
             y = est,
             group = CalendarYearGp,
             col = CalendarYearGp )) +
  scale_y_continuous( labels = label_percent() ) +
  scale_colour_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) + #blue, #red, #lightblue, #green, purple, peach, dark red, gry
  scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#0099B4FF", "#42B540FF", "#925E9FFF", "#FDAF91FF", "#AD002AFF", "grey")) +
  geom_line(aes(linetype = CalendarYearGp),size = 0.5) +
  scale_linetype_manual(values = c("dotted","dashed", "dotdash", "twodash","solid", "longdash")) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl, 
                  fill = CalendarYearGp), alpha = .15, color = NA, show.legend = FALSE) +
  labs(x = "Time (Years)",
       y = "Survival Probability",
       col = "Calendar Year Group",
       linetype = "Calendar Year Group") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6), 
        strip.background = element_rect(color = "black", size = 0.6) ,
        panel.background = element_blank() ,
        #axis.line = element_line(colour = "black", size = 0.6) ,
        panel.grid.major = element_line(color = "grey", size = 0.2, linetype = "dashed"),
        legend.box.spacing = unit(0, "pt") ,
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.position='bottom') +
  guides(col = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  ggh4x::facet_grid2(cols = vars(Database), scales="free", independent = "y") 

survival_risk_table_cy <- read.csv(here("3_Plots", "/survival_risk_table_by_calendaryr_prostate.csv"))


survival_risk_table_cy_prostate <- survival_risk_table_cy %>% 
  select(-c(Calendar.Year))
  
colnames(survival_risk_table_cy_prostate) <- NULL

rownames(survival_risk_table_cy_prostate) <- c("", 
                                            " ",
                                            "  ",
                                            "   ",
                                            "    "
                                            )

table_theme <- ttheme_minimal(
  core = list(fg_params = list(cex = 0.8, fontface = "plain")),  # Change core text size and fontface
  colhead = list(fg_params = list(cex = 0.8, fontface = "plain")),  # Change column header text size and fontface
  rowhead = list(fg_params = list(cex = 0.8, fontface = "plain")) # ,
  #,  # Change row header text size and fontface
  # padding = unit(c(2, 5), "mm")
  
  
)

# # Define column widths based on plot x-axis ranges
# # Adjust these ranges according to your actual data
# table_grob <- gridExtra::tableGrob(survival_risk_table_cy_prostate, theme = table_theme)
# 
# #table_grob$widths <- rep(max(table_grob$widths), (length(table_grob$widths)))
# 
# # 
# # # Define relative column widths based on plot width
# 
# # Combine plot and table using ggarrange
# combined_plot <- ggarrange(
#   plotlist = list(survivalFigureData, table_grob),
#   ncol = 1, nrow = 2,
#   heights = c(2, 0.6) 
# )

ncol_table <- ncol(survival_risk_table_cy_prostate)
table_grob$widths <- rep(grid::unit(1, "cm"), length(table_grob$widths)) 



# Combine the plot and table
combined_plot <- plot_grid(
  survivalFigureData,
  table_grob,
  ncol = 1,
  rel_heights = c(2, 0.6)
)


combined_plot



plotname <- paste0("FIGURE6_KMSurvival_prostate_cy.pdf")

pdf(paste0(datapath ,"/", plotname), width = 6, height = 6)

print(combined_plot, newpage = FALSE)
dev.off()



# combined_plot <- combined_plot +
#   annotation_custom(grob = grid::textGrob("Number at Risk", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.28, xmax = 0, ymin = 0.18, ymax = 0) +
#   annotation_custom(grob = grid::textGrob("2000 to 2004", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.125, ymax = 0) +
#   annotation_custom(grob = grid::textGrob("2005 to 2009", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.125, ymax = 0) +
#   annotation_custom(grob = grid::textGrob("2010 to 2014", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.125, ymax = 0) +
#   annotation_custom(grob = grid::textGrob("2015 to 2019", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.125, ymax = 0) +
#   annotation_custom(grob = grid::textGrob("2020 to 2021", gp = grid::gpar(fontsize = 10, fontface = "bold")),  xmin = 0.10, xmax = 0, ymin = 0.06, ymax = 0)
# 
# # Print the combined plot
# print(combined_plot)
