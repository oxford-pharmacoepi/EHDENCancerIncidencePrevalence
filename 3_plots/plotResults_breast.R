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
library(patchwork)

pathResults <- "C:/Users/dnewby/Desktop/"

ncases <- read.csv(here("3_Plots", "Breast cancer_cases.csv"))

# ncases <- ncases %>%
#   mutate(percentage = ifelse(sex == "Female", -percentage, percentage)) %>% 
#   mutate(shade = ifelse(percentage < 0, "Dark", "Light")) 
# 
# 
# ggplot(ncases, aes(x = Age.group, y = percentage, fill = sex)) +
#   geom_bar(stat = "identity", width = 1, color = "black") +
#   coord_flip() +
#   scale_y_continuous(labels = abs, expand = c(0.1, 0.1)) +
#   labs( x = NULL, y = NULL) +
#   theme_minimal() +
#   theme(
#     panel.grid.major = element_blank(),  # Remove major grid lines
#     panel.grid.minor = element_blank() ,
#     legend.position = "none"   ,
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     plot.margin = margin(7, 20, 7, 20)
#   ) +
#   # Labels for the left side (negative values)
#   geom_text(data = subset(ncases, percentage < 0), 
#             aes(label = Age.group),
#             hjust = 1.2,  # Place text to the left of the bar
#             color = "black",
#             size = 3.5) +
#   # Labels for the right side (positive values)
#   geom_text(data = subset(ncases, percentage >= 0), 
#             aes(label = Age.group),
#             hjust = -0.2,  # Place text to the right of the bar
#             color = "black",
#             size = 3.5) 




ncases <- ncases %>%
  mutate(percentage = ifelse(sex == "Female", -percentage, percentage)) %>% 
  mutate(shade = rep(c("Dark", "Light"), length.out = n())) %>%
  mutate(fill = interaction(sex, shade)) 

ncases$Age.group <- factor(ncases$Age.group, levels = rev(levels(factor(ncases$Age.group))))


# Define colors for each sex and shade
colors <- c("Male.Dark" = "#569DD1", "Male.Light" = "#86B8DE", 
            "Female.Dark" = "#DE5F9F", "Female.Light" = "#E789B8")


# Create the butterfly plot with different shades for each bar
ncases1 <- ggplot(ncases, aes(x = Age.group, y = percentage, fill = interaction(sex, shade))) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_flip() +
  scale_y_continuous(labels = abs, expand = c(0.1, 0.1)) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none",  # Remove the legend
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    plot.margin = margin(5, 30, 5, 30)  # Increase margins further if needed
  ) +
  geom_text(data = subset(ncases, percentage < 0), 
            aes(label = Age.group),
            hjust = 1.2,  # Place text to the left of the bar
            color = "black",
            size = 3.5) +
  # Labels for the right side (positive values)
  geom_text(data = subset(ncases, percentage >= 0), 
            aes(label = Age.group),
            hjust = -0.2,  # Place text to the right of the bar
            color = "black",
            size = 3.5) 


plotname <- paste0("breast_numbers.png")

#png(paste0(pathResults ,"/ExtraPlots/", plotname), width = 8, height = 5, units = "in", res = 1200)
png(paste0("C:/Users/dnewby/OneDrive - Nexus365/Desktop/" , plotname), width = 6, height = 5, units = "in", res = 1200)

print(ncases1, newpage = FALSE)
dev.off()

