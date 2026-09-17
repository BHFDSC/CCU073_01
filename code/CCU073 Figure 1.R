######################################################################################
## TITLE: Create figure 1 for models with categorical exposure
##
## Author: Genevieve Cezard
######################################################################################
library(readr)
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("2. Projects/2. CVD-COVID-UK/CCU073 project/Results/2. Cox results")


results_dir <- "2. Projects/2. CVD-COVID-UK/CCU073 project/Results/2. Cox results"
output_dir <- "2. Projects/2. CVD-COVID-UK/CCU073 project/Results/2. Cox results/Figures/"


############################
# Set 1 - Adjusted for age #
############################

#-----------#
# Read data #
#-----------#
df_covid <- readr::read_csv("Models - T2DM age cat/HRs_covid_expo_DDSC_T2DM_noCVDhist_adj_set1.csv", show_col_types = FALSE)
df_non_covid <- readr::read_csv("Models - T2DM age cat/HRs_non_covid_expo_DDSC_T2DM_noCVDhist_adj_set1.csv", show_col_types = FALSE)

df_covid$Covid_history <- "Yes"
df_non_covid$Covid_history <- "No"

df <- rbind(df_covid,df_non_covid)


#---------------------------------------#
# Select results of interest and format #
#---------------------------------------#
estimates <- subset(df, variable == "expo_T2DM_DDSC_10y_age_cat")

estimates$ucl[is.na(estimates$ucl) & estimates$category =="No diabetes"] <- 1
estimates$lcl[is.na(estimates$lcl) & estimates$category =="No diabetes"] <- 1

# Add formatting specifications to the dataset
estimates$colour <- ""
estimates$colour <- factor(estimates$colour, levels=c("#FF0000","#4169E1"))
estimates$category <- factor(estimates$category, levels = c("No diabetes","0-29","30-39","40-49","50-59","60-69","70-79","80+"))
estimates$category <- relevel(estimates$category, ref="No diabetes")
estimates$Covid_history <- factor(estimates$Covid_history, levels=c("Yes","No"))
estimates$outcome <- factor(estimates$outcome, levels = c('Angina','CHD','ARR','HF','stroke','PVD','DVT','s1','s2','s3'))
estimates$outcome_labels <- fct_recode(estimates$outcome,
                                   "Angina (53,230 cases)"="Angina",
                                   "CHD (181,735 cases)"="CHD",
                                   "Arrhythmia (380,960 cases)"="ARR",
                                   "Heart Failure (232,175 cases)"="HF",
                                   "Stroke (168,885 cases)"="stroke",
                                   "PVD (112,650 cases)"="PVD",
                                   "DVT (75,440 cases)"="DVT",
                                   "CVD 1 (936,380 cases)"="s1",
                                   "CVD 2 (904,620 cases)"="s2",
                                   "CVD 3 (524,040 cases)"="s3" )
#--------------#
# Plot results #
#--------------#
y_lim <- c(0.2,16)
y_lim_breaks <- c(0.25,0.5,1,2,4,8,16)
ggplot2::ggplot(data = estimates, 
                mapping = ggplot2::aes(x = category, y = HR, color = Covid_history, shape = Covid_history, fill = Covid_history)) +
  ggplot2::labs(color = "COVID-19 history", shape = "COVID-19 history", fill = "COVID-19 history") +
  ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = 1), colour = "#A9A9A9") +
  ggplot2::geom_point(size=2)+
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lcl, ymax = ucl),width = 0.2)+
  ggplot2::labs(x = "Age at diagnosis of T2DM diabetes group", y = "HR (95% CI) vs. No Diabetes") +
  ggplot2::scale_y_continuous(lim = y_lim, breaks = y_lim_breaks, trans = "log") +
  ggplot2::scale_fill_manual(values = levels(estimates$colour), labels = levels(estimates$Covid_history)) +
  ggplot2::scale_color_manual(values = levels(estimates$colour), labels = levels(estimates$Covid_history)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                 legend.position = "bottom",legend.direction = "horizontal",
                 axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))+
  ggplot2::facet_wrap(outcome_labels~., ncol=5)

ggplot2::ggsave(paste0(output_dir,"Figure_T2DMagecat_set1_v2.pdf"), height = 175, width = 250, unit = "mm", dpi = 300, scale = 1)
ggplot2::ggsave(paste0(output_dir,"Figure_T2DMagecat_set1_v2.jpeg"), height = 175, width = 250, unit = "mm", dpi = 300, scale = 1)
dev.off()
