# Script to run the AK skate projection scenarios
# 2025 by C Tribuzio

# Setup ----
# load required packages and code
libs <- c("data.table",
          'tidyverse',
          "magrittr", 
          "r4ss",
          'foreach')

if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}

lapply(libs, library, character.only = TRUE)

# Current assessment year
new_year <- 2023

# load functions
source(here::here(new_year, "R", 'AK_Scenarios_For_SS', "Running_scenarios_Update.r"))

# Model names and directory
#Model_name_14 <- "M14_2d1"
#M14_2d_path <- here::here('2023/rsch/M14_2d1')

# Run management scenarios----
M14_2d_projrun <- Do_AK_TIER_3_Scenarios(DIR = here::here('2023/rsch/M14_2d1/PROJ'), 
                                         CYR = 2023, 
                                         SYR = 1977,  
                                         SEXES = 1, 
                                         FLEETS = c(1:3), 
                                         Scenario2 = 1, 
                                         S2_F = 0.4, 
                                         s4_F = 0.75, 
                                         do_fig = TRUE, 
                                         do_mark=FALSE, 
                                         URL='https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2024/BSAIskate.pdf')
M14_2d_projrun$FIGS


Do_AK_TIER_3_Scenarios(DIR = "M14_2d1/PROJ", 
                       CYR = 2023, 
                       SYR = 1977,  
                       SEXES = 1, 
                       FLEETS = c(1:2), 
                       Scenario2 = 1, 
                       S2_F = 0.4, 
                       s4_F = 0.75, 
                       do_fig = TRUE, 
                       do_mark=FALSE, 
                       URL='https://www.npfmc.org/wp-content/PDFdocuments/SAFE/2024/BSAIskate.pdf')