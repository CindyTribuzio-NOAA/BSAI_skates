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
new_year <- as.numeric(format(Sys.Date(), format = "%Y"))

# load functions
source(here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', "R", "mdl_anlys_fcns.r"))

# Model names and directory
Model_name_14 <- "M14_2d"
M14_2d_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/M14_2d')

# Run management scenarios----
M12_2d_projrun <- run_ss3_mgmnt_scen(dir = M14_2d_path,
                   cyr = new_year,
                   syr = 1950,
                   sexes = 1,
                   fleets = c(1:3), 
                   Scenario2 = 3, 
                   #s2_F = 0.4, if Scenario2 = 3 this should not be used
                   s4_F = 0.5, #F = 50% of maxFABC, as per the text in the last assessment. NOTE: las proj dat file says 0.6
                   output_name = 'M14_2d_proj')


