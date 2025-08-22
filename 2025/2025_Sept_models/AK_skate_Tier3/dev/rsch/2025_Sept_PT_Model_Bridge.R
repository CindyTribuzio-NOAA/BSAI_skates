# Bridging Model 14_2d to the proposed model for 2025
# developed by C Tribuzio August 2025

# setup ----
libs <- c("r4ss", "here", "tidyverse")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

# Default with no version downloads the latest release
# r4ss::get_ss3_exe()

exe_loc <- here::here('2025/2025_Sept_models/AK_skate_Tier3/ss.exe')

##########
#base model (14_2d) brought over without re-running
#version update ran in M14_2d_ss3version

# Model 14_2d----
# Previously accepted model, but didn't converge
M14_2d_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/M14_2d')

# Model 14_2d1----
# base model with version update and recruitment ramp corrections
M14_2d1_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/M14_2d1')

r4ss::run(dir = M14_2d1_path, skipfinished = FALSE, exe = exe_loc)

M14_2d1_out <- SS_output(M14_2d1_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M14_2d1_out)

# Model 25_0----
# fixed growth
M25_0_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/M25_0')

r4ss::run(dir = M25_0_path, skipfinished = FALSE, exe = exe_loc)

M25_0_out <- SS_output(M25_0_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_0_out)



# Compare bridging models ----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3/2025_final_Sept_Models/")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("M14_2d", 
                                     "M14_2d_ss3version",
                                     "M25_0_bias"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/2025_final_Sept_Models/plots_compare')
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c('base', 
                                   'vers',
                                   'bias'
                                   ))
