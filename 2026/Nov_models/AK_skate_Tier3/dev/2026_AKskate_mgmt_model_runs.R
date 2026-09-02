# Model 14_2d1 to the proposed model for 2026
# developed by C Tribuzio October 2026

#TODO
# update LOO
# run M sensitivies and growth option, not for doc, but for presentation

# setup ----
libs <- c("r4ss", "here", "tidyverse")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

# Default with no version downloads the latest release
# r4ss::get_ss3_exe()

AYR <- 2026

exe_loc <- here::here(paste0(AYR, '/Nov_models/AK_skate_Tier3/ss.exe'))

##########
# Model 14_2d1----
# base model
M14_2d1_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M14_2d1')

r4ss::run(dir = M14_2d1_path, skipfinished = FALSE, exe = exe_loc)

M14_2d1_out <- SS_output(M14_2d1_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M14_2d1_out)

##########
# Model 25_3----
# From Sept 2025, but this version fixed the mis-specified CV young/old parameters, so not a new model
# Q = 1
M25_3_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_3')

r4ss::run(dir = M25_3_path, skipfinished = FALSE, exe = exe_loc)

M25_3_out <- SS_output(M25_3_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_3_out)

# Compare models ----
# bridge
datapath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt')
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("M14_2d1",
                                     "M25_3"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'model_comparison')
if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c("M14_2d1",
                                   "M25_3"
                                   ))

