# Bridging Model 14_2d to the proposed model for 2025
# developed by C Tribuzio October 2026

#TODO
# switch "a" models to main models
# run 25_3 with growth estimated
# data peels on 25_3a
# M sensitivities

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
# Model 14_2d1_fixedslx----
# fixed all selectivities to see if that changed outcome
# note: did this because when I fixed the selectivies in 25_2 to be the same as 14_d1, fit was different
M14_2d1_fs_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M14_2d1_fixedslx')

r4ss::run(dir = M14_2d1_fs_path, skipfinished = FALSE, exe = exe_loc)

M14_2d1_fs_out <- SS_output(M14_2d1_fs_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M14_2d1_fs_out)

##########
# Model 25_2a----
# growth from 14_2d1_fixedslx_cvgr
# selectivities from 25_2
M25_2a_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_2a')

r4ss::run(dir = M25_2a_path, skipfinished = FALSE, exe = exe_loc)

M25_2a_out <- SS_output(M25_2a_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_2a_out)

##########
# Model 25_2----
# original 25_2 from Sept 2025
# Q = 0.836 (or something close to that)
M25_2_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_2')

r4ss::run(dir = M25_2_path, skipfinished = FALSE, exe = exe_loc)

M25_2_out <- SS_output(M25_2_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_2_out)

##########
# Model 25_3a----
# growth from 14_2d1_fixedslx_cvgr
# selectivities from 25_2
M25_3a_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_3a')

r4ss::run(dir = M25_3a_path, skipfinished = FALSE, exe = exe_loc)

M25_3a_out <- SS_output(M25_3a_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_3a_out)

##########
# Model 25_3----
# original from Sept 2025
# Q = 1
M25_3_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_3')

r4ss::run(dir = M25_3_path, skipfinished = FALSE, exe = exe_loc)

M25_3_out <- SS_output(M25_3_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_3_out)

##########
# Model 25_3a_growth----
# growth from 14_2d1_fixedslx_cvgr
# selectivities from 25_2
# letting model estimate growth
M25_3ag_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_3a_growth')

r4ss::run(dir = M25_3ag_path, skipfinished = FALSE, exe = exe_loc)

M25_3ag_out <- SS_output(M25_3ag_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_3ag_out)




# Compare models ----
# bridge
datapath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt')
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("M14_2d1",
                                     "M25_3a",
                                     "M25_3a_growth"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'model_comparison')
if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c("M14_2d1",
                                   "M25_3a",
                                   "M25_3a_growth"
                                   ))

#catchability
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("M25_2",
                                     "M25_2a",
                                     "M25_3a",
                                     "M25_3"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'model_comparison', 'catchability')
if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c("M25_2",
                                   "M25_2a",
                                   "M25_3a",
                                   "M25_3"))

#mortality
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("M25_2a_M",
                                     "M25_2a",
                                     "M25_3a",
                                     "M25_3a_M"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'model_comparison', 'mortality')
if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c("M25_2a_M",
                                   "M25_2a",
                                   "M25_3a",
                                   "M25_3a_M"))

