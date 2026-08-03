# Bridging Model 14_2d to the proposed model for 2025
# developed by C Tribuzio October 2026

#TODO
# confirm 25_2a
# run 25_3 and 25_3a
# run both with high M

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
M14_2d1_fs_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'M14_2d1_fixedslx')

r4ss::run(dir = M14_2d1_fs_path, skipfinished = FALSE, exe = exe_loc)

M14_2d1_fs_out <- SS_output(M14_2d1_fs_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M14_2d1_fs_out)

# Model 14_2d1_fixedslxgrowth----
# fixed all selectivities like 14_2d1
# fixing growth and CV on age like 25_2
M14_2d1_fsg_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'M14_2d1_fixedslxgrowth')

r4ss::run(dir = M14_2d1_fsg_path, skipfinished = FALSE, exe = exe_loc)

M14_2d1_fsg_out <- SS_output(M14_2d1_fsg_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M14_2d1_fsg_out)

# Model 14_2d1_fixedslx_cvgr----
# fixed all selectivities like 14_2d1
# fixed CV on age like 14_2d1
# fixed growth parameters from outside the model
M14_2d1_fscv_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'M14_2d1_fixedslx_cvgr')

r4ss::run(dir = M14_2d1_fscv_path, skipfinished = FALSE, exe = exe_loc)

M14_2d1_fscv_out <- SS_output(M14_2d1_fscv_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M14_2d1_fscv_out)


##########
# Model 25_2a----
# growth from 14_2d1_fixedslx_cvgr
# selectivities from 25_2
M25_2a_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'M25_2a')

r4ss::run(dir = M25_2a_path, skipfinished = FALSE, exe = exe_loc)

M25_2a_out <- SS_output(M25_2a_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_2a_out)

# Model 25_2M----
# only changed NatM ~0.2 to see if that fixed anything
M25_2M_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'M25_2M')

r4ss::run(dir = M25_2M_path, skipfinished = FALSE, exe = exe_loc)

M25_2M_out <- SS_output(M25_2M_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_2M_out)

# Model 25_3M----
# only changed NatM ~0.2 to see if that fixed anything
M25_3M_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'M25_3M')

r4ss::run(dir = M25_3M_path, skipfinished = FALSE, exe = exe_loc)

M25_3M_out <- SS_output(M25_3M_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M25_3M_out)




# Compare models ----
datapath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch')
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("M14_2d1",
                                     "M25_2",
                                     "M25_3",
                                     "M25_2M",
                                     "M25_3M"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'rsch', 'model_comparison')
if (!dir.exists(plotpath)) dir.create(plotpath, recursive = TRUE)
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c("M14_2d1",
                                   "M25_2",
                                   "M25_3",
                                   "M25_2M",
                                   "M25_3M"
                                   ))
