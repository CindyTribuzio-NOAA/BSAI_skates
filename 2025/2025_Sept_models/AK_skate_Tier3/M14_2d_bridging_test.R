# Code to test run a few model changes to bridge to RTMB
# developed by C Tribuzio July 2025

# setup ----
libs <- c("r4ss", "here", "tidyverse")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

# Default with no version downloads the latest release
# r4ss::get_ss3_exe()

exe_loc <- here::here('2025/2025_Sept_models/AK_skate_Tier3/ss.exe')

# run old model with updated SS3 version----
vers_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run1_M14_2d_ss3version')

r4ss::run(dir = vers_mode_path, skipfinished = FALSE, exe = exe_loc)

mnew <- SS_output(vers_mode_path, printstats = FALSE, verbose = FALSE)
mnew_selexpars <- mnew$SelSizeAdj %>% 
  filter(Yr == 2023) %>% 
  mutate(model = 'updatedSS3')
mnew_selex <- mnew$sizeselex %>% 
  filter(Yr == 2023,
         Factor == 'Lsel') %>% 
  mutate(model = 'updatedSS3')

# compare old vs. new SS3 versions
oldv_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/base_M14_2d_fixedcatch')

mold <- SS_output(oldv_path, printstats = FALSE, verbose = FALSE)
mold_selexpars <- mold$SelSizeAdj %>% 
  filter(Yr == 2023) %>% 
  mutate(model = 'base')
mold_selex <- mold$sizeselex %>% 
  filter(Yr == 2023,
         Factor == 'Lsel') %>% 
  mutate(model = 'base')

mnew_out <- SS_output(dir = vers_mode_path, verbose = TRUE)

mnew_SS <- SSsummarize(mnew_out)

# plots the results
SS_plots(mnew_out)

# Arun1 fixed growth parameters based on model estimated with Amax = 20----
Arun1_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/agerun1_fixed_params_Amx20')
run(dir = Arun1_mode_path, skipfinished = FALSE, exe = exe_loc)

mArun1 <- SS_output(Arun1_mode_path, printstats = FALSE, verbose = FALSE)
mArun1_out <- SS_output(dir = Arun1_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mArun1_out)

# Arun2 fixed growth parameters based on model estimated with Amax = 26----
Arun2_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/agerun2_fixed_params_Amx26')
run(dir = Arun2_mode_path, skipfinished = FALSE, exe = exe_loc)

mArun2 <- SS_output(Arun2_mode_path, printstats = FALSE, verbose = FALSE)
mArun2_out <- SS_output(dir = Arun2_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mArun2_out)

# run4 Arun2 with wider length selectivity ----
run4_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run4_Arun2_slxbnds')
run(dir = run4_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun4 <- SS_output(run4_mode_path, printstats = FALSE, verbose = FALSE)
mrun4_out <- SS_output(dir = run4_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mrun4_out)

#overall comparison of model runs----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch", 
                                     "agerun1_fixed_params_Amx20", 
                                     "agerun2_fixed_params_Amx26",
                                     "run4_Arun2_slxbnds"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")

model_comp <- SSsummarize(bridge_out)

SSplotComparisons(model_comp,
                        print = TRUE,
                        plotdir = here::here(datapath),
                        legendlabels = c('base', 'Amax20', 'Amax26', 'slxbnds'))




























