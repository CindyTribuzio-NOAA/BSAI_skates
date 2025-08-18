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


# run5 base model with wider bounds on growth----
run5_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run5_growth_bnds')
run(dir = run5_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun5 <- SS_output(run5_mode_path, printstats = FALSE, verbose = FALSE)
mrun5_out <- SS_output(dir = run5_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mrun5_out)

#run7 estimating CV_young/CV_old-----
run7_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run7_A26_cvbnds')
run(dir = run7_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun7 <- SS_output(run7_mode_path, printstats = FALSE, verbose = FALSE)
mrun7_out <- SS_output(dir = run7_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mrun7_out)

#run8 estimating CV_young/CV_old with wide bounds-----
run8_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run8_A26_cvwidebnds')
run(dir = run8_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun8 <- SS_output(run8_mode_path, printstats = FALSE, verbose = FALSE)
mrun8_out <- SS_output(dir = run8_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun8_out)

#run9 fixing CV_young/CV_old-----
# values empirically derived from data for A0. 
# No CV for Amax 26, so set at 1
run9_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run9_A26_cvfixed')
run(dir = run9_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun9 <- SS_output(run9_mode_path, printstats = FALSE, verbose = FALSE)
mrun9_out <- SS_output(dir = run9_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun9_out)

#fixed growth comparison of model runs----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch", 
                                     "run5_growth_bnds",
                                     "agerun1_fixed_params_Amx20", 
                                     "agerun2_fixed_params_Amx26",
                                     "run7_A26_cvbnds",
                                     "run8_A26_cvwidebnds",
                                     "run9_A26_cvfixed"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)

SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(datapath),
                  legendlabels = c('base',  'grbnds', 'Amax20', 'Amax26', 'estCV', 'CVwidebnds', 'CVfixed'))

# Selectivity explorations----
# run4 Arun2 with wider length selectivity ----
run4_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run4_Arun2_slxbnds')
run(dir = run4_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun4 <- SS_output(run4_mode_path, printstats = FALSE, verbose = FALSE)
mrun4_out <- SS_output(dir = run4_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mrun4_out)

# run6 with new selex paramters-----
# and fixed descending limb at 1 for survey, but with old growth parameters
run6_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run6_arun2_slxparams')
run(dir = run6_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun6 <- SS_output(run6_mode_path, printstats = FALSE, verbose = FALSE)
mrun6_out <- SS_output(dir = run6_mode_path, verbose = TRUE)

#mArun1_SS <- SSsummarize(mArun1_out)

# plots the results
SS_plots(mrun6_out)

# run10 with slx parameters-----
# and fixed descending limb at 1 for survey and trawl fishery
# used S Barbeaux's control file
run10_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run10_run9slx')
run(dir = run10_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun10 <- SS_output(run10_mode_path, printstats = FALSE, verbose = FALSE)
mrun10_out <- SS_output(dir = run10_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun10_out)

# run11 with slx parameters-----
# and fixed descending limb at 1 for survey, allowed trawl fishery to be dome
# used S Barbeaux's control file
run11_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run11_run9slx')
run(dir = run11_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun11 <- SS_output(run11_mode_path, printstats = FALSE, verbose = FALSE)
mrun11_out <- SS_output(dir = run11_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun11_out)

# run12 with S-R parameters-----
# used S Barbeaux's control file
run12_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run12_run10rec')
run(dir = run12_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun12 <- SS_output(run12_mode_path, printstats = FALSE, verbose = FALSE)
mrun12_out <- SS_output(dir = run12_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun12_out)

# run13 with estimating q-----
# 
run13_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run13_run12estq')
run(dir = run13_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun13 <- SS_output(run13_mode_path, printstats = FALSE, verbose = FALSE)
mrun13_out <- SS_output(dir = run13_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun13_out)

# run13 allowing selectivity and fixq at values-----
# 
run14_mode_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run14_wdslx_estq')
run(dir = run14_mode_path, skipfinished = FALSE, exe = exe_loc)

mrun14 <- SS_output(run14_mode_path, printstats = FALSE, verbose = FALSE)
mrun14_out <- SS_output(dir = run14_mode_path, verbose = TRUE)

# plots the results
SS_plots(mrun14_out)


# Compare selected models ----
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch", 
                                     "run9_A26_cvfixed",
                                     "run10_run9slx",
                                     "run11_run9slx",
                                     "run12_run10rec",
                                     "run13_run12estq",
                                     "run14_wdslx_estq"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/plots_compare')
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c('base', 'fixgA26', 'SBslx', 'SBslx_dometwl', 'SR', 'q', 'slx+q=0.75'))


# kitchen sink model----
# see notes in 2025_AK_skate google doc
KS_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_run')
run(dir = KS_mod_path, skipfinished = FALSE, exe = exe_loc)

mrun14 <- SS_output(KS_mod_path, printstats = FALSE, verbose = FALSE)
mrun14_out <- SS_output(dir = KS_mod_path, verbose = TRUE)
mrun14_out$breakpoints_for_bias_adjustment_ramp

# plots the results
SS_plots(mrun14_out)

# kitchen sink model with rec ramp corrections----
# see notes in 2025_AK_skate google doc
KSbias_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_bias')
run(dir = KSbias_mod_path, skipfinished = FALSE, exe = exe_loc)

mrun15 <- SS_output(KSbias_mod_path, printstats = FALSE, verbose = FALSE)
mrun15_out <- SS_output(dir = KSbias_mod_path, verbose = TRUE)
mrun15_out$breakpoints_for_bias_adjustment_ramp

# plots the results
SS_plots(mrun15_out)

# kitchen sink rec ramp model ----
KSbias_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_bias')
run(dir = KSbias_mod_path, skipfinished = FALSE, exe = exe_loc)

mrun15 <- SS_output(KSbias_mod_path, printstats = FALSE, verbose = FALSE)
mrun15_out <- SS_output(dir = KSbias_mod_path, verbose = TRUE)

# plots the results
SS_plots(mrun15_out)

# kitchen sink changing slx ----
# estimating the logistic survey selectivity
KSslx_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_slx')
run(dir = KSslx_mod_path, skipfinished = FALSE, exe = exe_loc)

#mKSslx <- SS_output(KSslx_mod_path, printstats = FALSE, verbose = FALSE)
mKSslx_out <- SS_output(dir = KSslx_mod_path, verbose = TRUE)

# plots the results
SS_plots(mKSslx_out)

# kitchen sink changing survey slx back to double normal ----
KSsurvslx_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_survdnslx')
run(dir = KSsurvslx_mod_path, skipfinished = FALSE, exe = exe_loc)

mKSsurvslx_out <- SS_output(dir = KSsurvslx_mod_path, verbose = TRUE)

# plots the results
SS_plots(mKSsurvslx_out)

# kitchen sink short model ----
# start year 1980 with rec ramp correction
KSshort_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_short')
run(dir = KSshort_mod_path, skipfinished = FALSE, exe = exe_loc)

mrun16 <- SS_output(KSshort_mod_path, printstats = FALSE, verbose = FALSE)
mrun16_out <- SS_output(dir = KSshort_mod_path, verbose = TRUE)

# plots the results
SS_plots(mrun16_out)

# Compare kitchen sinks ----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch", 
                                     "kitchen_sink_run",
                                     "kitchen_sink_bias"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/plot_kitchen_sink_compare')
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c('base', 'kitchen sink', 'kitchen sink rec ramp'))




























