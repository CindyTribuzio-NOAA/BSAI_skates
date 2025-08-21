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

# kitchen sink rec ramp model ----
KSrecr_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/Model_14_2d_fixes/kitchen_sink_recramp')
run(dir = KSrecr_mod_path, skipfinished = FALSE, exe = exe_loc)

mKSrecr_out <- SS_output(dir = KSrecr_mod_path, verbose = TRUE)

# plots the results
SS_plots(mKSrecr_out)

# kitchen sink changing slx ----
# estimating the logistic survey selectivity
KSslx_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/Model_14_2d_fixes/kitchen_sink_slx')
run(dir = KSslx_mod_path, skipfinished = FALSE, exe = exe_loc)

#mKSslx <- SS_output(KSslx_mod_path, printstats = FALSE, verbose = FALSE)
mKSslx_out <- SS_output(dir = KSslx_mod_path, verbose = TRUE)

# plots the results
SS_plots(mKSslx_out)

# kitchen sink changing survey slx back to double normal ----
KSsurvslx_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/Model_14_2d_fixes/kitchen_sink_survdnslx')
run(dir = KSsurvslx_mod_path, skipfinished = FALSE, exe = exe_loc)

mKSsurvslx_out <- SS_output(dir = KSsurvslx_mod_path, verbose = TRUE)

# plots the results
SS_plots(mKSsurvslx_out)

# kitchen sink with dn surv and wider bounds on fishery slx ----
KSallslx_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/rsch/Model_14_2d_fixes/kitchen_sink_allslx')
run(dir = KSallslx_mod_path, skipfinished = FALSE, exe = exe_loc)

mKSallslx_out <- SS_output(dir = KSallslx_mod_path, verbose = TRUE)

# plots the results
SS_plots(mKKSallslx_out)


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




























