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

# compare old vs. new SS3 versions applied to 2017 inputs
oldv_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/base_M14_2d_fixedcatch')

mold <- SS_output(oldv_path, printstats = FALSE, verbose = FALSE)
mnew <- SS_output(vers_mode_path, printstats = FALSE, verbose = FALSE)
SStableComparisons(
  SSsummarize(list(mold, mnew)),
  names = c("NatM", "SSB_Virgin", "SSB_2017", "Bratio_2017"),
  likenames = NULL,
  modelnames = c("3.30.21.00", "3.30.23.2"),
) |>
  dplyr::mutate(across(-1, ~ round(.x, 3)))

# change fishery selectivity to logistic----


