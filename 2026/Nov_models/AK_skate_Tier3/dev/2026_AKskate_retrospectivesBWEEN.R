#BETWEEN model retrospectives
# Attempt at mining historic files
# 2014 model 14_2: C:\Users\cindy.tribuzio\Work\SAFE\Assessments\AFSC_BSAI_SKATE_Assessment\Archived_Ormseth\Dropbox\BSAI skates\2014 BSAI skates\oct14_pref_model
# 2016 model 14_2: C:\Users\cindy.tribuzio\Work\SAFE\Assessments\AFSC_BSAI_SKATE_Assessment\Archived_Ormseth\Dropbox\BSAI skates\2016 BSAI skates\model_14_2_Oct_2016
# 2018 model 14_2: C:\Users\cindy.tribuzio\Work\SAFE\Assessments\AFSC_BSAI_SKATE_Assessment\Archived_Ormseth\Dropbox\BSAI skates\AK skate SS\AKSK_model_14_2_Oct2018
# 2020 model 14_2: C:\Users\cindy.tribuzio\Work\SAFE\Assessments\AFSC_BSAI_SKATE_Assessment\Archived_Ormseth\Dropbox\BSAI skates\AK skate SS\AKSK_model_14_2_OCT2020
# 2023 model 14_2d: C:\Users\cindy.tribuzio\Work\SAFE\Assessments\AFSC_BSAI_SKATE_Assessment\2023\Tier3\Model_Runs\M14_2d_fixedcatch
# 2026 author preferred model

# setup ----
libs <- c("r4ss", "here", "tidyverse", 'viridis', 'patchwork')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

# Default with no version downloads the latest release
# r4ss::get_ss3_exe()

AYR <- 2026

