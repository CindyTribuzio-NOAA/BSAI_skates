# Breaking down the various component of the Kitchen Sink model
# individual changes, not sequentially building
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

# Model 14_2d + bias correction and updated SS3 version----
M142_v_bias_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/M14_2d_vers_bias')

r4ss::run(dir = M142_v_bias_path, skipfinished = FALSE, exe = exe_loc)

M142_v_bias_out <- SS_output(M142_v_bias_path, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M142_v_bias_out)

# Model 14_2d + vers/bias and jitter/lamda----
# jitter changed from 0.01 to 0.1 in starter file
# #_maxlambdaphase from 4 to 1 in control file
M142_jitter <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/M14_2d_vb_jitter_lambda')

r4ss::run(dir = M142_jitter, skipfinished = FALSE, exe = exe_loc)

M142_jitter_out <- SS_output(M142_jitter, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M142_jitter_out)

# Model 14_2d + vers/bias and fixed growth----
# Schnute growth params from AK_skate_Richards_growth analysis
# cv_young and CV_old from AK_skate_Richards_growth analysis
# sd_offset from 1 to 0 in control file
M142_growth <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/M14_2d_vb_growth')

r4ss::run(dir = M142_growth, skipfinished = FALSE, exe = exe_loc)

M142_growth_out <- SS_output(M142_growth, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M142_growth_out)

# Model 14_2d + vers/bias and fixed catchability----
# Q = ln(0.836) as per Kotwicki and Weinberg 2005
M142_q <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/M14_2d_vb_q')

r4ss::run(dir = M142_q, skipfinished = FALSE, exe = exe_loc)

M142_q_out <- SS_output(M142_q, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M142_q_out)

# Model 14_2d + vers/bias and fixed survey slx----
# fixed at values from 2007 which were empirically derived from Kotwiki and Weinberg 2005
M142_survslx <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/M14_2d_vb_fixsurveyslx')

r4ss::run(dir = M142_survslx, skipfinished = FALSE, exe = exe_loc)

M142_survslx_out <- SS_output(M142_survslx, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M142_survslx_out)

# Model 14_2d + vers/bias and fixed top param fishery slx----
# widened bounds on fishery selectivity
# fixed top parameter at -6, when it's negative, the value doesn't matter?
M142_fshslx <- here::here('2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/M14_2d_vb_fisheryslx')

r4ss::run(dir = M142_fshslx, skipfinished = FALSE, exe = exe_loc)

M142_fshslx_out <- SS_output(M142_fshslx, printstats = FALSE, verbose = FALSE)

# plots the results
SS_plots(M142_fshslx_out)

# Compare kitchen sinks ----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch", 
                                     "M14_2d_vers_bias",
                                     "M14_2d_vb_jitter_lambda",
                                     "M14_2d_vb_growth",
                                     "M14_2d_vb_q",
                                     "M14_2d_vb_fixsurveyslx",
                                     "M14_2d_vb_fisheryslx"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")
model_comp <- SSsummarize(bridge_out)
plotpath <- paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/kitchen_sink_decompose/plots_compare')
SSplotComparisons(model_comp,
                  print = TRUE,
                  plotdir = here::here(plotpath),
                  legendlabels = c('base', 
                                   'vers_bias',
                                   'jitter_lambda',
                                   'growth',
                                   'q',
                                   'survey_slx',
                                   'fishery_slx'))
