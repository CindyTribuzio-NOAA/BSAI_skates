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

# change widen selectivity bounds----
#changing bounds directly in the control file for now
bnd_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run3_M14_2d_widebounds_lselex')

r4ss::run(dir = bnd_mod_path, skipfinished = FALSE, exe = exe_loc)

SSplotSelex(mbnd)

#look at selectivity parameters output
mbnd <- SS_output(bnd_mod_path, printstats = FALSE, verbose = FALSE)
mbnd$parameters

mbnd_selex <- mbnd$sizeselex %>% 
  filter(Yr == 2023,
         Factor == 'Lsel') %>% 
  mutate(model = 'widebounds')

mbnd_selexpars <- mbnd$SelSizeAdj %>% 
  filter(Yr == 2023) %>% 
  mutate(model = 'widebounds')

selex_results <- mnew_selex %>% 
  bind_rows(mbnd_selex, mold_selex) %>% 
  select(-c(Factor, Yr, Sex, Label)) %>% 
  pivot_longer(!c(Fleet, model), names_to = 'Lbin', values_to = 'Selectivity')

ggplot(selex_results, aes(x = as.numeric(Lbin), y = Selectivity, color = model, shape = model))+
  geom_point()+
  geom_line()+
  facet_grid(Fleet~.)+
  theme_bw()

#overall comparison of model runs----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch","run1_M14_2d_ss3version", "run3_M14_2d_widebounds_lselex"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")

model_comp <- SSsummarize(bridge_out)

SSplotComparisons(model_comp,
                        print = TRUE,
                        plotdir = here::here(datapath),
                        legendlabels = c('base', 'updated SS3', 'bounds'))


























# compare old vs. new SS3 versions applied to 2017 inputs
oldv_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/base_M14_2d_fixedcatch')

mold <- SS_output(oldv_path, printstats = FALSE, verbose = FALSE)

mbnd <- SS_output(bnd_mod_path, printstats = FALSE, verbose = FALSE)
SStableComparisons(
  SSsummarize(list(mold, mnew, mbnd)),
  names = c("NatM", "SSB_Virgin", "SSB_2017", "Bratio_2017"),
  likenames = NULL,
  modelnames = c("3.30.21.00", "3.30.23.2", "Bounds"),
) |>
  dplyr::mutate(across(-1, ~ round(.x, 3)))

SSplotComparisons(list(mold, mnew, mbnd), subplots = "selectivity")

summarymodels <- SSsummarize(list(mold, mnew, mbnd))
SSplotComparisons(summarymodels, subplots = "selectivity")

