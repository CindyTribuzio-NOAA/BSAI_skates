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

# logistic selectivity for comparison to RTMB----
# same as above, changing selectivity directly in control fil
log_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run2_M14_2d_logistic_lselex')

r4ss::run(dir = log_mod_path, skipfinished = FALSE, exe = exe_loc)

mlog <- SS_output(log_mod_path, printstats = FALSE, verbose = FALSE)

SSplotSelex(mlog)

mlog_selexpars <- mlog$SelSizeAdj %>% 
  filter(Yr == 2023) %>% 
  mutate(model = 'logistic')

# fixing growth outside of the model----
grow_mod_path <- here::here('2025/2025_Sept_models/AK_skate_Tier3/run4_M14_2d_fixedgrowth')

waa <- r4ss::SS_readwtatage(file = paste0(grow_mod_path, '/wtatage.ss_new'))
out <- r4ss::SS_output(grow_mod_path, verbose = TRUE)

# Extract Data ------------------------------------------------------------

# define dimensions
n_regions <- 1 # number of regions
n_sexes <- 1 # number of sexes
n_fish_fleets <- 2 # number of fishery fleets
n_srv_fleets <- 1 # number of survey fleets
years <- 1950:2023 # years
n_yrs <- length(years) # number of years
ages <- 0:25 # vector of ages
n_ages <- length(ages) # number of ages
lens <- out$lbins # vector of lengths
n_lens <- length(lens) # number of length bins
spwn_month <- 1 # spawning month (start of year)

### Demographics ------------------------------------------------------------

# weight at age
waa <- out$ageselex %>% filter(Label == '1952_1_bodywt')
waa_arr <- array(0, dim = c(n_regions, n_yrs, n_ages, n_sexes))
waa_arr[] <- rep(out$endgrowth$Wt_Mid, each = n_yrs)

# get length at age and sd to construct ALK
laa <- out$endgrowth$Len_Mid # length at age
sd <- out$endgrowth$SD_Mid # sd of length at age

# construct alk
get_al_trans_matrix = function(age_bins, len_bins, mean_length, sd) {
  
  # Construct age length matrix
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))
  
  for(age_ndx in 1:length(age_bins)) {
    for(len_ndx in 1:length(len_bins)) {
      if (len_ndx == 1) {
        age_length[age_ndx, len_ndx] = pnorm(len_bins[2], mean_length[age_ndx], sd[age_ndx])
      } else if (len_ndx == length(len_bins)) {
        age_length[age_ndx, len_ndx] = 1 - pnorm(len_bins[length(len_bins)], mean_length[age_ndx], sd[age_ndx])
      } else {
        age_length[age_ndx, len_ndx] = pnorm(len_bins[len_ndx+1], mean_length[age_ndx], sd[age_ndx]) -  
          pnorm(len_bins[len_ndx], mean_length[age_ndx], sd[age_ndx])
      }
    }
  }
  return(age_length)
} # end function

alk <- get_al_trans_matrix(age_bins = ages,
                           len_bins  = lens,
                           mean_length =  as.numeric(unlist(laa)),
                           sd = as.numeric(unlist(sd))
)

# normalize ALK in case it doesnt sum to 1
alk <- t(apply(alk, 1, function(x) x / sum(x)))

# populate size age transition matrix
sizeage <- array(0, dim = c(n_regions, n_yrs, n_lens, n_ages, n_sexes))
for(i in 1:n_yrs) sizeage[1,i,,,1] <- t(alk)










#overall comparison of model runs----
datapath <- paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3")
setwd(datapath)
bridge_out <- SSgetoutput(dirvec = c("base_M14_2d_fixedcatch", 
                                     "run1_M14_2d_ss3version", 
                                     "run2_M14_2d_logistic_lselex", 
                                     "run3_M14_2d_widebounds_lselex"))
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/BSAI_skates")

model_comp <- SSsummarize(bridge_out)

SSplotComparisons(model_comp,
                        print = TRUE,
                        plotdir = here::here(datapath),
                        legendlabels = c('base', 'updated SS3', 'logistic', 'bounds'))


























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

