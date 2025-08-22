## Script to play around with bsai skate assessment


# user-defined function arguments ----

# data file name
dat_name <- "data_aksk14_2_2023.ss"

# model folder name
mdl <- "kitchen_sink_allslx_rebin"


# load necessary packages ----
library(tidytable)

## github packages ----
pkg_git <- c("r4ss")

# if not installed, then install
if(!isTRUE("r4ss" %in% rownames(installed.packages()))) {
  devtools::install_github("r4ss/r4ss", force = TRUE)
}

# load packages
lapply(pkg_git, library, character.only = TRUE)

## load functions ----
source_files <- c(list.files(here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'R'), pattern = "*.R$"))
purrr::map(here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'R', source_files), source)

# read in data file ----
dat <- r4ss::SS_readdat(here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', 'kitchen_sink_allslx', dat_name), verbose = TRUE)

# adjust length bins ----

## change pop'n length structure to match input for data ----
dat$lbin_method <- 1

## rebin length comp ----
data.frame(dat$lencomp %>% 
             tidytable::select(-c(year, month, fleet, sex, part, Nsamp)) %>% 
             # start at 20 cm
             tidytable::mutate(l20 = l4 + l8 + l12 + l16 + l20) %>%
             tidytable::select(-c(l4, l8, l12, l16)) %>% 
             # end at 112 cm
             tidytable::mutate(l112 = l112 + l116 + l120 + l124 + l128 + l132) %>%
             tidytable::select(-c(l116, l120, l124, l128, l132))) -> new_comp
new_bins <- as.numeric(gsub("l", "", colnames(new_comp)))

## redefine length comp data ----
dat$lencomp <- data.frame(dat$lencomp %>% 
                            tidytable::select(c(year, month, fleet, sex, part, Nsamp)) %>% 
                            tidytable::bind_cols(new_comp))

## redefine number of length bins ----
dat$N_lbins <- length(new_bins)

## redefine length bins ----
dat$lbin_vector <- new_bins

# clean up errors from dat file ----
## fishing fleet survey timing (set to -1 rather than 1) ----
dat$fleetinfo$surveytiming[1:2] <- -1 # set to -1 for all fleets

# write out new data file ----
r4ss::SS_writedat_3.30(dat,
                       here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl, dat_name), 
                       overwrite = TRUE)


# clean up errors from control file ----
## read in control file ----
ctl <- r4ss::SS_readctl_3.30(here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', 'kitchen_sink_allslx', 'control.ss'))

## change rec_dist_method to 4, and remove from MG_parms ----
ctl$recr_dist_method <- 4
ctl$MG_parms <- ctl$MG_parms[-which(grepl('RecrDist', rownames(ctl$MG_parms)) == TRUE),]

## change parameter bounds to match new length bins ----
ctl$size_selex_parms$LO[which(ctl$size_selex_parms$LO == 7.6)] <- (new_bins[2] + new_bins[3]) / 2
ctl$size_selex_parms$HI[which(ctl$size_selex_parms$HI == 126.2)] <- (new_bins[length(new_bins)] + new_bins[length(new_bins) - 1]) / 2

## write new ctl file ----
r4ss::SS_writectl_3.30(ctl,
                       here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl, 'control.ss'), 
                       overwrite = TRUE)

# run base model ----
r4ss::get_ss3_exe(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', 'kitchen_sink_allslx'))
r4ss::run(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', 'kitchen_sink_allslx'),
          skipfinished = FALSE,
          show_in_console = TRUE)


# run new model ----
r4ss::get_ss3_exe(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl))
r4ss::run(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl),
          skipfinished = FALSE,
          show_in_console = TRUE)

# plot results ----
model_run <- r4ss::SS_output(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl))
unlink(here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl, 'plots'), recursive = TRUE)
r4ss::SS_plots(model_run)




# testing stuff ----



# try a model that opens up all selex params ----
sel_turnedoff <- which(ctl$size_selex_parms$PHASE < 0)
ctl$size_selex_parms$PHASE[sel_turnedoff] <- ctl$size_selex_parms$PHASE[sel_turnedoff] * -1

# rest initial values
ctl$size_selex_parms$INIT <- (ctl$size_selex_parms$LO + ctl$size_selex_parms$HI) / 2



r4ss::SS_writectl_3.30(ctl,
                       here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl, 'control.ss'), 
                       overwrite = TRUE)

# run new model ----
r4ss::get_ss3_exe(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl))
r4ss::run(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl),
          skipfinished = FALSE,
          show_in_console = TRUE)

# plot results ----
model_run <- r4ss::SS_output(dir = here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl))
unlink(here::here('2025', '2025_Sept_models', 'AK_skate_Tier3', 'rsch', 'Model_14_2d_fixes', mdl, 'plots'), recursive = TRUE)
r4ss::SS_plots(model_run)
