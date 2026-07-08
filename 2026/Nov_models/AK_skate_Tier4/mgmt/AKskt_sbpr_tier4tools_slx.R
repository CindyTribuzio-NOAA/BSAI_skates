# Tier 4 calculations developed by K Omori and J Sullivan
# adapted for AK skates April 2026 by C Tribuzio

# Setup----
#devtools::install_github("noaa-afsc/tier4tools", build_vignettes = TRUE)

libs <- c("dplyr", "tidyr", 'tidyverse', "ggplot2", 'tier4tools', 'r4ss')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

#ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
#                     cowplot::background_grid() +
#                     cowplot::panel_border())

# Setting up inputs----
AYR <- 2026

#proportion of catch by gear type----
#mean proportion over the most recent 5 years
#from 2_catch_gear_target_summary BSAI_m5gear
#update each assessment cycle, things change
pHAL <- 0.834 
pTRW <- 0.166

#age and growth parameters----
ages <- 0:25 # max age - 1 from AGP data, allows for matching to selectivity vector length

# ignoring the built in function and making curve based on Matta MS Thesis
ages <- 0:25
mata_fem <- -15.280
matb_fem <- 1.574
maa <- 1/(1+exp(-(mata_fem+matb_fem*ages)))

# selectivity----
# using Model 25_2 fixed for this analysis
# Load your SS3 model output
m25_2_out <- SS_output(dir = here::here(2025, 'Nov_models', 'AK_skate_Tier3', 'mgmt', 'M25_2'))

# 1. Inspect your fleet names and IDs to find your exact gears
print(m25_2_out$FleetNames) 
# Note down the exact index or string name (e.g., Fleet 1 = "Trawl", Fleet 2 = "Longline")

# 2. Get the maximum model age used (to match vector length)
m25_2_max_age <- m25_2_out$accuage  # or use max(ss_out$ageselex$Age)

# 3. Pull the calculated "Asel" (Age-selectivity) rows for the terminal year
#    Note: Factor == "Asel" represents the combined age-based selectivity
m25_2_terminal_year <- m25_2_out$endyr

# EXTRACTION FOR TRAWL
m25_2_twl_selex <- m25_2_out$ageselex %>%
  filter(Factor == "Asel2", 
         Yr == m25_2_terminal_year, 
         Fleet == 2,          # Swap with your actual Trawl Fleet ID/Name
         Sex == 1) %>%        # 1 = Females/unisex, 2 = Males
  select(as.character(0:m25_2_max_age)) %>%
  as.numeric()

# EXTRACTION FOR LONGLINE
m25_2_lgl_selex <- m25_2_out$ageselex %>%
  filter(Factor == "Asel2", 
         Yr == m25_2_terminal_year, 
         Fleet == 1,          # Swap with your actual Longline Fleet ID/Name
         Sex == 1) %>%
  select(as.character(0:m25_2_max_age)) %>%
  as.numeric()

#SPR inputs----
inp_25_2 <- spr_input(
  ages = ages,
  species = list(
    AK_skate = list(
      len_at_age = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), #Matta and Gunderson 2007
      wt_at_age  = list(type = "wl", alpha = 0.000009, beta = 2.9617), #Matta and Gunderson 2007
      maturity  = list(type = "vector", m = maa),
      M = 0.13, #current assessment value
      selectivity = list(
        type = "fleets",
        fleets = list(
          twl = list(
            propF = pTRW,
            selex = list(type = "vector", s = m25_2_twl_selex)
          ),
          ll = list(
            propF = pHAL,
            selex = list(type = "vector", s= m25_2_lgl_selex)
          )
        )
      )
    )
  ),
  use_plus_group = FALSE
)

# take a look at the inputs
str(inp_14_2d1, max.level = 3)
str(inp_25_2, max.level = 3)
str(inp_logistic, max.level = 3)

# checking input assumptions
inp_14_2d1_plots <- plot_spr_inputs(inp_14_2d1, panels = c("mat_selex", "wt", "len", "M"))
grid_spr_inputs(inp_14_2d1_plots, order = c("mat_selex", "wt","len","M"), ncol = 2, guides = "keep")

inp_25_2_plots <- plot_spr_inputs(inp_25_2, panels = c("mat_selex", "wt", "len", "M"))
grid_spr_inputs(inp_25_2_plots, order = c("mat_selex", "wt","len","M"), ncol = 2, guides = "keep")

inp_logistic_plots <- plot_spr_inputs(inp_logistic, panels = c("mat_selex", "wt", "len", "M"))
grid_spr_inputs(inp_logistic_plots, order = c("mat_selex", "wt","len","M"), ncol = 2, guides = "keep")

# build selex compare plots with composite----
#extract the composite combined selectivity curves
slxdat_14_2d1 <- as_tibble(inp_14_2d1$species[['AK_skate']]$selex_at_age) %>%
  mutate(select_fun = 'M14_2d1',
         type = 'composite') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)
slxdat_25_2 <- as_tibble(inp_25_2$species[['AK_skate']]$selex_at_age) %>% 
  mutate(select_fun = 'M25_2',
         type = 'composite') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)
slxdat_logistic <- as_tibble(inp_logistic$species[['AK_skate']]$selex_at_age) %>% 
  mutate(select_fun = 'logistic',
         type = 'composite') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)

#build the separate logistic selectivity curves
#TWL
slx_log_twl <- as_tibble(1 / (1 + exp((ages - a50_twl) / del_twl))) %>% 
  mutate(select_fun = 'logistic',
         type = 'TWL') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)
#LGL
slx_log_lgl <- as_tibble(1 / (1 + exp((ages - a50_lgl) / del_lgl))) %>% 
  mutate(select_fun = 'logistic',
         type = 'LGL') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)

#format the T3 selectivity vectors
slx_25_2_lgl <- as_tibble(m25_2_lgl_selex) %>% 
  mutate(select_fun = 'M25_2',
         type = 'LGL') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)
slx_25_2_twl <- as_tibble(m25_2_twl_selex) %>% 
  mutate(select_fun = 'M25_2',
         type = 'TWL') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)
slx_14_2d1_lgl <- as_tibble(m14_2d1_lgl_selex) %>% 
  mutate(select_fun = 'M14_2d1',
         type = 'LGL') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)
slx_14_2d1_twl <- as_tibble(m14_2d1_twl_selex) %>% 
  mutate(select_fun = 'M14_2d1',
         type = 'TWL') %>% 
  bind_cols(ages) %>% 
  rename(selex = value,
         age = ...4)


# make graph data
slx_fig_dat <- slxdat_14_2d1 %>% 
  bind_rows(slxdat_25_2,
            slxdat_logistic,
            slx_log_lgl,
            slx_log_twl,
            slx_25_2_lgl,
            slx_25_2_twl,
            slx_14_2d1_lgl,
            slx_14_2d1_twl)
  

slx_fig_14_2d1 <- ggplot(slx_fig_dat, aes(x = age, y = selex, color = type))+
  geom_line()+
  facet_grid(select_fun~.)+
  theme_tier4()

composite_dat <- slx_fig_dat %>% 
  filter(type == 'composite')

comp_comp <- ggplot(composite_dat, aes(x = age, y = selex, color = select_fun))+
  geom_line()+
  theme_tier4()


# SPR and YPR----
spr_14_2d1_out <- run_spr(
  inp_14_2d1,
  diagnostics = TRUE
)
spr_14_2d1_out$F_spr_total
plot_spr_curves(spr_14_2d1_out)

spr_25_2_out <- run_spr(
  inp_25_2,
  diagnostics = TRUE
)
spr_25_2_out$F_spr_total
plot_spr_curves(spr_25_2_out)

spr_logistic_out <- run_spr(
  inp_logistic,
  diagnostics = TRUE
)
spr_logistic_out$F_spr_total
plot_spr_curves(spr_logistic_out)

# Diagnostics----
plots_14_2d1 <- plot_spr_decomp(spr_14_2d1_out,
                         which_panels = c("contrib", "removed", "survivorship"),
                         drop_plus_from_plot = TRUE)
grid_spr_decomp(plots_14_2d1)

plots_25_2 <- plot_spr_decomp(spr_25_2_out,
                                which_panels = c("contrib", "removed", "survivorship"),
                                drop_plus_from_plot = TRUE)
grid_spr_decomp(plots_25_2)

plots_logistic <- plot_spr_decomp(spr_logistic_out,
                                which_panels = c("contrib", "removed", "survivorship"),
                                drop_plus_from_plot = TRUE)
grid_spr_decomp(plots_logistic)


# Tier 4 harvest recommendations----
AKsk_rema <- read_csv(here::here(AYR, 'Nov_models', 'AK_skate_Tier5', 'mgmt', 'Tier5_AKskt_alloutput.csv')) %>% 
  filter(year == AYR,
         strata == 'EBS') %>% 
  rename(rbiom = pred) %>% 
  select(year, rbiom, strata) %>% 
  mutate(rbiom = round(rbiom, 0))

T4_14_2d1 <- AKsk_rema %>% 
  mutate(OFL = round(rbiom * spr_14_2d1_out$F_spr_total[2], 0),
         ABC = round(rbiom * spr_14_2d1_out$F_spr_total[1], 0),
         FOFL = spr_14_2d1_out$F_spr_total[2],
         FABC = spr_14_2d1_out$F_spr_total[1],
         selectivity = "Model14_2_d1",
         q = 1,
         Tier = 4) %>% 
  select(!strata)

T4_25_2 <- AKsk_rema %>% 
  mutate(OFL = round(rbiom * spr_25_2_out$F_spr_total[2], 0),
         ABC = round(rbiom * spr_25_2_out$F_spr_total[1], 0),
         FOFL = spr_25_2_out$F_spr_total[2],
         FABC = spr_25_2_out$F_spr_total[1],
         selectivity = "Model25_2",
         q = 1,
         Tier = 4) %>% 
  select(!strata)

T4_logistic <- AKsk_rema %>% 
  mutate(OFL = round(rbiom * spr_logistic_out$F_spr_total[2], 0),
         ABC = round(rbiom * spr_logistic_out$F_spr_total[1], 0),
         FOFL = spr_logistic_out$F_spr_total[2],
         FABC = spr_logistic_out$F_spr_total[1],
         selectivity = "Logistic",
         q = 1,
         Tier = 4) %>% 
  select(!strata)

T3_base_output <- tibble::tibble(
  year = 2026,
  rbiom = 372517, #find this #
  OFL = 29508,
  ABC = 24242,
  FOFL = 0.092,
  FABC = 0.079,
  selectivity = "Model14_2_d1",
  q = 1,
  Tier = 3
)

T3_m25_2_output <- tibble::tibble(
  year = 2026,
  rbiom = 464159,
  OFL = 30123,
  ABC = 26026,
  FOFL = 0.113,
  FABC = 0.097,
  selectivity = "Model25_2",
  q = 0.836,
  Tier = 3
)

T3_m25_3_output <- tibble::tibble(
  year = 2026,
  rbiom = 390230,
  OFL = 25692,
  ABC = 22188,
  FOFL = 0.114,
  FABC = 0.098,
  selectivity = "Model25_2",
  q = 1,
  Tier = 3
)

T5_output <- tibble::tibble(
  year = 2026,
  rbiom = 385402,
  OFL = round(0.13*385402, 0),
  ABC = round(0.0975*385402, 0),
  FOFL = 0.13,
  FABC = 0.0975,
  selectivity = NA,
  q = 1,
  Tier = 5
)


#add in Tier 5 and Tier 3 M25_3
T4_selex_compare <- T4_14_2d1 %>% 
  bind_rows(T4_25_2) %>% 
  bind_rows(T4_logistic) %>% 
  bind_rows(T3_base_output) %>% 
  bind_rows(T3_m25_2_output) %>% 
  bind_rows(T3_m25_3_output) %>% 
  bind_rows(T5_output)
write_csv(T4_selex_compare, here::here(AYR, 'Nov_models', 'AK_skate_Tier4', 'mgmt', 'Tier4_AKskt_slxcompare.csv'))


