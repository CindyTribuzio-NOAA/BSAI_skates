# Tier 4 calculations developed by K Omori and J Sullivan
# adapted for AK skates April 2026 by C Tribuzio

# Setup----
#devtools::install_github("noaa-afsc/tier4tools", build_vignettes = TRUE)

libs <- c("dplyr", "tidyr", "ggplot2", 'tier4tools')
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
ages <- 0:26 # max age from AGP data

# ignoring the built in function and making curve based on Matta MS Thesis
ages <- 0:26
mata_fem <- -15.280
matb_fem <- 1.574
maa <- 1/(1+exp(-(mata_fem+matb_fem*ages)))

#SPR inputs----
inp <- spr_input(
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
            selex = list(type = "logistic", a50 = 3, delta = -0.5)#need to make these vectors of selex @ age
          ),
          ll = list(
            propF = pHAL,
            selex = list(type = "logistic", a50 = 7, delta = -2)
          )
        )
      )
    )
  ),
  use_plus_group = FALSE
)

# take a look at the inputs
str(inp, max.level = 3)

# checking input assumptions
inp_plots <- plot_spr_inputs(inp, panels = c("mat_selex", "wt", "len", "M"))

grid_spr_inputs(inp_plots, order = c("mat_selex", "wt","len","M"), ncol = 2, guides = "keep")


# SPR and YPR----
spr_out <- run_spr(
  inp,
  diagnostics = TRUE
)

spr_out$F_spr_total

plot_spr_curves(spr_out)

# Diagnostics----
plots <- plot_spr_decomp(spr_out,
                         which_panels = c("contrib", "removed", "survivorship"),
                         drop_plus_from_plot = TRUE)

## individual plots
# plots$contrib
# plots$removed
# plots$survivorship

grid_spr_decomp(plots)


# Tier 4 harvest recommendations----
AKsk_rema <- read_csv(here::here(AYR, 'Nov_models', 'AK_skate_Tier5', 'mgmt', 'Tier5_AKskt_alloutput.csv')) %>% 
  filter(year == AYR) %>% 
  rename(rbiom = pred) %>% 
  select(year, rbiom, strata) %>% 
  mutate(rbiom = round(rbiom, 0))

T4regionrecs <- AKsk_rema %>% 
  mutate(OFL = round(rbiom * res$F35, 0),
         ABC = round(rbiom * res$F40, 0),
         FOFL = res$F35,
         FABC = res$F40,
         model_name = c("Tier 4 EBS", "Tier 4 NBS")) %>% 
  select(!strata)

T4combrecs <- read_csv(here::here(AYR, 'Nov_models', 'AK_skate_Tier5', 'mgmt', 'Tier5_AKskt_alloutput.csv')) %>% 
  filter(year == AYR) %>% 
  group_by(year) %>% 
  summarise(rbiom = round(sum(pred), 0),
            OFL = round(rbiom * res$F35, 0),
            ABC = round(rbiom * res$F40, 0)) %>% 
  filter(year == AYR) %>% 
  mutate(FOFL = res$F35,
         FABC = res$F40,
         model_name = "Tier 4 combined")

# all T4 options----
T4recs_combined <- T4regionrecs %>% 
  bind_rows(T4combrecs) 
write_csv(T4recs_combined, here::here(AYR, 'Nov_models', 'AK_skate_Tier4', 'mgmt', 'Tier4_AKskt_options.csv'))


