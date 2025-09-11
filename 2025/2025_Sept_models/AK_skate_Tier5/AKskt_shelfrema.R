# Tier 5 Alaska Skate model alternative

# Setup ----
libs <- c("tidyverse", "janitor", "Hmisc", "RColorBrewer", "gridExtra", "gtable", 
          "grid", "flextable", "officer", "lubridate", "RODBC", "DBI", "gtable", 
          "patchwork", 'cowplot', 'rema', 'r4ss', 'viridis')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

# set up ----
# folder set up
dat_path <- paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier5'); dir.create(dat_path)
out_path <- paste0(getwd(), "/", AYR, "/Tier5/Output"); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 13) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

pYEAR <- 2024

# Get Data ----
biomass_dat <- read_csv(paste0(dat_path, '/RACE_biomass_skates', SYR, ".csv"))

#EBS shelf Survey for Alaska skate
EBS_dat <- biomass_dat %>% 
  filter(RACE_name == "Alaska skate",
         year >= 1999,
         survey == "EBS_SHELF") %>% 
  group_by(year) %>% 
  summarise(biomass = sum(biomass),
            var = sum(biomass_var)) %>% 
  ungroup() %>% 
  mutate(cv = sqrt(var)/biomass,
         strata = 'AKskate') %>% 
  select(strata, year, biomass, cv)

# prep data for rema  ----
EBS_input <- prepare_rema_input(model_name = 'AKskate',
                            biomass_dat = EBS_dat,
                            end_year = pYEAR,
                            zeros = list(assumption = "NA"),
                            PE_options = list(pointer_PE_biomass = c(1)))

# fit rema models ----
m_AKskt <- fit_rema(EBS_input)
m_AKskt_out <- tidy_rema(m_AKskt)

# plot and clean up output ----
m_AKskt_plots <- plot_rema(tidy_rema = m_AKskt_out)

m_AKskt_tot <- tidy_rema(m_AKskt)$biomass_by_strata %>% 
  select(strata, model_name, strata, year, variable, pred, pred_lci, pred_uci, obs, obs_cv, obs_lci, obs_uci) 

write_csv(m_AKskt_tot, paste0(dat_path, "/Tier5_AKskt_output.csv"))

# Tier 5 harvest recommendations ----
T5recs <- m_AKskt_tot %>% 
  group_by(year) %>% 
  summarise(rbiom = sum(pred),
            OFL = rbiom * 0.1,
            ABC = OFL * 0.75) %>% 
  filter(year == 2023) %>% 
  mutate(FOFL = 0.13,
         FABC = 0.75,
         model_name = "Tier 5")

# compare to Tier 3 models----
T3_out <- read_csv(paste0(dat_path, '/Tier3_AKskt_survfits.csv')) %>% 
  clean_names() %>% 
  filter(yr >= 1999) %>% 
  select(yr, obs, exp, model) %>% 
  rename(year = yr,
         pred = exp,
         model_name = model)

T5_out <- m_AKskt_tot %>% 
  select(year, obs, pred, pred_lci, pred_uci, obs_lci, obs_uci) %>% 
  mutate(model_name = 'rema')

model_comp <- T3_out %>% 
  bind_rows(T5_out)

ggplot(data = model_comp,
       aes(x = year, y = pred,
           col = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                  fill = model_name), col = NA,
              alpha = 0.25) +
  geom_line() +
  geom_point(aes(x = year, y = obs, col = model_name, shape = model_name)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci, col = model_name)) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = 'Year', y = 'Biomass',
       fill = NULL, colour = NULL, shape = NULL)
  #ggplot2::scale_fill_viridis_d(direction = 1) +
  #ggplot2::scale_colour_viridis_d(direction = 1)
