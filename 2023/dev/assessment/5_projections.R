# Script to run the AK skate projection scenarios
# 2025 by C Tribuzio

# Setup ----
# load required packages and code
libs <- c("data.table",
          'tidyverse',
          "magrittr", 
          "r4ss",
          'foreach')

if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}

lapply(libs, library, character.only = TRUE)

# Current assessment year
new_year <- 2023

# load functions
source(here::here(new_year, "R", "mdl_anlys_fcns.r"))

# Model names and directory
Model_name_14 <- "M14_2d1"
M14_2d_path <- here::here('2023/rsch/M14_2d1')

# Run management scenarios----
M12_2d_projrun <- run_ss3_mgmnt_scen(dir = M14_2d_path,
                   cyr = new_year,
                   syr = 1950,
                   sexes = 1,
                   fleets = c(1:3), 
                   Scenario2 = 2, 
                   s2_F = 0.67, #if Scenario2 = 3 this should not be used
                   s4_F = 0.5, #F = 50% of maxFABC, as per the text in the last assessment. NOTE: las proj dat file says 0.6
                   output_name = 'M14_2d_proj')

# plot scenarios----
# projection figure ----
setwd("C:/Users/cindy.Tribuzio/Work/SAFE/Assessments/AFSC_BSAI_SKATE_Assessment")
#note, saved bigfile.out as .csv in excel, need to figure out how to streamline that
projout <- read_csv(paste0(getwd(), "/", AYR, "/Tier3/Projections/prg_AKSK_14_2_2023/bigfile.csv")) %>% 
  select(Alternative, Yr, SSB) %>% 
  group_by(Alternative, Yr) %>% 
  summarise(SSB = mean(SSB))

B40 <- 71370
B35 <- 62449

plot_proj <- ggplot(projout, aes(x = Yr, y = SSB/1000, color = as.factor(Alternative), shape = as.factor(Alternative)))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = B40/1000, linetype = "dashed")+
  geom_hline(yintercept = B35/1000)+
  labs(y = "Spawning Biomass (1,000s t)", x = "Year", color = "")+
  scale_y_continuous(labels = scales::comma_format())+
  scale_color_viridis(discrete = T) +
  #coord_cartesian(ylim = c(75, 350), xlim = c(1950, 2025))+
  theme_bw()

ggsave(path = paste0(getwd(), "/", AYR, "/Tier3/Output"),
       "Projections.png",plot=plot_proj,dpi=600,width = 6, height = 4)


