# explorations of length data to inform on selectivity ----
# Contact: cindy.tribuzio@noaa.gov
# Last Updated: July 2025

# Setup ----
libs <- c("tidyverse", "RODBC", 'stats', 'janitor', 'fishgrowth', 'readxl', 'nlstools')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function
round_any = function(x, accuracy, f=floor){f(x/ accuracy) * accuracy}

# db connection----
dbname <- "akfin"
db <- read_csv('database.csv')
database_akfin=db %>% filter(database == dbname) %>% select(database) #need to add filter for AKFIN user/pass only
username_akfin=db %>% filter(database == dbname) %>% select(username)
password_akfin=db %>% filter(database == dbname) %>% select(password)

channel_akfin <- odbcConnect(dbname, uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

# survey and fishery length distribution, overall----
#survey length data

gapLdat <- sqlQuery(channel_akfin, query = ("
                select b.cruisejoin, a.hauljoin, c.year, c.survey_name, a.species_code, a.sex, a.frequency, a.length_mm, b.date_time_start, b.depth_m, b.gear_temperature_c
                from gap_products.akfin_length a
                left join gap_products.akfin_haul b on a.hauljoin = b.hauljoin
                left join gap_products.akfin_cruise c on b.cruisejoin = c.cruisejoin
                where a.species_code = 471")) %>% 
  clean_names() %>% 
  mutate(length_cm = length_mm/10) 

gapLsumm <- gapLdat %>% 
  mutate(lbin = round_any(length_cm, 10)) %>% 
  group_by(survey_name, lbin) %>% 
  summarise(n_lengths = sum(frequency)) %>% 
  mutate(survey = if_else(survey_name == 'Aleutian Islands Bottom Trawl Survey', 'AI',
                          if_else(survey_name == 'Gulf of Alaska Bottom Trawl Survey', 'GOA',
                                  if_else(survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey', 'EBS Slope',
                                          if_else(survey_name == 'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey', 'EBS Shelf',
                                                  'NBS'))))) %>% 
  ungroup() %>% 
  select(!survey_name)

ggplot(gapLsumm, aes(x = lbin, y = n_lengths))+
  geom_bar(stat = 'identity')+
  facet_grid(survey~.)+
  theme_bw()

lrgskt <- gapLsumm %>% 
  filter(lbin > 110)

ggplot(lrgskt, aes(x = lbin, y = n_lengths))+
  geom_bar(stat = 'identity')+
  facet_grid(survey~.)+
  theme_bw()

#fishery length data
fishLdat <- sqlQuery(channel_akfin, query = ("
                SELECT * from NORPAC.debriefed_length_mv
                WHERE species = 88"))
fishLdat <- fishLdat %>% 
  clean_names() %>% 
  mutate(lbin = round_any(length, 10)) %>% 
  group_by(fmp_subarea, lbin) %>% 
  summarise(n_lengths = sum(frequency))

ggplot(fishLdat, aes(x = lbin, y = n_lengths))+
  geom_bar(stat = 'identity')+
  facet_grid(fmp_subarea~.)+
  theme_bw()

lrgskt_fsh <- fishLdat %>% 
  filter(lbin > 110)

ggplot(lrgskt_fsh, aes(x = lbin, y = n_lengths))+
  geom_bar(stat = 'identity')+
  facet_grid(fmp_subarea~.)+
  theme_bw()

# Lengths by depth----
gapDsum <- gapLdat %>% 
  mutate(lbin = round_any(length_cm, 10),
         dbin = round_any(depth_m, 50))
ggplot(gapDsum, aes(x = lbin, y = depth_m, color = gear_temperature_c))+
  geom_jitter()+
  scale_color_gradient2(low = "blue", high = "red")

totlengths <- sum(gapDsum$frequency)
gapD2 <- gapDsum %>% 
  group_by(dbin) %>% 
  summarise(n_lengths = sum(frequency),
    pdepth = round((n_lengths/totlengths)*100,2))

mdepthyr <- gapDsum %>% 
  mutate(survey = if_else(survey_name == 'Aleutian Islands Bottom Trawl Survey', 'AI',
                          if_else(survey_name == 'Gulf of Alaska Bottom Trawl Survey', 'GOA',
                                  if_else(survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey', 'EBS Slope',
                                          if_else(survey_name == 'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey', 'EBS Shelf',
                                                  'NBS'))))) %>% 
  filter(survey == 'EBS Shelf') %>% 
  group_by(year, lbin) %>% 
  summarise(mdepth = mean(depth_m, na.rm = T),
            mtemp = mean(gear_temperature_c, na.rm = T))

ggplot(mdepthyr, aes(x = lbin, y = mdepth, color = mtemp))+
  geom_point()+
  scale_color_gradient2(low = "blue", high = "red")+
  facet_wrap(.~year, ncol = 3)

# length by temp----
gapD3 <- gapDsum %>% 
  mutate(survey = if_else(survey_name == 'Aleutian Islands Bottom Trawl Survey', 'AI',
                          if_else(survey_name == 'Gulf of Alaska Bottom Trawl Survey', 'GOA',
                                  if_else(survey_name == 'Eastern Bering Sea Slope Bottom Trawl Survey', 'EBS Slope',
                                          if_else(survey_name == 'Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey', 'EBS Shelf',
                                                  'NBS'))))) %>% 
  group_by(lbin, dbin, survey, year) %>% 
  summarise(mtemp = mean(gear_temperature_c, na.rm = T), n_lengths = sum(frequency)) %>% 
  filter(survey == "EBS Shelf")

ggplot(gapD3, aes(x = dbin))
