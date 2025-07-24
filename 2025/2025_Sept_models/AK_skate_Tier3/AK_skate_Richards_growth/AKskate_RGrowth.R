# Query length at age data and develop Richards Growth model parameters ----
# Contact: cindy.tribuzio@noaa.gov
# Last Updated: July 2025

# TODO: #need to check size at age vectors too

# Setup ----
libs <- c("tidyverse", "RODBC", 'stats', 'janitor', 'fishgrowth', 'readxl', 'nlstools')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

# db connection----
dbname <- "akfin"
db <- read_csv('database.csv')
database_akfin=db %>% filter(database == dbname) %>% select(database) #need to add filter for AKFIN user/pass only
username_akfin=db %>% filter(database == dbname) %>% select(username)
password_akfin=db %>% filter(database == dbname) %>% select(password)

channel_akfin <- odbcConnect(dbname, uid = username_akfin, pwd = password_akfin, believeNRows=FALSE)

#query gap_products----
AKskt_dat <- sqlQuery(channel_akfin, query = ("
                select * from gap_products.akfin_specimen a
                left join gap_products.akfin_haul b
                on a.hauljoin = b.hauljoin
                where a.species_code = 471
                and a.age >= 0")) %>% #not sure about including zeros
  clean_names() %>% 
  filter(!is.na(final_age) & readability %in% c(1,2)) %>% 
  mutate(haul_year = year(date_time_start),
         yrmday = ymd(date_time_start),
         length_cm = length_mm/10)

#data were compared to a direct pull from the AGP database and everything matches
write_csv(AKskt_dat, paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/AK_skate_Richards_growth/AKskate_gap_ages.csv'))

AKskt_dat <- AKskt_dat %>% 
  select(c('cruisejoin', 'hauljoin', 'haul_year', 'yrmday', 'specimen_id', 'species_code', 'length_cm', 'sex', 'weight_g', 'age', 'maturity'))

nyr <- AKskt_dat %>% 
  group_by(haul_year) %>% 
  summarise(naged = length(age))

# add in B. Matta thesis samples not included in AGP database----
#NOTE: there are about 350 fishery samples available, however, they are from different months as the survey samples, excluded for now
# see folder: AGP_Matta_data for fishery data

surv2003 <- read_excel(paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/AK_skate_Richards_growth/AGP_Matta_data/2003survey_AKskate_specimen data.xlsx')) %>% 
  clean_names() %>% 
  mutate(specimen_id = as.numeric(gsub( "*.-", "", specimen))) %>% 
  filter(final_age != -9)
surv2004 <- read_excel(paste0(getwd(), '/2025/2025_Sept_models/AK_skate_Tier3/AK_skate_Richards_growth/AGP_Matta_data/2004survey_AKskate_specimen data.xlsx')) %>% 
  clean_names() %>% 
  mutate(specimen_id = as.numeric(gsub( "*.-", "", specimen))) %>% 
  filter(final_age != -9)
surv_ages <- surv2003 %>% 
  bind_rows(surv2004) %>% 
  mutate(haul_year = year(date),
         yrmday = ymd(date),
         species_code = 471) %>% 
  rename(cruisejoin = cruise,
         age = final_age,
         length_cm = total_length,
         weight_g = weight) %>% 
  select(c('cruisejoin', 'haul_year', 'yrmday', 'specimen_id', 'species_code', 'length_cm', 'sex', 'weight_g', 'age', 'maturity' ))

akdat2 <- AKskt_dat %>% 
  bind_rows(surv_ages)

# Schnute/Richards growth model----
# matches formulation in SS3, see user manuals: https://github.com/nmfs-ost/ss3-doc/releases

start <- list(L0=15, Linf=102, k=0.2, g=-1)
m1 <- akdat2$length_cm ~ (L0^g + (Linf^g - L0^g)* ((1-exp(-k*(akdat2$age - 0)))/(1-exp(-k*(26-0)))))^(1/g)
fitm1 <- nls(m1, data = akdat2, start = start)
summary(fitm1)
confint(fitm1)
fitpm1 <- as.data.frame(coef(fitm1))

ages <- seq(0, 26)
estLm1 <- (fitpm1[1,1]^fitpm1[4,1] + (fitpm1[2,1]^fitpm1[4,1] - fitpm1[1,1]^fitpm1[4,1])* ((1-exp(-fitpm1[3,1]*(ages - 0)))/(1-exp(-fitpm1[3,1]*(26-0)))))^(1/fitpm1[4,1])
modestm1 <- bind_cols(ages, estLm1)
names(modestm1) <- c('age', 'length_cm')


AKskt_growthfit <- ggplot(akdat2, aes(x=age, y=length_cm))+
  geom_point(alpha=0.5, aes(color=factor(haul_year))) +
  geom_smooth(lwd=1.3, method = "nls",
              se = FALSE,
              method.args = list(formula = y ~ (L0^g + (Linf^g - L0^g)* ((1-exp(-k*(x - 0)))/(1-exp(-k*(26-0)))))^(1/g),
                                 start = start)) +
  geom_line(data = modestm1, aes(x = age, y = length_cm), color = 'red')+ #just tests that the geom_smooth is the same as the converged nls model
  scale_color_viridis_d() +
  labs(x = "Age (yr)", y= "Total length (mm)", color = "Year")+
  theme_bw()

ggsave(path = paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3/AK_skate_Richards_growth"),
       "AKskt_fitgrowth.png",plot = AKskt_growthfit,dpi=600,width = 8, height = 8)


ggplot(akdat2, aes(x= age, y=length_cm))+
  #geom_boxplot(alpha=0.5) +
  geom_jitter(alpha = 0.25)+
  #geom_line(modest, aes(x = age, y = length_cm))+
  geom_smooth(lwd=1.3, method = "nls",
              se = FALSE,
              method.args = list(formula = y ~ (L0^g + (Linf^g - L0^g)* ((1-exp(-k*(x - 0)))/(1-exp(-k*(20-0)))))^(1/g),
                                 start = start))


# model fit by year ----
AKskt_growthyr <- ggplot(akdat2, aes(x=age, y=length_cm, colour = as.factor(haul_year)))+
  geom_point(alpha=0.5) +
  geom_smooth(lwd=1.3, method = "nls",
              se = FALSE,
              method.args = list(formula = y ~ (L0^g + (Linf^g - L0^g)* (1-exp(-k*(x - 0))/(1-exp(-k*(26-0)))))^(1/g),
                                 start = start)) +
  scale_color_viridis_d() +
  labs(x = "Age (yr)", y= "Total length (mm)", color = "Year")+
  theme_bw()
ggsave(path = paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3/AK_skate_Richards_growth"),
       "AKskt_growthyr.png",plot = AKskt_growthyr, dpi=600,width = 8, height = 8)

AKskt_growthyr2 <- ggplot(akdat2, aes(x=age, y=length_cm, colour = as.factor(haul_year)))+
  geom_point(alpha=0.5) +
  geom_smooth(lwd=1.3, method = "nls",
              se = FALSE,
              method.args = list(formula = y ~ (L0^g + (Linf^g - L0^g)* (1-exp(-k*(x - 0))/(1-exp(-k*(26-0)))))^(1/g),
                                 start = start)) +
  scale_color_viridis_d() +
  facet_wrap(.~haul_year, ncol = 2)+
  labs(x = "Age (yr)", y= "Total length (mm)", color = "Year")+
  theme_bw()

ggsave(path = paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier3/AK_skate_Richards_growth"),
       "AKskt_growthyr2.png",plot = AKskt_growthyr2, dpi=600,width = 8, height = 8)

# Amax = 20 instead of full dataset----
akdat3 <- akdat2 %>% 
  filter(age < 26)

m2 <- akdat3$length_cm ~ (L0^g + (Linf^g - L0^g)* ((1-exp(-k*(akdat3$age - 0)))/(1-exp(-k*(20-0)))))^(1/g)
fitm2 <- nls(m2, data = akdat3, start = start)
summary(fitm2)
confint(fitm2)

fitpm2 <- as.data.frame(coef(fitm2))

estLm2 <- (fitpm2[1,1]^fitpm2[4,1] + (fitpm2[2,1]^fitpm2[4,1] - fitpm2[1,1]^fitpm2[4,1])* ((1-exp(-fitpm2[3,1]*(ages - 0)))/(1-exp(-fitpm2[3,1]*(26-0)))))^(1/fitpm2[4,1])
modestm2 <- bind_cols(ages, estLm2)
names(modestm2) <- c('age', 'length_cm')

ggplot(akdat3, aes(x=age, y=length_cm))+
  geom_point(alpha=0.5, aes(color=factor(haul_year))) +
  geom_smooth(lwd=1.3, method = "nls",
              se = FALSE,
              method.args = list(formula = y ~ (L0^g + (Linf^g - L0^g)* ((1-exp(-k*(x - 0)))/(1-exp(-k*(20-0)))))^(1/g),
                                 start = start)) +
  scale_color_viridis_d() +
  labs(x = "Age (yr)", y= "Total length (mm)", color = "Year")+
  theme_bw()

# compare m1 and m2 (Amax = 26 vs Amax = 20) to what's in current assessment
estL14_2 <- (14.9558^-1 + (102.12^-1 - 14.9558^-1)* ((1-exp(-0.3669*(ages - 0)))/(1-exp(-0.3669*(26-0)))))^(1/-1)
modest14_2 <- bind_cols(ages, estL14_2)
names(modest14_2) <- c('age', 'length_cm')

ggplot(akdat2, aes(x= age, y=length_cm))+
  geom_point(alpha = 0.25)+
  geom_line(data = modestm1, aes(x = age, y = length_cm), color = 'red', lwd = 2)+
  geom_line(data = modestm2, aes(x = age, y = length_cm), color = 'green', lwd = 1.5)+
  geom_line(data = modest14_2, aes(x = age, y = length_cm), color = 'blue', lwd = 1.5)+
  labs(x = "Age (yr)", y= "Total length (cm)")+
  theme_bw()

# mean length at age-----
# direct comparison of growth parameters to data going into the model

# matches the data in the model data file
empLAA <- akdat2 %>%
  filter(haul_year != 2004) %>% #2004 not in model, probably due to small sample size
  group_by(age, haul_year) %>% 
  summarise(mLAA = mean(length_cm), nper = length(length_cm))

ggplot(empLAA, aes(x = age, y = mLAA))+
  geom_point()+
  geom_line(data = modestm1, aes(x = age, y = length_cm), color = 'red', lwd = 2)+
  geom_line(data = modestm2, aes(x = age, y = length_cm), color = 'green', lwd = 1.5)+
  geom_line(data = modest14_2, aes(x = age, y = length_cm), color = 'blue', lwd = 1.5)+
  labs(x = "Age (yr)", y= "Total length (cm)")+
  facet_grid(haul_year~.)+
  theme_bw()








