# script to compare gap products length comps to base model length comps

# setup ----
libs <- c("here", "tidyverse", 'janitor')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

new_year <- as.numeric(format(Sys.Date(), format = "%Y"))

# 2023 length comp file----
lcomp_old <- read_csv(here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', 'data', 'raw', 'AKskate_Lcomps2023.csv')) %>% 
  clean_names() %>% 
  filter(source == "survey") %>% 
  rename(racebase = lprop,
         totl_old = tot_l,
         yrtot_old = yrtot_l) %>% 
  select(-source)

# new gap products data----
sLcomp_dat <- read_csv(here::here(new_year, '2025_Sept_models', 'AK_skate_Tier3', 'data', 'raw', 'twl_srvy_lpop.csv'))

lbins <- seq(4, 132, by = 4) #set up length bins
llabels <- seq(0, 132, by = 4)

Lseq <- as.data.frame(seq(min(sLcomp_dat$length), max(sLcomp_dat$length)))
names(Lseq) <- "length"
sLcomp_dat <- sLcomp_dat %>% 
  full_join(as.data.frame(Lseq)) #add in blanks for lengths with no observations
#test <- sLcomp_dat %>% 
#  filter(is.na(YEAR))

sLcomp_dat_short <- sLcomp_dat %>% 
  filter(length < 132) %>% 
  mutate(lengbin = as.numeric(as.character(cut(length, breaks = lbins, right = F, labels = lbins[1:32]))))

sLcomp_dat <- sLcomp_dat %>% 
  filter(length >= 132) %>% 
  mutate(lengbin = 132) %>% 
  bind_rows(sLcomp_dat_short)

yrtot <- sLcomp_dat %>% 
  group_by(year) %>% 
  summarise(totL = sum(num))

sLcomps <- sLcomp_dat %>% 
  group_by(year, lengbin) %>% 
  summarise(yrtotL = sum(num)) %>% 
  left_join(yrtot) %>% 
  mutate(gap = yrtotL/totL) %>% 
  clean_names()

#compare RACE to GAP PRODUCTS----
comp_comp <- sLcomps %>% 
  left_join(lcomp_old) %>% 
  mutate(pdiff = (gap - racebase)/racebase)

dat2 <- comp_comp %>% 
  select(year, lengbin, gap, racebase) %>% 
  pivot_longer(!c(year, lengbin), names_to = 'version', values_to = 'prop')


ggplot(dat2, aes(x = lengbin, y = prop, fill = version))+
  geom_bar(stat = 'identity', position = position_dodge())+
  theme_bw()
         