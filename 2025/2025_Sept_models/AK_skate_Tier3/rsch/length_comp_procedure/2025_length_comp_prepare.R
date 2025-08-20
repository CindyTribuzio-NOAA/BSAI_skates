# Preparation of AK skate length comp data for SS  ----
# Updated Aug 18 2025 by C. Tribuzio
# initially built off AKSK_survey_sizecomp_thru2019.xlsx authored O. Ormseth
# updated to use GAP products 

# Setup ----
libs <- c("tidyverse", "RODBC", 'stats', 'janitor', 'readxl', 'nlstools', 'lubridate')
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

#query gap_products----
#don't run if data have already been queried
AKskt_dat <- sqlQuery(channel_akfin, query = ("
                select * from gap_products.akfin_specimen a
                left join gap_products.akfin_haul b
                on a.hauljoin = b.hauljoin
                where a.species_code = 471")) %>% #not sure about including zeros
  clean_names() %>% 
  mutate(haul_year = year(date_time_start),
         yrmday = ymd(date_time_start),
         length_cm = length_mm/10)








# Setup ----
datapath<-paste(getwd(),"/Data/",AYR,sep="")

lbins <- seq(4, 132, by = 4) #set up length bins
llabels <- seq(0, 132, by = 4)

# EBS Shelf Survey ----
sLcomp_dat <- read_csv(paste0(datapath, "/RACE_Shelflcomp_AKskt", AYR, ".csv")) %>% 
  filter(SUBAREA == 999) %>% 
  mutate(leng_cm = LENGTH/10)

Lseq <- as.data.frame(seq(min(sLcomp_dat$leng_cm), max(sLcomp_dat$leng_cm)))
names(Lseq) <- "leng_cm"
sLcomp_dat <- sLcomp_dat %>% 
  full_join(as.data.frame(Lseq)) #add in blanks for lengths with no observations
#test <- sLcomp_dat %>% 
#  filter(is.na(YEAR))

sLcomp_dat_short <- fLcomp_dat %>% 
  filter(leng_cm < 132) %>% 
  mutate(lengbin = as.numeric(as.character(cut(leng_cm, breaks = lbins, right = F, labels = lbins[1:32]))))

sLcomp_dat <- sLcomp_dat %>% 
  filter(leng_cm >= 132) %>% 
  mutate(lengbin = 132) %>% 
  bind_rows(sLcomp_dat_short)

yrtot <- sLcomp_dat %>% 
  group_by(YEAR) %>% 
  summarise(totL = sum(TOTAL))

sLcomps <- sLcomp_dat %>% 
  group_by(YEAR, lengbin) %>% 
  summarise(yrtotL = sum(TOTAL)) %>% 
  left_join(yrtot) %>% 
  mutate(Lprop = yrtotL/totL,
         source = "survey")

# Trawl Fishery length comps----
fLcomp_dat <- read_csv(paste0(datapath, "/fisherylcompall_AKskt", AYR, ".csv")) %>%
  filter(GEAR <= 4,
         YEAR >= 2013, 
         FMP_AREA == "BSAI") %>% 
  select(YEAR, leng_cm = LENGTH, TOTAL = FREQUENCY)
fLseq <- as.data.frame(seq(min(fLcomp_dat$leng_cm), max(fLcomp_dat$leng_cm)))
names(fLseq) <- "leng_cm"
fLcomp_dat <- fLcomp_dat %>% 
  full_join(as.data.frame(fLseq))

fNsamp <- read_csv(paste0(datapath, "/fisherylcompall_AKskt", AYR, ".csv")) %>%
  filter(GEAR <= 4,
         YEAR >= 2013) %>% 
  group_by(YEAR, HAUL_JOIN) %>% 
  summarise(somval = sum(FREQUENCY)) %>% 
  group_by(YEAR) %>% 
  summarise(Nhauls = length(HAUL_JOIN)) %>% 
  mutate(Nsamp = sqrt(Nhauls))

fLcomp_dat_short <- fLcomp_dat %>% 
  filter(leng_cm < 132) %>% 
  mutate(lengbin = as.numeric(as.character(cut(leng_cm, breaks = lbins, right = F, labels = lbins[1:32]))))

fLcomp_dat <- fLcomp_dat %>% 
  filter(leng_cm >= 132) %>% 
  mutate(lengbin = 132) %>% 
  bind_rows(fLcomp_dat_short)
  
# LL Fishery length comps----
llLcomp_dat <- read_csv(paste0(datapath, "/fisherylcompall_AKskt", AYR, ".csv")) %>%
  filter(GEAR <= 8,
         NMFS_AREA < 541,
         YEAR >= 2013,
         PERFORMANCE >= 1) %>% 
  select(YEAR, leng_cm = LENGTH, TOTAL = FREQUENCY)
llLseq <- as.data.frame(seq(min(llLcomp_dat$leng_cm), max(llLcomp_dat$leng_cm)))
names(llLseq) <- "leng_cm"
llLcomp_dat <- llLcomp_dat %>% 
  full_join(as.data.frame(llLseq))

LLNsamp <- read_csv(paste0(datapath, "/fisherylcompall_AKskt", AYR, ".csv")) %>%
  filter(YEAR >= 2013) %>% 
  group_by(YEAR, HAUL_JOIN, FMP_GEAR) %>% 
  summarise(somval = sum(FREQUENCY)) %>% 
  group_by(YEAR) %>% 
  summarise(Nhauls = length(HAUL_JOIN)) %>% 
  mutate(Nsamp = sqrt(Nhauls))

llLcomp_dat_short <- llLcomp_dat %>% 
  filter(leng_cm < 132) %>% 
  mutate(lengbin = as.numeric(as.character(cut(leng_cm, breaks = lbins, right = F, labels = lbins[1:32]))))

llLcomp_dat <- llLcomp_dat %>% 
  filter(leng_cm >= 132) %>% 
  mutate(lengbin = 132) %>% 
  bind_rows(llLcomp_dat_short)

llyrtot <- llLcomp_dat %>% 
  group_by(YEAR) %>% 
  summarise(totL = sum(TOTAL))

llLcomps <- llLcomp_dat %>% 
  group_by(YEAR, lengbin) %>% 
  summarise(yrtotL = sum(TOTAL)) %>% 
  left_join(llyrtot) %>% 
  mutate(Lprop = yrtotL/totL,
         source = "LL")

# Output results ----
outdat <- sLcomps %>% 
  bind_rows(fLcomps, llLcomps)
write_csv(outdat, paste0(datapath, "/AKskate_Lcomps", AYR, ".csv"))
