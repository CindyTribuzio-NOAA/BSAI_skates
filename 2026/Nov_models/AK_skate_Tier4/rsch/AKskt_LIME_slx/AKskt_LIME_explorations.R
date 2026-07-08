# Using LIME to generate selectivity curves from length comps
# Cindy Tribuzio
# June 26 2026

# 1. Install devtools if you don't have it
#install.packages("devtools")

# 2. Install the core LIME modeling package
#devtools::install_github("merrillrudd/LIME")

# 3. Launch the dedicated Shiny testing interface
#shiny::runGitHub("merrillrudd/LIME_shiny")

# setup----
libs <- c("tidyverse", "janitor")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

AYR <- 2026
raw_dat_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier3', 'data')
mod_path <- here::here(AYR, 'Nov_models', 'AK_skate_Tier4', 'rsch', 'AKskt_LIME_slx')

lbins <- seq(4, 132, by = 4) #set up length bins
llabels <- seq(0, 132, by = 4)

# read in and format data for LIME----
# Catch----
catch_dat <- read_csv(here::here(raw_dat_path, paste0('confidential_CAS_SKpart_', AYR, '.csv'))) %>% 
  clean_names() %>% 
  group_by(year, agency_gear_code) %>% 
  summarise(catch = sum(catch_weight)) %>% 
  mutate(fleet = if_else(agency_gear_code == 'HAL', 1, 2)) %>% #same order as in SS3 files
  select(!agency_gear_code) %>% 
  select(fleet, year, catch)
write_csv(catch_dat, here::here(mod_path, paste0('catch_multifleet_', AYR, '.csv')))

# HAL lengths----
#keeping all of the steps in here, even though for the first pass only using raw numbers
llLcomp_dat <- read_csv(here::here(raw_dat_path, paste0('confidential_AKskate_fisherysizecomp_', AYR, ".csv"))) %>%
  filter(GEAR <= 8,
         NMFS_AREA < 541,
         YEAR >= 2013,
         PERFORMANCE >= 1) %>% 
  select(YEAR, leng_cm = LENGTH, TOTAL = FREQUENCY)
llLseq <- as.data.frame(seq(min(llLcomp_dat$leng_cm), max(llLcomp_dat$leng_cm)))
names(llLseq) <- "leng_cm"
llLcomp_dat <- llLcomp_dat %>% 
  full_join(as.data.frame(llLseq))

LLNsamp <- read_csv(here::here(raw_dat_path, paste0('confidential_AKskate_fisherysizecomp_', AYR, ".csv"))) %>%
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

# TWL lengths----
#note: here TWL is called "f" for "fishery", it's a hold over from other code that I didn't want to change all of them
fLcomp_dat <- read_csv(here::here(raw_dat_path, paste0('confidential_AKskate_fisherysizecomp_', AYR, ".csv"))) %>%
  filter(GEAR <= 4,
         YEAR >= 2013, 
         FMP_AREA == "BSAI") %>% 
  select(YEAR, leng_cm = LENGTH, TOTAL = FREQUENCY)
fLseq <- as.data.frame(seq(min(fLcomp_dat$leng_cm), max(fLcomp_dat$leng_cm)))
names(fLseq) <- "leng_cm"
fLcomp_dat <- fLcomp_dat %>% 
  full_join(as.data.frame(fLseq))

fNsamp <- read_csv(here::here(raw_dat_path, paste0('confidential_AKskate_fisherysizecomp_', AYR, ".csv"))) %>%
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

fyrtot <- fLcomp_dat %>% 
  group_by(YEAR) %>% 
  summarise(totL = sum(TOTAL))

fLcomps <- fLcomp_dat %>% 
  group_by(YEAR, lengbin) %>% 
  summarise(yrtotL = sum(TOTAL)) %>% 
  left_join(fyrtot) %>% 
  mutate(Lprop = yrtotL/totL,
         source = "TWL")

#combine L freqs----
lfreq_dat <- llLcomps %>% 
  bind_rows(fLcomps) %>% 
  mutate(fleet = if_else(source == 'LL', 1, 2)) %>% 
  select(fleet, YEAR, lengbin, yrtotL) %>% 
  clean_names() %>% 
  pivot_wider(names_from = lengbin, values_from = yrtot_l, names_sort = T) %>% 
  filter(year != 2026) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)))

write_csv(lfreq_dat, here::here(mod_path, paste0('length_multifleet_', AYR, '.csv')))
  

##---------------------------------------------------------------
## *** LIME, Multiple years, single-fleet
##----------------------------------------------------------------

##----------------------------------------------------------------
## Step 1: Weight length data
##----------------------------------------------------------------
catch <- read_csv(here::here(mod_path, paste0("catch_multifleet_", AYR, ".csv")))

data_mymf <- read_csv(here::here(mod_path, paste0("length_multifleet_", AYR, ".csv")))

years <- unique(data_mymf$year)

fleets <- unique(data_mymf$fleet)

cldata <- full_join(data_mymf, catch)

## length comps
lcomps <- cldata %>% select(grep("X", colnames(cldata)))

## observed samples by bin by year
nsamps <- sum(cldata %>% select(grep("X", colnames(cldata))))

## multiply observed number of fish in bins by catch, keep fleet information
lweight1 <- ((cldata %>% select(grep("X", colnames(cldata)))) * cldata$catch) %>% mutate(fleet=cldata$fleet)

## add weighted lengths from each fleet
lweight2 <- ((lweight1 %>% filter(Fleet==1)) + (lweight1 %>% filter(Fleet==2))) %>% select(-Fleet)

## calculate proportions weighted by catch
lweight3 <- matrix(t(sapply(1:length(years), function(x) as.numeric(lweight2[x,]/sum(lweight2[x,])))), nrow=length(years), ncol=length(bins))

LF_new <- lweight3 * nsamps
rownames(LF_new) <- years
colnames(LF_new) <- bins

##----------------------------------------------------------------
## Step 3: Run LIME
##----------------------------------------------------------------

data_list <- list("years"=years, "LF"=LF_new)

## use lh - single fleet
input_data <- create_inputs(lh=lh, input_data=data_list)

res <- run_LIME(modpath=NULL, input=input_data, data_avail="LC")

## check TMB inputs
Inputs <- res$Inputs

## Report file
Report <- res$Report

## Standard error report
Sdreport <- res$Sdreport

## check convergence
hessian <- Sdreport$pdHess
gradient <- res$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE

check <- Check_Identifiable2(res$obj)


##----------------------------------------------------------------
## Step 4: Examine results
##----------------------------------------------------------------

spr <- Report$SPR_t[length(Report$SPR_t)]

## standard error
sd_spr <- summary(Sdreport)[which(rownames(summary(Sdreport))=="SPR_t"),2][length(Report$SPR_t)]

## lower confidence limit
lcl_spr <- max(0,spr - 1.96 * sd_spr)

## upper confidence limit
ucl_spr <- min(1,spr + 1.96 * sd_spr)


save <- rbind.data.frame(save, data.frame("model"="LIME", "run"="multiyear_singlefleet", "spr"=spr[length(spr)], "lcl_spr"=lcl_spr[length(spr)], "ucl_spr"=ucl_spr[length(spr)]))

ggplot(save) +
  geom_point(aes(x=run, y=spr, color=model), cex=5) +
  geom_linerange(aes(x=run, ymin=lcl_spr, ymax=ucl_spr, color=model)) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Run") + ylab("SPR")

plot_LCfits(Inputs=Inputs, Report=Report)

## plot model output
plot_output(Inputs=Inputs, 
            Report=Report,
            Sdreport=Sdreport, 
            lh=lh,
            true_years=years)
abline(v=lh$linf, lwd=2, lty=2, col="red")
abline(v=Report$ML_ft[length(Report$ML_ft)], lwd=2, lty=2, col="blue")

