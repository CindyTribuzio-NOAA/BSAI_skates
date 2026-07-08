# Tier 4 sbpr derivation of selectivity curve
# Cindy Tribuzio
# based on Kotwicki and Weinberg 2005 underbag study
# see email from Stan Kotwicki 9 Sept 2025
# https://www.adfg.alaska.gov/static/home/library/PDFs/afrb/kotwv11n2.pdf

# Setup----
libs <- c("tidyverse", "janitor")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

AYR <- 2026

# get data----
# run AFSC_OBS_length_comp_data.R in the Tier 3 folder structure first
L_dat <- read_csv(here::here(AYR, 'Nov_models', 'AK_skate_Tier4', 'data', 'fisherylcompall_AKskt2023.csv')) %>% 
  clean_names() %>% 
  select(c(sex, length, frequency, fmp_gear)) %>% 
  filter(fmp_gear != 'POT')

# get proportions by size----
gear_tots <- L_dat %>% 
  group_by(fmp_gear) %>% 
  summarise(gtot = sum(frequency))

alllength <- sum(gear_tots$gtot)

tot_lengths <- L_dat %>% 
  group_by(length, fmp_gear) %>% 
  summarise(ntotal = sum(frequency),
            lprop = ntotal/alllength)

ggplot(tot_lengths, aes(x = length, y = lprop, color = fmp_gear))+
  geom_point()+
  geom_line()

# convert data to ages----
#convert to ages using same parameters as Tier 3 model, sexes combined
L0 <- 23.95
Linf <- 104.77
g <- -1.90
k <- 0.36

#A parameters
A1 <- 0 # age 0's are caught in the data
A2 <- 15 # mean age at which Linf is achieved

lage_dat <- L_dat %>%
  mutate(
    expr = 1 - ((length^g - L0^g)/(Linf^g - L0^g)) * (1 - exp(-k*(A2 - A1))),
    age = ifelse(expr > 0,
                    round(A1 - (1/k) * log(expr), 0),
                    NA_real_)
  )

# change the -1 ages to be zero, some small females in the data
# all animals with TL > Linf return na, convert them to 20 for a plus group for this analysis
lage_dat <- lage_dat %>% 
  mutate(age = if_else(age < 0, 0, age),
         age = if_else(is.na(age), 20, age))

# mean length at age----
mlage <- lage_dat %>% 
  group_by(age) %>% 
  summarise(mlength = mean(length))


#proportions by ages ----
tot_ages <- lage_dat %>% 
  group_by(fmp_gear) %>% 
  summarise(totage = sum(frequency))
age_dat <- lage_dat %>% 
  group_by(age, fmp_gear) %>% 
  summarise(nages = sum(frequency),
            .groups = 'drop') %>% 
  left_join(tot_ages, by = 'fmp_gear') %>% 
  mutate(aprop = nages/totage)

ll_scaled_age <- age_dat %>%
  filter(fmp_gear == 'HAL') %>% 
  mutate(scprop = (aprop - min(aprop))/(max(aprop)-min(aprop)))

t_scaled_age <- age_dat %>%
  filter(fmp_gear == 'TRW') %>% 
  mutate(scprop = (aprop - min(aprop))/(max(aprop)-min(aprop)))

scaled_age <- ll_scaled_age %>% 
  bind_rows(t_scaled_age)

ggplot(scaled_age, aes(x = age, y = scprop, color = fmp_gear))+
  geom_point()+
  geom_line()

# base LL selectivity curve----
# define the linear ramp up
LLsrupval <- 0.05
LLsrupage <- 0
LLerupage <- 7
LLerupval <- 1
LLstartramp <- seq(LLsrupval, LLerupval, (LLerupval - LLsrupval)/(LLerupage-LLsrupage))

# define the descending limb
LLsrdnval <- 1
LLsrdnage <- 13
LLerdnage <- 26
LLerdnval <- 0.75
LLendramp <- seq(LLsrdnval, LLerdnval, (LLerdnval - LLsrdnval)/(LLerdnage-LLsrdnage))

LLtopwidth <- (LLsrdnage - LLerupage)-1
LLselexF <- c(LLstartramp, rep(1, LLtopwidth), LLendramp)
length(LLselexF)
plot(LLselexF) 

plotLLselexF <- LLselexF %>% 
  bind_cols(seq(0,26)) %>% 
  mutate(fmp_gear = 'HAL') %>% 
  rename(age = `...2`,
         selex = `...1`)

# base TWL selectivity curve----
# define the linear ramp up
TWLsrupval <- 0.05
TWLsrupage <- 0
TWLerupage <- 7
TWLerupval <- 1
TWLstartramp <- seq(TWLsrupval, TWLerupval, (TWLerupval - TWLsrupval)/(TWLerupage-TWLsrupage))

# define the descending limb
TWLsrdnval <- 1
TWLsrdnage <- 11
TWLerdnage <- 26
TWLerdnval <- 0.75
TWLendramp <- seq(TWLsrdnval, TWLerdnval, (TWLerdnval - TWLsrdnval)/(TWLerdnage-TWLsrdnage))

TWLtopwidth <- (TWLsrdnage - TWLerupage)-1
TWLselexF <- c(TWLstartramp, rep(1, TWLtopwidth), TWLendramp)
length(TWLselexF)
plot(TWLselexF) 

plotTWLselexF <- TWLselexF %>% 
  bind_cols(seq(0,26)) %>% 
  mutate(fmp_gear = 'TRW') %>% 
  rename(age = `...2`,
         selex = `...1`)

selexF_dat <- plotLLselexF %>% 
  bind_rows(plotTWLselexF)

ggplot(scaled_age, aes(x = age, y = scprop, color = fmp_gear, fill = fmp_gear))+
  geom_bar(stat = 'identity')+
  facet_grid(fmp_gear~.)+
  scale_fill_viridis_d(option = 'magma')+
  scale_color_viridis_d(option = 'magma')+
  geom_line(data = selexF_dat, aes(x = age, y = selex))
