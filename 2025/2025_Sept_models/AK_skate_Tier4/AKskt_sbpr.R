# Tier 4 calculations developed by K Omori and J Sullivan
# adapted for AK skates Sept 2025 by C Tribuzio

# Setup----
libs <- c("dplyr", "tidyr", "ggplot2")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# Functions----
# SPR -----
spr_fun <- function(
    R0 = 1,
    ages = 0:60,
    sex_ratio = 0.5,
    M = 0.106,
    Fspr = 0.4,
    
    # Length-at-age
    Linf = 35.02, # cm
    k = 0.122,
    t0 = -0.75,
    
    # Weight-at-length
    alpha =  0.0000119,
    beta = 3.06,
    
    # Maturity
    a50 = 10.06,
    delta = 1.96,
    maa = NULL,
    
    # Fishery selectivity
    Fdelta = 0.69,
    Fa50 = 13.50,
    selex = selexF
){
  
  # Life history
  nages <- length(ages)
  if(is.null(maa)) {
  maturity = 1/(1+exp(-delta*((ages) - a50)))
  } else {
    maturity = maa
  }
  if(length(maturity) != nages) stop("maturity vector doesn't match nages")
  
  # growth
  length = Linf * (1 - exp(-k*((ages) - t0))) # cm
  wt = alpha * length^beta
  
  # fishery selectivity
  if(is.null(selex)) {
    slx = 1/(1+exp(-Fdelta*((ages) - Fa50)))
  } else {
    slx = selex
  }
  if(length(slx) != nages) stop("selectivity vector doesn't match nages")
  
  
  # Run through ages
  N <- rep(NA, nages)
  Nf <- rep(NA, nages)
  N[1] <- R0
  Nf[1] <- R0
  for(i in 2:nages){
    N[i] = N[i-1] * exp(-M)
    Nf[i] = Nf[i-1] * exp(-M - Fspr*slx[i-1])
  }
  
  # Geometric series for plus group
  N[nages] = N[nages-1]*exp(-M)/(1-exp(-M))
  Nf[nages] = Nf[nages-1]*exp(-M - Fspr*slx[nages])/(1-exp(-M - Fspr*slx[nages]))
  
  # Yield
  Nc = 0
  for(age in 1:nages){
    # Baranov
    Nc[age] = Fspr*slx[age] / (Fspr*slx[age] + M) * (1 - exp(-M - Fspr*slx[nages])) * Nf[age]
  }
  
  # Calculate spawning biomass
  SB0 = sum( N * wt * maturity * sex_ratio)
  SBFspr = sum( Nf * wt * maturity * sex_ratio)
  Yield = sum(Nc * wt)
  
  # Return
  return(list(SB0 = SB0, SBFspr = SBFspr, Yield = Yield))
}

run_spr_fun <- function(F_vec = seq(0, 1, by = 0.001),
                        R0 = 1,
                        ages = 0:60,
                        sex_ratio = 0.5,
                        M = 0.106,

                        # Length-at-age
                        Linf = 35.02, # cm
                        k = 0.122,
                        t0 = -0.75,
                        
                        # Weight-at-length
                        alpha =  0.0000119,
                        beta = 3.06,
                        
                        # Maturity
                        a50 = 10.06,
                        delta = 1.96,
                        maa = NULL,

                        # Fishery selectivity
                        Fdelta = 0.69,
                        Fa50 = 13.50,
                        selex = selexF){
  
  # Calculated YPR/SBPR ----

  results <- lapply(F_vec, function(Fspr) {
    spr_fun(
      R0 = R0, ages = ages,
      sex_ratio = sex_ratio,
      M = M, Fspr = Fspr,
      Linf = Linf, k = k, t0 = t0,
      alpha = alpha, beta = beta,
      a50 = a50, delta = delta, maa = maa,
      Fdelta = Fdelta, Fa50 = Fa50,
      selex = selex
    )
  })
  
  SBFspr <- sapply(results, function(x) x$SBFspr)
  YPR <- sapply(results, function(x) x$Yield)
  SPRx <- SBFspr / SBFspr[1] # Relative depletion
  
  F40 <- F_vec[which.min(abs(SPRx - 0.40))]
  F35 <- F_vec[which.min(abs(SPRx - 0.35))]
  
  return(list(
    results = data.frame(F = F_vec, SPRx = SPRx, YPR = YPR),
    F40 = F40,
    F35 = F35
  ))
}
  
# Plot function
plot_spr_ypr <- function(res) {

    df <- res$results
    # YPR doesn't make much sense to me without an informed R0 so I made all of the YPRs relative
    df <- df %>% mutate(Rel_YPR = YPR / max(YPR))

    p1 <- ggplot(df, aes(x = SPRx, y = Rel_YPR)) +
      geom_line(size = 1.2) +
      geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue", size = 1) +
      geom_vline(xintercept = 0.35, color = "red", linetype = "dashed", size = 1) +
      labs(x = "SPR", y = "Relative YPR", title = "Relative YPR vs. SPR")
    
    p2 <- ggplot(df, aes(x = F, y = SPRx)) +
      geom_line(size = 1.2) +
      geom_vline(xintercept = res$F40, color = "blue", linetype = "dashed", size = 1) +
      geom_vline(xintercept = res$F35, color = "red", linetype = "dashed", size = 1) +
      labs(x = "Fishing Mortality (F)", y = "SPR", title = "SPR vs. F")
    
    p3 <- ggplot(df, aes(x = F, y = Rel_YPR)) +
      geom_line(size = 1.2) +
      geom_vline(xintercept = res$F40, color = "blue", linetype = "dashed", size = 1) +
      geom_vline(xintercept = res$F35, color = "red", linetype = "dashed", size = 1) +
      labs(x = "F", y = "Relative YPR", title = "Relative YPR vs. F")
    
    p4 <- ggplot(tidyr::pivot_longer(df, cols = c(SPRx, Rel_YPR), names_to = "Metric", values_to = "Value"), 
                 aes(x = F, y = Value, color = Metric)) +
      geom_line(size = 1.2) +
      geom_vline(xintercept = res$F40, color = "blue", linetype = "dashed", size = 1) +
      geom_vline(xintercept = res$F35, color = "red", linetype = "dashed", size = 1) +
      labs(x = "Fishing Mortality (F)", y = "Value", title = "SPR and Relative YPR vs. F")
    
    return(list(
      ypr_spr = p1,
      spr_F = p2,
      ypr_F = p3,
      spr_ypr_F = p4
    ))
  }

# Calculate values for AK skates----
# TODO need to reparameterize function for Gompertz growth, big diff in Linf
data_inputs <- list(#mat_a50 = 15.280, #Matta MS thesis
                    #mat_sigmaa = 1.574, #Matta MS thesis
                    vb_Linf = 135.39, #Matta and Gunderson 2007
                    vb_t0 = -1.60, #Matta and Gunderson 2007
                    vb_k = 0.10, #Matta and Gunderson 2007
                    natM = 0.13, #current assessment value
                    lw_a = 0.000009, #current assessment value
                    lw_b = 2.9617, #current assessment value
                    Amax = 26, #AGP data
                    selex = 1, # 1= knife edge, 2= logistic #TBD for now
                    selex_dat1 = 5, #TBD
                    selex_dat2 = NA #tBD
)

#Selectivity set up
#AK skates are caught at age 0, generally have a dome shaped selectivity, best approximation would be from Kotwicki and Weinberg 2005
# See Figure 5
# https://www.adfg.alaska.gov/static/home/library/PDFs/afrb/kotwv11n2.pdf

#placeholder values for now
# this curve highly influences results!
#selexF <- c(seq(0.5, 1, 0.25), rep(1, 13), rev(seq(0.5, 1, 0.05))) # results in low F rates
#length(selexF) #nages = 27 (max age plus 1 for age 0)
#plot(selexF) #Goal is for a reasonable dome shaped curve

selexF <- c(seq(0.1, 1, 0.2), rep(1, 16), rev(seq(0.5, 1, 0.1)))
length(selexF)
plot(selexF) #Goal is for a reasonable dome shaped curve

#selexF <- c(seq(0.1, 1, 0.2), rep(1, 22)) #ramp up and fully selected for whole age range, not a huge diff in F rates from the above, not a lot of large animals
#length(selexF)
#plot(selexF)

#selexF <- c(rep(0.1, data_inputs$selex_dat1), rep(1, data_inputs$Amax-data_inputs$selex_dat1+1)) #knife edge
#length(selexF)
#plot(selexF)

#selexF <- c(rep(1, data_inputs$Amax+1)) #fully selected at all age classes
#length(selexF)
#plot(selexF)

#selexF <- c(rep(0.1, data_inputs$selex_dat1), rep(1, 15), rep(0.75, 7)) #knife edge
#length(selexF)
#plot(selexF)

# maturity curve set up
# ignoring the built in function and making curve based on Matta MS Thesis
ages <- 0:26
mata_fem <- -15.280
matb_fem <- 1.574
maa <- 1/(1+exp(-(mata_fem+matb_fem*ages)))

# run the functions to calc F rates
res <- run_spr_fun(F_vec = seq(0, 1, by = 0.001),
                   R0 = 1,
                   ages = 0:26,
                   sex_ratio = 0.5,
                   M = data_inputs$natM,
                   # Length-at-age
                   Linf = data_inputs$vb_Linf, # cm
                   k = data_inputs$vb_k,
                   t0 = data_inputs$vb_t0,
                   # Weight-at-length
                   alpha =  data_inputs$lw_a,
                   beta = data_inputs$lw_b,
                   # Maturity
                   #a50 = data_inputs$mat_a50,
                   #delta = data_inputs$mat_sigma,
                   maa = maa,
                   # Fishery selectivity
                   #Fdelta = 0.69,
                   #Fa50 = 13.50,
                   selex = selexF)
res$F40
res$F35

plots <- plot_spr_ypr(res)
plots$ypr_spr
plots$spr_F
plots$ypr_F
plots$spr_ypr_F

# Tier 4 harvest recommendations----
AKsk_rema <- read_csv(paste0(getwd(), "/2025/2025_Sept_models/AK_skate_Tier5/Tier5_AKskt_output.csv"))
T4recs <- AKsk_rema %>% 
  group_by(year) %>% 
  summarise(rbiom = sum(pred),
            OFL = rbiom * res$F35,
            ABC = rbiom * res$F40) %>% 
  filter(year == 2023) %>% 
  mutate(FOFL = res$F35,
         FABC = res$F40,
         model_name = "Tier 4")

