
libs <- c("dplyr", "tidyr", "ggplot2")
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# inputs from sharpchin from Kristen Omori, this code has been cross-checked
# with kristen for our spp, we get the same values
data_inputs <- list(mat_a50 = 10.06,
                    mat_sigmaa = 1.96,
                    vb_Linf = 35.02,
                    vb_t0 = -0.75,
                    vb_k = 0.122,
                    natM = 0.106,
                    lw_a = 0.0000119,
                    lw_b =3.06,
                    Amax = 60,
                    selex = 1, # 1= knife edge, 2= logistic
                    selex_dat1 = 10,
                    selex_dat2 = NA
)

# knife edge
selexF <- c(rep(0, data_inputs$selex_dat1), rep(1, data_inputs$Amax-data_inputs$selex_dat1+1))
length(selexF)

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

res <- run_spr_fun() # can put in different life history inputs here, default to sharpcin
res$F40
res$F35

plots <- plot_spr_ypr(res)
plots$ypr_spr
plots$spr_F
plots$ypr_F
plots$spr_ypr_F

# check RE ----

mat_input_re <- read.csv(here::here('results', 'maturity_results_3_52.csv')) %>%
  dplyr::filter(species == 'RE'& Maturity == 'Biological') %>%
  dplyr::pull(p)

fish_selex_input <- c(0.000561359, 0.00129505, 0.00295123, 0.00658641, 0.0139612,
                      0.0275503, 0.0480452, 0.0659263, 0.063593, 0.0569547, 0.0522039,
                      0.0757491, 0.203113, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1)
length(mat_input_re)
length(fish_selex_input)
res <- run_spr_fun(ages = 3:52, M = 0.042, 
                   Linf = 54.4, k = 0.102, t0 = -0.36,
                   alpha =  0.00001143, beta = 3.098,
                   maa = mat_input_re,
                   selex = fish_selex_input)
res$F40
res$F35
names(res)
res$results %>% head
plots <- plot_spr_ypr(res)
plots$ypr_spr
plots$spr_F
plots$ypr_F
plots$spr_ypr_F

# Check BS ----

mat_input_bs <- read.csv(here::here('results', 'maturity_results_3_52.csv')) %>%
  dplyr::filter(species == 'BS'& Maturity == 'Functional') %>%
  dplyr::pull(p)

fish_selex_input <- c(0.000561359, 0.00129505, 0.00295123, 0.00658641, 0.0139612,
                      0.0275503, 0.0480452, 0.0659263, 0.063593, 0.0569547, 0.0522039,
                      0.0757491, 0.203113, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1)
res <- run_spr_fun(ages = 3:52, M = 0.042, 
                    Linf = 52.6, k = 0.062, t0 = 0.2,
                    alpha = 0.00000813,
                    beta = 3.177,
                    maa = mat_input_bs,
                    selex = fish_selex_input)
res$F40
res$F35

plots <- plot_spr_ypr(res)
plots$ypr_spr
plots$spr_F
plots$ypr_F
plots$spr_ypr_F
