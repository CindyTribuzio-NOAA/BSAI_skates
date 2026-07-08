# Tier 4 calculations developed by K Omori and J Sullivan
# adapted for AK skates April 2026 by C Tribuzio
# can not use the built in sensitivity functions because only look at M and knife edge
# this code breaks it up into chapters: M, logistic, dome-shaped

# Setup----
#devtools::install_github("noaa-afsc/tier4tools", build_vignettes = TRUE)

libs <- c("dplyr", "tidyr", "ggplot2", 'tier4tools', 'purrr')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# --- THE AGE-BASED SS3 DOUBLE NORMAL FUNCTION ---
ss3_double_normal_age <- function(ages, peak_age, top_logit, ascend_se, descend_se, start_logit, end_logit, max_age = 25) {
  
  # Step A: Convert logit/log values to normal space constants
  # Inverse-logit for top width proportion relative to max age
  p_top <- 1 / (1 + exp(-top_logit))
  top_width_age <- p_top * (max_age - peak_age)
  
  # Exponentiate widths to get standard deviations in age space
  sigma_asc  <- exp(ascend_se)
  sigma_desc <- exp(descend_se)
  
  # Inverse-logit for tail heights (handling SS3 special codes)
  p_start <- if (start_logit == -999) 0 else 1 / (1 + exp(-start_logit))
  p_end   <- if (end_logit == 99) 1 else 1 / (1 + exp(-end_logit))
  
  # Step B: Pre-calculate scaling constraints for the tail boundaries
  # t1 isolates the curve's baseline value at Age 0
  t1 <- exp(-(0 - peak_age)^2 / sigma_asc^2) 
  
  # t2 isolates the curve's baseline value at the maximum age
  end_age <- max(ages)
  t2 <- exp(-(end_age - (peak_age + top_width_age))^2 / sigma_desc^2)
  
  # Step C: Compute selectivity for every single age class
  sel <- sapply(ages, function(A) {
    if (A < peak_age) {
      # --- ASCENDING LIMB ---
      s <- exp(-(A - peak_age)^2 / sigma_asc^2)
      
      if (start_logit == -999) {
        # FIX: Rescale the ascending limb so it hits exactly 0 at Age = 0
        s <- (s - t1) / (1 - t1)
      } else {
        # Standard SS3 logit-bounded tail height math
        s <- p_start + (1 - p_start) * (s - t1) / (1 - t1)
      }
      return(s)
      
    } else if (A >= peak_age && A <= (peak_age + top_width_age)) {
      # --- FLAT TOP PLATEAU ---
      return(1.0)
      
    } else {
      # --- DESCENDING LIMB ---
      if (end_logit == 99) {
        return(1.0) # Special switch: Bypasses dome completely (Asymptotic/Logistic)
      } else {
        s <- exp(-(A - (peak_age + top_width_age))^2 / sigma_desc^2)
        s <- 1.0 + (p_end - 1.0) * (1.0 - s) / (1.0 - t2)
        return(s)
      }
    }
  })
  
  # Prevent tiny negative decimal artifacts from floating-point rounding
  return(pmax(0, sel)) 
}



# Setting up inputs----
AYR <- 2026

#proportion of catch by gear type----
#mean proportion over the most recent 5 years
#from 2_catch_gear_target_summary BSAI_m5gear
#update each assessment cycle, things change
pHAL <- 0.834 
pTRW <- 0.166

#age and growth parameters----
ages <- 0:26 # max age from AGP data

# ignoring the built in function and making curve based on Matta MS Thesis
ages <- 0:26
mata_fem <- -15.280
matb_fem <- 1.574
maa <- 1/(1+exp(-(mata_fem+matb_fem*ages)))

# Sens M----
# using logistic selectivity for simplicity in this case
inp_M <- spr_input(
  ages = ages,
  species = list(
    AK_skate = list(
      len_at_age = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), #Matta and Gunderson 2007
      wt_at_age  = list(type = "wl", alpha = 0.000009, beta = 2.9617), #Matta and Gunderson 2007
      maturity  = list(type = "vector", m = maa),
      M = 0.13, #current assessment value
      selectivity = list(
        type = "fleets",
        fleets = list(
          twl = list(
            propF = pTRW,
            selex = list(type = "logistic", a50 = 3, delta = -0.5)#need to make these vectors of selex @ age
          ),
          ll = list(
            propF = pHAL,
            selex = list(type = "logistic", a50 = 7, delta = -2)
          )
        )
      )
    )
  ),
  use_plus_group = FALSE
)

inp_M_plots <- plot_spr_inputs(inp_M, panels = c("mat_selex", "wt", "len", "M"))
grid_spr_inputs(inp_M_plots, order = c("mat_selex", "wt","len","M"), ncol = 2, guides = "keep")
spr_M_out <- run_spr(
  inp_M,
  diagnostics = TRUE
)
spr_M_out$F_spr_total
plot_spr_curves(spr_M_out)

# setting up the M scenarios
M_vec <- seq(0.07, 0.23, by = 0.01)

M_out <- as.data.frame(matrix(nrow = length(M_vec), ncol = 3, dimnames = list(NULL, c('F40', 'F35', 'M'))))

for(i in 1:length(M_vec)){
  #loop_inp <- inp_M
  #loop_inp$species$AK_skate$spec$M <- M_vec[i] #not sure why this doesn't work, have to do the whole inp over again each loop
  inp_loop <- spr_input(
    ages = ages,
    species = list(
      AK_skate = list(
        len_at_age = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), #Matta and Gunderson 2007
        wt_at_age  = list(type = "wl", alpha = 0.000009, beta = 2.9617), #Matta and Gunderson 2007
        maturity  = list(type = "vector", m = maa),
        M = M_vec[i],
        selectivity = list(
          type = "fleets",
          fleets = list(
            twl = list(
              propF = pTRW,
              selex = list(type = "logistic", a50 = 3, delta = -0.5)#need to make these vectors of selex @ age
            ),
            ll = list(
              propF = pHAL,
              selex = list(type = "logistic", a50 = 7, delta = -2)
            )
          )
        )
      )
    ),
    use_plus_group = FALSE
  )
  spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
  loop_out <-t(as.data.frame(spr_loop_out$F_spr_total)) %>% 
    bind_cols(M_vec[i])
  M_out[i,] <- loop_out
}

M_out <- M_out %>% 
  pivot_longer(!M, names_to = "type", values_to = "F_rate") 
h_lines <- M_out %>% 
  filter(M == 0.13)
ph_lines <- M_out %>% 
  filter(round(M, 2) == 0.21) #had to do the round due to some weird R precision floating intenger thing
T3_lines <- M_out %>% 
  filter(round(F_rate, 3) == 0.098)
M_plot <- ggplot(M_out, aes(x = M, y = F_rate, color = type))+
  geom_line(linewidth = 1)+
  geom_segment(data = h_lines, aes(x = M, xend = M, y = -Inf, yend = max(F_rate)), linetype = 'dashed', color = "black")+
  geom_segment(data = h_lines, aes(x = M, y = F_rate, xend = min(M_out$M), yend = F_rate), linetype = "dashed", color = 'black')+
  geom_segment(data = ph_lines, aes(x = M, xend = M, y = -Inf, yend = max(F_rate)), color = 'grey70', linetype = 'dashed')+
  geom_segment(data = ph_lines, aes(x = M, y = F_rate, xend = min(M_out$M), yend = F_rate), color = 'grey70', linetype = "dashed")+
  geom_segment(data = T3_lines, aes(x = M, y = F_rate, xend = min(M_out$M), yend = F_rate), color = 'grey70', linetype = "dashed")+
  geom_segment(data = T3_lines, aes(x = M, y = F_rate, xend = M, yend = -Inf), color = 'grey70', linetype = "dashed")+
  #geom_hline(yintercept = 0.098, linetype = "dashed", color = "grey70")+
  geom_text(
    data = h_lines,
    aes(x = 0.101, y = F_rate, label = round(F_rate, 3)), # Rounds to 3 decimal places
    color = "grey30",
    vjust = -0.5,     # Nudges text slightly above the dashed line
    hjust = 0,        # Left-aligns the text
    show.legend = FALSE # Prevents "a" from appearing in your legend box
  )+ 
  annotate(
    "text",
    x = 0.13, 
    y = 0.09, 
    label = "Current M=0.13",
    family = "sans",   # Matches ggplot2's default theme font
    hjust = -0.1
  ) + 
  geom_text(
    data = ph_lines,
    aes(x = M, y = F_rate, label = round(F_rate, 3)), # Rounds to 3 decimal places
    color = "grey70",
    vjust = -0.5,     # Nudges text slightly above the dashed line
    hjust = 1,        # Left-aligns the text
    show.legend = FALSE # Prevents "a" from appearing in your legend box
  )+ 
  annotate(
    "text",
    x = 0.21, 
    y = 0.09, 
    color = 'grey70',
    label = "Possible M=0.21",
    family = "sans",   # Matches ggplot2's default theme font
    hjust = -0.1
  ) + 
  annotate("text", x = 0.08, y = 0.099, color = 'grey70', label = "Tier 3 F40=0.098", family = 'sans')+
  annotate("text", x = 0.095, y = 0.09, color = 'grey70', label = "M=0.09", family = 'sans')+
  labs(y = "F rate", x = "Natural Mortality (M)")+
  #scale_color_viridis_d(option = 'turbo')+
  scale_color_brewer(palette = "Dark2", guide = guide_legend(title = NULL))+
  scale_x_continuous(expand = c(0, 0))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

# Sens logistic----
ages <- seq(0, 26, by = 1)

# set up for sensitivity runs
rana50 <- seq(0, 10, by = 1)
randelt <- seq(0, 5, by = 0.1)
ranlpars <- expand_grid(rana50 = rana50, randelt = randelt)

lslx_out <- as.data.frame(matrix(nrow = nrow(ranlpars), ncol = 4, dimnames = list(NULL, c('F40', 'F35', 'a50', 'delta'))))

#run spr and make a50 vs delt heat map with color by F40
# Pre-allocate lslx_out with NA before the loop to match expected final dimensions
# (Ensure columns match the total width of your loop_out dataframe)
for(i in 1:nrow(ranlpars)){
  loop_out <- tryCatch({
    inp_loop <- spr_input(
      ages = ages,
      species = list(
        AK_skate = list(
          len_at_age = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), 
          wt_at_age  = list(type = "wl", alpha = 0.000009, beta = 2.9617), 
          maturity  = list(type = "vector", m = maa),
          M = 0.13,
          selectivity = list(type = "logistic", a50 = ranlpars[i,1], delta = ranlpars[i,2])
        )
      ),
      use_plus_group = FALSE
    )
    spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
    # Combine results successfully
    t(as.data.frame(spr_loop_out$F_spr_total)) %>% 
      bind_cols(ranlpars[i, ])
  }, error = function(e) {
    # This prints when an error occurs, showing you exactly which row failed
    message(paste("Error encountered at row", i, ":", e$message))
    # Return NULL so we know this iteration failed
    return(NULL)
  })
  # Check if an error occurred
  if (is.null(loop_out)) {
    lslx_out[i,] <- c(rep(NA, 2), ranlpars[i, ])
    next # Skip immediately to the next row in ranlpars
  }
  # Save successful output
  lslx_out[i,] <- loop_out
}

# plot heatmap
lslx_plot <- ggplot(lslx_out, aes(x = a50, y = delta, fill = F40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  labs(
    title = "Logistic Selectivity Parameters",
    #subtitle = "AK Skate Sensitivity Analysis",
    x = "a50",
    y = "delta"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

# flip: logistic slx----
# what logistic curves will result in F40 similar to Tier 3 models?
#why don't all combos that have F40 values also have curves??? Rounding error somewhere??

lslx_thin <- lslx_out %>% 
  filter(between(F40, 0.097, 0.13)) %>% 
  expand_grid(age = ages) %>%
  mutate(
    slx = 1 / (1 + exp(-(age - a50) / delta)))

F40_lslx_plot <- ggplot(lslx_thin, aes(x = age, y = slx, color = F40)) +
  # Draw a line for every single distinct parameter combination
  geom_line(
    aes(group = interaction(a50, delta)), 
    alpha = 0.5,
    linewidth = 0.4
  ) +
  scale_color_viridis_c(option = "viridis") + # Continuous Viridis scale
  labs(
    title = "Logistic Curves with similar F40 results",
    x = "Age",
    y = "Selectivity"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )


# Sens M vs Log ----
ranlpars <- expand.grid(
  M   = seq(0.08, 0.23, by = 0.01),
  a50 = seq(0, 10, by = 1),
  delta = seq(0, 5, by = 0.1) 
) %>%
  mutate(
    scenario    = paste0("Scenario_", row_number())
  )

unique_scenarios <- unique(ranlpars$scenario)
total_scenarios  <- length(unique_scenarios)

#run spr and make a50 vs delt heat map with color by F40
results_collector <- tibble(
  scenario = unique_scenarios,
  F_40     = rep(NA_real_, total_scenarios),
  F_35     = rep(NA_real_, total_scenarios),
  status   = rep("Not Run", total_scenarios)
)

# 3. Execute the sensitivity loop
for(i in 1:total_scenarios) {
  
  current_scenario <- unique_scenarios[i]
  
  # A single, clean tryCatch to isolate the modeling process
  loop_out <- tryCatch({
    # Build your reference point engine input
    inp_loop <- spr_input(
      ages = ages,
      species = list(
        AK_skate = list(
          len_at_age  = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), 
          wt_at_age   = list(type = "wl", alpha = 0.000009, beta = 2.9617), 
          maturity    = list(type = "vector", m = maa),
          M           = ranlpars[i,1],
          selectivity = list(type = "logistic", a50 = ranlpars[i,2], delta = ranlpars[i,3]) 
        )
      ),
      use_plus_group = FALSE
    )
    
    # Calculate target reference points
    spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
    
    # Pass the calculated outputs forward to the collector assignment step
    # Note: Replace 'F40_total'/'F35_total' with the exact names your package uses
    list(
      F40 = spr_loop_out$F_spr_total["F40_total"], 
      F35 = spr_loop_out$F_spr_total["F35_total"], 
      status = "Success"
    )
    
  }, error = function(e) {
    # Catches crashes gracefully and prints a notification
    message(paste("Skipping", current_scenario, "due to evaluation error:", e$message))
    return(list(F40 = NA_real_, F35 = NA_real_, status = paste("Error:", e$message)))
  })
  
  # --- 3. Assign the results back to the pre-allocated data frame rows ---
  results_collector$F_40[i]   <- loop_out$F40
  results_collector$F_35[i]   <- loop_out$F35
  results_collector$status[i] <- loop_out$status
}

final_ML_results <- ranlpars %>%
  left_join(results_collector, by = "scenario") %>% 
  filter(!is.na(F_40))

Ml_plot <- ggplot(final_ML_results, aes(x = M, y = delta, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_wrap(~a50, ncol = 4,
             labeller = labeller(
               a50  = label_both,
             ))+
  labs(
    title = "M and logistic sensitivity",
    #subtitle = "Fixed top width and start value",
    x = "Natural Mortality (M)",
    y = "delta"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

# flip: M/logistic slx----
# what logistic curves will result in F40 similar to Tier 3 models?
#why don't all combos that have F40 values also have curves??? Rounding error somewhere??

Mlslx_thin <- final_ML_results %>% 
  filter(between(F_40, 0.097, 0.13)) %>% 
  expand_grid(age = ages) %>%
  mutate(
    slx = 1 / (1 + exp(-(age - a50) / delta)))

F40_Mlslx_plot <- ggplot(Mlslx_thin, aes(x = age, y = slx, color = F_40)) +
  # Draw a line for every single distinct parameter combination
  geom_line(
    aes(group = interaction(a50, delta)), 
    alpha = 0.5,
    linewidth = 0.4
  ) +
  facet_wrap(~M, ncol = 4,
             labeller = labeller(
               M  = label_both,
             ))+
  scale_color_viridis_c(option = "viridis") + # Continuous Viridis scale
  labs(
    title = "M impacted Logistic Curves with similar F40 results",
    x = "Age",
    y = "Selectivity"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )



# Sens dome----
# take a faceted heat map approach

# Define lengths and ages matching your stock assessment specifications
lengths_grid <- seq(10, 160, by = 1)
ages_grid <- ages # Maps to your existing 'ages' object

#fixed top and start----
# Create a grid varying Peak, Ascending Width, and Descending Width
# Keep top, start, and end parameters fixed to the Longline baseline defaults
ranlpars <- expand.grid(
  peak_age   = seq(6, 14, by = 2),        # Shifting peak selectivity age (Years)
  ascend_se  = seq(1.0, 3.0, by = 0.5),   # Entry scale parameter in age-variance space
  descend_se = seq(1.0, 3.0, by = 0.5),   # Dome decay rate parameter in age-variance space
  end_logit  = seq(-5.0, 5, by = 1.0)  
) %>%
  mutate(
    top_logit   = -9.52666, 
    start_logit = -999,
    scenario    = paste0("Scenario_", row_number())
  )

# look at the curves, do they do what you want in the sensitivity?
# Generate curve values using the SS3 math function
all_curves_age_df <- ranlpars %>%
  mutate(
    Curve_Data = pmap(
      list(peak_age, top_logit, ascend_se, descend_se, start_logit, end_logit), 
      function(p, tl, ase, dse, sl, el) {
        tibble(
          Age = ages_grid,
          Selectivity = ss3_double_normal_age(
            ages        = ages_grid,
            peak_age    = p,
            top_logit   = tl,
            ascend_se   = ase,
            descend_se  = dse,
            start_logit = sl,
            end_logit   = el,
            max_age     = max(ages_grid)
          )
        )
      }
    )
  ) %>%
  unnest(cols = c(Curve_Data))

# Generate the faceted structural diagnostic plot
ggplot(all_curves_age_df, aes(x = Age, y = Selectivity)) +
  geom_line(aes(group = scenario), color = "darkgreen", alpha = 0.04, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, max(ages_grid), by = 5)) +
  labs(
    title = "fixed top and start",
    subtitle = paste("Ensemble envelope across", length(unique(all_curves_age_df$scenario)), "age-structured scenarios"),
    x = "Age (Years)",
    y = "Selectivity Proportion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Get SPR outputs
# 1. Get the list of all unique scenarios to iterate over
unique_scenarios <- unique(all_curves_age_df$scenario)
total_scenarios  <- length(unique_scenarios)

# 2. Pre-allocate a data frame to store reference point outcomes cleanly
# (Using a data frame instead of a matrix handles strings and numbers safely)
results_collector <- tibble(
  scenario = unique_scenarios,
  F_40     = rep(NA_real_, total_scenarios),
  F_35     = rep(NA_real_, total_scenarios),
  status   = rep("Not Run", total_scenarios)
)

# 3. Execute the sensitivity loop
for(i in 1:total_scenarios) {
  
  current_scenario <- unique_scenarios[i]
  
  # A single, clean tryCatch to isolate the modeling process
  loop_out <- tryCatch({
    
    # Extract the pre-calculated age selectivity vector for this scenario
    sel_age_vector <- all_curves_age_df %>% 
      filter(scenario == current_scenario) %>% 
      arrange(Age) %>% 
      pull(Selectivity)
    
    # Build your reference point engine input
    inp_loop <- spr_input(
      ages = ages,
      species = list(
        AK_skate = list(
          len_at_age  = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), 
          wt_at_age   = list(type = "wl", alpha = 0.000009, beta = 2.9617), 
          maturity    = list(type = "vector", m = maa),
          M           = 0.13,
          selectivity = list(type = "vector", sel = sel_age_vector) 
        )
      ),
      use_plus_group = FALSE
    )
    
    # Calculate target reference points
    spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
    
    # Pass the calculated outputs forward to the collector assignment step
    # Note: Replace 'F40_total'/'F35_total' with the exact names your package uses
    list(
      F40 = spr_loop_out$F_spr_total["F40_total"], 
      F35 = spr_loop_out$F_spr_total["F35_total"], 
      status = "Success"
    )
    
  }, error = function(e) {
    # Catches crashes gracefully and prints a notification
    message(paste("Skipping", current_scenario, "due to evaluation error:", e$message))
    return(list(F40 = NA_real_, F35 = NA_real_, status = paste("Error:", e$message)))
  })
  
  # --- 3. Assign the results back to the pre-allocated data frame rows ---
  results_collector$F_40[i]   <- loop_out$F40
  results_collector$F_35[i]   <- loop_out$F35
  results_collector$status[i] <- loop_out$status
}

# 4. Bind the results directly back to your structural parameter grid
final_sensitivity_results <- ranlpars %>%
  left_join(results_collector, by = "scenario")

# make heatmap 
dome1 <- ggplot(final_sensitivity_results, aes(x = ascend_se, y = descend_se, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~end_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               end_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed top width and start value",
    x = "ascending se",
    y = "descending se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

dome1_zoom_dat <- filter(final_sensitivity_results) %>% 
  filter(between(F_40, 0.09, 0.13))
dome1_zoom <- ggplot(dome1_zoom_dat, aes(x = ascend_se, y = descend_se, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~end_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               end_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed top width and start value",
    x = "ascending se",
    y = "descending se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

#fixed ascend and descent----
# Create a grid varying Peak, Ascending Width, and Descending Width
# Keep top, start, and end parameters fixed to the Longline baseline defaults
ranlpars <- expand.grid(
  peak_age   = seq(6, 20, by = 2),        # Shifting peak selectivity age (Years)
  top_logit   = seq(-5.0, 5, by = 1.0),  
  start_logit = seq(-5.0, 5, by = 1.0),  
  end_logit  = seq(-5.0, 5, by = 1.0)  
) %>%
  mutate(
    ascend_se  = 2,   # Entry scale parameter in age-variance space
    descend_se = 2,   # Dome decay rate parameter in age-variance space
    scenario    = paste0("Scenario_", row_number())
  )

# look at the curves, do they do what you want in the sensitivity?
# Generate curve values using the SS3 math function
all_curves_age_df <- ranlpars %>%
  mutate(
    Curve_Data = pmap(
      list(peak_age, top_logit, ascend_se, descend_se, start_logit, end_logit), 
      function(p, tl, ase, dse, sl, el) {
        tibble(
          Age = ages_grid,
          Selectivity = ss3_double_normal_age(
            ages        = ages_grid,
            peak_age    = p,
            top_logit   = tl,
            ascend_se   = ase,
            descend_se  = dse,
            start_logit = sl,
            end_logit   = el,
            max_age     = max(ages_grid)
          )
        )
      }
    )
  ) %>%
  unnest(cols = c(Curve_Data))

# Generate the faceted structural diagnostic plot
ggplot(all_curves_age_df, aes(x = Age, y = Selectivity)) +
  geom_line(aes(group = scenario), color = "darkgreen", alpha = 0.04, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, max(ages_grid), by = 5)) +
  labs(
    title = "fixed ascent and descent",
    subtitle = paste("Ensemble envelope across", length(unique(all_curves_age_df$scenario)), "age-structured scenarios"),
    x = "Age (Years)",
    y = "Selectivity Proportion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Get SPR outputs
# 1. Get the list of all unique scenarios to iterate over
unique_scenarios <- unique(all_curves_age_df$scenario)
total_scenarios  <- length(unique_scenarios)

# 2. Pre-allocate a data frame to store reference point outcomes cleanly
# (Using a data frame instead of a matrix handles strings and numbers safely)
results_collector <- tibble(
  scenario = unique_scenarios,
  F_40     = rep(NA_real_, total_scenarios),
  F_35     = rep(NA_real_, total_scenarios),
  status   = rep("Not Run", total_scenarios)
)

# 3. Execute the sensitivity loop
for(i in 1:total_scenarios) {
  
  current_scenario <- unique_scenarios[i]
  
  # A single, clean tryCatch to isolate the modeling process
  loop_out <- tryCatch({
    
    # Extract the pre-calculated age selectivity vector for this scenario
    sel_age_vector <- all_curves_age_df %>% 
      filter(scenario == current_scenario) %>% 
      arrange(Age) %>% 
      pull(Selectivity)
    
    # Build your reference point engine input
    inp_loop <- spr_input(
      ages = ages,
      species = list(
        AK_skate = list(
          len_at_age  = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), 
          wt_at_age   = list(type = "wl", alpha = 0.000009, beta = 2.9617), 
          maturity    = list(type = "vector", m = maa),
          M           = 0.13,
          selectivity = list(type = "vector", sel = sel_age_vector) 
        )
      ),
      use_plus_group = FALSE
    )
    
    # Calculate target reference points
    spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
    
    # Pass the calculated outputs forward to the collector assignment step
    # Note: Replace 'F40_total'/'F35_total' with the exact names your package uses
    list(
      F40 = spr_loop_out$F_spr_total["F40_total"], 
      F35 = spr_loop_out$F_spr_total["F35_total"], 
      status = "Success"
    )
    
  }, error = function(e) {
    # Catches crashes gracefully and prints a notification
    message(paste("Skipping", current_scenario, "due to evaluation error:", e$message))
    return(list(F40 = NA_real_, F35 = NA_real_, status = paste("Error:", e$message)))
  })
  
  # --- 3. Assign the results back to the pre-allocated data frame rows ---
  results_collector$F_40[i]   <- loop_out$F40
  results_collector$F_35[i]   <- loop_out$F35
  results_collector$status[i] <- loop_out$status
}

# 4. Bind the results directly back to your structural parameter grid
final_sensitivity_results <- ranlpars %>%
  left_join(results_collector, by = "scenario")

# make heatmap 
dome2 <- ggplot(final_sensitivity_results, aes(x = start_logit, y = end_logit, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~top_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               top_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed ascending and descending se",
    x = "start logit",
    y = "end logit"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

dome2_zoom_dat <- filter(final_sensitivity_results) %>% 
  filter(between(F_40, 0.09, 0.13))
dome2_zoom <- ggplot(dome2_zoom_dat, aes(x = start_logit, y = end_logit, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~top_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               top_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed ascending and descending se",
    x = "ascending se",
    y = "descending se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

# fixed top and end logit----
ranlpars <- expand.grid(
  peak_age   = seq(6, 20, by = 2),        # Shifting peak selectivity age (Years)
  ascend_se  = seq(1.0, 3.0, by = 0.5),   # Entry scale parameter in age-variance space
  descend_se = seq(1.0, 3.0, by = 0.5),   # Dome decay rate parameter in age-variance space
  start_logit = seq(-5.0, 5, by = 1.0)
) %>%
  mutate(
    top_logit   = -10, 
    end_logit  = 0,
    scenario    = paste0("Scenario_", row_number())
  )

# look at the curves, do they do what you want in the sensitivity?
# Generate curve values using the SS3 math function
all_curves_age_df <- ranlpars %>%
  mutate(
    Curve_Data = pmap(
      list(peak_age, top_logit, ascend_se, descend_se, start_logit, end_logit), 
      function(p, tl, ase, dse, sl, el) {
        tibble(
          Age = ages_grid,
          Selectivity = ss3_double_normal_age(
            ages        = ages_grid,
            peak_age    = p,
            top_logit   = tl,
            ascend_se   = ase,
            descend_se  = dse,
            start_logit = sl,
            end_logit   = el,
            max_age     = max(ages_grid)
          )
        )
      }
    )
  ) %>%
  unnest(cols = c(Curve_Data))

# Generate the faceted structural diagnostic plot
ggplot(all_curves_age_df, aes(x = Age, y = Selectivity)) +
  geom_line(aes(group = scenario), color = "darkgreen", alpha = 0.04, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, max(ages_grid), by = 5)) +
  labs(
    title = "fixed top and end",
    subtitle = paste("Ensemble envelope across", length(unique(all_curves_age_df$scenario)), "age-structured scenarios"),
    x = "Age (Years)",
    y = "Selectivity Proportion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Get SPR outputs
# 1. Get the list of all unique scenarios to iterate over
unique_scenarios <- unique(all_curves_age_df$scenario)
total_scenarios  <- length(unique_scenarios)

# 2. Pre-allocate a data frame to store reference point outcomes cleanly
# (Using a data frame instead of a matrix handles strings and numbers safely)
results_collector <- tibble(
  scenario = unique_scenarios,
  F_40     = rep(NA_real_, total_scenarios),
  F_35     = rep(NA_real_, total_scenarios),
  status   = rep("Not Run", total_scenarios)
)

# 3. Execute the sensitivity loop
for(i in 1:total_scenarios) {
  
  current_scenario <- unique_scenarios[i]
  
  # A single, clean tryCatch to isolate the modeling process
  loop_out <- tryCatch({
    
    # Extract the pre-calculated age selectivity vector for this scenario
    sel_age_vector <- all_curves_age_df %>% 
      filter(scenario == current_scenario) %>% 
      arrange(Age) %>% 
      pull(Selectivity)
    
    # Build your reference point engine input
    inp_loop <- spr_input(
      ages = ages,
      species = list(
        AK_skate = list(
          len_at_age  = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), 
          wt_at_age   = list(type = "wl", alpha = 0.000009, beta = 2.9617), 
          maturity    = list(type = "vector", m = maa),
          M           = 0.13,
          selectivity = list(type = "vector", sel = sel_age_vector) 
        )
      ),
      use_plus_group = FALSE
    )
    
    # Calculate target reference points
    spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
    
    # Pass the calculated outputs forward to the collector assignment step
    # Note: Replace 'F40_total'/'F35_total' with the exact names your package uses
    list(
      F40 = spr_loop_out$F_spr_total["F40_total"], 
      F35 = spr_loop_out$F_spr_total["F35_total"], 
      status = "Success"
    )
    
  }, error = function(e) {
    # Catches crashes gracefully and prints a notification
    message(paste("Skipping", current_scenario, "due to evaluation error:", e$message))
    return(list(F40 = NA_real_, F35 = NA_real_, status = paste("Error:", e$message)))
  })
  
  # --- 3. Assign the results back to the pre-allocated data frame rows ---
  results_collector$F_40[i]   <- loop_out$F40
  results_collector$F_35[i]   <- loop_out$F35
  results_collector$status[i] <- loop_out$status
}

# 4. Bind the results directly back to your structural parameter grid
final_sensitivity_results <- ranlpars %>%
  left_join(results_collector, by = "scenario")

# make heatmap 
dome3 <- ggplot(final_sensitivity_results, aes(x = ascend_se, y = descend_se, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~start_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               start_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed top width and end logit",
    x = "ascend_se",
    y = "descend_se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

dome3_zoom_dat <- filter(final_sensitivity_results) %>% 
  filter(between(F_40, 0.09, 0.13))
dome3_zoom <- ggplot(dome3_zoom_dat, aes(x = ascend_se, y = descend_se, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~start_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               start_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed top width and ending value",
    x = "ascending se",
    y = "descending se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

# fixed top and start logit----
ranlpars <- expand.grid(
  peak_age   = seq(7, 11, by = 1),        # Shifting peak selectivity age (Years)
  ascend_se  = seq(1.0, 3.0, by = 0.5),   # Entry scale parameter in age-variance space
  descend_se = seq(1.0, 3.0, by = 0.5),   # Dome decay rate parameter in age-variance space
  end_logit = seq(-2, 2, by = 1.0)
) %>%
  mutate(
    top_logit   = -10, 
    start_logit  = 0,
    scenario    = paste0("Scenario_", row_number())
  )

# look at the curves, do they do what you want in the sensitivity?
# Generate curve values using the SS3 math function
all_curves_age_df <- ranlpars %>%
  mutate(
    Curve_Data = pmap(
      list(peak_age, top_logit, ascend_se, descend_se, start_logit, end_logit), 
      function(p, tl, ase, dse, sl, el) {
        tibble(
          Age = ages_grid,
          Selectivity = ss3_double_normal_age(
            ages        = ages_grid,
            peak_age    = p,
            top_logit   = tl,
            ascend_se   = ase,
            descend_se  = dse,
            start_logit = sl,
            end_logit   = el,
            max_age     = max(ages_grid)
          )
        )
      }
    )
  ) %>%
  unnest(cols = c(Curve_Data))

# Generate the faceted structural diagnostic plot
ggplot(all_curves_age_df, aes(x = Age, y = Selectivity)) +
  geom_line(aes(group = scenario), color = "darkgreen", alpha = 0.04, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(0, max(ages_grid), by = 5)) +
  labs(
    title = "fixed top and start",
    subtitle = paste("Ensemble envelope across", length(unique(all_curves_age_df$scenario)), "age-structured scenarios"),
    x = "Age (Years)",
    y = "Selectivity Proportion"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Get SPR outputs
# 1. Get the list of all unique scenarios to iterate over
unique_scenarios <- unique(all_curves_age_df$scenario)
total_scenarios  <- length(unique_scenarios)

# 2. Pre-allocate a data frame to store reference point outcomes cleanly
# (Using a data frame instead of a matrix handles strings and numbers safely)
results_collector <- tibble(
  scenario = unique_scenarios,
  F_40     = rep(NA_real_, total_scenarios),
  F_35     = rep(NA_real_, total_scenarios),
  status   = rep("Not Run", total_scenarios)
)

# 3. Execute the sensitivity loop
for(i in 1:total_scenarios) {
  
  current_scenario <- unique_scenarios[i]
  
  # A single, clean tryCatch to isolate the modeling process
  loop_out <- tryCatch({
    
    # Extract the pre-calculated age selectivity vector for this scenario
    sel_age_vector <- all_curves_age_df %>% 
      filter(scenario == current_scenario) %>% 
      arrange(Age) %>% 
      pull(Selectivity)
    
    # Build your reference point engine input
    inp_loop <- spr_input(
      ages = ages,
      species = list(
        AK_skate = list(
          len_at_age  = list(type = "vb", Linf = 135.39, k = 0.10, t0 = -1.60), 
          wt_at_age   = list(type = "wl", alpha = 0.000009, beta = 2.9617), 
          maturity    = list(type = "vector", m = maa),
          M           = 0.13,
          selectivity = list(type = "vector", sel = sel_age_vector) 
        )
      ),
      use_plus_group = FALSE
    )
    
    # Calculate target reference points
    spr_loop_out <- run_spr(inp_loop, diagnostics = TRUE)
    
    # Pass the calculated outputs forward to the collector assignment step
    # Note: Replace 'F40_total'/'F35_total' with the exact names your package uses
    list(
      F40 = spr_loop_out$F_spr_total["F40_total"], 
      F35 = spr_loop_out$F_spr_total["F35_total"], 
      status = "Success"
    )
    
  }, error = function(e) {
    # Catches crashes gracefully and prints a notification
    message(paste("Skipping", current_scenario, "due to evaluation error:", e$message))
    return(list(F40 = NA_real_, F35 = NA_real_, status = paste("Error:", e$message)))
  })
  
  # --- 3. Assign the results back to the pre-allocated data frame rows ---
  results_collector$F_40[i]   <- loop_out$F40
  results_collector$F_35[i]   <- loop_out$F35
  results_collector$status[i] <- loop_out$status
}

# 4. Bind the results directly back to your structural parameter grid
final_sensitivity_results <- ranlpars %>%
  left_join(results_collector, by = "scenario")

# make heatmap 
dome4 <- ggplot(final_sensitivity_results, aes(x = ascend_se, y = descend_se, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~end_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               end_logit = label_both    # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed top width and start logit",
    x = "ascend_se",
    y = "descend_se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

dome4_zoom_dat <- filter(final_sensitivity_results) %>% 
  filter(between(F_40, 0.09, 0.13))
dome4_zoom <- ggplot(dome4_zoom_dat, aes(x = ascend_se, y = descend_se, fill = F_40)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "F40") +
  facet_grid(peak_age~end_logit,
             labeller = labeller(
               peak_age  = label_both,   # Displays as "peak_age: 6", "peak_age: 8", etc.
               end_logit = label_both   # Displays as "end_logit: -5", "end_logit: -3", etc.
             ))+
  labs(
    title = "Dome sensitivity",
    subtitle = "Fixed top width and start value",
    x = "ascending se",
    y = "descending se"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )


# dome notes----
# based on those 4 runs, Peak, ascending_se and start_logit appear to have the most impact
# does that mean dome doesn't matter and we can just use logistic to be simpler?
# how do we determine the optimal logistic curve?
# flip the question, what kind of selectivity curves would it take to get F rates similar to Tier 3?




# Summary Plots----
M_plot / lslx_plot
