s <- 1 / (1 + exp((ages - a50) / delta))

ra50 <- round(runif(10000, 5, 15), 0)
rdelt <- round(runif(10000, -2.5, -1), 1)

rvals <- ra50 %>% 
  bind_cols(rdelt)

names(rvals) <- c("ra50", 'rdelt')

#need to make a loop to run the logistic curve for all ages
ages <- 5

rvals <- rvals %>% 
  mutate(s = 1 / (1 + exp((ages - ra50) / rdelt)),
         ratios = ra50/rdelt)


theme_tier4 <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.position = "right"#,
      # Increase legend key size - opted to use linetype scale instead
      # legend.key.size = grid::unit(1.6, "lines"),
    )
}


ggplot(rvals, aes(x = as.factor(ra50), y = as.factor(rdelt), fill = s))+
  geom_tile()+
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_discrete(expand = expansion(mult = 0, add = 0)) +
  theme_tier4()

ggplot(rvals, aes(x = as.factor(ra50), y = as.factor(rdelt), fill = ratios))+
  geom_tile()+
  scale_x_discrete(expand = expansion(mult = 0, add = 0)) +
  scale_y_discrete(expand = expansion(mult = 0, add = 0)) +
  theme_tier4()
