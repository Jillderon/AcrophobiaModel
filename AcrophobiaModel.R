###############################################################################
#                             Acrophobia Model                                #
###############################################################################

# Load Libraries -----
library(pheatmap)
library(proxy)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggplotify)
library(tidyverse)
library(smacpod)

# Load in supporting functions -----
source("SupportingFunctionsRL.R")

# Run simPhobia -----
simPhobia <- function(control,
                      gamma = .95, 
                      reward = 1, 
                      meters = 10,
                      size = c(11, 22), 
                      grid_space = c("abyss", "abyss2", 
                                     "bridge", "bridge2","bridge3", "bridge4", 
                                     "cliff", 
                                     "ground"), 
                      rate_somatic = 0.5,
                      exponent_somatic = .25,
                      exponent_reward = 1,
                      alpha = 10,
                      beta = .5,
                      punishment_h = .2,
                      max_iter = 25, # for debugging purposes with size 11, normally I would do 100
                      tol = 0.0001, 
                      plot = FALSE) {
  
  punishment <- -1 / (1 + (exp(1)^(-5 * (log10(meters) - 2 * punishment_h))))
  grid <- determine_grid(size, grid_space, reward, punishment, plot = FALSE)
  
  # Define MPD environment: 
  MDP_info <- define_MDP(size = size, 
                         punishment = punishment,
                         punishment_locs = grid$punishment_locs,
                         reward = reward,
                         reward_locs = grid$reward_locs)
  
  # Define MPD environment with only punishment: 
  MDP_info_p <- define_MDP(size = size, 
                           punishment = punishment,
                           punishment_locs = grid$punishment_locs,
                           reward = NA,
                           reward_locs = NA)
  
  # Initialize state variables: 
  threat <- somatic <- fear <- behavior <- c <-c()
  behavior[1] <- grid$start
  threat[1] <- somatic[1] <- fear[1] <- 0
  c[1] <- control
  
  # Make a list for the outcomes of the MDP
  out <- list()
  out_p <- list()
  grids_v <- list()
  
  for (i in 1:max_iter) {
    
    print(i) 
    
    # Solve MDP via weighted Bellmans equation: 
    out[[i]] <- solve_MDP(MDP = MDP_info$MDP, 
                          control = c[i], 
                          gamma = gamma, 
                          terminal_states = MDP_info$terminal_states, 
                          start_state = behavior[i], 
                          max_iter = max_iter, 
                          tol = tol, 
                          size = size)
  
    # Solve MDP with only punishment via weighted Bellmans equation: 
    out_p[[i]] <- solve_MDP(MDP = MDP_info_p$MDP,
                            control = c[i],
                            gamma = gamma,
                            terminal_states = MDP_info_p$terminal_states,
                            start_state = behavior[i],
                            max_iter = max_iter,
                            tol = tol,
                            size = size)

    behavior[i+1]    <- out[[i]]$policy[2]
    threat[i+1]      <- abs(out_p[[i]]$V_values[behavior[i+1]]) 
    somatic[i+1]     <- somatic[i] + rate_somatic*(threat[i]^exponent_somatic - somatic[i])
    fear[i+1]        <- sqrt(somatic[i+1]*threat[i+1])
    c[i+1]           <- c[i] + (control / (1 + (exp(1)^(10 * (somatic[i+1] - .5)))) - c[i]) # CHANGE TO SOMATIC AGAIN

    # Plotting
    if (plot) {
      
      # Place the value function in a grid format:
      grid_v <- matrix(out[[i]]$V_values, ncol = size)
      grid_v[grid$punishment_locs] <- punishment
      grid_v[grid$reward_locs] <- reward
      grid_v[behavior] <- NA

      # Plot the value function and policy in a heatmap:
      ggplotify::as.grob(pheatmap::pheatmap(grid_v, legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
                                            display_numbers = FALSE,
                                            color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(99),
                                            breaks=c(seq(-1,0,length.out=50)[1:49],0,seq(0,1,length.out=50)[2:50]),
                                            na_col = "black",  width = 3, height = 3))


      grids_v[[i]] <- grid_v
    }
    
    # Check if the agent is at a terminal state: 
    if (behavior[i+1] %in% MDP_info$terminal_states) { 
      break 
    } 
  }
  
  simulation_name <- paste0("PCwSomatic_W", control,
                            "_M", meters,
                            "_h", punishment_h,
                            "_r", reward,
                            ".pdf")
  
  if (plot) {
  
    plot1 <- pheatmap::pheatmap(grids_v[[i]], legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
                                display_numbers = FALSE, 
                                color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
                                na_col = "black",  width = 3, height = 3)
    
    # Add a plot of the arousal/fear response 
    affect <- data.frame(somatic = somatic[1:length(somatic)-1], # do we need a value at step zero?
                         fear = fear[1:length(fear)-1],
                         control_steps = c[1:length(c)-1],
                         x = 1:(length(somatic)-1))
    
    plot2 <- ggplot(affect, aes(x = x, y = somatic)) + 
                    geom_point(size = 4) +
                    geom_point(aes(x = x, y = control_steps), col = "blue", size=3, shape = 17)+
                    geom_point(aes(x = x, y = fear), col = "orange", size=3, shape = 15)+
                    ylim(0,1) +
                    scale_x_discrete(name = "Sequence", limits = factor(1:c(length(somatic)-1))) +
                    labs(title = "Somatic Response and Perceived Control", y = "")+
                    theme_minimal() 
    
    pdf(simulation_name, width = 10, height = 5)
    grid.arrange(plot1[[4]], plot2, ncol = 2)
    dev.off()
    
  }
  
  return(list(simulation_name = simulation_name, 
              grid_reward_locs = grid$reward_locs,
              grid_punishment_locs = grid$punishment_locs,
              V_values = out[[i]]$V_values, 
              fear = fear, 
              behavior = behavior, 
              somatic = somatic, 
              threat = threat,
              grids_v = grids_v,
              c = c))
  
}

####################
# Running Simulation
###################

# Define simulation conditions
conditions <- expand.grid(PC     = c(1, 0.8, 0.6, 0.4, 0.2), 
                          grid   = "abyss2", #c("abyss", "bridge", "cliff")
                          meter = c(0, 0.5, 5, 10, 20, 40, 70, 100), 
                          h     = c(0.2, 0.5))
# Save in data frame: 
data <- data.frame()

# Loop over conditions and run simulation:
for (i in 1:nrow(conditions)) {
  result <- simPhobia(control = conditions$PC[i],
                      punishment_h = conditions$h[i],
                      gamma = .95,
                      size = 11,
                      meters = conditions$meter[i],
                      reward = 1,
                      grid_space = conditions$grid[i],
                      plot = TRUE)
  
  row <- data.frame(result$simulation_name, 
                    paste0(result$behavior, collapse = " "), # path
                    length(result$behavior), # number of steps
                    paste0(result$c, collapse = " "), # perceived control 
                    paste0(result$V_values[result$behavior], collapse = " "), # speed (vector with v-values of the path)
                    paste0(result$somatic, collapse = " "), # somatic sensation 
                    paste0(result$threat, collapse = " "), # threat 
                    paste0(result$fear, collapse = " ")) # fear 
  
  write.table(row, file = "output.csv", append = TRUE, col.names = FALSE)
}



