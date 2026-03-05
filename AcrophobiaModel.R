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
library(colorspace)

# Load baseline supporting functions (grid + MDP constructors) -----
source("SupportingFunctionsRL.R")

# Run simPhobia -----
simPhobia <- function(control,
                      gamma = .95,
                      reward = 1,
                      meters = 10,
                      size = c(11, 22),
                      grid_space = c("cove", "cove2",
                                     "bridge", "bridge2",
                                     "holes", "cliff", "ground"),
                      rate_somatic = 0.5,
                      exponent_somatic = .25,
                      exponent_reward = 1,
                      alpha = 10,
                      beta = .5,
                      punishment_h = .2,
                      max_iter = 25,
                      tol = 0.0001,
                      pdf_trajectory = FALSE,
                      pdf_grid = FALSE) {

  # Convert objective elevation into the perceived punishment level:
  punishment <- -1 / (1 + exp(-5 * (log10(meters) - 2 * punishment_h)))

  # Build the environment (reward/punishment locations and start state): 
  grid <- determine_grid(size, grid_space, reward, punishment, plot = pdf_grid)

  # MDP with reward and punishment (to determine behavior): 
  MDP_info <- define_MDP(size = size, 
                         punishment = punishment, punishment_locs = grid$punishment_locs,
                         reward = reward, reward_locs = grid$reward_locs)

  # MDP with punishment only (to estimate threat):
  MDP_info_p <- define_MDP(size = size,
                           punishment = punishment, punishment_locs = grid$punishment_locs,
                           reward = NA, reward_locs = NA)

  # Define state variables: 
  threat <- numeric(max_iter + 1)
  somatic <- numeric(max_iter + 1)
  fear <- numeric(max_iter + 1)
  behavior <- numeric(max_iter + 1)
  control_trace <- numeric(max_iter + 1)
  
  # Initialize state variables: 
  behavior[1] <- grid$start
  control_trace[1] <- control

  out <- list()
  out_p <- list()
  grids_v <- list()

  state_index <- build_state_index(MDP_info$MDP$S, size * size)
  state_index_p <- build_state_index(MDP_info_p$MDP$S, size * size)
  q_prev <- NULL
  q_prev_p <- NULL

  for (i in seq_len(max_iter)) {
    out[[i]] <- solve_MDP(MDP_info = MDP_info$MDP,
                          control = control_trace[i],
                          gamma = gamma,
                          terminal_states = MDP_info$terminal_states,
                          start_state = behavior[i],
                          max_iter = max_iter,
                          tol = tol,
                          size = size,
                          init_Q = q_prev,
                          state_index = state_index)

    out_p[[i]] <- solve_MDP(MDP_info = MDP_info_p$MDP,
                            control = control_trace[i],
                            gamma = gamma,
                            terminal_states = MDP_info_p$terminal_states,
                            start_state = behavior[i],
                            max_iter = max_iter,
                            tol = tol,
                            size = size,
                            init_Q = q_prev_p,
                            state_index = state_index_p)

    q_prev <- out[[i]]$Q_values
    q_prev_p <- out_p[[i]]$Q_values

    behavior[i+1] <- out[[i]]$policy[2]
    threat[i+1]   <- abs(out_p[[i]]$V_values[behavior[i + 1]])
    somatic[i+1]  <- somatic[i] + rate_somatic * (threat[i]^exponent_somatic - somatic[i])
    
    # Fear is modeled as the geometric mean of somatic arousal and threat:
    fear[i+1] <- sqrt(somatic[i + 1] * threat[i + 1])
    
    # Perceived control decreases as somatic arousal increases:
    control_trace[i+1] <- control_trace[i] + (control / (1 + exp(10 * (somatic[i + 1] - .5))) - control_trace[i])

    if (pdf_trajectory) {
      grid_v <- matrix(out[[i]]$V_values, ncol = size)
      grid_v[grid$punishment_locs] <- punishment
      grid_v[grid$reward_locs] <- reward
      grid_v[behavior] <- NA

      ggplotify::as.grob(pheatmap::pheatmap(
        grid_v,
        legend = FALSE,
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        display_numbers = FALSE,
        color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))(99),
        breaks = c(seq(-1, 0, length.out = 50)[1:49], 0, seq(0, 1, length.out = 50)[2:50]),
        na_col = "black",
        width = 3,
        height = 3
      ))

      grids_v[[i]] <- grid_v
    }

    if (behavior[i + 1] %in% MDP_info$terminal_states) {
      break
    }
  }

  simulation_name <- paste0(
    "PC", control,
    "_M", meters,
    "_h", punishment_h,
    "_r", reward,
    ".pdf"
  )

  if (pdf_trajectory) {
    plot1 <- pheatmap::pheatmap(
      grids_v[[i]],
      legend = FALSE,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      display_numbers = FALSE,
      color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
      na_col = "black",
      width = 3,
      height = 3
    )

    n_points <- length(somatic) - 1
    affect <- data.frame(
      somatic = somatic[seq_len(n_points)],
      fear = fear[seq_len(n_points)],
      control_steps = control_trace[seq_len(n_points)],
      x = seq_len(n_points)
    )

    plot2 <- ggplot(affect, aes(x = x, y = somatic)) +
      geom_point(size = 4) +
      geom_point(aes(x = x, y = control_steps), col = "blue", size = 3, shape = 17) +
      geom_point(aes(x = x, y = fear), col = "orange", size = 3, shape = 15) +
      ylim(0, 1) +
      scale_x_discrete(name = "Sequence", limits = factor(1:c(length(somatic)-1))) +
      labs(title = "Somatic Response and Perceived Control", y = "") +
      theme_minimal()

    pdf(simulation_name, width = 10, height = 5)
    grid.arrange(plot1[[4]], plot2, ncol = 2)
    dev.off()
  }

  list(
    simulation_name = simulation_name,
    grid_reward_locs = grid$reward_locs,
    grid_punishment_locs = grid$punishment_locs,
    V_values = out[[i]]$V_values,
    fear = fear,
    behavior = behavior,
    somatic = somatic,
    threat = threat,
    grids_v = grids_v,
    c = control_trace,
    q_iterations = vapply(out, function(x) x$q_iterations, numeric(1)),
    q_iterations_p = vapply(out_p, function(x) x$q_iterations, numeric(1))
  )
}

########
# Running the output
########

simPhobia(control = .2, 
          size = 11, 
          grid_space = "bridge2", 
          pdf_trajectory = FALSE,
          pdf_grid = TRUE)


c("cove", "cove2",
  "bridge", "bridge2",
  "holes", "cliff", "ground")



# Define simulation conditions
conditions <- expand.grid(PC = c(1, 0.8, 0.6, 0.4, 0.2),
                          grid = "cove2", # c("cove", "bridge", "cliff")
                          meter = c(0, 0.5, 5, 10, 20, 40, 70, 100),
                          h = c(0.2, 0.5, 0.7))

output_file <- "output.csv"
if (file.exists(output_file)) {
  file.remove(output_file)
}

# Loop over conditions and write summary output:
for (i in 1:nrow(conditions)) {
  result <- simPhobia(control = conditions$PC[i],
                          punishment_h = conditions$h[i],
                          gamma = .95,
                          size = 11,
                          meters = conditions$meter[i],
                          reward = 1,
                          grid_space = conditions$grid[i],
                          pdf_trajectory = FALSE,
                          pdf_grid = FALSE)

  # One row per condition summarizing end-state and averaged affective measures.
  row <- data.frame(
    Control = conditions$PC[i],
    Punishment_h = conditions$h[i],
    Meter = conditions$meter[i],
    Gridspace = conditions$grid[i],
    Goal = tail(result$behavior, 1) %in% result$grid_reward_locs,
    Steps = length(result$behavior),
    MeanFear = mean(result$fear, na.rm = TRUE),
    MeanBodily = mean(result$somatic, na.rm = TRUE),
    MeanThreat = mean(result$threat, na.rm = TRUE),
    MeanDistAbyss = mean_distance_to_abyss(result$behavior, result$grid_punishment_locs, size = 11)
  )

  write.table(row, file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)
}
