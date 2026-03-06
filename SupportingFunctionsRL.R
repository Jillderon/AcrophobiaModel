###############################################################################
# 
# Supporting function for the Acrophobia Model
#
# Note: Significant portions of this code were adapted from Zorowitz et al. (2020)
# See https://github.com/ndawlab/seqanx for the original code for Zorowitz model
# 
###############################################################################

define_MDP <- function(size, 
                       punishment,
                       punishment_locs,
                       reward,
                       reward_locs,
                       standing_still = TRUE,
                       terminal = NA) {
  
  # Validate indices against the flattened grid indexing scheme.
  if (any(!na.omit(c(punishment_locs, reward_locs)) %in% 1:(size*size))) {
    stop("The locations of the punishments and rewards are incorrect given the size of the grid.")
  }
  
  # Each state is one cell in a size x size grid.
  viable_coordinates = data.frame(x = rep(1:size), y = rep(1:size, each = size))
  
  # Transition matrix marks legal moves (distance 1) and optionally staying put.
  transition_matrix <- proxy::dist(viable_coordinates, viable_coordinates, method = "Euclidean")
  if (standing_still) {
    transition_matrix <- ifelse(transition_matrix == 1 | transition_matrix == 0, 1, NA) 
  } else {
    transition_matrix <- ifelse(transition_matrix == 1, 1, NA) 
  }
  
  # By default, treat reward and punishment states as terminal.
  if (any(is.na(terminal))) terminal <- na.omit(c(reward_locs, punishment_locs))
  for (t in terminal) {
    transition_matrix[t,] <- NA
    transition_matrix[t,t] <- 1
  }
  
  # Define rewards/punishments (R) that can be obtained FROM a given state (s):
  R <- matrix(0, nrow = size*size, ncol = size*size)
  for (r in reward_locs){
    R[,r] <- reward; R[r,r] <- 0 # could be a problem when standing still is TRUE
  }
  for (p in punishment_locs){
    R[,p] <- punishment; R[p,p] <- 0 # could be a problem when standing still is TRUE
  }
  R <- ifelse(transition_matrix == 1, R, NA)
  
  # Iteratively define MDP information:
  info <- data.frame(S = 0, S_prime = 0, R = 0)
  index = 1
  n_states <- size*size
  for (s in 1:n_states) {
    
    # What are the possible next states?
    s_prime = which(transition_matrix[s,] %in% 1)  # check on which position are the ones
    
    # What are the rewards for these next states?
    r = R[s, s_prime]
    
    for (i in seq_along(s_prime)) {
      info[index,1] = s; info[index,2] <- s_prime[i]; info[index,3] <- r[i]
      index = index + 1
    }
  }

  list(MDP = info, terminal_states = terminal)
}

Bellman <- function(arr, control) { 
  # Weighted blend between optimistic (max) and pessimistic (min) valuation.
  control * max(arr) + (1 - control) * min(arr)
}

# Build row indices for each source state once, then reuse in value/policy loops.
build_state_index <- function(source_states, n_states) {
  split_idx <- split(seq_along(source_states), source_states)
  state_index <- vector("list", n_states)
  for (name in names(split_idx)) {
    state_index[[as.integer(name)]] <- split_idx[[name]]
  }
  state_index
}

solve_Qvalues <- function(MDP_info, control, max_iter, gamma, tol, size,
                          state_index = NULL, init_Q = NULL) {
  
  n_transitions <- nrow(MDP_info)

  # Warm-start from previous Q-values when available; otherwise zeros.
  if (!is.null(init_Q) && length(init_Q) == n_transitions) {
    Q <- as.numeric(init_Q)
  } else {
    Q <- numeric(n_transitions)
  }

  rewards <- MDP_info$R
  successor_states <- MDP_info$S_prime
  n_states <- size * size
  if (is.null(state_index)) {
    state_index <- build_state_index(MDP_info$S, n_states)
  }
  
  iter_used <- 0

  # Solve for Q-values:
  for (iter in 1:max_iter) {
    q_prev <- Q
    
    # Precompute successor value: 
    V_prime <- numeric(n_states)
    for (s in seq_len(n_states)) {
      idx <- state_index[[s]]
      if (length(idx) == 0) {
        next
      }
      arr <- q_prev[idx]
      V_prime[s] <- Bellman(arr, control) # V_prime is the value function for every possible state given the policy
    }
    
    # Compute Q-values:
    Q <- rewards + gamma * V_prime[successor_states]
    iter_used <- iter
    
    # Check for termination (i.e., all values of delta below threshold):
    if (max(abs(Q - q_prev)) < tol) {
      break
    }
  }
  
  return(list(Q = Q, iterations = iter_used))
}

solve_MDP <- function(MDP_info, control, terminal_states, start_state, max_iter, gamma, tol, size, init_Q = NULL, state_index = NULL) {
  
  n_states <- size * size
  if (is.null(state_index)) {
    state_index <- build_state_index(MDP_info$S, n_states)
  }
  q_out <- solve_Qvalues(MDP_info, control, max_iter, gamma, tol, size, state_index, init_Q)
  q_values <- q_out$Q
  
  # Identify max by state:
  V <- numeric(n_states)
  for (s in seq_len(n_states)) {
    idx <- state_index[[s]]
    if (length(idx) == 0) {
      next
    }
    arr <- q_values[idx]
    V[s] <- Bellman(arr, control) # Before this was max(arr$Q), but this leaves a "cross" around the reward
  }
  
  # Compute policy 
  policy <- c(start_state)
  state <- start_state
  i = 1
  is_terminal <- logical(n_states)
  is_terminal[terminal_states] <- TRUE
  while (!is_terminal[state]) {
    
    # Compute optimal q(s, a):
    idx <- state_index[[state]]
    state <- MDP_info$S_prime[idx][which.max(q_values[idx])]
    
    # Append it to policy
    policy <- c(policy, state)
    
    # Make sure it doesn't get stuck: 
    i = i + 1
    if (i > 100 || is_terminal[state]) {
      break
    }
  }
  
  return(list(policy   = policy, 
              V_values = V, 
              Q_values = q_values,
              q_iterations = q_out$iterations))
}

# Mean Manhattan distance from the visited path to nearest punishment state.
mean_distance_to_abyss <- function(path_states, punishment_locs, size) {
  if (length(punishment_locs) == 0 || length(path_states) == 0) {
    return(NA_real_)
  }
  x_path <- ((path_states - 1) %% size) + 1
  y_path <- ((path_states - 1) %/% size) + 1
  x_pun <- ((punishment_locs - 1) %% size) + 1
  y_pun <- ((punishment_locs - 1) %/% size) + 1

  min_dist <- vapply(seq_along(path_states), function(i) {
    min(abs(x_path[i] - x_pun) + abs(y_path[i] - y_pun))
  }, numeric(1))

  mean(min_dist)
}

determine_grid <- function(size,
                           grid_space,
                           reward,
                           punishment,
                           plot) {
  size <- as.integer(size[1])
  grid_space <- as.character(grid_space[1])
  valid_spaces <- c("cove", "cove2", "cove3", "abyss2", "bridge", "bridge2", "holes", "cliff", "ground")

  if (is.na(size) || size <= 10 || size >= 50) {
    stop("`size` must be an integer between 11 and 49.")
  }
  if (!grid_space %in% valid_spaces) {
    stop(paste0("Unknown grid_space: ", grid_space))
  }

  total_states <- size * size
  all_states <- seq_len(total_states)

  # Convert (row, col) coordinates to flattened state index.
  idx <- function(row, col) row + (col - 1) * size
  idx_rect <- function(rows, cols) as.vector(outer(rows, cols, Vectorize(idx)))

  center_col <- ceiling(size / 2)
  start <- idx(size - 2, center_col)
  r_locs <- idx(3, center_col)

  edge_top <- idx(1, seq_len(size))
  edge_bottom <- idx(size, seq_len(size))
  edge_left <- idx(seq_len(size), 1)
  edge_right <- idx(seq_len(size), size)

  # Shared "abyss-like" interior threat block.
  plus_rows <- seq(max(2, floor(size * 0.25)), min(size - 1, ceiling(size * 0.7)))
  plus_cols <- seq(max(2, ceiling(size * 0.55)), size - 1)
  plus <- idx_rect(plus_rows, plus_cols)

  # Core bridge corridors (safe sets used to define punishment complement).
  corridor_half_width <- max(0, floor(size * 0.06))
  bridge_cols <- seq(max(2, center_col - corridor_half_width), min(size - 1, center_col + corridor_half_width))
  bridge_safe <- idx_rect(2:(size - 1), bridge_cols)

  # Base zig-zag safe path (used for bridge2).
  path_rows <- 2:(size - 1)
  max_shift <- max(1, floor(size * 0.1))
  zigzag_shift <- round(sin(seq(0, pi, length.out = length(path_rows))) * max_shift)
  zigzag_cols <- pmin(size - 1, pmax(2, center_col + zigzag_shift))
  bridge2_safe <- unique(c(mapply(idx, path_rows, zigzag_cols),
                           mapply(idx, path_rows, pmin(size - 1, zigzag_cols + 1))))

  if (grid_space == "cove") {
    # Keep original start/reward positions, but move the cove threat block left
    # so the straight center path passes punishment on its right side.
    cove_cols <- seq(min(size - 1, center_col + 1), size - 1)
    cove_rows <- 3:(size - 2)
    cove_plus <- idx_rect(cove_rows, cove_cols)
    p_locs <- c(edge_top, edge_bottom, edge_left, cove_plus)
  } else if (grid_space == "cove2" || grid_space == "cove3") {
    # Scaled cove templates:
    # cove3: top+bottom+left+right borders + mid-right inward block
    # cove2: same as cove3, but without left and bottom border punishments
    if (grid_space == "cove2") {
      start <- idx(size, center_col)
    } else {
      start <- idx(size - 1, center_col)
    }
    goal_col <- min(size - 2, center_col + max(1, floor(size * 0.18)))
    r_locs <- idx(3, goal_col)

    if (size == 11) {
      # Match the reference figure exactly for the 11x11 case.
      block_rows <- 5:7
      block_cols <- 7:10
    } else {
      block_rows <- seq(max(4, floor(size * 0.45)), min(size - 3, ceiling(size * 0.65)))
      block_cols <- seq(min(size - 1, center_col + 1), size - 1)
    }
    mid_right_block <- idx_rect(block_rows, block_cols)
    if (grid_space == "cove3") {
      p_locs <- c(edge_top, edge_bottom, edge_left, edge_right, mid_right_block)
    } else {
      p_locs <- c(edge_top, edge_right, mid_right_block)
    }
  } else if (grid_space == "abyss2") {
    # Preserve legacy abyss2 coordinates for the original 11x11 setup.
    if (size == 11) {
      start <- 66
      r_locs <- 80
      plus_legacy <- c(104:106, 93:95, 82:84, 71:73)
      p_locs <- c(edge_right, edge_top, plus_legacy)
    } else {
      # Scaled fallback for other sizes.
      start <- idx(size - 2, center_col)
      r_locs <- idx(3, min(size - 2, center_col + max(1, floor(size * 0.15))))
      abyss2_rows <- seq(max(2, floor(size * 0.25)), min(size - 1, ceiling(size * 0.7)))
      abyss2_cols <- seq(max(2, ceiling(size * 0.55)), size - 1)
      abyss2_block <- idx_rect(abyss2_rows, abyss2_cols)
      p_locs <- c(edge_right, edge_top, abyss2_block)
    }
  } else if (grid_space == "ground") {
    p_locs <- integer(0)
  } else if (grid_space == "bridge") {
    p_locs <- setdiff(all_states, c(bridge_safe, start, r_locs))
  } else if (grid_space == "bridge2") {
    r_locs <- idx(3, max(2, center_col - 1))
    p_locs <- setdiff(all_states, c(bridge2_safe, start, r_locs))
  } else if (grid_space == "holes") {
    start <- idx(size - 1, min(size - 1, center_col + 1))
    r_locs <- idx(2, center_col)
    stripe_rows <- seq(2, size - 1, by = 2)
    stripe_cols <- seq(3, size - 2, by = 3)
    p_locs <- idx_rect(stripe_rows, stripe_cols)
  } else if (grid_space == "cliff") {
    start <- idx(size, 1)
    r_locs <- idx(size, size)
    p_locs <- idx(size, 2:(size - 1))
  }

  p_locs <- sort(unique(setdiff(p_locs, c(start, r_locs))))

  if (plot) {
    # Plot grid as four classes so each tile type has a stable color.
    # -1 = punishment, 0 = neutral, 1 = reward, 2 = start
    grid_plot <- rep(0, total_states)
    grid_plot[p_locs] <- -1
    grid_plot[r_locs] <- 1
    grid_plot[start] <- 2
    grid_matrix <- matrix(grid_plot, nrow = size, ncol = size)

    pheatmap::pheatmap(grid_matrix,
                       legend = FALSE,
                       cluster_col = FALSE,
                       cluster_row = FALSE,
                       display_numbers = FALSE,
                       color = c("gray70", "white", "dodgerblue3", "forestgreen"),
                       breaks = c(-1.5, -0.5, 0.5, 1.5, 2.5),
                       na_col = "white",
                       width = 3,
                       height = 3,
                       filename = paste0("Grid_", grid_space, "_size", size, ".pdf"))

  }

  list(punishment_locs = p_locs, reward_locs = r_locs, start = start)
}
