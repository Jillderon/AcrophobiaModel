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
  
  # Check if the punishment are correct given the size of grid: 
  if (any(!na.omit(c(punishment_locs, reward_locs)) %in% 1:(size*size))) {
    stop("The locations of the punishments and rewards are incorrect given the size of the grid.")
  }
  
  # Identify coordinates of viable states:
  viable_coordinates = data.frame(x = rep(1:size), y = rep(1:size, each = size))
  
  # Make transition matrix with NA's and 1's:
  transition_matrix <- proxy::dist(viable_coordinates, viable_coordinates, method = "Euclidean")
  if (standing_still) {
    transition_matrix <- ifelse(transition_matrix == 1 | transition_matrix == 0, 1, NA) 
  } else {
    transition_matrix <- ifelse(transition_matrix == 1, 1, NA) 
  }
  
  # Identify terminal points
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
    
    for (i in 1:length(s_prime)) {
      info[index,1] = s; info[index,2] <- s_prime[i]; info[index,3] <- r[i]
      index = index + 1
    }
  }

  return(list(MDP = info, terminal_states = terminal))
}

Bellman <- function(arr, control) { 
  return(control * max(arr) + (1 - control) * min(arr))
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

solve_Qvalues <- function(MDP_info, control, max_iter, gamma, tol, size, state_index = NULL, init_Q = NULL) {
  
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

determine_grid <- function(size, grid_space, reward, punishment, plot) {
  
  # Determine grid parameters --------
  length <- size*size
  if (size == 22) {
    
    # Determine starting state:
    start <- 242
    
    # Reward Location 
    r_locs <- size * (size/2 -1 ) + 3  
    
    # Punishment Locations
    edge <- c(1:size,(length-size+1):length,seq(1,(length-size+1),size)) # creates a border of threats
    plus <- c(247:258, 269:280, 291:302, 313:324, 335:346, 357:368, 379:390, 401:412, 423:434, 445:456) # ugly hard-coded, change!
    p_locs <- c(edge,plus) 
    
    if (grid_space == "abyss2") {
      edge <- c((length-size+1):length,seq(1,(length-size+1),size))
      p_locs <- c(edge,plus) 
    }
    
    if (grid_space == "ground") {
      p_locs <- c()
    }
    
    if (grid_space == "bridge") {
      
      # add more negative states:
      p_locs <- c(p_locs, 26:44, 48:66, 70:88, 92:110, 114:132, 136:154, 
                           158:176, 180:198, 202:220, 246:264, 268:286, 290:308, 
                           312:330, 334:352, 356:374, 378:396, 400:418, 422:440, 444:462)
      
    }
  } 
  
  if (size == 11) {
    
    # Determine starting state
    grid <- matrix(seq(1:(size*size)), ncol = size)
    start <- ceiling(length(grid)*.5)+floor(size/2)
    
    # Reward Location
    r_locs <- ceiling(length(grid)*.5)-floor(size/2)+2
    
    # Punishment Location
    edge <- c(1:size,(length(grid)-size+1):length(grid),seq(1,(length(grid)-size+1),size)) # creates a border of threats
    plus <- c((size*floor(size/3)+floor(size/2)):(size*floor(size/3)+floor(size/2)+2)) # specifies additional punishment spots
    plus <- c(104:106,93:95,82:84,71:73) # specifies additional punishment spots
    p_locs <- c(edge,plus) # specifies punishments
    
    
    if (grid_space == "abyss2") {
      edge <- c((length(grid)-size+1):length(grid),seq(1,(length(grid)-size+1),size))
      p_locs <- c(edge,plus)
      r_locs<-80
    }
    
    
    if (grid_space == "ground") {
      start <- 64 
      r_locs <- 58
      p_locs <- c()
      
    }
    
    
    if (grid_space == "bridge") {
      
      start <- 64 
      r_locs <- 58
      
      # add more negative states:
      p_locs <- c(p_locs, 13:33, 35:42, 46:53, 68:75, 79:86, 90:99, 101:110)
    }
    
    if (grid_space == "bridge2") {
      start <- 64
      r_locs <- 57
      p_locs <- c(1:33, 34:42, 48:53, 70:75, 78:86, 89:121)
    }
    
    if (grid_space == "bridge3") {
      start <- 64
      r_locs <- 57
      p_locs <- c(1:33, 34:42, 78:86, 89:121)
    }
    
    if (grid_space == "bridge4") {
      start <- 66
      r_locs <- 56
      p_locs <- c(3:6,14:17,25:28, 36:39,47:50,69:72, 80:83, 91:94,102:105,113:116)
    }
  }
  
  # Determine grid parameters --------
  if (grid_space == "cliff") {
    start <- size 
    r_locs <- length
    p_locs <- c(size*2, size*3, size*4, size*5, size*6, size*7, size*8, size*9, size*10)
    
    if (size == 22) {
      p_locs <- c(p_locs, size*11, size*12, size*13, size*14, size*15, size*16, 
                          size*17, size*18, size* 19, size*20, size*21)
    }
     
  }

  
  if (plot) {
    # Plot grid: 
    grid <- rep(NA, size*size); grid[r_locs] <- reward; grid[p_locs] <- punishment
    
    # Visualize environment:
    pdf(paste0(grid_space,".pdf"))
    pheatmap::pheatmap(matrix(grid, nrow = size, ncol = size), legend = FALSE, cluster_col = FALSE, cluster_row = FALSE,
                       display_numbers = FALSE, color = c("gray", "white"), 
                       na_col = "white", width = 3, height = 3)
    dev.off()
  }
  
  return(list(punishment_locs = p_locs, reward_locs = r_locs, start = start))
}

solve_one_step <- function(MDP, Qvalues, weight, terminal_states, current_state) {
  
  

  
}
