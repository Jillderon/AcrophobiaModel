################################################################################
#                      Computational Model of Acrophobia                       #
#                   Visualizing different movement patterns                    #
################################################################################ 

# Load functions: 
source("AcrophobiaModel.R")

# =========================
# Fleeing and Freezing
# =========================

fleeing <- simPhobia(control = 0.4,
                     punishment_h = 0.5,
                     gamma = .95,
                     size = 11,
                     meters = 20,
                     reward = 1,
                     grid_space = "abyss2",
                     plot = FALSE)

# Set up standard grid (also used for other figures)
grid <- matrix(0, nrow = 11, ncol = 11)
grid[fleeing$grid_punishment_locs] <- -1
grid[fleeing$grid_reward_locs] <- 1

# Set up heatmap: 
grid_fleeing <- grid; grid_fleeing[fleeing$behavior] <- NA
plot1_fleeing <- pheatmap::pheatmap(grid_fleeing, legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
                            display_numbers = FALSE, 
                            color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
                            na_col = "black",  width = 3, height = 3)

# Set up figure 2: 
data_fleeing <- data.frame(fear = fleeing$fear[1:length(fleeing$fear)-1],
                           control_steps = fleeing$c[1:length(fleeing$c)-1],
                           x = 1:(length(fleeing$somatic)-1))

plot2_fleeing <- ggplot(data_fleeing, aes(x = x)) + 
  geom_point(aes(y = control_steps, color = "Perceived Control", shape = "Perceived Control"), size = 3) +
  geom_point(aes(y = fear, color = "Fear", shape = "Fear"), size = 3) +
  ylim(0, 1) +
  scale_y_continuous(name = "", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(name = "Steps", limits = factor(1:c(length(data_fleeing$fear)-1))) +
  scale_color_manual(name = "", 
                     values = c("Perceived Control" = "blue", "Fear" = "orange")) +
  scale_shape_manual(name = "", 
                     values = c("Perceived Control" = 17, "Fear" = 15)) +  # Corrected here
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),   # Increase legend text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.text.x = element_text(size = 10),   # Increase x-axis tick labels
        axis.text.y = element_text(size = 18),   # Increase y-axis tick labels
        plot.title = element_text(size = 18, face = "bold")) # Increase plot title size

pdf("Fleeing.pdf", width = 10, height = 5)
grid.arrange(plot1_fleeing[[4]], plot2_fleeing, ncol = 2)
dev.off()

# =========================
# Approach and Retract  
# =========================

retract <- simPhobia(control = 0.8,
                     punishment_h = 0.2,
                     gamma = .95,
                     size = 11,
                     meters = 20,
                     reward = 1,
                     grid_space = "abyss2",
                     plot = FALSE)

# Set up heatmap: 
grid_retract <- grid; grid_retract[retract$behavior] <- NA
plot1_retract <- pheatmap::pheatmap(grid_retract, legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
                                    display_numbers = FALSE, 
                                    color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
                                    na_col = "black",  width = 3, height = 3)

# Set up figure 2: 
data_retract <- data.frame(fear = retract$fear[1:length(retract$fear)-1],
                           control_steps = retract$c[1:length(retract$c)-1],
                           x = 1:(length(retract$fear)-1))

plot2_retract <- ggplot(data_retract, aes(x = x)) + 
  geom_point(aes(y = control_steps, color = "Perceived Control", shape = "Perceived Control"), size = 3) +
  geom_point(aes(y = fear, color = "Fear", shape = "Fear"), size = 3) +
  ylim(0, 1) +
  scale_y_continuous(name = "", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(name = "Steps", limits = factor(1:c(length(data_retract$fear)-1))) +
  scale_color_manual(name = "", 
                     values = c("Perceived Control" = "blue", "Fear" = "orange")) +
  scale_shape_manual(name = "", 
                     values = c("Perceived Control" = 17, "Fear" = 15)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),   # Increase legend text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.text.x = element_text(size = 10),   # Increase x-axis tick labels
        axis.text.y = element_text(size = 18)) # Increase plot title size

pdf("Retract.pdf", width = 10, height = 5)
grid.arrange(plot1_retract[[4]], plot2_retract, ncol = 2)
dev.off()

# =========================
# Courageous Behavior / Detour 
# =========================

detour <- simPhobia(control = 0.8,
                    punishment_h = 0.5,
                    gamma = .95,
                    size = 11,
                    meters = 20,
                    reward = 1,
                    grid_space = "abyss2",
                    plot = FALSE)

# Set up heatmap: 
grid_detour <- grid; grid_detour[detour$behavior] <- NA
plot1_detour <- pheatmap::pheatmap(grid_detour, legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
                                    display_numbers = FALSE, 
                                    color = c(colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(1), "white"),
                                    na_col = "black",  width = 3, height = 3)

# Set up figure 2: 
data_detour <- data.frame(fear = detour$fear[1:length(detour$fear)-1],
                          control_steps = detour$c[1:length(detour$c)-1],
                          x = 1:(length(detour$fear)-1))

plot2_detour <- ggplot(data_detour, aes(x = x)) + 
  geom_point(aes(y = control_steps, color = "Perceived Control", shape = "Perceived Control"), size = 3) +
  geom_point(aes(y = fear, color = "Fear", shape = "Fear"), size = 3) +
  ylim(0, 1) +
  scale_y_continuous(name = "", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(name = "Steps", limits = factor(1:c(length(data_detour$fear)-1))) +
  scale_color_manual(name = "", 
                     values = c("Perceived Control" = "blue", "Fear" = "orange")) +
  scale_shape_manual(name = "", 
                     values = c("Perceived Control" = 17, "Fear" = 15)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),   # Increase legend text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.text.x = element_text(size = 10),   # Increase x-axis tick labels
        axis.text.y = element_text(size = 18)) # Increase plot title size

pdf("Detour.pdf", width = 10, height = 5)
grid.arrange(plot1_detour[[4]], plot2_detour, ncol = 2)
dev.off()

 # =========================
# Nerves of Steal (direct route)
# =========================

direct <- simPhobia(control = 1,
                    punishment_h = 0.7,
                    gamma = .95,
                    size = 11,
                    meters = 20,
                    reward = 1,
                    grid_space = "abyss2",
                    plot = FALSE)

# Set up heatmap: 
grid_direct <- grid; grid_direct[direct$behavior] <- NA
plot1_direct <- pheatmap::pheatmap(grid_direct, legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
                                   display_numbers = FALSE, 
                                   color = c(colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(1), "white"),
                                   na_col = "black",  width = 3, height = 3)

# Set up figure 2: 
data_direct <- data.frame(fear = direct$fear[1:length(direct$fear)-1],
                          control_steps = direct$c[1:length(direct$c)-1],
                          x = 1:(length(direct$fear)-1))

plot2_direct <- ggplot(data_direct, aes(x = x)) + 
  geom_point(aes(y = control_steps, color = "Perceived Control", shape = "Perceived Control"), size = 3) +
  geom_point(aes(y = fear, color = "Fear", shape = "Fear"), size = 3) +
  ylim(0, 1) +
  scale_y_continuous(name = "", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(name = "Steps", limits = factor(1:c(length(data_direct$fear)-1))) +
  scale_color_manual(name = "", 
                     values = c("Perceived Control" = "blue", "Fear" = "orange")) +
  scale_shape_manual(name = "", 
                     values = c("Perceived Control" = 17, "Fear" = 15)) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.text = element_text(size = 18),   # Increase legend text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.text.x = element_text(size = 10),   # Increase x-axis tick labels
        axis.text.y = element_text(size = 18)) # Increase plot title size


pdf("Direct.pdf", width = 10, height = 5)
grid.arrange(plot1_direct[[4]], plot2_direct, ncol = 2)
dev.off()
 
