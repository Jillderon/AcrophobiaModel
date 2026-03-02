################################################################################
#                       Visualizing the results from the                       #
#                       computational model of Acrophobia                      #
#                      (phenomena as statistical patterns)                     #
################################################################################ 

# Load libraries: 
library("ggplot2")
library("dplyr")
library("patchwork")

# Load in data: 
data <- read.table("output.csv", header = FALSE, stringsAsFactors = FALSE)
data <- data[, -1, drop = FALSE]

colnames(data) <- c("Control", "Punishment_h", "Meter", "Gridspace", "Goal",
                    "Steps", "MeanFear", "MeanBodily", "MeanThreat", "MeanDistAbyss")


# Ensure numeric columns are parsed consistently.
data$Control <- as.numeric(data$Control)
data$Punishment_h <- as.numeric(data$Punishment_h)
data$Meter <- as.numeric(data$Meter)
data$MeanFear <- as.numeric(data$MeanFear)
data$MeanDistAbyss <- as.numeric(data$MeanDistAbyss)

######################################
# Create Figure
######################################

for (punishment_h in c(0.2, 0.5, 0.7)) {
  dataTemp <- data %>% filter(Punishment_h == punishment_h)
  
  # Ensure that 'control' and 'elevation' are treated as factors (ordinal variables)
  dataTemp$Control <- factor(dataTemp$Control, ordered = TRUE)
  dataTemp$Meter <- factor(dataTemp$Meter, ordered = TRUE)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Mean Fear per Elevation level and Perceived Control  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Fear <- ggplot(dataTemp, aes(x = Meter, y = MeanFear, group = Control, color = Control)) +
    geom_line() +   # Use lines to represent different levels of control
    geom_point() +  # Optionally, add points to show data points
    ylim(0,.8) +
    labs(x = "Elevation (in meters)", y = "Mean Fear Response", color = "Perceived Control") +  # Label axes and legend
    theme_minimal() +
    theme(
      axis.line = element_line(),       # Add x and y axis lines
      axis.ticks = element_line()       # Add ticks on both axes
    ) 
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Mean Distance to Threat per Elevation level and Perceived Control 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Distance <- ggplot(dataTemp, aes(x = Meter, y = MeanDistAbyss, group = Control, color = Control)) +
    geom_line() +   # Use lines to represent different levels of control
    geom_point() +  # Optionally, add points to show data points
    labs(x = "Elevation (in meters)", y = "Mean Distance to Abyss", color = "Perceived Control") +  # Label axes and legend
    ylim(2,10) +
    theme_minimal() +
    theme(
      axis.line = element_line(),       # Add x and y axis lines
      axis.ticks = element_line()       # Add ticks on both axes
    ) 
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Save as one figure
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #pdf(paste0("Phenomena_h", punishment_h ,".pdf"), width = 8, height = 5)
  Fear + Distance + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom") # or "top"
  #dev.off()
}
