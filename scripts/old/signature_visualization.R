# Visualization of mutation signatures
# signature_visualization.R

library(ggplot2)
library(reshape2)

# Create clustering visualization for JRIU01000031.1
# We'll visualize the position distribution using a density plot
for(treatment in c("WT", "STC", "CAS", "WTA")) {
  positions <- read.table(paste0("results/analysis/position_clustering/", treatment, "_JRIU01000031.1_positions.txt"))
  
  # Plot position density
  pdf(paste0("results/analysis/position_clustering/", treatment, "_position_density.pdf"), width=10, height=6)
  ggplot(positions, aes(x=V1)) + 
    geom_density(fill="blue", alpha=0.5) +
    labs(title=paste(treatment, "- Variant Position Density in JRIU01000031.1"),
         x="Position", y="Density") +
    theme_minimal()
  dev.off()
}