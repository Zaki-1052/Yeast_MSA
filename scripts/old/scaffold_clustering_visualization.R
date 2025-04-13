# scaffold_clustering_visualization.R

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create output directory
dir.create("results/visualizations/clustering", recursive = TRUE, showWarnings = FALSE)

# 1. Load position data for JRIU01000031.1 for all treatments
treatments <- c("WT", "STC", "CAS", "WTA")
position_data <- data.frame()

for (treatment in treatments) {
  if (file.exists(paste0("results/analysis/position_clustering/", treatment, "_positions.txt"))) {
    positions <- read.table(paste0("results/analysis/position_clustering/", treatment, "_positions.txt"))
    positions$Treatment <- treatment
    position_data <- rbind(position_data, positions)
  }
}

# Rename columns for clarity
colnames(position_data) <- c("Position", "Treatment")

# Get range of positions
min_pos <- min(position_data$Position)
max_pos <- max(position_data$Position)
scaffold_length <- max_pos - min_pos + 1

# 2. Create full scaffold position distribution plot
p1 <- ggplot(position_data, aes(x = Position, y = Treatment, color = Treatment)) +
  geom_jitter(height = 0.2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Variant Positions in JRIU01000031.1",
       subtitle = paste0("Scaffold span: ", scaffold_length, " bp"),
       x = "Position (bp)",
       y = "") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# 3. Create density distribution to better visualize clusters
p2 <- ggplot(position_data, aes(x = Position, fill = Treatment)) +
  geom_density(alpha = 0.5) +
  facet_grid(Treatment ~ ., scales = "free_y") +
  theme_minimal() +
  labs(title = "Variant Density Distribution",
       x = "Position (bp)",
       y = "Density") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"))

# 4. Calculate & visualize clustering metrics
cluster_summary <- data.frame(
  Treatment = treatments,
  Total_Variants = numeric(length(treatments)),
  Clustered_Variants = numeric(length(treatments)),
  Percent_Clustered = numeric(length(treatments))
)

for (i in 1:length(treatments)) {
  treatment <- treatments[i]
  
  # Get statistics from the text file
  stats_file <- paste0("results/analysis/position_clustering/", treatment, "_stats.txt")
  if (file.exists(stats_file)) {
    stats_text <- readLines(stats_file)
    
    # Extract values using regex
    total_match <- regexpr("Total variants:\\s+(\\d+)", stats_text[2])
    if (total_match > 0) {
      total_str <- regmatches(stats_text[2], total_match)
      total_variants <- as.numeric(gsub("Total variants:\\s+", "", total_str))
      cluster_summary$Total_Variants[i] <- total_variants
    }
    
    clustered_match <- regexpr("Variants <10bp apart:\\s+(\\d+)", stats_text[5])
    if (clustered_match > 0) {
      clustered_str <- regmatches(stats_text[5], clustered_match)
      clustered_variants <- as.numeric(gsub("Variants <10bp apart:\\s+", "", clustered_str))
      cluster_summary$Clustered_Variants[i] <- clustered_variants
      cluster_summary$Percent_Clustered[i] <- round(clustered_variants / total_variants * 100, 1)
    }
  }
}

# Create clustering metrics visualization
p3 <- ggplot(cluster_summary, aes(x = Treatment, y = Percent_Clustered, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percent_Clustered, "%")), vjust = -0.5) +
  ylim(0, 100) +
  theme_minimal() +
  labs(title = "Percentage of Variants in Clusters (<10bp apart)",
       x = "",
       y = "Percent") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"))

# 5. Save plots
ggsave("results/visualizations/clustering/position_distribution.png", p1, width = 10, height = 6)
ggsave("results/visualizations/clustering/density_distribution.png", p2, width = 10, height = 8)
ggsave("results/visualizations/clustering/clustering_metrics.png", p3, width = 8, height = 6)

# 6. Save summary table
write.table(cluster_summary, 
            "results/visualizations/clustering/clustering_summary.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Print summary to console
print("Clustering summary:")
print(cluster_summary)
print("Visualization complete. Check results/visualizations/clustering/ directory.")