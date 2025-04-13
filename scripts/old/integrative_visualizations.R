# simplified_integrative_viz_fixed.R

# Load required libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(igraph)
library(patchwork)
library(viridis)

# Create output directory
dir.create("results/visualizations/integrative", recursive = TRUE, showWarnings = FALSE)

#-----------------------------------------
# 1. Load all our analysis data
#-----------------------------------------

# Load clustering data
clustering_summary <- read.table("results/visualizations/clustering/clustering_summary.tsv", 
                                 header = TRUE, sep = "\t")

# Load transition/transversion data
trans_data <- read.table("results/visualizations/signatures/transition_transversion.tsv", 
                         header = TRUE, sep = "\t")

# Load mutation frequencies
mutation_freqs <- read.table("results/visualizations/signatures/mutation_frequencies.tsv", 
                             header = TRUE, sep = "\t")

# Get treatment scaffold information
treatments <- c("WT", "STC", "CAS", "WTA")
hotspot_scaffolds <- c()

for (treatment in treatments) {
  file_path <- paste0("results/analysis/hotspots/", treatment, "_top_hotspots.txt")
  if (file.exists(file_path)) {
    lines <- readLines(file_path)
    if (length(lines) > 1) {
      scaffold <- strsplit(lines[2], "\\s+")[[1]][1]
      hotspot_scaffolds <- c(hotspot_scaffolds, scaffold)
    }
  }
}

names(hotspot_scaffolds) <- treatments

#-----------------------------------------
# 2. Calculate treatment similarities
#-----------------------------------------

# Create a similarity matrix based on multiple metrics
similarity_matrix <- matrix(0, nrow = length(treatments), ncol = length(treatments))
rownames(similarity_matrix) <- treatments
colnames(similarity_matrix) <- treatments

# Compare treatments pairwise
for (i in 1:length(treatments)) {
  for (j in 1:length(treatments)) {
    if (i != j) {
      # Metric 1: Similarity in clustering percentage (normalized difference)
      clust_i <- clustering_summary$Percent_Clustered[clustering_summary$Treatment == treatments[i]]
      clust_j <- clustering_summary$Percent_Clustered[clustering_summary$Treatment == treatments[j]]
      clust_sim <- 1 - abs(clust_i - clust_j) / 100
      
      # Metric 2: Similarity in transition percentage
      trans_i <- trans_data$Percentage[trans_data$Treatment == treatments[i] & 
                                         trans_data$Mutation_Class == "Transition"]
      trans_j <- trans_data$Percentage[trans_data$Treatment == treatments[j] & 
                                         trans_data$Mutation_Class == "Transition"]
      if(length(trans_i) == 0) trans_i <- 0
      if(length(trans_j) == 0) trans_j <- 0
      trans_sim <- 1 - abs(trans_i - trans_j) / 100
      
      # Metric 3: Similarity in mutation spectrum
      # Extract mutation percentages for each treatment
      mut_i <- mutation_freqs %>% 
        filter(Treatment == treatments[i]) %>%
        select(Mutation, Percentage)
      
      mut_j <- mutation_freqs %>% 
        filter(Treatment == treatments[j]) %>%
        select(Mutation, Percentage)
      
      # Calculate cosine similarity between mutation spectra
      if(nrow(mut_i) > 0 && nrow(mut_j) > 0) {
        mut_sim_df <- merge(mut_i, mut_j, by = "Mutation")
        if(nrow(mut_sim_df) > 0) {
          # Cosine similarity
          cosine_sim <- sum(mut_sim_df$Percentage.x * mut_sim_df$Percentage.y) / 
            (sqrt(sum(mut_sim_df$Percentage.x^2)) * sqrt(sum(mut_sim_df$Percentage.y^2)))
          mut_sim <- cosine_sim
        } else {
          mut_sim <- 0
        }
      } else {
        mut_sim <- 0
      }
      
      # Combine metrics (equal weights for now)
      similarity_matrix[i, j] <- (clust_sim + trans_sim + mut_sim) / 3
    }
  }
}

#-----------------------------------------
# 3. Create treatment relationship network (SIMPLIFIED APPROACH)
#-----------------------------------------

# Create a data frame for the treatment nodes
nodes_df <- data.frame(
  Treatment = treatments,
  Clustering = clustering_summary$Percent_Clustered,
  stringsAsFactors = FALSE
)

# Add transition percentages
nodes_df$Transitions <- 0
for (i in 1:nrow(nodes_df)) {
  trans_pct <- trans_data$Percentage[trans_data$Treatment == nodes_df$Treatment[i] & 
                                       trans_data$Mutation_Class == "Transition"]
  if(length(trans_pct) > 0) {
    nodes_df$Transitions[i] <- trans_pct
  }
}

# Create a data frame for the edges (connections between treatments)
edges_df <- data.frame(
  from = character(),
  to = character(),
  weight = numeric(),
  stringsAsFactors = FALSE
)

# Fill the edges data frame
for (i in 1:length(treatments)) {
  for (j in 1:length(treatments)) {
    if (i != j && similarity_matrix[i, j] > 0) {
      edges_df <- rbind(edges_df, 
                        data.frame(from = treatments[i], 
                                   to = treatments[j], 
                                   weight = similarity_matrix[i, j],
                                   stringsAsFactors = FALSE))
    }
  }
}

# Calculate node positions manually for a simple circular layout
n_nodes <- length(treatments)
angle <- seq(0, 2*pi, length.out = n_nodes + 1)[1:n_nodes]
nodes_df$x <- 5 + 3 * cos(angle)
nodes_df$y <- 5 + 3 * sin(angle)

# Create edges data with positions
edges_plot_df <- data.frame()
for (i in 1:nrow(edges_df)) {
  # Get positions for source node
  from_pos <- nodes_df[nodes_df$Treatment == edges_df$from[i], c("x", "y")]
  # Get positions for target node
  to_pos <- nodes_df[nodes_df$Treatment == edges_df$to[i], c("x", "y")]
  
  # Add to plotting data frame
  edges_plot_df <- rbind(edges_plot_df,
                         data.frame(x = from_pos$x, 
                                    y = from_pos$y,
                                    xend = to_pos$x,
                                    yend = to_pos$y,
                                    weight = edges_df$weight[i],
                                    stringsAsFactors = FALSE))
}

# Create treatment relationship network plot manually
p1 <- ggplot() +
  # Add edges
  geom_segment(data = edges_plot_df, 
               aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight),
               color = "gray50", alpha = 0.8) +
  # Add nodes
  geom_point(data = nodes_df, 
             aes(x = x, y = y, size = Clustering, fill = Transitions),
             shape = 21, color = "black") +
  # Add labels
  geom_text(data = nodes_df,
            aes(x = x, y = y, label = Treatment),
            size = 5, nudge_y = 0.4) +
  # Styling
  scale_fill_viridis(option = "plasma", name = "% Transitions") +
  scale_size_continuous(name = "% Clustered", range = c(5, 12)) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5),
         fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  theme_void() +
  labs(title = "Treatment Relationship Network",
       subtitle = "Node size = clustering %, Node color = transition %\nEdge thickness = overall similarity") +
  theme(legend.position = "right",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

#-----------------------------------------
# 4. Create comprehensive treatment summary
#-----------------------------------------

# Extract dominant mutations for each treatment
dominant_mutations <- data.frame(
  Treatment = character(0),
  DominantMutation = character(0),
  Percentage = numeric(0)
)

for (treatment in treatments) {
  treatment_data <- mutation_freqs %>% 
    filter(Treatment == treatment) %>%
    arrange(desc(Percentage))
  
  if(nrow(treatment_data) > 0) {
    top_mutation <- treatment_data[1, ]
    dominant_mutations <- rbind(dominant_mutations, 
                                data.frame(Treatment = treatment,
                                           DominantMutation = top_mutation$Mutation,
                                           Percentage = top_mutation$Percentage))
  }
}

# Format data for summary visualization
summary_data <- clustering_summary %>%
  left_join(dominant_mutations, by = "Treatment") %>%
  left_join(trans_data %>% 
              filter(Mutation_Class == "Transition") %>%
              select(Treatment, TransitionPct = Percentage), 
            by = "Treatment")

# Create a data frame with top hotspot information
hotspot_data <- data.frame(
  Treatment = treatments,
  TopHotspot = hotspot_scaffolds,
  stringsAsFactors = FALSE
)

summary_data <- summary_data %>%
  left_join(hotspot_data, by = "Treatment")

# Create summary table visualization - FIXED ANNOTATION ISSUE
p2 <- ggplot(summary_data) +
  geom_segment(aes(x = 0, xend = Percent_Clustered, y = Treatment, yend = Treatment),
               color = "darkblue", size = 10, alpha = 0.5) +
  geom_text(aes(x = Percent_Clustered/2, y = Treatment, 
                label = paste0(Percent_Clustered, "% clustered")), 
            color = "white", fontface = "bold") +
  geom_text(aes(x = Percent_Clustered + 5, y = Treatment,
                label = paste0(TransitionPct, "% transitions")),
            hjust = 0) +
  geom_text(aes(x = Percent_Clustered + 35, y = Treatment,
                label = paste0(DominantMutation, " (", round(Percentage, 1), "%)")),
            hjust = 0) +
  geom_text(aes(x = Percent_Clustered + 65, y = Treatment,
                label = TopHotspot),
            hjust = 0) +
  labs(title = "Treatment Effects Summary",
       subtitle = "Key genetic signatures of each treatment") +
  # Fix annotation to use fixed x-positions
  annotate("text", x = 100, y = length(treatments) + 0.5, 
           label = "Dominant Mutation", fontface = "bold", hjust = 0) +
  annotate("text", x = 130, y = length(treatments) + 0.5, 
           label = "Top Hotspot", fontface = "bold", hjust = 0) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 12))

#-----------------------------------------
# 5. Save visualizations
#-----------------------------------------

# Save individual plots
ggsave("results/visualizations/integrative/treatment_network.png", p1, width = 10, height = 8)
ggsave("results/visualizations/integrative/treatment_summary.png", p2, width = 12, height = 6)

# Create a combined figure
p_combined <- p1 / p2 + patchwork::plot_layout(heights = c(1.5, 1))
ggsave("results/visualizations/integrative/combined_summary.png", p_combined, width = 12, height = 14)

# Create data summary file
write.table(summary_data, 
            "results/visualizations/integrative/treatment_summary.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Print summary
print("Integrative analysis summary:")
print(summary_data)
cat("\nRelationship strength between treatments:\n")
print(round(similarity_matrix, 2))