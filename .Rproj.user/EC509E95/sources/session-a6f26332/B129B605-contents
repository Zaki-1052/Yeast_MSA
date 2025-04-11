# mutation_signature_visualization.R

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)

# Create output directory
dir.create("results/visualizations/signatures", recursive = TRUE, showWarnings = FALSE)

# 1. Load mutation pattern data for each treatment
treatments <- c("WT", "STC", "CAS", "WTA")
mutation_data <- data.frame()

for (treatment in treatments) {
  # Get file path for treatment-specific variants
  file_path <- paste0("results/analysis/signatures/", treatment, "_specific_variants.tsv")
  
  if (file.exists(file_path)) {
    # Read file
    variants <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(variants) >= 4) {
      colnames(variants) <- c("Scaffold", "Position", "Ref", "Alt", "Quality", "Depth")
      
      # Filter for SNPs (single base substitutions)
      snps <- variants %>% 
        filter(nchar(Ref) == 1 & nchar(Alt) == 1)
      
      if (nrow(snps) > 0) {
        # Create mutation type column
        snps$Mutation <- paste0(snps$Ref, ">", snps$Alt)
        snps$Treatment <- treatment
        
        # Add to combined data
        mutation_data <- rbind(mutation_data, 
                               select(snps, Scaffold, Position, Ref, Alt, Mutation, Treatment))
      }
    }
  }
}

# 2. Create mutation spectrum visualization
# Get frequencies
mutation_counts <- mutation_data %>%
  group_by(Treatment, Mutation) %>%
  summarize(Count = n(), .groups = 'drop') %>%
  group_by(Treatment) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Define mutation categories and colors
mutation_types <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
mutation_colors <- c("A>C" = "#8dd3c7", "A>G" = "#ffffb3", "A>T" = "#bebada", 
                     "C>A" = "#fb8072", "C>G" = "#80b1d3", "C>T" = "#fdb462",
                     "G>A" = "#b3de69", "G>C" = "#fccde5", "G>T" = "#d9d9d9",
                     "T>A" = "#bc80bd", "T>C" = "#ccebc5", "T>G" = "#ffed6f")

# Ensure all mutation types are represented
all_combinations <- expand.grid(Treatment = treatments, Mutation = mutation_types, stringsAsFactors = FALSE)
mutation_counts <- merge(all_combinations, mutation_counts, by = c("Treatment", "Mutation"), all.x = TRUE)
mutation_counts[is.na(mutation_counts)] <- 0

# Plot mutation spectrum
p1 <- ggplot(mutation_counts, aes(x = Mutation, y = Percentage, fill = Mutation)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Treatment, ncol = 1) +
  scale_fill_manual(values = mutation_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold")) +
  labs(title = "Mutation Spectrum by Treatment",
       subtitle = "Percentage of each mutation type",
       x = "Mutation Type",
       y = "Percentage (%)")

# 3. Create comparative heatmap
# Reshape data for heatmap
mutation_matrix <- dcast(mutation_counts, Treatment ~ Mutation, value.var = "Percentage")
row.names(mutation_matrix) <- mutation_matrix$Treatment
mutation_matrix <- mutation_matrix[, -1]  # Remove Treatment column

# Convert to matrix
mutation_matrix <- as.matrix(mutation_matrix)

# Create heatmap using ggplot
heatmap_data <- melt(mutation_matrix)
colnames(heatmap_data) <- c("Treatment", "Mutation", "Percentage")

p2 <- ggplot(heatmap_data, aes(x = Mutation, y = Treatment, fill = Percentage)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "plasma") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), color = "black", size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Mutation Type Comparison Across Treatments",
       x = "Mutation Type",
       y = "Treatment",
       fill = "Percentage (%)")

# 4. Analyze transitions vs transversions
mutation_data <- mutation_data %>%
  mutate(Mutation_Class = case_when(
    Mutation %in% c("A>G", "G>A", "C>T", "T>C") ~ "Transition",
    TRUE ~ "Transversion"
  ))

transition_summary <- mutation_data %>%
  group_by(Treatment, Mutation_Class) %>%
  summarize(Count = n(), .groups = 'drop') %>%
  group_by(Treatment) %>%
  mutate(Percentage = Count / sum(Count) * 100)

p3 <- ggplot(transition_summary, aes(x = Treatment, y = Percentage, fill = Mutation_Class)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("Transition" = "#4daf4a", "Transversion" = "#e41a1c")) +
  theme_minimal() +
  labs(title = "Transitions vs Transversions by Treatment",
       x = "Treatment",
       y = "Percentage (%)",
       fill = "Mutation Class")

# 5. Save plots
ggsave("results/visualizations/signatures/mutation_spectrum.png", p1, width = 10, height = 10)
ggsave("results/visualizations/signatures/mutation_heatmap.png", p2, width = 12, height = 6)
ggsave("results/visualizations/signatures/transition_transversion.png", p3, width = 9, height = 6)

# 6. Save summary data
write.table(mutation_counts, 
            "results/visualizations/signatures/mutation_frequencies.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(transition_summary, 
            "results/visualizations/signatures/transition_transversion.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Print summary
print("Mutation signature analysis complete.")
print("Transition vs Transversion summary:")
print(transition_summary)