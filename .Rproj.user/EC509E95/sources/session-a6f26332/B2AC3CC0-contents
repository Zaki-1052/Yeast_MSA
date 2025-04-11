# hotspot_visualization.R - Create visualizations for variant hotspots

# Load libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)

# Set working directory if needed
# setwd("/path/to/yeast_analysis")

# Create output directory
dir.create("results/analysis/hotspots/plots", showWarnings = FALSE, recursive = TRUE)

# 1. Load and combine density data from all treatments
treatments <- c("WT", "STC", "CAS", "WTA")
density_data <- data.frame()

for (treatment in treatments) {
  file_path <- paste0("results/analysis/hotspots/", treatment, "_density.tsv")
  if (file.exists(file_path)) {
    data <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
    data$Treatment <- treatment
    density_data <- rbind(density_data, data)
  }
}

# 2. Identify top 10 scaffolds by total variants across all treatments
top_scaffolds <- density_data %>%
  group_by(Scaffold) %>%
  summarize(Total_Variants = sum(Variants)) %>%
  arrange(desc(Total_Variants)) %>%
  head(10)

# 3. Plot variant density for top scaffolds across treatments
top_scaffold_data <- density_data %>%
  filter(Scaffold %in% top_scaffolds$Scaffold)

# Create horizontal bar plot
p1 <- ggplot(top_scaffold_data, aes(x = reorder(Scaffold, -Density_per_kb), y = Density_per_kb, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Variant Density in Top Hotspot Scaffolds",
       x = "Scaffold", 
       y = "Variants per kb") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/analysis/hotspots/plots/top_hotspot_density.pdf", p1, width = 10, height = 7)

# 4. Analyze variant positions in the top hotspot scaffold
pos_data <- data.frame()

for (treatment in treatments) {
  file_path <- paste0("results/analysis/hotspots/position_analysis/", treatment, "_JRIU01000031.1_variants.tsv")
  if (file.exists(file_path)) {
    data <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(data) <- c("Scaffold", "Position", "Ref", "Alt", "Quality", "Depth")
    data$Treatment <- treatment
    pos_data <- rbind(pos_data, data)
  }
}

# Add mutation type
pos_data$Type <- ifelse(nchar(pos_data$Ref) == 1 & nchar(pos_data$Alt) == 1, "SNP", "Indel")

# Create position distribution plot
p2 <- ggplot(pos_data, aes(x = Position, y = Treatment, color = Type)) +
  geom_point(aes(size = Quality)) +
  theme_minimal() +
  labs(title = "Variant Positions in JRIU01000031.1 Hotspot",
       x = "Position", 
       y = "Treatment Group") +
  scale_color_manual(values = c("SNP" = "blue", "Indel" = "red"))

ggsave("results/analysis/hotspots/plots/hotspot_positions.pdf", p2, width = 12, height = 6)

# 5. Create mutation type proportion plot
type_summary <- pos_data %>%
  group_by(Treatment, Type) %>%
  summarize(Count = n()) %>%
  group_by(Treatment) %>%
  mutate(Proportion = Count / sum(Count))

p3 <- ggplot(type_summary, aes(x = Treatment, y = Proportion, fill = Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Mutation Type Proportions in JRIU01000031.1",
       x = "Treatment", 
       y = "Proportion") +
  scale_fill_manual(values = c("SNP" = "blue", "Indel" = "red"))

ggsave("results/analysis/hotspots/plots/hotspot_mutation_types.pdf", p3, width = 8, height = 6)

print("Hotspot visualizations complete.")