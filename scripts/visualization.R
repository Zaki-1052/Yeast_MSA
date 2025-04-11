# Install required packages if needed
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

# Load libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Create output directory
dir.create("results/visualization", showWarnings = FALSE)

# 1. Variant Count and Type Visualization
# Read the group comparison data
group_data <- read.table("results/functional/group_comparison.tsv", header=TRUE, sep="\t")

# Prepare data for stacked bar plot
variant_types <- data.frame(
  Treatment = group_data$Group,
  SNPs = group_data$SNPs,
  Indels = group_data$Indels
)

# Reshape data for ggplot
variant_long <- melt(variant_types, id.vars="Treatment", 
                     variable.name="Variant_Type", value.name="Count")

# Create stacked bar plot
p1 <- ggplot(variant_long, aes(x=Treatment, y=Count, fill=Variant_Type)) +
  geom_bar(stat="identity", position="stack") +
  theme_minimal() +
  labs(title="Variant Distribution by Treatment",
       x="Treatment Group", y="Number of Variants",
       fill="Variant Type") +
  scale_fill_brewer(palette="Set1")

# Save the plot
ggsave("results/visualization/variant_distribution.png", p1, width=8, height=6)

# 2. Variant Sharing Heatmap
# Create a sharing matrix
treatments <- c("WT", "STC", "CAS", "WTA")
sharing_matrix <- matrix(0, nrow=4, ncol=4)
rownames(sharing_matrix) <- treatments
colnames(sharing_matrix) <- treatments

# Diagonal - total variants per treatment
for (i in 1:4) {
  sharing_matrix[i,i] <- group_data$Variant_Count[i]
}

# Fill in shared variants
sharing_matrix["WT", "STC"] <- 25
sharing_matrix["WT", "CAS"] <- 37
sharing_matrix["WT", "WTA"] <- 86
sharing_matrix["STC", "WT"] <- 25
sharing_matrix["STC", "CAS"] <- 12
sharing_matrix["STC", "WTA"] <- 20
sharing_matrix["CAS", "WT"] <- 37
sharing_matrix["CAS", "STC"] <- 12
sharing_matrix["CAS", "WTA"] <- 32
sharing_matrix["WTA", "WT"] <- 86
sharing_matrix["WTA", "STC"] <- 20
sharing_matrix["WTA", "CAS"] <- 32

# Create heatmap
png("results/visualization/variant_sharing_heatmap.png", width=800, height=600)
pheatmap(sharing_matrix, 
         main="Variant Sharing Between Treatment Groups",
         color=colorRampPalette(c("white", "steelblue", "darkblue"))(100),
         display_numbers=TRUE,
         cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

# 3. Variant Quality Metrics
quality_data <- data.frame(
  Treatment = group_data$Group,
  Average_Quality = group_data$Avg_Quality,
  Median_Depth = group_data$Median_Depth
)

# Create quality plot
p3 <- ggplot(quality_data, aes(x=Treatment, y=Average_Quality, fill=Treatment)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(Average_Quality, 1)), vjust=-0.5) +
  theme_minimal() +
  labs(title="Average Variant Quality by Treatment",
       x="Treatment Group", y="Average Quality Score") +
  scale_fill_brewer(palette="Set2")

# Save quality plot
ggsave("results/visualization/variant_quality.png", p3, width=8, height=6)

# 4. Scaffold Distribution Visualization
# Function to read scaffold distribution
read_scaffold_dist <- function(group) {
  file_path <- paste0("results/functional/", group, "/scaffold_distribution.txt")
  scaffolds <- read.table(file_path, col.names=c("Count", "Scaffold"))
  scaffolds$Group <- group
  return(scaffolds)
}

# Read all scaffold distributions
wt_scaffolds <- read_scaffold_dist("WT")
stc_scaffolds <- read_scaffold_dist("STC")
cas_scaffolds <- read_scaffold_dist("CAS")
wta_scaffolds <- read_scaffold_dist("WTA")

# Combine and get top 10 scaffolds overall
all_scaffolds <- rbind(wt_scaffolds, stc_scaffolds, cas_scaffolds, wta_scaffolds)
top_scaffolds <- names(sort(table(all_scaffolds$Scaffold), decreasing=TRUE)[1:10])

# Filter for just the top scaffolds
top_data <- all_scaffolds[all_scaffolds$Scaffold %in% top_scaffolds,]

# Create scaffold distribution plot
p4 <- ggplot(top_data, aes(x=Scaffold, y=Count, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Distribution of Variants in Top 10 Scaffolds",
       x="Scaffold", y="Number of Variants") +
  scale_fill_brewer(palette="Set3")

# Save scaffold plot
ggsave("results/visualization/scaffold_distribution.png", p4, width=10, height=6)

# 5. Indel Percentage Visualization
indel_percent <- data.frame(
  Treatment = group_data$Group,
  Percent_Indels = 100 * group_data$Indels / group_data$Variant_Count
)

p5 <- ggplot(indel_percent, aes(x=Treatment, y=Percent_Indels, fill=Treatment)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=sprintf("%.1f%%", Percent_Indels)), vjust=-0.5) +
  theme_minimal() +
  labs(title="Percentage of Indels in Treatment-Specific Variants",
       x="Treatment Group", y="Indel Percentage") +
  scale_fill_brewer(palette="Set3")

# Save indel percentage plot
ggsave("results/visualization/indel_percentage.png", p5, width=8, height=6)

cat("Visualization complete. Results saved to results/visualization/ directory.\n")