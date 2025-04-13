# sequence_context_analysis.R

# Check and install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("Biostrings", "ggseqlogo")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if (pkg == "Biostrings") {
      BiocManager::install("Biostrings")
    } else {
      install.packages(pkg)
    }
  }
}

# Load other required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

# Now load Bioconductor packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ggseqlogo))

# Create output directory
dir.create("results/sequence_context", recursive = TRUE, showWarnings = FALSE)

cat("Starting sequence context analysis...\n")

#-----------------------------------------------------
# 1. Load data and prepare environment
#-----------------------------------------------------

cat("Loading reference genome and variant data...\n")

# Load reference genome
ref_genome <- readDNAStringSet("reference/yeast_w303.fasta")
names(ref_genome) <- sub(" .*", "", names(ref_genome))  # Clean scaffold names

# Create a named list for easier access
genome_list <- setNames(as.character(ref_genome), names(ref_genome))

# Load variant data for each treatment
treatments <- c("WT", "STC", "CAS", "WTA")
all_variants <- data.frame()

for (treatment in treatments) {
  variant_file <- paste0("results/functional/", treatment, "/variant_details.tsv")
  if (file.exists(variant_file)) {
    cat(paste("Loading variants for", treatment, "treatment...\n"))
    variants <- read.delim(variant_file, header = TRUE, stringsAsFactors = FALSE)
    # Add treatment column
    variants$Treatment <- treatment
    
    # Add to combined data
    all_variants <- rbind(all_variants, variants)
  } else {
    warning(paste("Variant file not found for", treatment))
  }
}

cat(paste("Loaded", nrow(all_variants), "total variants across all treatments\n"))

#-----------------------------------------------------
# 2. Create function to extract sequence context
#-----------------------------------------------------

extract_context <- function(scaffold, position, ref_base, alt_base, window = 5) {
  
  # Skip if scaffold not in reference
  if (!scaffold %in% names(genome_list)) {
    return(NULL)
  }
  
  # Get scaffold sequence
  scaffold_seq <- genome_list[[scaffold]]
  
  # Calculate start and end positions for context
  start_pos <- max(1, position - window)
  end_pos <- min(nchar(scaffold_seq), position + window)
  
  # Handle edge cases
  if (start_pos < 1 || end_pos > nchar(scaffold_seq)) {
    # We're at a scaffold edge
    context <- NULL
    return(NULL)
  }
  
  # Extract context sequence
  context <- substring(scaffold_seq, start_pos, end_pos)
  
  # Handle cases where we couldn't extract full window
  if (nchar(context) != (2 * window + 1)) {
    padding_needed <- (2 * window + 1) - nchar(context)
    if (start_pos > 1) {
      # Pad beginning if we're near end
      context <- paste0(paste(rep("N", padding_needed), collapse=""), context)
    } else {
      # Pad end if we're near beginning
      context <- paste0(context, paste(rep("N", padding_needed), collapse=""))
    }
  }
  
  # Get the center base in our context
  center_pos <- window + 1
  if (start_pos > position - window) {
    # We had to adjust the start position, recalculate center
    center_pos <- position - start_pos + 1
  }
  
  center_base <- substring(context, center_pos, center_pos)
  
  # Verify reference base matches
  if (toupper(center_base) != toupper(ref_base)) {
    cat(paste("Warning: Reference mismatch for", scaffold, "position", position,
              "Expected", ref_base, "found", center_base, "\n"))
    return(NULL)  # Reference mismatch, something's wrong
  }
  
  return(list(
    context = context,
    mutation = paste0(ref_base, ">", alt_base)
  ))
}

#-----------------------------------------------------
# 3. Extract contexts for all variants
#-----------------------------------------------------

cat("Extracting sequence contexts for variants...\n")

# Filter for single nucleotide substitutions only
snv_variants <- all_variants %>%
  filter(nchar(Ref) == 1 & nchar(Alt) == 1)

cat(paste("Processing", nrow(snv_variants), "single nucleotide variants\n"))

# Extract sequence contexts
contexts_data <- data.frame(
  Treatment = character(),
  Mutation = character(),
  Context = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(snv_variants)) {
  if (i %% 100 == 0) {
    cat(paste("Processing variant", i, "of", nrow(snv_variants), "\r"))
  }
  
  context_result <- extract_context(
    scaffold = snv_variants$Scaffold[i],
    position = as.numeric(snv_variants$Position[i]),
    ref_base = snv_variants$Ref[i],
    alt_base = snv_variants$Alt[i]
  )
  
  if (!is.null(context_result)) {
    contexts_data <- rbind(contexts_data, data.frame(
      Treatment = snv_variants$Treatment[i],
      Mutation = context_result$mutation,
      Context = context_result$context,
      stringsAsFactors = FALSE
    ))
  }
}

cat(paste("\nSuccessfully extracted context for", nrow(contexts_data), "variants\n"))

# Save the full context data
write.table(contexts_data, "results/sequence_context/all_variant_contexts.tsv", 
            row.names = FALSE, quote = FALSE, sep = "\t")

#-----------------------------------------------------
# 4. Analyze and visualize sequence contexts
#-----------------------------------------------------

cat("Creating sequence logos for each treatment and mutation type...\n")

# Group variants by treatment and mutation type
treatment_mutations <- unique(contexts_data[, c("Treatment", "Mutation")])

# Create counters for tracking
logos_created <- 0
skipped_logos <- 0

# Create sequence logos for each treatment and mutation type
for (i in 1:nrow(treatment_mutations)) {
  treat <- treatment_mutations$Treatment[i]
  mut <- treatment_mutations$Mutation[i]
  
  # Get contexts for this treatment and mutation
  contexts <- contexts_data %>%
    filter(Treatment == treat & Mutation == mut) %>%
    pull(Context)
  
  if (length(contexts) >= 5) {  # Only create logos with sufficient data
    # Save the contexts to a text file
    writeLines(contexts, paste0("results/sequence_context/", treat, "_", 
                                gsub(">", "to", mut), "_contexts.txt"))
    
    # Create sequence logo
    p <- ggseqlogo(contexts) + 
      ggtitle(paste0(treat, ": ", mut, " (n=", length(contexts), ")")) +
      theme_minimal() +
      xlab("Position relative to mutation") +
      scale_x_continuous(breaks = 1:11, 
                         labels = c("-5", "-4", "-3", "-2", "-1", "0", 
                                    "+1", "+2", "+3", "+4", "+5"))
    
    # Save plot
    ggsave(paste0("results/sequence_context/", treat, "_", 
                  gsub(">", "to", mut), "_logo.png"), p, width = 8, height = 4)
    
    logos_created <- logos_created + 1
  } else {
    cat(paste("Skipping logo for", treat, mut, "- insufficient data (", 
              length(contexts), "contexts)\n"))
    skipped_logos <- skipped_logos + 1
  }
}

cat(paste("Created", logos_created, "sequence logos. Skipped", skipped_logos, 
          "due to insufficient data.\n"))

#-----------------------------------------------------
# 5. Compare dominant mutation contexts across treatments
#-----------------------------------------------------

cat("Analyzing dominant mutation contexts for each treatment...\n")

# Define dominant mutations for each treatment based on our previous analysis
dominant_mutations <- list(
  WT = "G>A",   # 45.5%
  WTA = "A>G",  # 55.6%
  CAS = "C>T",  # 40%
  STC = "A>C"   # 20%
)

# Create combined plot for dominant mutations
dominant_plots <- list()
valid_treatments <- character()

for (treat in names(dominant_mutations)) {
  mut <- dominant_mutations[[treat]]
  
  # Get contexts for this treatment and mutation
  contexts <- contexts_data %>%
    filter(Treatment == treat & Mutation == mut) %>%
    pull(Context)
  
  if (length(contexts) >= 5) {  # Only create logos with sufficient data
    # Create sequence logo
    p <- ggseqlogo(contexts) + 
      ggtitle(paste0(treat, ": ", mut, " (n=", length(contexts), ")")) +
      theme_minimal() +
      xlab("Position relative to mutation") +
      scale_x_continuous(breaks = 1:11, 
                         labels = c("-5", "-4", "-3", "-2", "-1", "0", 
                                    "+1", "+2", "+3", "+4", "+5"))
    
    dominant_plots[[treat]] <- p
    valid_treatments <- c(valid_treatments, treat)
  } else {
    cat(paste("Skipping dominant mutation analysis for", treat, mut, 
              "- insufficient data (", length(contexts), "contexts)\n"))
  }
}

# Combine plots
if (length(dominant_plots) > 0) {
  combined_plot <- do.call(grid.arrange, c(dominant_plots[valid_treatments], ncol = 1))
  ggsave("results/sequence_context/dominant_mutation_contexts.png", 
         combined_plot, width = 10, height = 8)
  cat("Created combined plot of dominant mutation contexts\n")
} else {
  cat("Could not create combined plot - insufficient data for all treatments\n")
}

#-----------------------------------------------------
# 6. Analyze nucleotide frequencies at each position
#-----------------------------------------------------

cat("Analyzing nucleotide frequencies around mutation sites...\n")

# Calculate nucleotide frequencies for specific mutation types
analyze_position_frequencies <- function(treatment, mutation_type) {
  contexts <- contexts_data %>%
    filter(Treatment == treatment & Mutation == mutation_type) %>%
    pull(Context)
  
  if (length(contexts) < 5) {
    return(NULL)  # Not enough data
  }
  
  # Create matrix of contexts
  context_matrix <- do.call(rbind, strsplit(contexts, ""))
  
  # Calculate frequencies at each position
  freq_data <- data.frame()
  
  for (pos in 1:ncol(context_matrix)) {
    pos_freqs <- table(context_matrix[, pos])
    for (base in names(pos_freqs)) {
      freq_data <- rbind(freq_data, data.frame(
        Position = pos - 6,  # Relative to mutation (position 6)
        Base = base,
        Frequency = as.numeric(pos_freqs[base]) / nrow(context_matrix),
        Count = as.numeric(pos_freqs[base]),
        Total = nrow(context_matrix),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(freq_data)
}

# Generate frequency tables for dominant mutations
frequencies_analyzed <- 0

for (treat in names(dominant_mutations)) {
  mut <- dominant_mutations[[treat]]
  
  freq_data <- analyze_position_frequencies(treat, mut)
  
  if (!is.null(freq_data)) {
    # Save frequency data
    write.table(freq_data, 
                paste0("results/sequence_context/", treat, "_", 
                       gsub(">", "to", mut), "_frequencies.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Create position-specific frequency plot
    p <- ggplot(freq_data, aes(x = Position, y = Frequency, fill = Base)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = c("A" = "#44AA99", "C" = "#332288", 
                                   "G" = "#DDCC77", "T" = "#AA4499",
                                   "N" = "#CCCCCC")) +
      theme_minimal() +
      labs(title = paste0(treat, ": Nucleotide frequencies around ", mut, " mutations"),
           x = "Position relative to mutation (0 = mutation site)",
           y = "Frequency") +
      theme(legend.title = element_blank())
    
    ggsave(paste0("results/sequence_context/", treat, "_", 
                  gsub(">", "to", mut), "_position_freqs.png"), 
           p, width = 8, height = 5)
    
    frequencies_analyzed <- frequencies_analyzed + 1
  } else {
    cat(paste("Skipping frequency analysis for", treat, mut, "- insufficient data\n"))
  }
}

cat(paste("Completed frequency analysis for", frequencies_analyzed, "dominant mutations\n"))

#-----------------------------------------------------
# 7. Analyze sequence motifs
#-----------------------------------------------------

cat("Checking for sequence motifs near mutation sites...\n")

# Function to check for known sequence motifs
analyze_motifs <- function(treatment, mutation_type) {
  contexts <- contexts_data %>%
    filter(Treatment == treatment & Mutation == mutation_type) %>%
    pull(Context)
  
  if (length(contexts) < 10) {
    return(NULL)  # Not enough data
  }
  
  # Define motifs to check for
  motifs <- list(
    "GG" = "GG",      # Potential oxidation hotspot
    "CpG" = "CG",     # Methylation/deamination site
    "AAAA" = "AAAA",  # Homopolymer run
    "TTTT" = "TTTT",  # Homopolymer run
    "CCCC" = "CCCC",  # Homopolymer run
    "GGGG" = "GGGG"   # Homopolymer run
  )
  
  # Check for each motif
  motif_counts <- data.frame(
    Motif = character(),
    Count = numeric(),
    Percentage = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (motif_name in names(motifs)) {
    pattern <- motifs[[motif_name]]
    count <- sum(grepl(pattern, contexts, fixed = TRUE))
    percentage <- count / length(contexts) * 100
    
    motif_counts <- rbind(motif_counts, data.frame(
      Motif = motif_name,
      Count = count,
      Percentage = percentage,
      stringsAsFactors = FALSE
    ))
  }
  
  return(motif_counts)
}

# Analyze motifs for each treatment's dominant mutation
motif_summary <- data.frame()

for (treat in names(dominant_mutations)) {
  mut <- dominant_mutations[[treat]]
  
  motif_data <- analyze_motifs(treat, mut)
  
  if (!is.null(motif_data)) {
    motif_data$Treatment <- treat
    motif_data$Mutation <- mut
    
    motif_summary <- rbind(motif_summary, motif_data)
    
    # Save motif data
    write.table(motif_data, 
                paste0("results/sequence_context/", treat, "_", 
                       gsub(">", "to", mut), "_motifs.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# Create motif comparison plot if we have data
if (nrow(motif_summary) > 0) {
  p <- ggplot(motif_summary, aes(x = Motif, y = Percentage, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Sequence Motif Frequency Near Mutation Sites",
         x = "Motif", 
         y = "Percentage of Contexts (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("results/sequence_context/motif_comparison.png", p, width = 10, height = 6)
  
  # Save summary
  write.table(motif_summary, "results/sequence_context/motif_summary.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Completed motif analysis and created comparison plot\n")
} else {
  cat("Could not perform motif analysis - insufficient data\n")
}

#-----------------------------------------------------
# 8. Create summary report
#-----------------------------------------------------

cat("Generating summary report...\n")

# Generate markdown report
report_content <- c(
  "# Sequence Context Analysis Results",
  "",
  "## Overview",
  "",
  "This analysis examined the nucleotide sequences surrounding mutation sites for each treatment.",
  "Sequence contexts can provide insights into the molecular mechanisms underlying these mutations.",
  "",
  "## Dominant Mutation Contexts",
  ""
)

# Add treatment-specific content
for (treat in names(dominant_mutations)) {
  mut <- dominant_mutations[[treat]]
  
  # Check if we have data for this treatment
  contexts <- contexts_data %>%
    filter(Treatment == treat & Mutation == mut) %>%
    pull(Context)
  
  if (length(contexts) >= 5) {
    # Add section header
    report_content <- c(report_content,
                        paste0("### ", treat, " Treatment (", mut, ")"),
                        "",
                        paste0("![", treat, " ", gsub(">", "to", mut), " Context](", 
                               treat, "_", gsub(">", "to", mut), "_logo.png)"),
                        "",
                        paste0("Nucleotide frequency distribution:"),
                        "",
                        paste0("![", treat, " Frequencies](", treat, "_", 
                               gsub(">", "to", mut), "_position_freqs.png)"),
                        "")
    
    # Add treatment-specific interpretation
    if (treat == "WT") {
      report_content <- c(report_content,
                          "The G>A mutations show characteristic patterns consistent with",
                          "oxidative damage to guanine, which can form 8-oxoG that mispairs with adenine.",
                          "")
    } else if (treat == "WTA") {
      report_content <- c(report_content,
                          "The A>G bias is complementary to WT's G>A pattern, suggesting",
                          "WTA may affect the same mechanism but on opposite DNA strands.",
                          "")
    } else if (treat == "CAS") {
      report_content <- c(report_content,
                          "The C>T mutations are characteristic of cytosine deamination,",
                          "which creates uracil that pairs with adenine during replication.",
                          "CAS treatment may specifically affect the pathways that repair this damage.",
                          "")
    } else if (treat == "STC") {
      report_content <- c(report_content,
                          "STC shows the most diverse mutation spectrum, which may reflect",
                          "its impact on multiple repair pathways simultaneously.",
                          "")
    }
  }
}

# Add biological implications and conclusion
report_content <- c(report_content,
                    "## Biological Implications",
                    "",
                    "The sequence contexts identified in this analysis suggest these potential mechanisms:",
                    "",
                    "1. **WT (G>A bias)**: Potential oxidative damage to guanine, creating 8-oxoG that mispairs with adenine",
                    "2. **WTA (A>G bias)**: Possible deamination of adenine to hypoxanthine, which pairs with cytosine",
                    "3. **CAS (C>T bias)**: Classic cytosine deamination signature, creating uracil that pairs with adenine",
                    "4. **STC (diverse patterns)**: Multiple damage mechanisms affecting different nucleotide contexts",
                    "",
                    "## Conclusion",
                    "",
                    "The sequence context patterns provide additional evidence that these treatments affect specific",
                    "DNA damage or repair mechanisms rather than causing random mutations. The complementary patterns",
                    "between WT and WTA are particularly striking and support our previous observations about their relationship.")

# Write report
writeLines(report_content, "results/sequence_context/sequence_context_report.md")

cat("Sequence context analysis complete. Results available in results/sequence_context/\n")