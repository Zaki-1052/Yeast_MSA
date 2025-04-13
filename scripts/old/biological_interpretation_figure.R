# biological_interpretation_figure.R

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(viridis)

# Create output directory
dir.create("results/visualizations/interpretation", recursive = TRUE, showWarnings = FALSE)

# Set up a multi-panel figure on a single device
png("results/visualizations/interpretation/biological_mechanisms.png", 
    width = 1600, height = 1200, res = 150)

# Create a 2x2 layout for our panels
layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE))

#-----------------------------------------
# Panel 1: DNA Repair Pathways in Yeast
#-----------------------------------------
par(mar=c(2,2,3,2))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="DNA Repair Pathways Affected by Treatments", 
     axes=FALSE, cex.main=1.5)

# Draw central DNA
rect(4, 1, 6, 9, border=NA, col=rgb(0.95, 0.95, 0.95))

# Mismatch Repair Pathway (high transitions)
arrows(3, 7, 1, 9, lwd=2, col="darkred")
text(1, 9.3, "Mismatch Repair", col="darkred", font=2)
text(0.8, 8.5, "• CAS (100% transitions)", col="darkred", cex=0.8)
text(0.8, 8.0, "• WT/WTA (77-78% transitions)", col="darkred", cex=0.8)
text(0.8, 7.5, "→ Affects base pairing fidelity", col="darkred", cex=0.8)

# Base Excision Repair (G>A, C>T bias)
arrows(3, 5, 1, 4, lwd=2, col="darkblue")
text(1, 3.7, "Base Excision Repair", col="darkblue", font=2)
text(0.8, 3.2, "• WT (G>A bias)", col="darkblue", cex=0.8)
text(0.8, 2.7, "• CAS (C>T bias)", col="darkblue", cex=0.8)
text(0.8, 2.2, "→ Affects oxidative damage repair", col="darkblue", cex=0.8)

# Nucleotide Excision Repair (transversions)
arrows(7, 7, 9, 9, lwd=2, col="darkgreen")
text(9, 9.3, "Nucleotide Excision Repair", col="darkgreen", font=2)
text(9.2, 8.5, "• STC (40% transversions)", col="darkgreen", cex=0.8, adj=1)
text(9.2, 8.0, "• Diverse mutation types", col="darkgreen", cex=0.8, adj=1)
text(9.2, 7.5, "→ Affects bulky lesion repair", col="darkgreen", cex=0.8, adj=1)

# DNA Replication (clustered mutations)
arrows(7, 5, 9, 4, lwd=2, col="purple")
text(9, 3.7, "DNA Replication", col="purple", font=2)
text(9.2, 3.2, "• High clustering (71-84%)", col="purple", cex=0.8, adj=1)
text(9.2, 2.7, "• Hotspot regions", col="purple", cex=0.8, adj=1)
text(9.2, 2.2, "→ Affects replication fidelity", col="purple", cex=0.8, adj=1)

# Draw DNA base pairs inside the rectangle
dna_y <- seq(2, 8, by=0.5)
for(i in 1:length(dna_y)) {
  if(i %% 2 == 0) {
    text(4.7, dna_y[i], "A", col="red", cex=0.8)
    text(5.3, dna_y[i], "T", col="blue", cex=0.8)
  } else {
    text(4.7, dna_y[i], "G", col="green", cex=0.8)
    text(5.3, dna_y[i], "C", col="purple", cex=0.8)
  }
}

#-----------------------------------------
# Panel 2: Treatment-Specific Mechanisms
#-----------------------------------------
par(mar=c(2,2,3,2))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="Treatment-Specific Mutation Mechanisms", 
     axes=FALSE, cex.main=1.5)

# WT Treatment
rect(0.5, 7.5, 4.5, 9.5, border="darkred", lwd=2)
text(2.5, 9.2, "WT Treatment", font=2, col="darkred")
text(2.5, 8.8, "G>A bias (45.5%)", cex=0.9)
text(2.5, 8.4, "High clustering (83.1%)", cex=0.9)
arrows(2.5, 7.5, 2.5, 6.7, lwd=2)
text(2.5, 6.5, "G deamination or oxidation", font=3, cex=0.9)
text(2.5, 6.1, "Potential mechanism:", cex=0.8)
text(2.5, 5.7, "Oxidative stress → 8-oxoG", cex=0.8)

# WTA Treatment
rect(5.5, 7.5, 9.5, 9.5, border="darkblue", lwd=2)
text(7.5, 9.2, "WTA Treatment", font=2, col="darkblue")
text(7.5, 8.8, "A>G bias (55.6%)", cex=0.9)
text(7.5, 8.4, "High clustering (81.5%)", cex=0.9)
arrows(7.5, 7.5, 7.5, 6.7, lwd=2)
text(7.5, 6.5, "A deamination or alkylation", font=3, cex=0.9)
text(7.5, 6.1, "Potential mechanism:", cex=0.8)
text(7.5, 5.7, "Adenosine deaminase activity", cex=0.8)

# STC Treatment
rect(0.5, 2.5, 4.5, 4.5, border="darkgreen", lwd=2)
text(2.5, 4.2, "STC Treatment", font=2, col="darkgreen")
text(2.5, 3.8, "Diverse mutations", cex=0.9)
text(2.5, 3.4, "Lowest clustering (70.6%)", cex=0.9)
arrows(2.5, 2.5, 2.5, 1.7, lwd=2)
text(2.5, 1.5, "Multiple repair pathway defects", font=3, cex=0.9)
text(2.5, 1.1, "Potential mechanism:", cex=0.8)
text(2.5, 0.7, "Global genome instability", cex=0.8)

# CAS Treatment
rect(5.5, 2.5, 9.5, 4.5, border="purple", lwd=2)
text(7.5, 4.2, "CAS Treatment", font=2, col="purple")
text(7.5, 3.8, "100% transitions (C>T bias)", cex=0.9)
text(7.5, 3.4, "Highest clustering (83.6%)", cex=0.9)
arrows(7.5, 2.5, 7.5, 1.7, lwd=2)
text(7.5, 1.5, "Specific mismatch repair defect", font=3, cex=0.9)
text(7.5, 1.1, "Potential mechanism:", cex=0.8)
text(7.5, 0.7, "Cytosine deamination to uracil", cex=0.8)

#-----------------------------------------
# Panel 3: Hotspot Vulnerability
#-----------------------------------------
par(mar=c(2,2,3,2))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="Genomic Hotspot Vulnerability", 
     axes=FALSE, cex.main=1.5)

# Draw chromosome with hotspot regions
rect(2, 3, 8, 4, border="black", lwd=2)
text(5, 2.5, "Scaffold JRIU01000031.1", font=2)

# Highlight the hotspot regions
rect(3, 3, 4, 4, border=NA, col="red", lwd=2)
rect(5, 3, 5.5, 4, border=NA, col="red", lwd=2)
rect(6.5, 3, 7, 4, border=NA, col="red", lwd=2)

# Add explanation for hotspot vulnerability
text(5, 9, "Why are certain regions more vulnerable?", font=2, cex=1.1)

# Reason 1: Secondary structure
arrows(3.5, 4, 3.5, 5, lwd=2)
text(3.5, 5.5, "DNA Secondary Structure", font=2, cex=0.9)
text(3.5, 6.2, "• Hairpins or cruciforms", cex=0.8)
text(3.5, 5.8, "• Challenges replication", cex=0.8)

# Reason 2: Repeat regions
arrows(5.25, 4, 5.25, 5, lwd=2)
text(5.25, 5.5, "Repeat Sequences", font=2, cex=0.9)
text(5.25, 6.2, "• Homopolymer runs", cex=0.8)
text(5.25, 5.8, "• Replication slippage", cex=0.8)

# Reason 3: Late replication
arrows(6.75, 4, 6.75, 5, lwd=2)
text(6.75, 5.5, "Replication Timing", font=2, cex=0.9)
text(6.75, 6.2, "• Late-replicating regions", cex=0.8)
text(6.75, 5.8, "• Reduced repair efficiency", cex=0.8)

# Common vulnerability across treatments
text(5, 8, "All treatments show clustering in JRIU01000031.1", cex=0.9)
text(5, 7.5, paste0("(", 
                    paste(c("WT: 83.1%", "CAS: 83.6%", "WTA: 81.5%", "STC: 70.6%"), 
                          collapse=", "), ")"), cex=0.8)

#-----------------------------------------
# Panel 4: Treatment Relationships and Biological Significance
#-----------------------------------------
par(mar=c(2,2,3,2))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="Treatment Relationships & Biological Significance", 
     axes=FALSE, cex.main=1.5)

# Plot positions for treatments
positions <- list(
  WT = c(3, 8),
  CAS = c(7, 8),
  STC = c(3, 3),
  WTA = c(7, 3)
)

# Draw nodes
treatments <- c("WT", "STC", "CAS", "WTA")
colors <- c("darkred", "darkgreen", "purple", "darkblue")
names(colors) <- treatments

for(i in 1:length(treatments)) {
  treat <- treatments[i]
  x <- positions[[treat]][1]
  y <- positions[[treat]][2]
  
  # Draw node
  symbols(x, y, circles=0.8, inches=FALSE, add=TRUE, fg=colors[treat], 
          bg=adjustcolor(colors[treat], alpha.f=0.3))
  text(x, y, treat, font=2, cex=1.2)
}

# Draw edges with similarity values
# WT-CAS
segments(positions$WT[1], positions$WT[2], positions$CAS[1], positions$CAS[2], 
         lwd=2*0.82, col=adjustcolor("black", alpha.f=0.5))
text((positions$WT[1] + positions$CAS[1])/2, (positions$WT[2] + positions$CAS[2])/2, 
     "0.82", cex=0.8)

# WT-STC
segments(positions$WT[1], positions$WT[2], positions$STC[1], positions$STC[2], 
         lwd=2*0.76, col=adjustcolor("black", alpha.f=0.5))
text((positions$WT[1] + positions$STC[1])/2 - 0.3, (positions$WT[2] + positions$STC[2])/2, 
     "0.76", cex=0.8)

# STC-WTA
segments(positions$STC[1], positions$STC[2], positions$WTA[1], positions$WTA[2], 
         lwd=2*0.77, col=adjustcolor("black", alpha.f=0.5))
text((positions$STC[1] + positions$WTA[1])/2, (positions$STC[2] + positions$WTA[2])/2, 
     "0.77", cex=0.8)

# CAS-WTA
segments(positions$CAS[1], positions$CAS[2], positions$WTA[1], positions$WTA[2], 
         lwd=2*0.65, col=adjustcolor("black", alpha.f=0.5))
text((positions$CAS[1] + positions$WTA[1])/2, (positions$CAS[2] + positions$WTA[2])/2, 
     "0.65", cex=0.8)

# Add biological significance
text(5, 1.5, "Biological Significance", font=2, cex=1.1)
text(5, 0.7, "• WT & WTA show complementary mutation patterns (G>A vs A>G)", cex=0.8)
text(5, 0.3, "• CAS's 100% transition bias indicates a specific repair defect", cex=0.8)

# Close the PNG device
dev.off()

cat("Biological mechanism interpretation figure created at: results/visualizations/biological_mechanisms.png\n")