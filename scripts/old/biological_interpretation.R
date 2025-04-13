# biological_interpretation.R

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggforce)
library(gridExtra)

# Create output directory
dir.create("results/visualizations/interpretation", recursive = TRUE, showWarnings = FALSE)

# Set up canvas for the interpretation diagram
png("results/visualizations/interpretation/biological_interpretation.png", 
    width = 1200, height = 800, res = 120)

# Create a 2x2 layout
layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE))

# 1. Top left: Mutation types diagram
par(mar=c(1,1,3,1))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="DNA Mutation Types", axes=FALSE)

# Draw DNA base pairs
text(2, 9, "Normal DNA", cex=1.2, font=2)
text(1, 8, "A", cex=1.5, col="red")
text(1, 7, "T", cex=1.5, col="blue")
text(1, 6, "G", cex=1.5, col="green")
text(1, 5, "C", cex=1.5, col="purple")
text(3, 8, "T", cex=1.5, col="blue")
text(3, 7, "A", cex=1.5, col="red")
text(3, 6, "C", cex=1.5, col="purple")
text(3, 5, "G", cex=1.5, col="green")

# Draw arrows
arrows(4, 7, 6, 7, lwd=2)

# Draw mutated DNA
text(8, 9, "Mutated DNA", cex=1.2, font=2)
text(7, 8, "A", cex=1.5, col="red")
text(7, 7, "T", cex=1.5, col="blue")
text(7, 6, "G", cex=1.5, col="green", font=2)
text(7, 5, "C", cex=1.5, col="purple")
text(9, 8, "T", cex=1.5, col="blue")
text(9, 7, "A", cex=1.5, col="red")
text(9, 6, "A", cex=1.5, col="red", font=2)  # Mutation
text(9, 5, "G", cex=1.5, col="green")

# Add legend
text(5, 3, "Transitions (common in our data):", cex=1.1)
text(5, 2, "A↔G  or  C↔T", cex=1, font=2)
text(5, 1, "More likely due to replication errors", cex=0.9)

# 2. Top right: Clustered vs scattered mutations
par(mar=c(1,1,3,1))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="Mutation Distribution Patterns", axes=FALSE)

# Draw scattered pattern
text(2.5, 9, "Random Mutations", cex=1.2, font=2)
set.seed(123)
points(runif(15, 0.5, 4.5), runif(15, 4, 8), pch=16, col="red", cex=1)

# Draw clustered pattern
text(7.5, 9, "Clustered Mutations", cex=1.2, font=2)
set.seed(123)
cluster1_x <- rnorm(10, 7, 0.3)
cluster1_y <- rnorm(10, 7, 0.3)
cluster2_x <- rnorm(5, 8, 0.2) 
cluster2_y <- rnorm(5, 5, 0.2)
points(cluster1_x, cluster1_y, pch=16, col="red", cex=1)
points(cluster2_x, cluster2_y, pch=16, col="red", cex=1)

# Add explanation
text(5, 3, "Our findings:", cex=1.1, font=2)
text(5, 2, "71-84% of mutations occur in clusters", cex=0.9)
text(5, 1, "Suggests specific genomic targets", cex=0.9)

# 3. Bottom left: Treatment-specific effects
par(mar=c(1,1,3,1))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="Treatment-Specific Effects", axes=FALSE)

# CAS specifics
text(1, 9, "CAS:", cex=1.2, font=2, col="darkblue")
text(1, 8.3, "• 100% transitions", cex=0.9, adj=0)
text(1, 7.6, "• 84% clustered", cex=0.9, adj=0)
text(1, 6.9, "• Dominant: T→C & C→T", cex=0.9, adj=0)

# STC specifics
text(6, 9, "STC:", cex=1.2, font=2, col="darkgreen")
text(6, 8.3, "• Most diverse mutations", cex=0.9, adj=0)
text(6, 7.6, "• Only 71% clustered", cex=0.9, adj=0)
text(6, 6.9, "• 40% transversions", cex=0.9, adj=0)

# WT specifics
text(1, 5, "WT:", cex=1.2, font=2, col="darkred")
text(1, 4.3, "• Strong G→A bias", cex=0.9, adj=0)
text(1, 3.6, "• 83% clustered", cex=0.9, adj=0)
text(1, 2.9, "• 77% transitions", cex=0.9, adj=0)

# WTA specifics
text(6, 5, "WTA:", cex=1.2, font=2, col="purple")
text(6, 4.3, "• Strong A→G bias", cex=0.9, adj=0)
text(6, 3.6, "• 82% clustered", cex=0.9, adj=0)
text(6, 2.9, "• 78% transitions", cex=0.9, adj=0)

# Add central explanation
text(5, 1.5, "Each treatment has a distinct mutation signature", cex=1, font=2)

# 4. Bottom right: Biological implications
par(mar=c(1,1,3,1))
plot(0, 0, type="n", xlim=c(0,10), ylim=c(0,10), 
     xlab="", ylab="", main="Biological Implications", axes=FALSE)

# Draw circular DNA with hotspots
theta <- seq(0, 2*pi, length.out=100)
x <- 5 + 3.5*cos(theta)
y <- 5 + 3.5*sin(theta)
lines(x, y, lwd=2)

# Draw hotspot regions
arc1 <- seq(0.2*pi, 0.5*pi, length.out=20)
x1 <- 5 + 3.5*cos(arc1)
y1 <- 5 + 3.5*sin(arc1)
lines(x1, y1, lwd=4, col="red")
text(6.5, 7.5, "JRIU01000031.1", cex=0.9)
text(6.5, 7, "(All treatments)", cex=0.8)

arc2 <- seq(0.9*pi, 1.1*pi, length.out=20)
x2 <- 5 + 3.5*cos(arc2)
y2 <- 5 + 3.5*sin(arc2)
lines(x2, y2, lwd=4, col="blue")
text(2, 5, "Treatment-", cex=0.9)
text(2, 4.5, "specific", cex=0.9)
text(2, 4, "hotspots", cex=0.9)

# Add key conclusions
text(5, 1.5, "Key implications:", cex=1.1, font=2)
text(5, 0.8, "Treatments affect specific DNA repair or replication pathways", cex=0.9)

# Close the PNG device
dev.off()

cat("Biological interpretation diagram created at results/visualizations/interpretation/biological_interpretation.png\n")