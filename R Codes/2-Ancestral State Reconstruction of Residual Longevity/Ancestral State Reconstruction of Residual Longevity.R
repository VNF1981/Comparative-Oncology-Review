# ============================================================
#  Ancestral State Reconstruction of Residual Longevity
#  By Vahid N Fard
# ============================================================

# Clear workspace
rm(list=ls(all=TRUE))
ls()

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(viridis)
  library(readxl)
  library(dplyr)
  library(writexl)
  library(geiger)
})

# ------------------------------------------------------------
# Paths and input files
# ------------------------------------------------------------
setwd("C:/Path to directory containing your files")

tre <- read.tree("Pruned_mammals_longevity_body_mass_large_tree.nwk")
mydata <- read_excel("Curated_combined_MPA_with_residuals.xlsx")

# Standardize and match names (they are already matched)
mydata$Species <- gsub(" ", "_", mydata$Species)
tre$tip.label <- gsub(" ", "_", tre$tip.label)

# unique_tree    <- setdiff(tre$tip.label, mydata$Species)
# unique_tree
# unique_data    <- setdiff(mydata$Species, tre$tip.label)
# unique_data

# ------------------------------------------------------------
# 1. Prepare vector and run ASR
# ------------------------------------------------------------
# Align vector to tree tip order
resid_vector <- setNames(mydata$residLongevity, mydata$Species)
resid_vector <- resid_vector[tre$tip.label]

# # Double-check data >> FIne
# head(resid_vector)
# summary(resid_vector)
# sum(is.na(resid_vector))

# # Double check tre
# summary(tre$edge.length)
# sum(is.na(tre$edge.length))
# any(tre$edge.length == 0)
# is.binary(tre)
# is.rooted(tre)

# We have a branch with zero length > We need to add a small positive value to fix it 
# Replace zero-length branches with a small positive value
tre$edge.length[tre$edge.length == 0] <- 1e-6


###################################
# 1- General ASR (Brownian Motion only) > my estimated lambda is 0.94 so brownian motion works
# asr_resid <- fastAnc(tre, resid_vector, vars = TRUE, CI = TRUE)

# 2- Likelihood-based method from phytools
# We need to adjust the tree for estimated lambda first (phytools)
tre_lambda <- rescaleTree(tre, 0.94)
asr_resid <- anc.ML(tre_lambda, resid_vector)
str(asr_resid)

# 3- The slowest methof from ape
#asr_resid <- ace(resid_vector, tre, type = "continuous", method = "REML")


# ------------------------------------------------------------
# 2. Save node estimates
# ------------------------------------------------------------
# # for fastanc()
# asr_table <- data.frame(Node = names(asr_resid$ace),
#                         ResidualLongevity = asr_resid$ace,
#                         Variance = asr_resid$var)
# write_xlsx(asr_table, "ASR_residual_longevity_nodes_fastAnc.xlsx")

# for anc.ML
asr_table <- data.frame(
  Node = names(asr_resid$ace),
  ResidualLongevity = asr_resid$ace
)

# Here’s what each part means:
#   
#   $ace: the estimated ancestral states (this is what I want)
#   $sig2: estimated Brownian variance
#   $logLik: the log-likelihood of the model
#   $model: model type ("BM" here, consistent with λ = 0.94 scaling)
#   $convergence: 0 means optimization converged normally

write_xlsx(asr_table, "ASR_residual_longevity_nodes_ancML.xlsx")


# ------------------------------------------------------------
# 3. Visualize (simple branch coloring)
# ------------------------------------------------------------
contmap_resid <- contMap(tre, resid_vector, plot = FALSE)
pdf("ASR_ResidualLongevity_Tree.pdf", width = 48, height = 64)
plot(contmap_resid, legend = TRUE, fsize = 0.2)
title("Ancestral Reconstruction of Residual Longevity in Mammals")
dev.off()

cat("Saved:\n  • ASR_residual_longevity_nodes.xlsx\n  • ASR_ResidualLongevity_Tree.pdf\n")







