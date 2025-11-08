# ============================================================
# Residual Longevity Analysis (This Code is Conducted on Cluster)
# By VNF
# ============================================================

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(ape)
  library(readxl)
  library(dplyr)
  library(nlme)
  library(phytools)
  library(caper)
  library(geiger)
  library(writexl)
})

setwd("/Path to the directory containing your files")

# ------------------------------------------------------------
# Load custom function
# ------------------------------------------------------------
source("pglsSEyPagel.R")

# ------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------
tre <- read.tree("Pruned_mammals_longevity_body_mass_large_tree.nwk")
mydata <- read_excel("Curated_combined_MPA.xlsx")

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------
mydata <- as.data.frame(mydata)
mydata$Species <- gsub(" ", "_", mydata$Species)
tre$tip.label <- gsub(" ", "_", tre$tip.label)

# ------------------------------------------------------------
# SE vector and log transforms
# ------------------------------------------------------------
mydata$SE <- 1
rownames(mydata) <- mydata$Species
SE <- setNames(mydata$SE, mydata$Species)[rownames(mydata)]
mydata$logtenBodyMass <- log10(mydata$BodyMass_g)
mydata$logtenLongevity <- log10(mydata$MaxLongevity_yrs)

# ------------------------------------------------------------
# Run PGLS model
# ------------------------------------------------------------
cat("Running PGLS for", nrow(mydata), "species...\n")
start <- Sys.time()

modelRes <- pglsSEyPagel(
  logtenLongevity ~ logtenBodyMass,
  data = mydata,
  tree = tre,
  se = SE,
  method = "ML"
)

end <- Sys.time()
cat("Model completed in", round(as.numeric(end - start, units = "mins"), 2), "minutes.\n")

# ------------------------------------------------------------
# Add residuals and LQ
# ------------------------------------------------------------
mydata$residLongevity <- residuals(modelRes, type = "response")
mydata$LQ <- 10^(mydata$residLongevity)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

# 1. Save updated dataset (same Excel file name)
write_xlsx(mydata, "Curated_combined_MPA_with_residuals.xlsx")

# 2. Backup residuals only
backup <- mydata[, c("Species", "residLongevity", "LQ")]
write_xlsx(backup, "backup_residuals.xlsx")
write.table(backup, "backup_residuals.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 3. Save model object and parameters
saveRDS(modelRes, "PGLS_residual_longevity_model.rds")

sink("PGLS_model_summary.txt")
cat("=== PGLS Residual Longevity Model Summary ===\n")
print(summary(modelRes))
cat("\nEstimated Pagel lambda:\n")
print(coef(modelRes$modelStruct$corStruct, unconstrained = FALSE))
cat("\nComputation time (minutes):", round(as.numeric(end - start, units = "mins"), 2), "\n")
sink()

cat("All results saved\n")

