# ============================================================================
# Mammalian Longevity Phylogeny with Ancestral State Reconstruction
# Author: Vahid N Fard
# Date: 2025
# ============================================================================

# Clear workspace
rm(list = ls(all = TRUE))
ls()

# Load required packages
suppressPackageStartupMessages({
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(treeio)
  library(ape)
  library(readxl)
  library(ggtreeExtra)
  library(phytools)
  library(ggnewscale)
})

# Set working directory
setwd("C:/Path to Files' Directory")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

tree <- read.tree("Pruned_mammals_longevity_body_mass_large_tree.nwk")

# Adjust selected tip branch lengths manually
tips_to_fix <- c("Cavia_tschudii", "Aotus_azarai", "Galago_senegalensis", "Ursus_thibetanus")
tip_index <- match(tips_to_fix, tree$tip.label)
tree$edge.length[tree$edge[,2] == tip_index[1]] <- tree$edge.length[tree$edge[,2] == tip_index[1]] * 0.68
tree$edge.length[tree$edge[,2] == tip_index[2]] <- tree$edge.length[tree$edge[,2] == tip_index[2]] * 0.75
tree$edge.length[tree$edge[,2] == tip_index[3]] <- tree$edge.length[tree$edge[,2] == tip_index[3]] * 0.7
tree$edge.length[tree$edge[,2] == tip_index[4]] <- tree$edge.length[tree$edge[,2] == tip_index[4]] * 1.75

# ============================================================================
# Branch length transformation: gradient root compression
# ============================================================================

# Function to compress only the deepest few root branches
compress_basal_branches <- function(tree, n_levels = 1, factor = 0.3) {
  tree_plot <- tree
  
  # Get parent and child nodes
  parents <- tree_plot$edge[,1]
  children <- tree_plot$edge[,2]
  
  # Find the root
  root_node <- Ntip(tree_plot) + 1
  
  # Start with branches directly under the root
  basal_edges <- which(parents == root_node)
  affected_nodes <- children[basal_edges]
  
  # Optionally go one level deeper
  if (n_levels > 1) {
    for (i in 2:n_levels) {
      next_edges <- which(parents %in% affected_nodes)
      affected_nodes <- c(affected_nodes, children[next_edges])
      basal_edges <- unique(c(basal_edges, next_edges))
    }
  }
  
  # Apply compression only to selected edges
  tree_plot$edge.length[basal_edges] <- tree_plot$edge.length[basal_edges] * factor
  
  return(tree_plot)
}

# Compress only the first 2 levels below root by 70%
tree_plot <- compress_basal_branches(tree, n_levels = 2, factor = 0.3)

# Re-ultrametrize so all tips align
tree_plot <- force.ultrametric(tree_plot, method = "extend")

# ============================================================================
# ============================================================================
data <- read_excel("Curated_combined_MPA_with_residuals.xlsx")
asr_data <- read_excel("ASR_residual_longevity_nodes_fastAnc.xlsx")

data$Species <- gsub(" ", "_", data$Species)
data$Species_display <- gsub("_", " ", data$Species)

# ============================================================================
# 2. PREPARE ASR DATA
# ============================================================================

n_tips <- length(data$Species)
n_nodes <- n_tips - 1

tip_df <- data.frame(node = 1:n_tips, Species = tree$tip.label, stringsAsFactors = FALSE)
tip_df <- left_join(tip_df, data, by = "Species")

node_df <- data.frame(
  node = as.numeric(asr_data$Node),
  residLongevity = asr_data$ResidualLongevity,
  stringsAsFactors = FALSE
)

all_nodes_df <- bind_rows(tip_df, node_df)

# ============================================================================
# 3. COLOR PALETTES
# ============================================================================

asr_color_palette <- colorRampPalette(c("#000049", "#0000ab", "#8181FF", "#dadafe", "#e2fafe",
                                        "#c4feb9", "#33CD31", "#158a00", "#083900"))(13)

body_mass_color <- "#F781BF"
longevity_color  <- "#ACA4E2"

# ============================================================================
# 4. GLOBAL LAYOUT VARIABLES
# ============================================================================

# --- Horizontal layout offsets ---
bar_base_offset <- 0.06         # distance of Residual Longevity bar from tree edge
bar_gap          <- 0         # gap between consecutive bars
bar_offset_resid <- bar_base_offset         # horizontal offset of the Residual Longevity bar (closest to the tree)
bar_offset_long  <- bar_base_offset - 0.07        # horizontal offset of the log10 Longevity bar (middle position, spaced by bar_gap)
bar_offset_mass  <- bar_base_offset - 0.07      # horizontal offset of the log10 Body Mass bar (outermost position, spaced twice bar_gap)

bar_pwidth <- 0.05               # relative width of bar blocks
bar_width  <- 0.8                # thickness of individual bars
bar_alpha  <- 1               # transparency (1 = opaque)

# --- Vertical spacing and annotation coordinates ---
ylim_lower <- -10                # space below lowest tip
ylim_upper_extra <- 7           # space above highest tip
annotate_y_bottom <- -1         # vertical position of baseline axes
annotate_y_label1 <- -3         # vertical position of axis labels (numbers)
annotate_y_title  <- n_tips + 5.5    # vertical position of top titles

# --- Title styling and positioning ---
title_size   <- 4            # overall font size for top titles
title_color  <- "grey20"      # font color
# individual horizontal shifts (fine control)
title_shift_resid <- -2.5     # shift for Residual Longevity
title_shift_long  <- 2.6     # shift for log10 Longevity
title_shift_mass  <- 7.85     # shift for log10 Body Mass

# --- Axis styling and positioning ---
axis_font_size  <- 3      # size of the axis labels (numbers)
axis_font_color <- "grey20" # color of axis text and lines
axis_font_face  <- "bold"   # "plain", "italic", or "bold"
axis_line_size  <- 0.75      # thickness of axis lines

# --- Tick styling ---
tick_length   <- 1           # tick length (in y units)
tick_size     <- 0.75           # thickness of ticks
tick_color    <- axis_font_color  # color (usually same as axis line)

# Optional manual control for segment spans (x coordinates)
# axis_tree_start  <- 0.0     # custom start for tree axis
# axis_tree_end    <- 0.1     # custom end for tree axis
axis_resid_start <- 127    # custom start for Residual Longevity axis
axis_resid_end   <- 133
axis_long_start  <- 133.75  # custom start for log10 Longevity axis
axis_long_end    <- 140
axis_mass_start  <- 140.8     # custom start for log10 Body Mass axis
axis_mass_end    <- 147


# ============================================================================
# 5. BUILD TREE WITH ASR
# ============================================================================

p_rect <- ggtree(tree_plot, layout = "rectangular", size = 1.75, color = "grey20")

p_rect <- p_rect %<+% all_nodes_df +
  geom_tree(aes(color = residLongevity), size = 1) +
  scale_color_gradientn(
    colors = asr_color_palette,
    name = "Residual\nLongevity\n(ASR)",
    na.value = "grey50",
    limits = c(-2.5, 1.5)
  ) +
  theme_tree() +
  theme(
    legend.position = c(0.7, 0.955),   # coordinate
    legend.justification = c("left", "top"),
    legend.direction = "vertical",   # vertical color bar
    legend.title = element_text(size = 16, face = "plain", color = "grey20"),
    legend.text = element_text(size = 16, color = "grey20"),
    legend.key.height = unit(2, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  ylim(ylim_lower, n_tips + ylim_upper_extra)

p_rect <- p_rect +
  geom_tiplab(
    aes(label = Species_display),
    size = 2.25,
    fontface = "bold.italic",
    color = "grey20",
    offset = -0.3 
    #vjust = -0.1             # this is the vertical fine-tune knob
  )


# ============================================================================
# 6. ADD BAR PLOTS
# ============================================================================

p_rect <- p_rect + new_scale_color()

# Residual Longevity
p_rect <- p_rect +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = Species, x = residLongevity, fill = residLongevity),
    stat = "identity",
    orientation = "y",
    alpha = 0.9,
    pwidth = bar_pwidth,
    offset = bar_offset_resid,
    width = bar_width,
    axis.params = list(axis = "none")
  ) +
  scale_fill_gradientn(colors = asr_color_palette, guide = "none")

# log10 Longevity
p_rect <- p_rect +
  new_scale_fill() +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = Species, x = logtenLongevity),
    stat = "identity",
    orientation = "y",
    fill = longevity_color,
    alpha = bar_alpha,
    pwidth = bar_pwidth,
    offset = bar_offset_long,
    width = bar_width,
    axis.params = list(axis = "none")
  )

# log10 Body Mass
p_rect <- p_rect +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = Species, x = logtenBodyMass),
    stat = "identity",
    orientation = "y",
    fill = body_mass_color,
    alpha = bar_alpha,
    pwidth = bar_pwidth,
    offset = bar_offset_mass,
    width = bar_width,
    axis.params = list(axis = "none")
  )

# ============================================================================
# 7. ADD TITLES (TOP) AND AXES (BOTTOM)
# ============================================================================

# Define numeric ranges
resid_min <- min(data$residLongevity, na.rm = TRUE)
resid_max <- max(data$residLongevity, na.rm = TRUE)
logten_long_min <- min(data$logtenLongevity, na.rm = TRUE)
logten_long_max <- max(data$logtenLongevity, na.rm = TRUE)
logten_mass_min <- min(data$logtenBodyMass, na.rm = TRUE)
logten_mass_max <- max(data$logtenBodyMass, na.rm = TRUE)

phylo_end <- max(tree$edge.length, na.rm = TRUE)
bar1_pos <- phylo_end + bar_offset_resid
bar2_pos <- phylo_end + bar_offset_long
bar3_pos <- phylo_end + bar_offset_mass

#+++++++++++++++++++++++++++++  
# Top titles
#+++++++++++++++++++++++++++++  
p_rect <- p_rect +
  annotate("text",
           x = phylo_end + bar_offset_resid + title_shift_resid,
           y = annotate_y_title,
           label = "Residual\nLongevity",
           size = title_size,
           color = title_color,
           face = "bold") +
  annotate("text",
           x = phylo_end + bar_offset_long + title_shift_long,
           y = annotate_y_title,
           label = "log10\nLongevity (yrs)",
           size = title_size,
           color = title_color,
           face = "bold") +
  annotate("text",
           x = phylo_end + bar_offset_mass + title_shift_mass,
           y = annotate_y_title,
           label = "log10\nBody Mass (g)",
           size = title_size,
           color = title_color,
           face = "bold")

#+++++++++++++++++++++++++++++
# Bottom axes
#+++++++++++++++++++++++++++++  

# These are only for labeling, not for drawing coordinates
resid_label_min <- min(data$residLongevity, na.rm = TRUE)
resid_label_max <- max(data$residLongevity, na.rm = TRUE)
long_label_min  <- min(data$logtenLongevity, na.rm = TRUE)
long_label_max  <- max(data$logtenLongevity, na.rm = TRUE)
mass_label_min  <- min(data$logtenBodyMass, na.rm = TRUE)
mass_label_max  <- max(data$logtenBodyMass, na.rm = TRUE)

# ============================================================================
# 7B. BOTTOM AXES (COORDINATES from global variables, LABELS from data ranges)
# ============================================================================

# Draw all three bottom axes
p_rect <- p_rect +
  
  # ---------------- Residual Longevity ----------------
annotate("segment",
         x = axis_resid_start, xend = axis_resid_end,
         y = annotate_y_bottom, yend = annotate_y_bottom,
         color = axis_font_color, size = axis_line_size) +
  annotate("segment",
           x = axis_resid_start, xend = axis_resid_start,
           y = annotate_y_bottom, yend = annotate_y_bottom - tick_length,
           color = tick_color, size = tick_size) +
  annotate("segment",
           x = axis_resid_end, xend = axis_resid_end,
           y = annotate_y_bottom, yend = annotate_y_bottom - tick_length,
           color = tick_color, size = tick_size) +
  annotate("text",
           x = axis_resid_start, y = annotate_y_label1,
           label = sprintf("%.1f", resid_label_min),
           size = axis_font_size, color = axis_font_color, fontface = axis_font_face) +
  annotate("text",
           x = axis_resid_end, y = annotate_y_label1,
           label = sprintf("%.1f", resid_label_max),
           size = axis_font_size, color = axis_font_color, fontface = axis_font_face) +
  
  # ---------------- log10 Longevity ----------------
annotate("segment",
         x = axis_long_start, xend = axis_long_end,
         y = annotate_y_bottom, yend = annotate_y_bottom,
         color = axis_font_color, size = axis_line_size) +
  annotate("segment",
           x = axis_long_start, xend = axis_long_start,
           y = annotate_y_bottom, yend = annotate_y_bottom - tick_length,
           color = tick_color, size = tick_size) +
  annotate("segment",
           x = axis_long_end, xend = axis_long_end,
           y = annotate_y_bottom, yend = annotate_y_bottom - tick_length,
           color = tick_color, size = tick_size) +
  annotate("text",
           x = axis_long_start, y = annotate_y_label1,
           label = sprintf("%.1f", long_label_min),
           size = axis_font_size, color = axis_font_color, fontface = axis_font_face) +
  annotate("text",
           x = axis_long_end, y = annotate_y_label1,
           label = sprintf("%.1f", long_label_max),
           size = axis_font_size, color = axis_font_color, fontface = axis_font_face) +
  
  # ---------------- log10 Body Mass ----------------
annotate("segment",
         x = axis_mass_start, xend = axis_mass_end,
         y = annotate_y_bottom, yend = annotate_y_bottom,
         color = axis_font_color, size = axis_line_size) +
  annotate("segment",
           x = axis_mass_start, xend = axis_mass_start,
           y = annotate_y_bottom, yend = annotate_y_bottom - tick_length,
           color = tick_color, size = tick_size) +
  annotate("segment",
           x = axis_mass_end, xend = axis_mass_end,
           y = annotate_y_bottom, yend = annotate_y_bottom - tick_length,
           color = tick_color, size = tick_size) +
  annotate("text",
           x = axis_mass_start, y = annotate_y_label1,
           label = sprintf("%.1f", mass_label_min),
           size = axis_font_size, color = axis_font_color, fontface = axis_font_face) +
  annotate("text",
           x = axis_mass_end, y = annotate_y_label1,
           label = sprintf("%.1f", mass_label_max),
           size = axis_font_size, color = axis_font_color, fontface = axis_font_face)


# ============================================================================
# 8. SAVE OUTPUT
# ============================================================================

print(p_rect)

cairo_pdf("Supplementary Figure 1.pdf", width = 60, height = 200)
print(p_rect)
dev.off()

