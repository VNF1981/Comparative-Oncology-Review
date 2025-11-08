# ============================================================================
# Phylogenetic Tree with Cancer Data Visualization
# ============================================================================

# Clear workspace
rm(list=ls(all=TRUE))
ls()

# Load required packages
library(ggtree)
library(ggplot2)
library(dplyr)
library(treeio)
library(ape)
library(readxl)
library(ggtreeExtra)
library(cowplot)
library(scales)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

setwd("C:/Path to files' directory")

# Read tree
tree <- read.tree("Vertebrates_tree.nwk")

# Read data
data <- read_excel("Cancer_Data.xlsx")

# Ensure Species names match tree tip labels (underscores instead of spaces)
data$Species <- gsub(" ", "_", data$Species)
data$Species_display <- gsub("_", " ", data$Species)

# Handle missing values
data$NeoplasiaPrevalence[is.na(data$NeoplasiaPrevalence)] <- 0
data$MalignancyPrevalence[is.na(data$MalignancyPrevalence)] <- 0
# Keep NA for CMR (non-mammals)

# ============================================================================
# 2. COLOR PALETTES
# ============================================================================

class_colors <- c(
  "Mammalia" = "#E495A5", 
  "Aves" = "#39BEB1",
  "Reptilia" = "lightsalmon3",      
  "Amphibia" = "#ACA4E2"
)

neoplasia_color <- "grey70"
malignancy_color <- "grey40"
cmr_color <- "grey20"

bar_colors <- list(
  Neoplasia = neoplasia_color,
  Malignancy = malignancy_color,
  CMR = cmr_color
)

# ============================================================================
# 3. GROUP TREE BY CLASS (color entire branches)
# ============================================================================

groups <- split(data$Species, data$Class)
tree <- groupOTU(tree, groups)

# ============================================================================
# 4. CREATE CIRCULAR PHYLOGENY
# ============================================================================

p_circular <- ggtree(tree, layout = "fan", open.angle = 10, size = 0.3, color = "gray40") %<+% data

# Dummy data for legend
dummy_lines <- data.frame(
  x = 0, xend = 1,
  y = seq_along(class_colors),
  yend = seq_along(class_colors),
  Class = names(class_colors)
)

# Add dummy segment for class legend
p_circular <- p_circular +
  geom_segment(
    data = dummy_lines,
    aes(x = x, xend = xend, y = y, yend = yend, color = Class),
    inherit.aes = FALSE,
    show.legend = TRUE,
    linewidth = 1.5
  ) +
  scale_color_manual(
    values = class_colors,
    name = "Vertebrate Class",
    guide = guide_legend(
      override.aes = list(
        linetype = "solid",
        linewidth = 1.5
      )
    )
  )

# Add colored labels, hide label legend
p_circular <- p_circular +
  geom_tiplab(
    aes(label = Species_display, color = Class),
    size = 1.2,
    offset = 1.5,
    align = FALSE,
    fontface = "italic",
    show.legend = FALSE
  ) +
  theme(
    legend.position = c(0.2, 0.8),
    legend.title = element_text(size = 12, face = "bold", color = "gray30"),
    legend.text = element_text(size = 10, color = "gray30", face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# ============================================================================
# 5. ADD BAR LAYERS
# ============================================================================

# Neoplasia (inner)
p_circular <- p_circular +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = Species, x = NeoplasiaPrevalence),
    stat = "identity",
    orientation = "y",
    fill = neoplasia_color,
    alpha = 0.85,
    offset = 0.14,
    pwidth = 0.10
  )

# Malignancy (middle)
p_circular <- p_circular +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = Species, x = MalignancyPrevalence),
    stat = "identity",
    orientation = "y",
    fill = malignancy_color,
    alpha = 0.85,
    offset = 0,
    pwidth = 0.10
  )

# CMR (outermost)
p_circular <- p_circular +
  geom_fruit(
    geom = geom_bar,
    mapping = aes(y = Species, x = CMR),
    stat = "identity",
    orientation = "y",
    fill = cmr_color,
    alpha = 0.85,
    offset = 0,
    pwidth = 0.10,
    na.rm = TRUE
  )

# ============================================================================
# 6. MINI LEGEND FOR BAR VARIABLES (perfectly aligned, normalized)
# ============================================================================

# prep data for legend
legend_data <- data.frame(
  variable = c("Neoplasia", "Malignancy", "CMR"),
  max_value = c(
    max(data$NeoplasiaPrevalence,  na.rm = TRUE),
    max(data$MalignancyPrevalence, na.rm = TRUE),
    max(data$CMR,                  na.rm = TRUE)
  )
)

# compute global maximum across all three measures
global_max <- max(
  data$NeoplasiaPrevalence,
  data$MalignancyPrevalence,
  data$CMR,
  na.rm = TRUE
)

# function builds one compact row (no magic constants)
make_bar <- function(var_name, color, max_val, global_max) {
  scale_factor <- max_val / global_max
  
  ggplot() +
    # bar (scaled to its own max)
    geom_segment(aes(x = 0, xend = scale_factor, y = 0.6, yend = 0.6),
                 color = color, linewidth = 2, lineend = "butt") +
    
    # axis line (matches bar length)
    geom_segment(aes(x = 0, xend = scale_factor, y = 0.45, yend = 0.45),
                 color = "gray40", linewidth = 0.6) +
    
    # vertical ticks at both ends
    geom_segment(aes(x = 0, xend = 0, y = 0.45, yend = 0.40),
                 color = "gray40", linewidth = 0.6) +
    geom_segment(aes(x = scale_factor, xend = scale_factor, y = 0.45, yend = 0.40),
                 color = "gray40", linewidth = 0.6) +
    
    # tick labels (use actual max)
    annotate("text", x = 0, y = 0.35, label = "0.00",
             color = "gray30", size = 2.3, hjust = 0.5, vjust = 1) +
    annotate("text", x = scale_factor, y = 0.35, label = sprintf("%.2f", max_val),
             color = "gray30", size = 2.3, hjust = 1, vjust = 1) +
    
    # left label
    annotate("text", x = -0.05, y = 0.6, label = var_name,
             hjust = 1, vjust = 0.5, fontface = "bold",
             color = "gray30", size = 3) +
    
    coord_cartesian(xlim = c(-0.1, 1), ylim = c(0.2, 0.8), clip = "off") +
    theme_void()
}

# build each bar (pass global_max as argument)
p1 <- make_bar("Neoplasia",  bar_colors$Neoplasia,  max(data$NeoplasiaPrevalence,  na.rm = TRUE), global_max)
p2 <- make_bar("Malignancy", bar_colors$Malignancy, max(data$MalignancyPrevalence, na.rm = TRUE), global_max)
p3 <- make_bar("Cancer Motality Rate",        bar_colors$CMR,        max(data$CMR,                  na.rm = TRUE), global_max)

# stack with consistent spacing
bar_legend <- plot_grid(p1, NULL, p2, NULL, p3,
                        ncol = 1, rel_heights = c(1, 0.25, 1, 0.25, 1),
                        align = "v")

# ============================================================================
# 7. SAVE AND DISPLAY
# ============================================================================

final_plot <- ggdraw() +
  draw_plot(p_circular, 0, 0, 1, 1) +
  # Align mini legend roughly with the main legend on the left
  draw_plot(bar_legend, 0.8, 0.75, 0.10, 0.10)

print(final_plot)

cairo_pdf("Figure 1.pdf", width = 20, height = 20)
print(final_plot)
dev.off()
