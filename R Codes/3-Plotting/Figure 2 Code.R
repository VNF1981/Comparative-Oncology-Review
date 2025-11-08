# ============================================================================
# SIMPLE ORDER-LEVEL RECTANGULAR TREE + INTEGRATED HORIZONTAL BOXPLOT PANEL
# WITH CUSTOMIZABLE VERTICAL REFERENCE LINE AT X=0
# By Vahid N Fard
# ============================================================================

# Clear workspace
rm(list = ls(all = TRUE))

# Load packages
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggtreeExtra)
  library(gridExtra)
  library(ggrepel)
  library(gridExtra)
  library(gtable)
  library(grid)
  library(patchwork)
})

# ============================================================================
# 1. READ INPUTS
# ============================================================================

setwd("C:/Path to the directory containing your files/")


data <- read_excel("Curated_combined_MPA_with_residuals.xlsx")
tree <- read.tree("Pruned_mammals_longevity_body_mass_large_tree.nwk")

# Standardize species names
data$Species <- gsub(" ", "_", data$Species)

# Replace the outdated order name
data$Order[data$Order == "Soricomorpha"] <- "Eulipotyphla"

# ============================================================================
# 2. COMPUTE ORDER-LEVEL AVERAGES (ALL ORDERS - NO FILTERING)
# ============================================================================

order_means <- data %>%
  group_by(Order) %>%
  summarise(
    mean_resid_longevity = mean(residLongevity, na.rm = TRUE),
    n_species = n()
  )

# ============================================================================
# 3. BUILD A SIMPLE TREE WHERE TIPS = ORDERS
# ============================================================================

# Replace species labels with their Orders
tree_order <- tree
species_to_order <- data %>% select(Species, Order)
tree_order$tip.label <- species_to_order$Order[match(tree_order$tip.label, species_to_order$Species)]

# Drop duplicates (keep one representative per order)
tree_order <- drop.tip(tree_order, which(duplicated(tree_order$tip.label)))

# Keep only Orders that exist in order_means (now includes all orders)
tree_order <- keep.tip(tree_order, intersect(tree_order$tip.label, order_means$Order))
order_means <- order_means[match(tree_order$tip.label, order_means$Order), ]

# Store the tree tip order for consistency
tree_tip_order <- tree_order$tip.label

# ============================================================================
# 4. USER-DEFINED COLOR PALETTE & FONT SIZES
# ============================================================================

# # Define custom color palette for Orders that will be PLOTTED (by mean_resid_longevity)
# order_colors <- c(
#   "Rodentia" = "#E41A1C",           
#   "Carnivora" = "#8181FF",          
#   "Primates" = "#4DAF4A",
#   "Chiroptera" = "#984EA3",
#   "Artiodactyla" = "#FF7F00",
#   "Diprotodontia" = "#A65628",
#   "Cetacea" = "#F781BF",
#   "Lagomorpha" = "#999999",
#   "Eulipotyphla" = "#66C2A5",
#   "Didelphimorphia" = "#FC8D62",
#   "Dasyuromorphia" = "#8DA0CB",
#   "Scandentia" = "#c4feb9"
#   # Add/remove orders as needed - colors will apply only to those in order_means
# )

# USER OPTIONS FOR SPECIES TEXT AND FONTS
species_text_size_tree    <- 2.5    # font size for species names in tree columns
species_text_size_boxplot <- 2.5      # font size for species names in boxplot
species_text_color        <- "grey20"  # color for species names text

tree_tip_label_size       <- 3      # font size for order names + mean values in tree
tree_branch_size          <- 0.5    # thickness of tree branches
boxplot_point_size        <- 1      # size of jittered points in boxplot
boxplot_label_color       <- "darkred"  # color for labeled species in boxplot

# Column header options for tree columns
tree_column_header_font_size  <- 3   # font size for column headers
tree_column_header_color      <- "grey20"  # color for column headers
tree_column_header_face       <- "bold"  # font face for headers
tree_column_header_bg_color   <- "#FFFF99"  # highlight color (yellow like PDF highlighter)
tree_column_header_bg_alpha   <- 0.6   # transparency (0-1)

# ============================================================================
# TREE TIP LABEL COLOR OPTION
# ============================================================================
# Choose one: Set to TRUE to use order colors, FALSE to use stable color
use_tree_tip_label_color <- FALSE

# If using stable color, define it here:
stable_tree_tip_label_color <- "#8b0000"  # dark red

# ============================================================================
# VERTICAL REFERENCE LINE CUSTOMIZATION (NEW)
# ============================================================================

ref_line_color       <- "#8b0000"        # dark red
ref_line_thickness   <- 0.8              # line width (0.5-2 recommended)
ref_line_type        <- "twodash"          # "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"
ref_line_alpha       <- 0.4             # transparency: 0 = invisible, 1 = opaque (75% transparent = 0.25)

# ============================================================================
# 5. IDENTIFY TOP 3 AND BOTTOM 3 SPECIES PER ORDER
# ============================================================================

# FILTER: Only orders in order_means
top_bottom <- data %>%
  filter(Order %in% order_means$Order) %>%
  group_by(Order) %>%
  arrange(desc(residLongevity)) %>%
  summarise(
    top3    = paste(paste0(head(Species, 3), " (", round(head(residLongevity, 3), 2), ")"), collapse = "\n"),
    bottom3 = paste(paste0(tail(Species, 3), " (", round(tail(residLongevity, 3), 2), ")"), collapse = "\n")
  )

# Merge with your order_means (to align with the same orders)
order_table <- left_join(order_means, top_bottom, by = "Order")

# Ensure order_colors only contains orders that will be plotted
order_colors_filtered <- order_colors[names(order_colors) %in% order_table$Order]

# Auto-generate colors for any missing orders (not in order_colors dictionary)
missing_orders <- setdiff(order_table$Order, names(order_colors_filtered))
if (length(missing_orders) > 0) {
  # Create automatic colors using a nice palette
  auto_colors <- rainbow(length(missing_orders), s = 0.7, v = 0.8)
  names(auto_colors) <- missing_orders
  order_colors_filtered <- c(order_colors_filtered, auto_colors)
}

# Ensure colors are in the same order as order_table for consistency
order_colors_filtered <- order_colors_filtered[order_table$Order]

# ============================================================================
# CREATE CONDITIONAL COLOR MAPPING FOR TREE TIP LABELS
# ============================================================================

if (use_tree_tip_label_color) {
  # Use the order colors (current behavior)
  tree_tip_label_colors <- order_colors_filtered
} else {
  # Use stable color for all tree tip labels
  tree_tip_label_colors <- setNames(rep(stable_tree_tip_label_color, length(order_colors_filtered)), 
                                    names(order_colors_filtered))
}

# ============================================================================
# 6. BUILD THE TREE WITH FIXED X-AXIS ALIGNMENT
# ============================================================================

# Get the tree's actual branching times/positions
tree_max_height <- max(branching.times(tree_order))

p_rect <- ggtree(tree_order, layout = "rectangular", size = tree_branch_size, color = "grey30")

p_rect <- p_rect %<+% order_table +
  geom_tiplab(
    aes(label = paste0(label, " (", n_species, ")"),
        color = label),
    size = tree_tip_label_size,
    hjust = -0.05,
    fontface = "bold"
  ) +
  scale_color_manual(values = tree_tip_label_colors, guide = "none") +
  theme_tree2() +
  coord_cartesian(clip = "off") +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(20, 100, 20, 20),
    text = element_text(size = 12)
  )

# ============================================================================
# FIXED X-AXIS: INDEPENDENT FROM TREE (COMMENTED OUT as we have orders at tips)
# ============================================================================

# # Manually set axis range (independent from tree data)
# x_axis_min <- 0
# x_axis_max <- 200  # Adjust this value to control overall axis extent
# 
# # Define manual ticks: at 0, 1/3, 2/3, and tree_max_height (end point)
# x_axis_breaks <- c(0, x_axis_max/3, 2*x_axis_max/3, tree_max_height)
# x_axis_labels <- c("180", "120", "60", "0")
# 
# # Add manual x-axis line using annotation_custom
# # This creates an independent axis line that extends exactly to tree_max_height
# p_rect <- p_rect +
#   annotate("segment", x = x_axis_min, xend = tree_max_height, 
#            y = -Inf, yend = -Inf, color = "grey30", linewidth = 0.5) +
#   coord_cartesian(clip = "off", xlim = c(x_axis_min, x_axis_max)) +
#   theme_tree2() +
#   theme(
#     plot.margin = margin(10, 100, 10, 10),
#     axis.line.x = element_blank(),
#     axis.ticks.x = element_line(color = "grey30", linewidth = 0.5),
#     axis.text.x = element_text(size = 10, color = "grey30"),
#     axis.title.x = element_text(size = 11, face = "bold", color = "grey30")
#   ) +
#   scale_x_continuous(
#     breaks = x_axis_breaks,
#     labels = x_axis_labels,
#     limits = c(x_axis_min, x_axis_max),
#     expand = c(0, 0),
#     position = "bottom"
#   ) +
#   labs(x = "Millions of Years Ago")

# ============================================================================
# FUNCTION: ADD PANEL LABEL
# ============================================================================

# Function to add panel label with adjustable position (i.e., (a) and (b))
# x_position: 0 = left, 1 = right (default 0.02 is near left)
# y_position: 0 = bottom, 1 = top (default 0.98 is near top)
add_panel_label <- function(p, label, x_position = 0, y_position = 1) {
  p + annotation_custom(
    grob = textGrob(label, x = unit(x_position, "npc"), y = unit(y_position, "npc"),
                    hjust = 0, vjust = 1,
                    gp = gpar(fontsize = 12, fontface = "bold")),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
}

p_rect_labeled <- add_panel_label(p_rect, "(a)", y_position = 1.02)

# ============================================================================
# 7. CREATE HORIZONTAL BOXPLOT MATCHED TO TREE ORDER
# ============================================================================

# Extract the plotted order of tips from your tree plot (bottom to top)
tip_order <- p_rect$data %>%
  dplyr::filter(isTip) %>%
  dplyr::arrange(y) %>%
  dplyr::pull(label)

# Filter data to only include orders that are in the tree
data_filtered <- data %>%
  filter(Order %in% order_means$Order)

# Convert Order to factor with levels matching the exact tree display order (bottom to top)
data_filtered$Order <- factor(data_filtered$Order, 
                              levels = tip_order,
                              ordered = TRUE)

# Identify top 3 and bottom 3 per order for labeling
species_to_label <- data_filtered %>%
  group_by(Order) %>%
  arrange(desc(residLongevity)) %>%
  mutate(
    rank_desc = row_number(),
    rank_asc = n() - row_number() + 1,
    is_top3 = rank_desc <= 1,
    is_bottom3 = rank_asc <= 1,
    label_this = is_top3 | is_bottom3
  ) %>%
  filter(label_this) %>%
  ungroup() %>%
  mutate(
    # Determine if label should be above or below
    label_position = if_else(rank_desc <= 1, "above", "below")
  )

# Create HORIZONTAL boxplot with ORDER matching tree
p_boxplot_horizontal <- ggplot(data_filtered, aes(
  x = residLongevity,
  y = Order
)) +
  # ADD VERTICAL REFERENCE LINE AT X=0 WITH CUSTOMIZABLE PROPERTIES
  geom_vline(
    xintercept = 0,
    color = ref_line_color,
    linetype = ref_line_type,
    linewidth = ref_line_thickness,
    alpha = ref_line_alpha
  ) +
  geom_boxplot(
    color = "grey30",
    fill = "grey50",
    alpha = 0.7,
    outlier.shape = NA,
    linewidth = 0.4,
    width = 0.6,
  ) +
  stat_boxplot(geom = "errorbar", width = 0.5, linewidth = 0.4, color = "grey30") +
  geom_point(
    data = species_to_label,
    aes(x = residLongevity, y = Order),
    color = "grey30",
    size = boxplot_point_size + 0.3,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = species_to_label,
    aes(label = Species, x = residLongevity),
    color = "grey30",
    size = species_text_size_boxplot,
    fontface = "bold",
    segment.color = NA,
    box.padding = unit(0.2, "lines"),
    point.padding = unit(0.8, "lines"),
    force = 2,
    force_pull = 0,
    direction = "x",
    max.overlaps = Inf,
    show.legend = FALSE,
    nudge_x = ifelse(species_to_label$rank_desc <= 1, 0.5, -0.5)
  ) +
  scale_x_continuous(
    limits = c(-3, 2),          # x-axis starts at -3 and ends at 2
    breaks = seq(-3, 2, 1),     # show ticks at -3, -2, -1, 0, 1, 2
    expand = c(0, 0)            # remove any extra padding space
  ) +
  labs(y = "Order", x = "Residual Longevity") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "grey30"),
    axis.title.x = element_text(size = 11, face = "bold", color = "grey30"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm")
  )

p_boxplot_labeled <- add_panel_label(p_boxplot_horizontal, "(b)", y_position = 1.02)

# ============================================================================
# 8. COMBINE PLOTS SIDE-BY-SIDE USING PATCHWORK
# ============================================================================

combined_plot <- p_rect_labeled + p_boxplot_labeled + 
  plot_layout(ncol = 2, widths = c(0.6, 1.4))

print(combined_plot)

# ============================================================================
# 9. SAVE TO PDF
# ============================================================================

ggsave(
  filename = "Figure 2.pdf",
  plot = combined_plot,
  width = 12,
  height = 8,
  units = "in",
  device = cairo_pdf
)