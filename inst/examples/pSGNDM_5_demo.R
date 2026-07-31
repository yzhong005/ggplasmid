library(ggplasmidZY)

gbk_file <- system.file("extdata", "pSGNDM_5.gbk", package = "ggplasmidZY")

gbk_data <- read_plasmid_annotation(gbk = gbk_file)
if (interactive()) {
  View(gbk_data)
}

p_circular <- ggplasmid(
  annotation = gbk_data,
  name = "pSGNDM-5",
  layout = "circular",
  palette = "npg",
  label_exclude_categories = "Other functions",
  legend_position = "right",
  max_labels = 24,
  label_wrap_width = 18,
  label_text_size = 3.0,
  label_text_colour = "black",
  label_line_colour = "grey70"
)
p_circular

p_linear <- ggplasmid(
  annotation = gbk_data,
  name = "pSGNDM-5",
  layout = "linear",
  palette = "npg",
  plot_line_num = 4,
  max_labels = 36,
  linear_label_wrap_width = 18,
  linear_label_max_lines = 2,
  linear_row_spacing = 4.0,
  linear_label_allow_gene_line_crossing = FALSE,
  label_text_angle = 0,
  gc_skew_height = 0.28,
  gc_content_height = 0.14,
  gc_content_linewidth = 0.8,
  gc_legend_columns = 1,
  label_text_colour = "black",
  label_line_colour = "grey70",
  legend_position = "right"
)
p_linear

# The same four-line map with diagonal connector lines.
p_linear_45 <- ggplasmid(
  annotation = gbk_data,
  name = "pSGNDM-5",
  layout = "linear",
  palette = "npg",
  plot_line_num = 4,
  max_labels = 36,
  linear_label_wrap_width = 18,
  linear_label_max_lines = 2,
  linear_row_spacing = 4.0,
  linear_label_allow_gene_line_crossing = FALSE,
  label_line_angle = 45,
  label_text_angle = 0,
  gc_skew_height = 0.28,
  gc_content_height = 0.14,
  gc_content_linewidth = 0.8,
  gc_legend_columns = 1,
  label_text_colour = "black",
  label_line_colour = "grey70",
  legend_position = "right"
)
p_linear_45

# Data used by either plot are also available from the plot object.
used_data <- plasmid_data(p_circular)
names(used_data)
used_data$fasta[c("name", "description", "length")]
