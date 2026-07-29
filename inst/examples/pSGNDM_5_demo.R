library(ggplasmidZY)

gbk_file <- system.file("extdata", "pSGNDM_5.gbk", package = "ggplasmidZY")
fasta_file <- system.file("extdata", "pSGNDM_5.fasta", package = "ggplasmidZY")

gbk_data <- read_plasmid_annotation(gbk = gbk_file)
fasta_data <- read_plasmid_fasta(fasta_file)
if (interactive()) {
  View(gbk_data)
  View(fasta_data)
}

p_circular <- ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  layout = "circular",
  label_exclude_categories = "Other functions",
  legend_position = "right",
  max_labels = 24,
  label_wrap_width = 18,
  label_text_size = 3.0,
  label_text_colour = "category",
  label_line_colour = "category"
)
p_circular

p_linear <- ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  layout = "linear",
  rows = 5,
  max_labels = 36
)
p_linear

# Data used by either plot are also available from the plot object.
used_data <- plasmid_data(p_circular)
names(used_data)
used_data$fasta[c("name", "description", "length")]
