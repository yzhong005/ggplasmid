# ggplasmid

`ggplasmid` is a plasmid-first mapping package built on `ggplot2`. It draws
linear or circular feature maps from annotation tables or GenBank files, with an
optional GC-skew track from FASTA/GenBank sequence.

The main function is `ggplasmid()`. It returns a normal `ggplot` object, so you
can add ggplot layers, themes, titles, and scales.

## Install During Development

```r
install.packages(c("ggplot2", "ggrepel"))
install.packages("devtools")
devtools::load_all("C:/Users/xinmatrix/OneDrive/R_package_development/ggplasmid_package_development")
```

## pSGNDM-5 Demo

The package includes a real 84,257 bp pSGNDM-5 example plasmid:

- BioSample: SAMN18579051
- SRA project: SRP313016
- Plasmid: pSGNDM-5

```r
library(ggplasmid)

gbk_file <- system.file("extdata", "pSGNDM_5.gbk", package = "ggplasmid")
fasta_file <- system.file("extdata", "pSGNDM_5.fasta", package = "ggplasmid")

gbk_data <- read_gbk(gbk_file)
fasta_data <- read_fasta(fasta_file)

p <- ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  layout = "circular"
)

p

View(gbk_data)
View(fasta_data)
```

`gbk_data` is the parsed GenBank feature data frame. `fasta_data` is a data
frame with the FASTA record name, description, sequence, and length. The exact
data used by an existing plot can also be retrieved with `plasmid_data(p)`.

Use more labels while exploring:

```r
ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  layout = "circular",
  label_unknown = TRUE
)
```

Draw the same plasmid as a linear map:

```r
ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  layout = "linear",
  rows = 5,
  max_labels = 36
)
```

## Your Own GenBank Input

```r
ggplasmid(
  gbk = "plasmid.gbk",
  fasta = "plasmid.fasta",
  name = "pExample",
  layout = "circular",
  output = "pExample_circular.png"
)
```

## Annotation Table Input

The table needs at least start and end coordinates. Common column names are
detected automatically: `start`, `end`/`stop`, `strand`/`frame`, `product`/
`annot`/`annotation`/`function`, `category`/`group`, and `type`/`feature`.

```r
features <- data.frame(
  start = c(100, 1700, 3300),
  end = c(900, 2700, 4200),
  strand = c("+", "-", "+"),
  product = c("replication protein", "mobilization protein", "beta-lactamase")
)

ggplasmid(features, genome_length = 5000, name = "pExample")
```

## Label Control

Circular plasmids can contain too many features for every label to be readable.
The default circular map considers all eligible non-hypothetical labels, sorts
them from shorter to longer text, and places them on radial outside lanes.
Crowded labels are moved to farther outer lanes with leader lines. Very dense
maps may still need a larger output size, smaller `label_text_size`, wider
`label_wrap_width`, or a finite `max_labels`.
`max_labels` limits only the number of labels; it does not affect how many gene
arrows are shown. Hypothetical/unknown genes still have their own grey color
category, but their text labels are omitted unless `label_unknown = TRUE`.

```r
ggplasmid(annotation = gbk_data, fasta = fasta_data,
          name = "pSGNDM-5", max_labels = Inf)
ggplasmid(annotation = gbk_data, fasta = fasta_data,
          name = "pSGNDM-5",
          label_pattern = "NDM|CTX|OXA|TEM|Tra|Rep")
ggplasmid(annotation = gbk_data, fasta = fasta_data,
          name = "pSGNDM-5",
          max_labels = 20,
          label_min_gap_deg = 8)
ggplasmid(annotation = gbk_data, fasta = fasta_data,
          name = "pSGNDM-5", show_labels = FALSE)
```

If Windows reports a font database warning, keep the default portable font or
choose an installed family:

```r
ggplasmid(annotation = gbk_data, fasta = fasta_data,
          name = "pSGNDM-5", font_family = "sans")
ggplasmid(annotation = gbk_data, fasta = fasta_data,
          name = "pSGNDM-5",
          font_family = "Times New Roman")
```

## Styling

Most plot geometry and text settings have defaults but can be tuned directly:

```r
ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  palette = "npg",
  gene_highlight = gene_highlight(
    "Antimicrobial resistance" = "#B2182B",
    "Replication" = "#2166AC"
  ),
  gene_border_linewidth = 0.35,
  gene_arrow_head_fraction = 0.45,
  inner_radius = 0.38,
  gc_skew_radius = 0.80,
  gc_content_radius = 0.68,
  gc_content_linewidth = 0.35,
  gc_legend_linewidth = 1.6,
  ruler_linewidth = 0.30,
  label_text_size = 3.4,
  label_text_colour = "category",
  legend_position = "bottom",
  legend_columns = 3,
  gc_legend_columns = 1
)
```

`gene_highlight` changes only the categories you name, and the remaining
categories keep the selected palette. You can also pass a data frame with
`category` and `color` columns. Circular `*_radius` settings are distances from
the center of the plot to that ring.

You can also add highlights after creating the plot:

```r
p <- ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  palette = "npg"
)

p + gene_highlight(
  "Antimicrobial resistance" = "#B2182B",
  "Replication" = "#2166AC"
)
```

For a right-side legend, keep it vertical or use two compact columns:

```r
ggplasmid(
  annotation = gbk_data,
  fasta = fasta_data,
  name = "pSGNDM-5",
  legend_position = "right",
  legend_columns = 2,
  gc_legend_columns = 1
)
```

`read_plasmid_fasta()` also includes a FASTA-formatted text column:

```r
fasta_data$fasta[[1]]
format_plasmid_fasta(fasta_data, width = 70)
```

For table input, use `read_annotation_table()` when you want the table reader
directly:

```r
annotation_data <- read_annotation_table("features.tsv")
```

## Phage Wrapper

For a phage map, use the phage color/classification scheme:

```r
plot_phage_map(
  gbk = "phage.gbk",
  layout = "circular",
  phage_topology = "linear"
)
```

`phage_topology = "linear"` only affects circular layout. It adds a terminal
boundary line at the genome break point. It does not change the map into a
linear plot.
