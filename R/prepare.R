# Prefer compact feature labels on maps, especially for GenBank records where
# product names are often too long for a circular layout.
choose_plot_label <- function(label, gene_id, type, label_mode) {
  label <- clean_text(label)
  gene_id <- clean_text(gene_id)
  type <- clean_text(type)

  short_gene <- nzchar(gene_id) &
    !grepl("^feature_[0-9]+$", gene_id) &
    !grepl("^[A-Z]{2,}[0-9]+[.][0-9]+$", gene_id) &
    nchar(gene_id, type = "chars") <= 16

  if (label_mode == "gene") {
    return(ifelse(short_gene, gene_id, label))
  }
  if (label_mode == "product") {
    return(label)
  }

  out <- label
  out[type == "rep_origin"] <- "ori"
  mobile <- type == "mobile_element"
  out[mobile] <- sub("^[^:]+:", "", out[mobile])
  out[mobile & !nzchar(out)] <- "mobile element"
  out[short_gene] <- gene_id[short_gene]
  out[!nzchar(out)] <- gene_id[!nzchar(out)]
  out
}

#' Prepare annotation features for ggplasmid plotting
#'
#' @param annotation A data frame or path to a CSV/TSV annotation table.
#' @param gbk Optional GenBank file.
#' @param genome_length Optional genome/plasmid length.
#' @param fasta Optional FASTA file used only to infer length when needed.
#' @param feature_types GenBank feature types to keep.
#' @param label_mode `"auto"` prefers short gene symbols when available,
#'   `"product"` uses product/annotation text, and `"gene"` prefers gene IDs.
#' @param category_scheme Either `"plasmid"` or `"phage"`.
#' @return A normalized feature data frame.
#' @export
prepare_plasmid_features <- function(annotation = NULL, gbk = NULL,
                                     genome_length = NULL, fasta = NULL,
                                     feature_types = c(
                                       "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                                       "rep_origin", "mobile_element"
                                     ),
                                     label_mode = c("auto", "product", "gene"),
                                     category_scheme = c("plasmid", "phage")) {
  category_scheme <- match.arg(category_scheme)
  label_mode <- match.arg(label_mode)
  raw <- read_plasmid_annotation(
    annotation = annotation,
    gbk = gbk,
    feature_types = feature_types
  )
  gbk_length <- attr(raw, "genome_length", exact = TRUE)
  gbk_name <- attr(raw, "name", exact = TRUE)
  gbk_sequence <- attr(raw, "sequence", exact = TRUE) %||% ""

  if (is.null(genome_length) || !is.finite(genome_length)) {
    genome_length <- gbk_length
  }
  if ((is.null(genome_length) || !is.finite(genome_length)) && !is.null(fasta)) {
    genome_length <- nchar(read_fasta_sequence(fasta))
  }

  start_col <- find_column(raw, c("start", "begin", "from", "left"), label = "start coordinate")
  end_col <- find_column(raw, c("end", "stop", "to", "right"), label = "end coordinate")
  strand_col <- find_column(raw, c("strand", "frame", "direction"), required = FALSE, label = "strand")
  label_col <- find_column(
    raw,
    c("label", "product", "annot", "annotation", "function", "gene", "name", "note"),
    required = FALSE,
    label = "feature label"
  )
  type_col <- find_column(raw, c("type", "feature", "region"), required = FALSE, label = "feature type")
  id_col <- find_column(
    raw,
    c("gene_id", "gene", "locus_tag", "locustag", "protein_id", "proteinid", "id"),
    required = FALSE,
    label = "feature id"
  )
  category_col <- find_column(
    raw,
    c("category", "group", "pharokka_category", "function_category", "class"),
    required = FALSE,
    label = "category"
  )
  wrap_col <- find_column(
    raw,
    c("wraps_origin", "wrap_origin", "wrap", "circular_wrap", "crosses_origin"),
    required = FALSE,
    label = "origin-wrap flag"
  )

  start <- as_number(raw[[start_col]])
  end <- as_number(raw[[end_col]])
  keep <- is.finite(start) & is.finite(end)
  if (!any(keep)) {
    stop("No valid feature coordinates found.", call. = FALSE)
  }
  raw <- raw[keep, , drop = FALSE]
  start <- start[keep]
  end <- end[keep]

  if (is.na(strand_col)) {
    strand <- rep("+", nrow(raw))
  } else {
    strand <- clean_text(raw[[strand_col]])
  }
  if (is.na(label_col)) {
    label <- rep("", nrow(raw))
  } else {
    label <- clean_text(raw[[label_col]])
  }
  if (is.na(type_col)) {
    type <- rep("CDS", nrow(raw))
  } else {
    type <- clean_text(raw[[type_col]])
    type[!nzchar(type)] <- "CDS"
  }
  if (is.na(id_col)) {
    gene_id <- paste0("feature_", seq_len(nrow(raw)))
  } else {
    gene_id <- clean_text(raw[[id_col]])
    gene_id[!nzchar(gene_id)] <- paste0("feature_", which(!nzchar(gene_id)))
  }

  label[!nzchar(label)] <- gene_id[!nzchar(label)]
  plot_label <- choose_plot_label(label, gene_id, type, label_mode)
  plot_label <- vapply(plot_label, shorten_feature_label, character(1))
  direction <- strand_direction(strand)
  if (is.na(wrap_col)) {
    wraps_origin <- start > end & direction >= 0L
  } else {
    wraps_text <- tolower(clean_text(raw[[wrap_col]]))
    wraps_origin <- wraps_text %in% c("true", "t", "yes", "y", "1", "wrap", "wrapped")
  }
  reverse_table_coords <- start > end & !wraps_origin
  if (any(reverse_table_coords)) {
    start2 <- pmin(start, end)
    end2 <- pmax(start, end)
    start[reverse_table_coords] <- start2[reverse_table_coords]
    end[reverse_table_coords] <- end2[reverse_table_coords]
  }

  if (is.null(genome_length) || !is.finite(genome_length)) {
    genome_length <- max(c(start, end), na.rm = TRUE)
  }

  features <- data.frame(
    feature_id = seq_len(nrow(raw)),
    gene_id = gene_id,
    type = type,
    start = start,
    end = end,
    strand = ifelse(direction < 0L, "-", "+"),
    direction = direction,
    label_raw = label,
    label_display = plot_label,
    label = vapply(plot_label, wrap_label, character(1)),
    stringsAsFactors = FALSE
  )

  features$wraps_origin <- wraps_origin
  features$xmid <- feature_midpoint(features$start, features$end, genome_length)
  features$category <- classify_features(
    data.frame(
      label = features$label_raw,
      type = features$type,
      category = if (!is.na(category_col)) raw[[category_col]] else rep("", nrow(raw)),
      stringsAsFactors = FALSE
    ),
    category_col = "category",
    scheme = category_scheme
  )

  attr(features, "genome_length") <- as.numeric(genome_length)
  attr(features, "name") <- gbk_name %||% NA_character_
  attr(features, "sequence") <- gbk_sequence
  features
}
