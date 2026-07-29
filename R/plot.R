utils::globalVariables(c("colour", "line_colour"))

feature_segments <- function(features, genome_length, circular = TRUE) {
  pieces <- vector("list", nrow(features))
  for (i in seq_len(nrow(features))) {
    row <- features[i, , drop = FALSE]
    wraps_origin <- if ("wraps_origin" %in% names(row)) {
      isTRUE(row$wraps_origin[[1]])
    } else {
      row$start > row$end
    }
    if (circular && wraps_origin) {
      pieces[[i]] <- rbind(
        transform(row, segment_start = row$start, segment_end = genome_length),
        transform(row, segment_start = 1, segment_end = row$end)
      )
    } else {
      pieces[[i]] <- transform(
        row,
        segment_start = pmin(row$start, row$end),
        segment_end = pmax(row$start, row$end)
      )
    }
  }
  out <- do.call(rbind, pieces)
  row.names(out) <- NULL
  out$segment_id <- seq_len(nrow(out))
  out
}

gene_polygon_one <- function(row, y, gene_height, arrow_head_bp,
                             arrow_head_fraction = 0.55) {
  x1 <- as.numeric(row$segment_start)
  x2 <- as.numeric(row$segment_end)
  width <- max(abs(x2 - x1), 1)
  head <- min(width * arrow_head_fraction, arrow_head_bp)
  y_low <- y - gene_height / 2
  y_high <- y + gene_height / 2

  if (row$direction >= 0L) {
    x <- c(x1, x2 - head, x2, x2 - head, x1)
  } else {
    x <- c(x2, x1 + head, x1, x1 + head, x2)
  }
  data.frame(
    feature_id = row$feature_id,
    polygon_id = row$segment_id %||% row$feature_id,
    x = x,
    y = c(y_low, y_low, y, y_high, y_high),
    category = row$category,
    stringsAsFactors = FALSE
  )
}

circular_gene_polygons <- function(features, genome_length,
                                    radius = 1, gene_height = 0.10,
                                    arrow_head_bp = NULL,
                                    arrow_head_fraction = 0.55,
                                    gene_gap_bp = 0) {
  if (is.null(arrow_head_bp)) {
    arrow_head_bp <- max(80, genome_length / 240)
  }
  segments <- feature_segments(features, genome_length, circular = TRUE)
  gene_gap_bp <- as.numeric(gene_gap_bp %||% 0)
  if (is.finite(gene_gap_bp) && gene_gap_bp > 0) {
    width <- pmax(segments$segment_end - segments$segment_start, 0)
    shrink <- pmin(gene_gap_bp / 2, pmax((width - 1) / 2, 0))
    segments$segment_start <- segments$segment_start + shrink
    segments$segment_end <- segments$segment_end - shrink
  }
  polys <- lapply(seq_len(nrow(segments)), function(i) {
    gene_polygon_one(
      segments[i, , drop = FALSE],
      radius,
      gene_height,
      arrow_head_bp,
      arrow_head_fraction = arrow_head_fraction
    )
  })
  do.call(rbind, polys)
}

linear_gene_polygons <- function(features, genome_length, rows = 4,
                                 gene_height = 0.13, arrow_head_bp = NULL,
                                 arrow_head_fraction = 0.55) {
  if (is.null(arrow_head_bp)) {
    arrow_head_bp <- max(80, genome_length / 240)
  }
  segment <- ceiling(genome_length / rows / 1000) * 1000
  if (!is.finite(segment) || segment <= 0) {
    segment <- genome_length
  }
  row_start <- seq(0, max(genome_length - 1, 0), by = segment)
  rows <- length(row_start)
  xmid <- feature_midpoint(features$start, features$end, genome_length)
  row_index <- pmin(floor((xmid - 1) / segment) + 1L, rows)

  segments <- feature_segments(features, genome_length, circular = TRUE)
  segment_mid <- (segments$segment_start + segments$segment_end) / 2
  segments$row_index <- pmin(floor((segment_mid - 1) / segment) + 1L, rows)
  segments$row_y <- rows - segments$row_index + 1
  segments$segment_start <- pmax(segments$segment_start - row_start[segments$row_index], 0)
  segments$segment_end <- pmin(segments$segment_end - row_start[segments$row_index], segment)

  polys <- lapply(seq_len(nrow(segments)), function(i) {
    gene_polygon_one(
      segments[i, , drop = FALSE],
      y = segments$row_y[[i]],
      gene_height = gene_height,
      arrow_head_bp = arrow_head_bp,
      arrow_head_fraction = arrow_head_fraction
    )
  })
  gene_polys <- do.call(rbind, polys)
  gene_polys$row_index <- segments$row_index[match(gene_polys$feature_id, segments$feature_id)]

  list(
    polygons = gene_polys,
    segment = segment,
    rows = rows,
    row_start = row_start,
    row_end = pmin(row_start + segment, genome_length),
    row_index = row_index
  )
}

prepare_gc_skew_for_plot <- function(skew_table = NULL, fasta = NULL, sequence = NULL,
                                     genome_length, circular = TRUE,
                                     window = 100L, step = 1L) {
  if (!is.null(skew_table)) {
    skew <- read_gc_skew_table(skew_table)
  } else if (!is.null(fasta) || (is.character(sequence) && nzchar(sequence))) {
    skew <- compute_gc_skew(
      fasta = fasta,
      sequence = sequence,
      genome_length = genome_length,
      window = window,
      step = step,
      circular = circular
    )
  } else {
    return(NULL)
  }
  skew <- skew[is.finite(skew$position) & is.finite(skew$gc_skew), , drop = FALSE]
  skew <- skew[skew$position >= 1 & skew$position <= genome_length, , drop = FALSE]
  if (!nrow(skew)) {
    return(NULL)
  }
  max_abs <- max(abs(skew$gc_skew), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) {
    max_abs <- 1
  }
  skew$scaled_score <- skew$gc_skew / max_abs
  if ("gc_content" %in% names(skew)) {
    gc_content <- as.numeric(skew$gc_content)
    mean_gc <- mean(gc_content, na.rm = TRUE)
    max_gc_delta <- max(abs(gc_content - mean_gc), na.rm = TRUE)
    if (!is.finite(max_gc_delta) || max_gc_delta <= 0) {
      max_gc_delta <- 1
    }
    skew$scaled_gc_content <- (gc_content - mean_gc) / max_gc_delta
  } else {
    skew$scaled_gc_content <- NA_real_
  }
  skew
}

nice_tick_step <- function(genome_length, target_ticks = 8) {
  raw <- genome_length / target_ticks
  base <- 10 ^ floor(log10(raw))
  candidates <- c(1, 2, 5, 10) * base
  candidates[[which.min(abs(candidates - raw))]]
}

circular_ruler_data <- function(genome_length, major_bp = NULL, minor_bp = NULL,
                                ruler_radius = 0.59,
                                ruler_minor_tick = 0.015,
                                ruler_major_tick = 0.035,
                                ruler_label_radius = 0.47) {
  if (is.null(major_bp) || !is.finite(major_bp)) {
    major_bp <- nice_tick_step(genome_length, target_ticks = 8)
  }
  if (is.null(minor_bp) || !is.finite(minor_bp)) {
    minor_bp <- major_bp / 5
  }

  major <- seq(0, genome_length - 1, by = major_bp)
  minor <- seq(0, genome_length - 1, by = minor_bp)
  minor <- setdiff(round(minor), round(major))

  major_labels <- if (major_bp >= 1000) {
    paste0(format(round(major / 1000, 1), trim = TRUE), " kb")
  } else {
    paste0(format(round(major), trim = TRUE), " bp")
  }
  major_labels[[1]] <- "0"

  list(
    ring = data.frame(
      x = seq(0, genome_length, length.out = 721),
      y = ruler_radius
    ),
    ticks = rbind(
      data.frame(
        x = minor,
        y = ruler_radius - ruler_minor_tick,
        yend = ruler_radius + ruler_minor_tick,
        type = "minor"
      ),
      data.frame(
        x = major,
        y = ruler_radius - ruler_major_tick,
        yend = ruler_radius + ruler_major_tick,
        type = "major"
      )
    ),
    labels = data.frame(
      x = major,
      y = ruler_label_radius,
      label = major_labels,
      theta = (major / genome_length) * 360
    )
  )
}

default_label_pattern <- function(category_scheme) {
  terms <- if (category_scheme == "phage") {
    c(
      "terminase", "portal", "capsid", "head", "tail", "fiber", "fibre",
      "baseplate", "holin", "lysin", "spanin", "integrase", "recombinase",
      "polymerase", "helicase", "nuclease", "tRNA"
    )
  } else {
    c(
      "NDM", "beta-lactamase", "CTX", "OXA", "TEM", "Tet\\(", "Qac",
      "AAC", "ANT", "AadA", "sulfonamide", "dihydrofolate", "Rep",
      "replication", "Tra", "conjugative transfer", "relaxase", "FinO",
      "ParA", "ParB", "partition", "integrase", "transposase",
      "\\bIS[0-9A-Za-z-]*\\b", "toxin", "antitoxin"
    )
  }
  paste(terms, collapse = "|")
}

select_label_data <- function(features, label_unknown = FALSE,
                              label_exclude_categories = NULL,
                              label_pattern = NULL, max_labels = Inf,
                              min_feature_bp = 1,
                              category_scheme = "plasmid",
                              genome_length = NULL,
                              label_min_gap_deg = 0) {
  label_data <- features[nzchar(features$label), , drop = FALSE]
  if (!label_unknown) {
    unknown_label <- label_data$category == "Unknown function, hypothetical protein" |
      is_unknown_feature_label(
        label_data$label_raw,
        label_data$label_display,
        label_data$gene_id
      )
    label_data <- label_data[
      !unknown_label,
      ,
      drop = FALSE
    ]
  }
  if (!is.null(label_exclude_categories) && length(label_exclude_categories)) {
    label_exclude_categories <- as.character(label_exclude_categories)
    label_exclude_categories <- label_exclude_categories[nzchar(label_exclude_categories)]
    if (length(label_exclude_categories)) {
      label_data <- label_data[
        !(label_data$category %in% label_exclude_categories),
        ,
        drop = FALSE
      ]
    }
  }
  if (!nrow(label_data)) {
    return(label_data)
  }

  label_data$bp <- abs(label_data$end - label_data$start) + 1
  label_data <- label_data[label_data$bp >= min_feature_bp, , drop = FALSE]
  if (!nrow(label_data)) {
    return(label_data)
  }

  pattern <- label_pattern %||% default_label_pattern(category_scheme)
  label_display <- label_data$label_display %||% label_data$label_raw
  label_data$pattern_hit <- if (is.null(pattern) || !nzchar(pattern)) {
    FALSE
  } else {
    grepl(pattern, label_data$label_raw, ignore.case = TRUE, perl = TRUE) |
      grepl(pattern, label_display, ignore.case = TRUE, perl = TRUE) |
      grepl(pattern, label_data$gene_id, ignore.case = TRUE, perl = TRUE)
  }

  priority <- c(
    "Antimicrobial resistance" = 1,
    "Replication" = 2,
    "Conjugation and transfer" = 3,
    "Partition and maintenance" = 4,
    "Mobile element" = 5,
    "Virulence" = 6,
    "Regulation" = 7,
    "Metabolism" = 8,
    "Other functions" = 9,
    "Unknown function, hypothetical protein" = 10,
    "DNA, RNA and nucleotide metabolism" = 3,
    "Head and packaging" = 3,
    "Tail" = 3,
    "Lysis" = 3,
    "Connectors" = 4,
    "Transcription regulation" = 7,
    "Moron, auxiliary metabolic gene and host takeover" = 8,
    "tRNA" = 8
  )
  label_data$priority <- unname(priority[label_data$category])
  label_data$priority[is.na(label_data$priority)] <- 9
  label_data$priority <- label_data$priority - ifelse(label_data$pattern_hit, 4, 0)
  label_data <- label_data[order(label_data$priority, -label_data$bp, label_data$xmid), , drop = FALSE]

  if (is.finite(label_min_gap_deg) && label_min_gap_deg > 0 &&
      is.finite(genome_length) && nrow(label_data) > 1L) {
    label_data$theta <- (label_data$xmid %% genome_length) / genome_length * 360
    keep <- rep(FALSE, nrow(label_data))
    kept_theta <- numeric()
    for (i in seq_len(nrow(label_data))) {
      distance <- if (!length(kept_theta)) {
        Inf
      } else {
        abs(((label_data$theta[[i]] - kept_theta + 180) %% 360) - 180)
      }
      if (all(distance >= label_min_gap_deg)) {
        keep[[i]] <- TRUE
        kept_theta <- c(kept_theta, label_data$theta[[i]])
      }
    }
    label_data <- label_data[keep, , drop = FALSE]
    label_data$theta <- NULL
  }

  if (is.finite(max_labels) && nrow(label_data) > max_labels) {
    label_data <- label_data[seq_len(max_labels), , drop = FALSE]
  }

  label_data[order(label_data$xmid), , drop = FALSE]
}

radial_npc <- function(x, y, genome_length, y_max = 1.95,
                       inner_radius = 0.38) {
  theta <- 2 * pi * as.numeric(x) / genome_length
  radius <- 0.4 * (
    inner_radius + (as.numeric(y) / y_max) * (1 - inner_radius)
  )
  data.frame(
    npcx = 0.5 + radius * sin(theta),
    npcy = 0.5 + radius * cos(theta)
  )
}

spread_label_positions <- function(y, lower = 0.08, upper = 0.92,
                                   min_gap = 0.055) {
  if (!length(y)) {
    return(numeric())
  }
  if (length(y) == 1L) {
    return(pmin(pmax(y, lower), upper))
  }

  ordering <- order(y)
  sorted <- pmin(pmax(y[ordering], lower), upper)
  gap <- min(min_gap, (upper - lower) / (length(sorted) - 1L))
  for (i in 2:length(sorted)) {
    sorted[[i]] <- max(sorted[[i]], sorted[[i - 1L]] + gap)
  }
  if (sorted[[length(sorted)]] > upper) {
    sorted <- sorted - (sorted[[length(sorted)]] - upper)
  }
  for (i in (length(sorted) - 1L):1L) {
    sorted[[i]] <- min(sorted[[i]], sorted[[i + 1L]] - gap)
  }
  if (sorted[[1]] < lower) {
    sorted <- sorted + (lower - sorted[[1]])
  }

  out <- numeric(length(y))
  out[ordering] <- sorted
  out
}

spread_label_boxes <- function(y, height, lower = 0.04, upper = 0.96,
                               padding = 0.008) {
  if (!length(y)) {
    return(numeric())
  }
  height <- pmax(as.numeric(height), 0.01)
  if (length(y) == 1L) {
    return(pmin(pmax(y, lower + height / 2), upper - height / 2))
  }

  ordering <- order(y)
  sorted <- y[ordering]
  sorted_height <- height[ordering]
  total_height <- sum(sorted_height) + padding * (length(sorted) - 1L)
  if (total_height >= upper - lower) {
    centers <- seq(lower, upper, length.out = length(sorted))
    out <- numeric(length(y))
    out[ordering] <- centers
    return(out)
  }

  sorted <- pmin(pmax(sorted, lower + sorted_height / 2), upper - sorted_height / 2)
  for (i in 2:length(sorted)) {
    sorted[[i]] <- max(
      sorted[[i]],
      sorted[[i - 1L]] + sorted_height[[i - 1L]] / 2 + sorted_height[[i]] / 2 + padding
    )
  }
  overflow <- sorted[[length(sorted)]] + sorted_height[[length(sorted)]] / 2 - upper
  if (overflow > 0) {
    sorted <- sorted - overflow
  }
  for (i in (length(sorted) - 1L):1L) {
    sorted[[i]] <- min(
      sorted[[i]],
      sorted[[i + 1L]] - sorted_height[[i + 1L]] / 2 - sorted_height[[i]] / 2 - padding
    )
  }
  underflow <- lower - (sorted[[1L]] - sorted_height[[1L]] / 2)
  if (underflow > 0) {
    sorted <- sorted + underflow
  }

  out <- numeric(length(y))
  out[ordering] <- sorted
  out
}

wrap_label_text <- function(text, width = 22L) {
  text <- clean_text(text)
  text <- gsub("[\r\n]+", " ", text)
  if (is.null(width) || !is.finite(width) || width <= 0) {
    return(text)
  }
  vapply(
    text,
    function(one) {
      paste(strwrap(one, width = width), collapse = "\n")
    },
    character(1)
  )
}

max_label_line_chars <- function(text) {
  vapply(
    strsplit(as.character(text), "\n", fixed = TRUE),
    function(lines) max(nchar(lines, type = "chars"), 1),
    numeric(1)
  )
}

label_box_dimensions <- function(text, text_size, lineheight = 0.9) {
  line_count <- lengths(strsplit(as.character(text), "\n", fixed = TRUE))
  data.frame(
    width = 0.012 + max_label_line_chars(text) * text_size * 0.0024,
    height = 0.006 + pmax(line_count, 1L) * text_size * 0.0065 * lineheight
  )
}

label_box_overlaps <- function(candidate, placed, padding = 0.006) {
  if (!nrow(placed)) {
    return(FALSE)
  }
  any(
    candidate$xmin < placed$xmax + padding &
      candidate$xmax > placed$xmin - padding &
      candidate$ymin < placed$ymax + padding &
      candidate$ymax > placed$ymin - padding
  )
}

label_box_frame <- function(x, y, width, height, right = TRUE, hjust = NULL) {
  if (is.null(hjust)) {
    hjust <- if (isTRUE(right)) 0 else 1
  }
  data.frame(
    xmin = x - width * hjust,
    xmax = x + width * (1 - hjust),
    ymin = y - height / 2,
    ymax = y + height / 2
  )
}

label_overlay_npc_bounds <- function(label_overlay, padding = 0.006) {
  if (is.null(label_overlay) || !nrow(label_overlay)) {
    return(NULL)
  }
  keep <- is.finite(label_overlay$npcx) & is.finite(label_overlay$npcy)
  if (!any(keep)) {
    return(NULL)
  }
  label_overlay <- label_overlay[keep, , drop = FALSE]
  dims <- label_box_dimensions(label_overlay$label, label_overlay$size)
  boxes <- label_box_frame(
    label_overlay$npcx,
    label_overlay$npcy,
    dims$width,
    dims$height,
    hjust = label_overlay$hjust
  )
  c(
    xmin = min(boxes$xmin, na.rm = TRUE) - padding,
    xmax = max(boxes$xmax, na.rm = TRUE) + padding,
    ymin = min(boxes$ymin, na.rm = TRUE) - padding,
    ymax = max(boxes$ymax, na.rm = TRUE) + padding
  )
}

auto_legend_plot_spacing <- function(spacing, position, label_bounds,
                                     cm_per_npc = 18, padding_cm = 0.8) {
  if (is.null(label_bounds) || !length(label_bounds)) {
    return(spacing)
  }
  overflow <- switch(
    position,
    right = max(0, unname(label_bounds[["xmax"]]) - 1),
    left = max(0, -unname(label_bounds[["xmin"]])),
    top = max(0, unname(label_bounds[["ymax"]]) - 1),
    bottom = max(0, -unname(label_bounds[["ymin"]])),
    0
  )
  if (!is.finite(overflow) || overflow <= 0) {
    return(spacing)
  }
  max(spacing, overflow * cm_per_npc + padding_cm)
}

circular_plot_margin <- function(label_bounds = NULL,
                                 base = c(top = 16, right = 28, bottom = 16, left = 28),
                                 overflow_scale = 420) {
  margin_values <- base
  if (!is.null(label_bounds) && length(label_bounds)) {
    overflow <- c(
      top = max(0, unname(label_bounds[["ymax"]]) - 1),
      right = max(0, unname(label_bounds[["xmax"]]) - 1),
      bottom = max(0, -unname(label_bounds[["ymin"]])),
      left = max(0, -unname(label_bounds[["xmin"]]))
    )
    margin_values <- margin_values + overflow * overflow_scale
  }

  ggplot2::margin(
    t = unname(margin_values[["top"]]),
    r = unname(margin_values[["right"]]),
    b = unname(margin_values[["bottom"]]),
    l = unname(margin_values[["left"]])
  )
}

pack_side_label_boxes <- function(y, height, lower = 0.015, upper = 0.985,
                                  padding = 0.006) {
  n <- length(y)
  if (!n) {
    return(list(y = numeric(), lane = integer()))
  }

  y <- as.numeric(y)
  y[!is.finite(y)] <- 0.5
  height <- pmax(as.numeric(height), 0.01)
  available <- upper - lower
  if (!is.finite(available) || available <= 0) {
    return(list(y = y, lane = integer(n)))
  }

  total_height <- function(idx) {
    if (!length(idx)) {
      return(0)
    }
    sum(height[idx]) + padding * max(length(idx) - 1L, 0L)
  }

  ordering <- order(y)
  column_count <- max(
    1L,
    ceiling(total_height(seq_len(n)) / available)
  )
  column_count <- min(column_count, n)

  repeat {
    lane_sorted <- (seq_along(ordering) - 1L) %% column_count
    overfull <- FALSE
    for (lane in 0:(column_count - 1L)) {
      idx <- ordering[lane_sorted == lane]
      if (total_height(idx) > available && column_count < n) {
        overfull <- TRUE
        break
      }
    }
    if (!overfull || column_count >= n) {
      break
    }
    column_count <- column_count + 1L
  }

  y_out <- numeric(n)
  lanes <- integer(n)
  lane_sorted <- (seq_along(ordering) - 1L) %% column_count
  for (lane in 0:(column_count - 1L)) {
    idx <- ordering[lane_sorted == lane]
    if (!length(idx)) {
      next
    }
    y_out[idx] <- spread_label_boxes(
      y[idx],
      height[idx],
      lower = lower,
      upper = upper,
      padding = padding
    )
    lanes[idx] <- lane
  }

  list(y = y_out, lane = lanes)
}

place_radial_label_boxes <- function(x, genome_length, text, right,
                                     inner_radius = 0.38,
                                     label_anchor_radius = 1.07,
                                     label_line_length = 0.02,
                                     label_lane_spacing = 0.24,
                                     label_lane_step = NULL,
                                     label_text_size = 3.4,
                                     priority = NULL,
                                     lower = 0.04, upper = 0.96,
                                     y_max = 1.95) {
  n <- length(x)
  if (!n) {
    return(data.frame(npcx = numeric(), npcy = numeric(), lane = integer()))
  }

  dims <- label_box_dimensions(text, label_text_size)
  data_per_npc_radius <- y_max / (0.4 * (1 - inner_radius))
  lane_height <- stats::quantile(
    dims$height[is.finite(dims$height)],
    0.75,
    names = FALSE,
    na.rm = TRUE
  )
  if (!is.finite(lane_height)) {
    lane_height <- 0.02
  }
  if (is.null(label_lane_step)) {
    lane_step <- max(
      0.026,
      lane_height * data_per_npc_radius * 0.55,
      label_lane_spacing * 0.12
    )
  } else {
    lane_step <- as.numeric(label_lane_step)
    if (length(lane_step) != 1L || !is.finite(lane_step) || lane_step <= 0) {
      stop("`label_lane_step` must be a positive number.", call. = FALSE)
    }
  }
  padding <- 0.006

  right <- as.logical(right)
  right[is.na(right)] <- TRUE
  theta <- 2 * pi * as.numeric(x) / genome_length
  label_hjust <- ifelse(abs(sin(theta)) < 0.20, 0.5, ifelse(right, 0, 1))
  vertical_label <- abs(sin(theta)) < 0.20
  vertical_direction <- sign(cos(theta))
  vertical_direction[!vertical_label] <- 0
  vertical_edge_gap <- 0.002
  npcx <- rep(NA_real_, n)
  npcy <- rep(NA_real_, n)
  lanes <- rep(NA_integer_, n)
  placed <- data.frame(
    xmin = numeric(), xmax = numeric(),
    ymin = numeric(), ymax = numeric()
  )

  label_order <- order(
    dims$width,
    dims$height,
    max_label_line_chars(text),
    x
  )
  remaining <- label_order
  lane <- 0L
  max_lanes <- max(32L, n * 4L)

  candidate_box <- function(i, lane) {
    position <- radial_npc(
      x[[i]],
      y = label_anchor_radius + label_line_length + lane * lane_step,
      genome_length = genome_length,
      y_max = y_max,
      inner_radius = inner_radius
    )
    text_x <- position$npcx[[1L]]
    text_y <- position$npcy[[1L]] +
      vertical_direction[[i]] * (dims$height[[i]] / 2 + vertical_edge_gap)
    list(
      x = text_x,
      y = text_y,
      box = label_box_frame(
        text_x,
        text_y,
        dims$width[[i]],
        dims$height[[i]],
        hjust = label_hjust[[i]]
      )
    )
  }

  while (length(remaining) && lane <= max_lanes) {
    next_remaining <- integer()
    for (i in remaining[order(dims$width[remaining], dims$height[remaining], x[remaining])]) {
      candidate <- candidate_box(i, lane)
      if (!label_box_overlaps(candidate$box, placed, padding = padding)) {
        npcx[[i]] <- candidate$x
        npcy[[i]] <- candidate$y
        lanes[[i]] <- lane
        placed <- rbind(placed, candidate$box)
      } else {
        next_remaining <- c(next_remaining, i)
      }
    }
    remaining <- next_remaining
    lane <- lane + 1L
  }

  for (i in remaining) {
    repeat {
      candidate <- candidate_box(i, lane)
      if (!label_box_overlaps(candidate$box, placed, padding = padding) ||
          lane > max_lanes + n * 8L) {
        npcx[[i]] <- candidate$x
        npcy[[i]] <- candidate$y
        lanes[[i]] <- lane
        placed <- rbind(placed, candidate$box)
        lane <- lane + 1L
        break
      }
      lane <- lane + 1L
    }
  }

  data.frame(npcx = npcx, npcy = npcy, lane = lanes, hjust = label_hjust)
}

radial_label_lanes <- function(theta, text, min_gap_deg = 5,
                               char_width_deg = 0.35, max_lanes = Inf) {
  if (!length(theta)) {
    return(integer())
  }
  if (length(theta) == 1L) {
    return(0L)
  }

  theta <- theta %% 360
  sorted <- sort(theta)
  gaps <- diff(c(sorted, sorted[[1L]] + 360))
  start_angle <- sorted[(which.max(gaps) %% length(sorted)) + 1L]
  unwrapped <- ifelse(theta < start_angle, theta + 360, theta)

  line_count <- lengths(strsplit(as.character(text), "\n", fixed = TRUE))
  span <- pmax(
    5,
    pmin(70, (max_label_line_chars(text) + pmax(line_count - 1L, 0L) * 6) * char_width_deg)
  )
  label_start <- unwrapped - span / 2
  label_end <- unwrapped + span / 2

  if (!is.finite(max_lanes)) {
    max_lanes <- length(theta)
  }
  max_lanes <- max(1L, as.integer(max_lanes))
  lanes <- integer(length(theta))
  lane_end <- rep(-Inf, max_lanes)
  for (i in order(unwrapped)) {
    available <- which(label_start[[i]] > lane_end + min_gap_deg)
    lane <- if (length(available)) {
      available[[1L]]
    } else {
      which.min(lane_end)
    }
    lanes[[i]] <- lane - 1L
    lane_end[[lane]] <- label_end[[i]]
  }
  lanes
}

screen_label_lanes <- function(y, height, right, padding = 0.006,
                               max_lanes = 3L) {
  if (!length(y)) {
    return(integer())
  }
  lanes <- integer(length(y))
  max_lanes <- max(1L, as.integer(max_lanes))

  for (side in c(FALSE, TRUE)) {
    index <- which(right == side)
    if (!length(index)) {
      next
    }
    ordering <- index[order(y[index])]
    lane_end <- rep(-Inf, max_lanes)
    for (i in ordering) {
      label_start <- y[[i]] - height[[i]] / 2
      label_end <- y[[i]] + height[[i]] / 2
      available <- which(label_start > lane_end + padding)
      lane <- if (length(available)) {
        available[[1L]]
      } else {
        which.min(lane_end)
      }
      lanes[[i]] <- lane - 1L
      lane_end[[lane]] <- max(lane_end[[lane]], label_end)
    }
  }

  lanes
}

GeomPlasmidTextOverlay <- ggplot2::ggproto(
  "GeomPlasmidTextOverlay",
  ggplot2::Geom,
  required_aes = c("label"),
  default_aes = ggplot2::aes(
    npcx = NA_real_,
    npcy = NA_real_,
    panel_x = NA_real_,
    panel_y = NA_real_,
    angle = 0,
    hjust = 0.5,
    vjust = 0.5,
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = NA_real_,
    y1 = NA_real_,
    x2 = NA_real_,
    y2 = NA_real_,
    line_colour = NA_character_,
    label_colour = NA_character_
  ),
  draw_key = function(...) grid::nullGrob(),
  draw_panel = function(data, panel_params, coord, family = "sans",
                        text_colour = "black", text_size = 3,
                        fontface = "plain", lineheight = 0.9,
                        segment_colour = "grey45",
                        segment_linewidth = 0.25,
                        segment_linetype = "solid", na.rm = FALSE) {
    grobs <- list()

    text_x <- data$npcx
    text_y <- data$npcy
    has_panel_text <- is.finite(data$panel_x) & is.finite(data$panel_y)
    if (any(has_panel_text)) {
      panel_text <- coord$transform(
        data.frame(
          x = data$panel_x[has_panel_text],
          y = data$panel_y[has_panel_text]
        ),
        panel_params
      )
      text_x[has_panel_text] <- panel_text$x
      text_y[has_panel_text] <- panel_text$y
    }

    line_x0 <- data$x0
    line_y0 <- data$y0
    has_anchor <- is.finite(data$anchor_x) & is.finite(data$anchor_y)
    if (any(has_anchor)) {
      anchor <- coord$transform(
        data.frame(
          x = data$anchor_x[has_anchor],
          y = data$anchor_y[has_anchor]
        ),
        panel_params
      )
      line_x0[has_anchor] <- anchor$x
      line_y0[has_anchor] <- anchor$y
    }

    line_x1 <- data$x1
    line_y1 <- data$y1
    line_x2 <- data$x2
    line_y2 <- data$y2
    has_leader <- is.finite(data$leader_x) & is.finite(data$leader_y)
    if (any(has_leader)) {
      leader <- coord$transform(
        data.frame(
          x = data$leader_x[has_leader],
          y = data$leader_y[has_leader]
        ),
        panel_params
      )
      line_x1[has_leader] <- leader$x
      line_y1[has_leader] <- leader$y
      line_x2[has_leader] <- leader$x
      line_y2[has_leader] <- leader$y
    }
    use_text_endpoint <- has_anchor & is.finite(text_x) & is.finite(text_y) &
      !has_leader &
      (!is.finite(line_x1) | !is.finite(line_y1) |
        !is.finite(line_x2) | !is.finite(line_y2))
    if (any(use_text_endpoint)) {
      gap <- 0.008
      side_gap <- ifelse(
        data$hjust <= 0.25,
        gap,
        ifelse(data$hjust >= 0.75, -gap, 0)
      )
      line_x1[use_text_endpoint] <- text_x[use_text_endpoint] - side_gap[use_text_endpoint]
      line_y1[use_text_endpoint] <- text_y[use_text_endpoint]
      line_x2[use_text_endpoint] <- line_x1[use_text_endpoint]
      line_y2[use_text_endpoint] <- text_y[use_text_endpoint]
    }

    has_line <- is.finite(line_x0) & is.finite(line_y0) &
      is.finite(line_x1) & is.finite(line_y1) &
      is.finite(line_x2) & is.finite(line_y2)
    if (any(has_line)) {
      line_data <- data[has_line, , drop = FALSE]
      line_data$x0 <- line_x0[has_line]
      line_data$y0 <- line_y0[has_line]
      line_data$x1 <- line_x1[has_line]
      line_data$y1 <- line_y1[has_line]
      line_data$x2 <- line_x2[has_line]
      line_data$y2 <- line_y2[has_line]
      for (j in seq_len(nrow(line_data))) {
        line_colour <- if ("line_colour" %in% names(line_data) &&
          is.character(line_data$line_colour) &&
          !is.na(line_data$line_colour[[j]]) &&
          nzchar(line_data$line_colour[[j]])) {
          line_data$line_colour[[j]]
        } else {
          segment_colour
        }
        grobs[[length(grobs) + 1L]] <- grid::polylineGrob(
          x = grid::unit(
            c(line_data$x0[[j]], line_data$x1[[j]], line_data$x2[[j]]),
            "npc"
          ),
          y = grid::unit(
            c(line_data$y0[[j]], line_data$y1[[j]], line_data$y2[[j]]),
            "npc"
          ),
          gp = grid::gpar(
            col = line_colour,
            lwd = segment_linewidth * 2.845276,
            lty = segment_linetype,
            lineend = "round",
            linejoin = "round"
          )
        )
      }
    }

    for (i in seq_len(nrow(data))) {
      label_colour <- if ("label_colour" %in% names(data) &&
        is.character(data$label_colour) &&
        !is.na(data$label_colour[[i]]) &&
        nzchar(data$label_colour[[i]])) {
        data$label_colour[[i]]
      } else {
        text_colour
      }
      grobs[[length(grobs) + 1L]] <- grid::textGrob(
        label = data$label[[i]],
        x = grid::unit(text_x[[i]], "npc"),
        y = grid::unit(text_y[[i]], "npc"),
        hjust = data$hjust[[i]],
        vjust = data$vjust[[i]],
        rot = data$angle[[i]],
        gp = grid::gpar(
          col = label_colour,
          fontsize = text_size * 2.845276,
          fontfamily = family,
          fontface = fontface,
          lineheight = lineheight
        )
      )
    }
    if (!length(grobs)) {
      return(grid::nullGrob())
    }
    do.call(grid::grobTree, grobs)
  }
)

geom_plasmid_text_overlay <- function(data, family = "sans",
                                      text_colour = "black",
                                      text_size = 3,
                                      fontface = "plain",
                                      lineheight = 0.9,
                                      segment_colour = "grey45",
                                      segment_linewidth = 0.25,
                                      segment_linetype = "solid") {
  if (!"line_colour" %in% names(data)) {
    data$line_colour <- NA_character_
  }
  ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(
      npcx = npcx,
      npcy = npcy,
      panel_x = panel_x,
      panel_y = panel_y,
      label = label,
      angle = angle,
      hjust = hjust,
      vjust = vjust,
      anchor_x = anchor_x,
      anchor_y = anchor_y,
      leader_x = leader_x,
      leader_y = leader_y,
      x0 = x0,
      y0 = y0,
      x1 = x1,
      y1 = y1,
      x2 = x2,
      y2 = y2,
      line_colour = line_colour,
      label_colour = colour
    ),
    stat = "identity",
    geom = GeomPlasmidTextOverlay,
    position = "identity",
    inherit.aes = FALSE,
    show.legend = FALSE,
    params = list(
      family = family,
      text_colour = text_colour,
      text_size = text_size,
      fontface = fontface,
      lineheight = lineheight,
      segment_colour = segment_colour,
      segment_linewidth = segment_linewidth,
      segment_linetype = segment_linetype,
      na.rm = TRUE
    )
  )
}

resolve_label_text_colours <- function(label_data, feature_colors,
                                       label_text_colour = "black") {
  n <- nrow(label_data)
  if (!n) {
    return(character())
  }
  if (is.null(label_text_colour)) {
    return(rep("black", n))
  }

  if (length(label_text_colour) == 1L &&
      identical(unname(label_text_colour), "category")) {
    colours <- unname(feature_colors[as.character(label_data$category)])
    colours[is.na(colours)] <- "black"
    return(colours)
  }

  if (!is.null(names(label_text_colour)) && any(nzchar(names(label_text_colour)))) {
    colours <- unname(label_text_colour[as.character(label_data$category)])
    colours[is.na(colours)] <- "black"
    return(colours)
  }

  if (length(label_text_colour) == n) {
    return(as.character(label_text_colour))
  }

  rep(as.character(label_text_colour[[1L]]), n)
}

circular_feature_label_overlay <- function(label_data, genome_length,
                                           inner_radius = 0.38,
                                           label_line_length = 0.02,
                                           label_wrap_width = 28L,
                                           label_anchor_radius = 1.07,
                                           label_line_anchor_radius = label_anchor_radius,
                                           label_lane_spacing = 0.24,
                                           label_lane_step = NULL,
                                           label_text_size = 3.4,
                                           label_text_colour = "black",
                                           label_line_colour = "grey70",
                                           label_y_max = 1.95,
                                           feature_colors = NULL) {
  theta <- (label_data$xmid %% genome_length) / genome_length * 360
  right <- theta < 180

  text <- wrap_label_text(
    label_data$label_display %||% label_data$label,
    width = label_wrap_width
  )
  label_colours <- resolve_label_text_colours(
    label_data,
    feature_colors %||% stats::setNames(character(), character()),
    label_text_colour = label_text_colour
  )
  line_colours <- resolve_label_text_colours(
    label_data,
    feature_colors %||% stats::setNames(character(), character()),
    label_text_colour = label_line_colour
  )

  label_position <- place_radial_label_boxes(
    x = label_data$xmid,
    genome_length = genome_length,
    text = text,
    right = right,
    inner_radius = inner_radius,
    label_anchor_radius = label_anchor_radius,
    label_line_length = label_line_length,
    label_lane_spacing = label_lane_spacing,
    label_lane_step = label_lane_step,
    label_text_size = label_text_size,
    priority = label_data$priority,
    y_max = label_y_max
  )
  keep <- is.finite(label_position$npcx) & is.finite(label_position$npcy)
  label_position <- label_position[keep, , drop = FALSE]
  label_data <- label_data[keep, , drop = FALSE]
  text <- text[keep]
  label_colours <- label_colours[keep]
  line_colours <- line_colours[keep]
  theta <- theta[keep]
  if (!nrow(label_position)) {
    return(data.frame(
      npcx = numeric(),
      npcy = numeric(),
      panel_x = numeric(),
      panel_y = numeric(),
      label = character(),
      angle = numeric(),
      hjust = numeric(),
      vjust = numeric(),
      colour = character(),
      size = numeric(),
      fontface = character(),
      lineheight = numeric(),
      anchor_x = numeric(),
      anchor_y = numeric(),
      leader_x = numeric(),
      leader_y = numeric(),
      x0 = numeric(),
      y0 = numeric(),
      x1 = numeric(),
      y1 = numeric(),
      x2 = numeric(),
      y2 = numeric(),
      line_colour = character(),
      stringsAsFactors = FALSE
    ))
  }

  label_dims <- label_box_dimensions(text, label_text_size)
  line_x <- label_position$npcx
  line_y <- label_position$npcy
  top_bottom <- abs(sin(theta * pi / 180)) < 0.20
  top_label <- top_bottom & cos(theta * pi / 180) > 0
  bottom_label <- top_bottom & cos(theta * pi / 180) < 0
  edge_gap <- 0.002
  line_y[top_label] <- label_position$npcy[top_label] -
    label_dims$height[top_label] / 2 - edge_gap
  line_y[bottom_label] <- label_position$npcy[bottom_label] +
    label_dims$height[bottom_label] / 2 + edge_gap

  data.frame(
    npcx = label_position$npcx,
    npcy = label_position$npcy,
    panel_x = NA_real_,
    panel_y = NA_real_,
    label = text,
    angle = 0,
    hjust = label_position$hjust,
    vjust = 0.5,
    colour = label_colours,
    size = label_text_size,
    fontface = "bold",
    lineheight = 0.9,
    anchor_x = label_data$xmid,
    anchor_y = label_line_anchor_radius,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = line_x,
    y1 = line_y,
    x2 = line_x,
    y2 = line_y,
    line_colour = line_colours,
    stringsAsFactors = FALSE
  )
}

circular_ruler_label_overlay <- function(ruler, genome_length,
                                         inner_radius = 0.38) {
  theta <- ruler$labels$theta
  angle <- -theta
  flip <- theta > 90 & theta < 270
  angle[flip] <- angle[flip] + 180
  data.frame(
    npcx = NA_real_,
    npcy = NA_real_,
    panel_x = ruler$labels$x,
    panel_y = ruler$labels$y,
    label = ruler$labels$label,
    angle = angle,
    hjust = 0.5,
    vjust = 0.5,
    colour = "grey25",
    size = 2.25,
    fontface = "plain",
    lineheight = 0.9,
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = NA_real_,
    y1 = NA_real_,
    x2 = NA_real_,
    y2 = NA_real_,
    stringsAsFactors = FALSE
  )
}

center_text_overlay <- function(name, genome_length) {
  data.frame(
    npcx = 0.5,
    npcy = 0.5,
    panel_x = NA_real_,
    panel_y = NA_real_,
    label = paste(c(name, paste0(format_bp(genome_length), " bp")), collapse = "\n"),
    angle = 0,
    hjust = 0.5,
    vjust = 0.5,
    colour = "black",
    size = 4.3,
    fontface = "bold",
    lineheight = 0.9,
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = NA_real_,
    y1 = NA_real_,
    x2 = NA_real_,
    y2 = NA_real_,
    stringsAsFactors = FALSE
  )
}

plot_circular_plasmid <- function(features, genome_length, name = NULL,
                                  phage_topology = c("circular", "linear"),
                                  category_scheme = "plasmid",
                                  gc_skew = NULL, show_labels = TRUE,
                                  label_unknown = FALSE,
                                  label_exclude_categories = NULL,
                                  label_pattern = NULL, max_labels = Inf,
                                  label_min_gap_deg = 0,
                                  label_line_length = 0.02,
                                  mini_label_line_length = NULL,
                                  label_wrap_width = 28L,
                                  label_anchor_radius = 1.07,
                                  label_lane_spacing = 0.24,
                                  label_lane_step = NULL,
                                  label_line_colour = "grey70",
                                  label_linewidth = 0.30,
                                  label_line_linetype = "dashed",
                                  min_feature_bp = 1,
                                  show_ruler = TRUE,
                                  ruler_major_bp = NULL,
                                  ruler_minor_bp = NULL,
                                  palette = "default",
                                  gene_highlight = NULL,
                                  gene_radius = 1,
                                  gene_height = 0.10,
                                  gene_linewidth = 0.22,
                                  gene_border_linewidth = NULL,
                                  gene_arrow_head_bp = NULL,
                                  gene_arrow_head_fraction = 0.55,
                                  gene_gap_bp = 0,
                                  gc_skew_radius = 0.78,
                                  gc_skew_height = 0.10,
                                  gc_content_radius = 0.67,
                                  gc_content_height = 0.045,
                                  gc_content_linewidth = 0.28,
                                  gc_legend_linewidth = 1.8,
                                  ruler_radius = 0.59,
                                  ruler_minor_tick = 0.015,
                                  ruler_major_tick = 0.035,
                                  ruler_label_radius = 0.47,
                                  ruler_linewidth = 0.28,
                                  ruler_minor_linewidth = 0.18,
                                  ruler_major_linewidth = 0.38,
                                  label_text_size = 3.4,
                                  label_text_colour = "black",
                                  ruler_text_size = 2.25,
                                  center_text_size = 4.3,
                                  legend_position = "bottom",
                                  legend_text_size = 8.5,
                                  legend_columns = NULL,
                                  gc_legend_columns = NULL,
                                  legend_plot_spacing = NULL,
                                  font_family = "sans",
                                  inner_radius = 0.38) {
  phage_topology <- match.arg(phage_topology)
  if (!is.null(mini_label_line_length)) {
    label_line_length <- mini_label_line_length
  }
  gene_linewidth <- gene_border_linewidth %||% gene_linewidth
  legend_columns <- validate_legend_columns(legend_columns, legend_position)
  gc_legend_columns <- validate_gc_legend_columns(gc_legend_columns, legend_position)
  legend_plot_spacing_auto <- is.null(legend_plot_spacing)
  legend_plot_spacing <- validate_legend_plot_spacing(legend_plot_spacing, legend_position)
  plot_colors <- ggplasmid_resolve_colors(
    category_scheme,
    palette = palette,
    gene_highlight = gene_highlight
  )
  feature_colors <- plot_colors[names(plot_colors) %in% unique(c(features$category, "GC skew+", "GC skew-"))]
  feature_breaks <- gene_legend_breaks(feature_colors)
  gc_colors <- gc_track_colors(plot_colors)
  gene_polys <- circular_gene_polygons(
    features,
    genome_length,
    radius = gene_radius,
    gene_height = gene_height,
    arrow_head_bp = gene_arrow_head_bp,
    arrow_head_fraction = gene_arrow_head_fraction,
    gene_gap_bp = gene_gap_bp
  )
  gene_polys$category <- factor(gene_polys$category, levels = names(plot_colors))
  ruler <- NULL

  p <- ggplot2::ggplot()

  if (!is.null(gc_skew)) {
    skew <- gc_skew
    skew$base_y <- gc_skew_radius
    skew$pos_y <- skew$base_y + pmax(skew$scaled_score, 0) * gc_skew_height
    skew$neg_y <- skew$base_y + pmin(skew$scaled_score, 0) * gc_skew_height
    skew$content_y <- gc_content_radius + skew$scaled_gc_content * gc_content_height
    p <- p +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = position, ymin = base_y, ymax = pos_y, fill = "GC skew+"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = position, ymin = base_y, ymax = neg_y, fill = "GC skew-"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_segment(
        data = gc_skew_legend_data(),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = track),
        inherit.aes = FALSE,
        alpha = 0,
        linewidth = gc_content_linewidth,
        show.legend = TRUE
      )
    if (any(is.finite(skew$scaled_gc_content))) {
      p <- p +
        ggplot2::geom_path(
          data = skew,
          ggplot2::aes(x = position, y = content_y, colour = "GC content"),
          inherit.aes = FALSE,
          linewidth = gc_content_linewidth
        )
    }
  }

  p <- p +
    ggplot2::geom_polygon(
      data = gene_polys,
      ggplot2::aes(x = x, y = y, group = polygon_id, fill = category),
      colour = "black",
      linewidth = gene_linewidth
    )

  if (isTRUE(show_ruler)) {
    ruler <- circular_ruler_data(
      genome_length,
      major_bp = ruler_major_bp,
      minor_bp = ruler_minor_bp,
      ruler_radius = ruler_radius,
      ruler_minor_tick = ruler_minor_tick,
      ruler_major_tick = ruler_major_tick,
      ruler_label_radius = ruler_label_radius
    )

    p <- p +
      ggplot2::geom_path(
        data = ruler$ring,
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        colour = "grey35",
        linewidth = ruler_linewidth
      ) +
      ggplot2::geom_segment(
        data = ruler$ticks,
        ggplot2::aes(x = x, xend = x, y = y, yend = yend, linewidth = type),
        inherit.aes = FALSE,
        colour = "grey35"
      ) +
      ggplot2::scale_linewidth_manual(
        values = c(minor = ruler_minor_linewidth, major = ruler_major_linewidth),
        guide = "none"
      )
  }

  if (phage_topology == "linear") {
    p <- p +
      ggplot2::geom_segment(
        data = data.frame(x = 1, xend = 1, y = 0.70, yend = 1.10),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        inherit.aes = FALSE,
        linewidth = 0.75,
        colour = "black"
      )
  }

  label_data <- select_label_data(
    features,
    label_unknown = label_unknown,
    label_exclude_categories = label_exclude_categories,
    label_pattern = label_pattern,
    max_labels = max_labels,
    label_min_gap_deg = label_min_gap_deg,
    min_feature_bp = min_feature_bp,
    category_scheme = category_scheme,
    genome_length = genome_length
  )
  y_limit <- max(
    1.50,
    gene_radius + gene_height / 2 + 0.45,
    gc_skew_radius + gc_skew_height / 2 + 0.24,
    gc_content_radius + gc_content_height / 2 + 0.22,
    ruler_radius + ruler_major_tick + 0.18
  )
  label_bounds <- NULL
  if (show_labels && nrow(label_data)) {
    label_count <- nrow(label_data)
    effective_label_text_size <- if (label_count > 40L) {
      label_size_floor <- 2.2
      max(label_size_floor, min(label_text_size, label_text_size * sqrt(40 / label_count)))
    } else {
      label_text_size
    }
    effective_label_wrap_width <- if (label_count > 18L) {
      if (is.null(label_wrap_width) || !is.finite(label_wrap_width)) {
        18L
      } else {
        min(as.integer(label_wrap_width), 18L)
      }
    } else {
      label_wrap_width
    }
    effective_label_lane_spacing <- label_lane_spacing
    effective_label_anchor_radius <- max(
      label_anchor_radius,
      gene_radius + gene_height / 2,
      na.rm = TRUE
    )
    label_overlay <- circular_feature_label_overlay(
      label_data,
      genome_length = genome_length,
      inner_radius = inner_radius,
      label_line_length = label_line_length,
      label_wrap_width = effective_label_wrap_width,
      label_anchor_radius = effective_label_anchor_radius,
      label_line_anchor_radius = gene_radius,
      label_lane_spacing = effective_label_lane_spacing,
      label_lane_step = label_lane_step,
      label_text_size = effective_label_text_size,
      label_text_colour = label_text_colour,
      label_line_colour = label_line_colour,
      label_y_max = y_limit,
      feature_colors = plot_colors
    )
    label_bounds <- label_overlay_npc_bounds(label_overlay)
    if (isTRUE(legend_plot_spacing_auto)) {
      legend_plot_spacing <- auto_legend_plot_spacing(
        legend_plot_spacing,
        legend_position,
        label_bounds
      )
    }
    if (any(is.finite(label_overlay$panel_y))) {
      y_limit <- max(y_limit, max(label_overlay$panel_y, na.rm = TRUE) + 0.12)
    }
    p <- p + geom_plasmid_text_overlay(
      label_overlay,
      family = font_family,
      text_colour = "black",
      text_size = effective_label_text_size,
      fontface = "bold",
      lineheight = 0.9,
      segment_colour = label_line_colour,
      segment_linewidth = label_linewidth,
      segment_linetype = label_line_linetype
    )
  }

  if (!is.null(ruler)) {
    p <- p + geom_plasmid_text_overlay(
      circular_ruler_label_overlay(
        ruler,
        genome_length = genome_length,
        inner_radius = inner_radius
      ),
      family = font_family,
      text_colour = "grey25",
      text_size = ruler_text_size
    )
  }

  p <- p +
    geom_plasmid_text_overlay(
      center_text_overlay(name, genome_length),
      family = font_family,
      text_colour = "black",
      text_size = center_text_size,
      fontface = "bold",
      lineheight = 0.9
    ) +
    ggplot2::scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = feature_colors,
      breaks = feature_breaks,
      drop = TRUE,
      na.value = "#999999"
    ) +
    ggplot2::scale_colour_manual(
      values = gc_colors,
      breaks = names(gc_colors),
      drop = TRUE
    ) +
    ggplot2::labs(fill = NULL, colour = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 2,
        ncol = legend_columns,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.45, "cm")
      ),
      colour = ggplot2::guide_legend(
        order = 1,
        ncol = gc_legend_columns,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.70, "cm"),
        override.aes = list(alpha = 1, linewidth = gc_legend_linewidth)
      )
    ) +
    radial_coord(theta = "x", start = 0, inner.radius = inner_radius, clip = "off", expand = FALSE) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.justification = "center",
      legend.box = "vertical",
      legend.box.spacing = grid::unit(legend_plot_spacing, "cm"),
      legend.spacing.y = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(face = "bold", size = legend_text_size, colour = "black"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = circular_plot_margin(label_bounds)
    )

  attr(p, "ggplasmid_layout") <- "circular"
  attr(p, "ggplasmid_style") <- list(
    category_scheme = category_scheme,
    palette = palette,
    fill_colors = plot_colors,
    fill_breaks = feature_breaks,
    legend_columns = legend_columns,
    gc_legend_columns = gc_legend_columns,
    legend_plot_spacing = legend_plot_spacing
  )
  p
}

plot_linear_plasmid <- function(features, genome_length, name = NULL, rows = 4,
                                category_scheme = "plasmid", gc_skew = NULL,
                                show_labels = TRUE, label_unknown = FALSE,
                                label_exclude_categories = NULL,
                                label_pattern = NULL, max_labels = Inf,
                                label_min_gap_deg = 0,
                                min_feature_bp = 1,
                                palette = "default",
                                gene_highlight = NULL,
                                gene_height = 0.13,
                                gene_linewidth = 0.24,
                                gene_border_linewidth = NULL,
                                gene_arrow_head_bp = NULL,
                                gene_arrow_head_fraction = 0.55,
                                gc_skew_height = 0.11,
                                gc_content_height = 0.045,
                                gc_content_linewidth = 0.28,
                                gc_legend_linewidth = 1.8,
                                sequence_linewidth = 0.55,
                                label_text_size = 3.3,
                                label_text_colour = "black",
                                row_label_text_size = 3.4,
                                legend_position = "right",
                                legend_text_size = 8.5,
                                legend_columns = NULL,
                                gc_legend_columns = NULL,
                                legend_plot_spacing = NULL,
                                font_family = "sans") {
  gene_linewidth <- gene_border_linewidth %||% gene_linewidth
  legend_columns <- validate_legend_columns(legend_columns, legend_position)
  gc_legend_columns <- validate_gc_legend_columns(gc_legend_columns, legend_position)
  legend_plot_spacing <- validate_legend_plot_spacing(legend_plot_spacing, legend_position)
  plot_colors <- ggplasmid_resolve_colors(
    category_scheme,
    palette = palette,
    gene_highlight = gene_highlight
  )
  feature_breaks <- gene_legend_breaks(plot_colors)
  gc_colors <- gc_track_colors(plot_colors)
  linear <- linear_gene_polygons(
    features,
    genome_length,
    rows = rows,
    gene_height = gene_height,
    arrow_head_bp = gene_arrow_head_bp,
    arrow_head_fraction = gene_arrow_head_fraction
  )
  rows <- linear$rows
  gene_polys <- linear$polygons
  gene_polys$category <- factor(gene_polys$category, levels = names(plot_colors))
  row_id <- sprintf(
    "%s: %g-%g kbp",
    name %||% "sequence",
    linear$row_start / 1000,
    linear$row_end / 1000
  )
  seq_lines <- data.frame(
    x = 0,
    xend = pmax(linear$row_end - linear$row_start, 1),
    y = rows - seq_len(rows) + 1,
    yend = rows - seq_len(rows) + 1,
    label = row_id
  )

  p <- ggplot2::ggplot()

  if (!is.null(gc_skew)) {
    skew <- gc_skew
    skew$row_index <- pmin(floor((skew$position - 1L) / linear$segment) + 1L, rows)
    skew$row_y <- rows - skew$row_index + 1
    skew$x <- skew$position - linear$row_start[skew$row_index]
    skew$base_y <- skew$row_y - 0.20
    skew$pos_y <- skew$base_y + pmax(skew$scaled_score, 0) * gc_skew_height
    skew$neg_y <- skew$base_y + pmin(skew$scaled_score, 0) * gc_skew_height
    skew$content_y <- skew$base_y - 0.13 + skew$scaled_gc_content * gc_content_height
    p <- p +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = x, ymin = base_y, ymax = pos_y, group = row_index, fill = "GC skew+"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = x, ymin = base_y, ymax = neg_y, group = row_index, fill = "GC skew-"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_segment(
        data = gc_skew_legend_data(),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = track),
        inherit.aes = FALSE,
        alpha = 0,
        linewidth = gc_content_linewidth,
        show.legend = TRUE
      )
    if (any(is.finite(skew$scaled_gc_content))) {
      p <- p +
        ggplot2::geom_path(
          data = skew,
          ggplot2::aes(x = x, y = content_y, group = row_index, colour = "GC content"),
          inherit.aes = FALSE,
          linewidth = gc_content_linewidth
        )
    }
  }

  p <- p +
    ggplot2::geom_segment(
      data = seq_lines,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = sequence_linewidth,
      colour = "grey35"
    ) +
    ggplot2::geom_polygon(
      data = gene_polys,
      ggplot2::aes(x = x, y = y, group = polygon_id, fill = category),
      colour = "black",
      linewidth = gene_linewidth
    )

  label_data <- select_label_data(
    features,
    label_unknown = label_unknown,
    label_exclude_categories = label_exclude_categories,
    label_pattern = label_pattern,
    max_labels = max_labels,
    label_min_gap_deg = label_min_gap_deg,
    min_feature_bp = min_feature_bp,
    category_scheme = category_scheme,
    genome_length = genome_length
  )
  if (show_labels && nrow(label_data)) {
    label_data$row_index <- linear$row_index[match(label_data$feature_id, features$feature_id)]
    label_data$row_y <- rows - label_data$row_index + 1
    label_data$xmid_linear <- label_data$xmid - linear$row_start[label_data$row_index]
    label_data$xmid_linear <- pmin(pmax(label_data$xmid_linear, 0), linear$segment)
    label_data$label_colour <- resolve_label_text_colours(
      label_data,
      plot_colors,
      label_text_colour = label_text_colour
    )

    p <- p +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(x = xmid_linear, y = row_y + 0.12, label = label),
        inherit.aes = FALSE,
        angle = 0,
        hjust = 0,
        vjust = 0,
        direction = "x",
        nudge_y = 0.08,
        min.segment.length = 0.35,
        segment.color = "grey45",
        segment.size = 0.30,
        box.padding = 0.35,
        point.padding = 0.02,
        max.overlaps = Inf,
        seed = 20260709,
        colour = label_data$label_colour,
        size = label_text_size,
        fontface = "bold",
        lineheight = 0.82
      )
  }

  p <- p +
    ggplot2::geom_text(
      data = seq_lines,
      ggplot2::aes(x = -250, y = y - 0.42, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = row_label_text_size
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-300, linear$segment),
      ylim = c(0.15, rows + 0.75),
      clip = "off"
    ) +
    ggplot2::scale_fill_manual(
      values = plot_colors,
      breaks = feature_breaks,
      drop = TRUE,
      na.value = "#999999"
    ) +
    ggplot2::scale_colour_manual(
      values = gc_colors,
      breaks = names(gc_colors),
      drop = TRUE
    ) +
    ggplot2::labs(fill = NULL, colour = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 2,
        ncol = legend_columns,
        keyheight = grid::unit(0.32, "cm"),
        keywidth = grid::unit(0.45, "cm")
      ),
      colour = ggplot2::guide_legend(
        order = 1,
        ncol = gc_legend_columns,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.70, "cm"),
        override.aes = list(alpha = 1, linewidth = gc_legend_linewidth)
      )
    ) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.box = "vertical",
      legend.box.spacing = grid::unit(legend_plot_spacing, "cm"),
      legend.spacing.y = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(face = "bold", size = legend_text_size, colour = "black"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )

  attr(p, "ggplasmid_layout") <- "linear"
  attr(p, "ggplasmid_style") <- list(
    category_scheme = category_scheme,
    palette = palette,
    fill_colors = plot_colors,
    fill_breaks = feature_breaks,
    legend_columns = legend_columns,
    gc_legend_columns = gc_legend_columns,
    legend_plot_spacing = legend_plot_spacing
  )
  p
}

#' Draw a plasmid map
#'
#' @param annotation A data frame or path to a CSV/TSV annotation table.
#' @param gbk Optional GenBank file.
#' @param fasta Optional FASTA file for GC-skew calculation and length inference.
#' @param skew_table Optional precomputed GC-skew table with `position` and `gc_skew`.
#' @param output Optional output path. When supplied, the plot is saved and
#'   invisibly returned.
#' @param layout Either `"circular"` or `"linear"`.
#' @param phage_topology Either `"circular"` or `"linear"`. Only used when
#'   `layout = "circular"`; `"linear"` adds a terminal boundary line.
#' @param name Display name.
#' @param genome_length Optional plasmid/genome length.
#' @param feature_types GenBank feature types to keep.
#' @param label_mode `"auto"` prefers compact gene labels, `"product"` uses
#'   product names, and `"gene"` prefers gene IDs.
#' @param category_scheme Either `"plasmid"` or `"phage"`.
#' @param palette Color palette. `"default"` uses ggplasmid colors; when
#'   `ggsci` is installed, common ggsci palette names such as `"npg"`,
#'   `"aaas"`, `"lancet"`, `"jco"`, `"ucscgb"`, `"d3"`, and `"igv"` can be used.
#' @param gene_highlight Optional named character vector, or data frame with
#'   category/color columns, used to override specific category colors.
#' @param rows Number of rows for linear layout.
#' @param show_gc_skew Whether to draw GC-skew and GC-content tracks when
#'   sequence/skew data are available.
#' @param show_labels Whether to draw feature labels.
#' @param label_unknown Whether to label unknown/hypothetical features.
#' @param label_exclude_categories Optional category names whose features should
#'   still be drawn but not labeled.
#' @param label_pattern Optional regular expression. Matching features are
#'   prioritized for labels.
#' @param max_labels Maximum number of labels to consider. The default is `Inf`,
#'   so all eligible labels are considered. Circular labels are placed from
#'   shorter to longer text on farther outer lanes when needed. Use a finite
#'   value to label only the top eligible features.
#' @param label_min_gap_deg Minimum angular gap between circular labels. Larger
#'   values omit more crowded labels. Set to `0` to disable angular thinning.
#' @param label_line_length Distance from the feature anchor to radial labels,
#'   in circular plot radius units.
#' @param mini_label_line_length Clearer alias for `label_line_length`. If
#'   supplied, it overrides `label_line_length`.
#' @param label_wrap_width Number of characters per label line before wrapping.
#'   Set to `Inf` or `NULL` to disable wrapping.
#' @param label_lane_spacing Extra spacing for circular labels when nearby
#'   labels need a small outward offset.
#' @param label_lane_step Optional direct radial step between circular label
#'   lanes, in circular plot radius units. Smaller values pack lanes tighter.
#'   When `NULL`, the step is calculated from label height.
#' @param label_line_colour,label_linewidth,label_line_linetype Circular label
#'   leader line colour, width, and linetype.
#' @param min_feature_bp Minimum feature length to label.
#' @param show_ruler Whether to draw the circular clock/ruler track. Only used
#'   when `layout = "circular"`.
#' @param ruler_major_bp,ruler_minor_bp Optional major/minor ruler tick spacing
#'   in bp. When omitted, sensible values are chosen from genome length.
#' @param gene_radius,gene_height Circular gene-ring radius and gene arrow
#'   height. For linear maps, `gene_height` controls arrow height.
#' @param gene_linewidth,gene_border_linewidth Outline width for gene arrows.
#'   `gene_border_linewidth` is a clearer alias; if supplied, it overrides
#'   `gene_linewidth`.
#' @param gene_arrow_head_bp Optional gene arrow-head length in bp. The default
#'   scales with genome length.
#' @param gene_arrow_head_fraction Maximum fraction of each feature that can be
#'   used as the arrow head.
#' @param gene_gap_bp Gap removed from adjacent circular feature arrows, in bp.
#'   Use `0` for no gap.
#' @param gc_skew_radius,gc_skew_height Circular GC-skew track radius and
#'   height.
#' @param gc_content_radius,gc_content_height Circular GC-content track radius
#'   and height. For linear maps, `gc_content_height` controls line amplitude.
#' @param gc_content_linewidth Line width for the GC-content track.
#' @param gc_legend_linewidth Line width for GC content/skew entries in the
#'   legend.
#' @param ruler_minor_tick,ruler_major_tick,ruler_label_radius Circular ruler
#'   tick lengths and label radius.
#' @param ruler_linewidth,ruler_minor_linewidth,ruler_major_linewidth Circular
#'   ruler line widths.
#' @param sequence_linewidth Linear layout baseline width.
#' @param label_text_size,ruler_text_size,center_text_size,legend_text_size Text
#'   sizes.
#' @param label_text_colour Label text colour. Use `"black"` for one fixed
#'   colour, `"category"` to match each label to its gene category colour, or a
#'   named vector such as `c("Antimicrobial resistance" = "#B2182B")`.
#' @param label_anchor_radius Starting radius for outside circular labels.
#'   Increase it to push labels farther from the gene ring.
#' @param row_label_text_size Text size for linear row labels.
#' @param legend_position Legend position passed to ggplot2, such as
#'   `"bottom"`, `"right"`, `"left"`, `"top"`, or `"none"`.
#' @param legend_columns Number of columns in the feature-category legend. When
#'   omitted, right/left legends use one column and top/bottom legends use
#'   three columns.
#' @param gc_legend_columns Number of columns in the GC content/skew legend.
#'   When omitted, right/left legends use one column and top/bottom legends use
#'   three columns.
#' @param legend_plot_spacing Space in cm between the plot panel and legend.
#'   When omitted, right/left legends are moved outward based on the measured
#'   outer label range.
#' @param inner_radius Inner radius passed to `coord_radial()` for circular
#'   maps. Larger values make a larger center hole and push tracks outward.
#'   Use `NULL` to choose a layout-aware default.
#' @param font_family Font family passed to ggplot2 text themes. The default
#'   `"sans"` avoids Windows font database warnings.
#' @param gc_window GC-skew window.
#' @param gc_step GC-skew step. Use `NULL` to choose an automatic step that is
#'   faster for long genomes.
#' @param width,height,dpi Output size for `output`.
#' @return A ggplot object. Use [plasmid_data()] to retrieve the parsed
#'   annotation, FASTA, normalized features, and GC-skew data used by the plot.
#' @export
ggplasmid <- function(annotation = NULL, gbk = NULL, fasta = NULL, skew_table = NULL,
                      output = NULL, layout = c("circular", "linear"),
                      phage_topology = c("circular", "linear"),
                      name = NULL, genome_length = NULL,
                      feature_types = c(
                        "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                        "rep_origin", "mobile_element"
                      ),
                      label_mode = c("auto", "product", "gene"),
                      category_scheme = c("plasmid", "phage"),
                      palette = "default", gene_highlight = NULL,
                      rows = 4, show_gc_skew = TRUE, show_labels = TRUE,
                      label_unknown = FALSE, label_pattern = NULL,
                      label_exclude_categories = NULL,
                      max_labels = NULL, label_min_gap_deg = NULL,
                      label_line_length = 0.02,
                      mini_label_line_length = NULL,
                      label_wrap_width = 28L,
                      label_anchor_radius = 1.07,
                      label_lane_spacing = 0.24,
                      label_lane_step = NULL,
                      label_line_colour = "grey70",
                      label_linewidth = 0.30,
                      label_line_linetype = "dashed",
                      min_feature_bp = 1,
                      show_ruler = TRUE,
                      ruler_major_bp = NULL, ruler_minor_bp = NULL,
                      gene_radius = 1,
                      gene_height = NULL,
                      gene_linewidth = NULL,
                      gene_border_linewidth = NULL,
                      gene_arrow_head_bp = NULL,
                      gene_arrow_head_fraction = 0.55,
                      gene_gap_bp = NULL,
                      gc_skew_radius = 0.78,
                      gc_skew_height = NULL,
                      gc_content_radius = 0.67,
                      gc_content_height = NULL,
                      gc_content_linewidth = 0.28,
                      gc_legend_linewidth = 1.8,
                      ruler_minor_tick = 0.015,
                      ruler_major_tick = 0.035,
                      ruler_label_radius = 0.47,
                      ruler_linewidth = 0.28,
                      ruler_minor_linewidth = 0.18,
                      ruler_major_linewidth = 0.38,
                      sequence_linewidth = 0.55,
                      label_text_size = NULL,
                      label_text_colour = "black",
                      row_label_text_size = 3.4,
                      ruler_text_size = 2.25,
                      center_text_size = 4.3,
                      legend_position = NULL,
                      legend_text_size = 8.5,
                      legend_columns = NULL,
                      gc_legend_columns = NULL,
                      legend_plot_spacing = NULL,
                      inner_radius = NULL,
                      font_family = "sans",
                      gc_window = 100L, gc_step = NULL,
                      width = NULL, height = NULL, dpi = 300) {
  layout <- match.arg(layout)
  phage_topology <- match.arg(phage_topology)
  label_mode <- match.arg(label_mode)
  category_scheme <- match.arg(category_scheme)
  if (!is.null(mini_label_line_length)) {
    label_line_length <- mini_label_line_length
  }

  annotation_data <- read_plasmid_annotation(
    annotation = annotation,
    gbk = gbk,
    feature_types = feature_types
  )
  fasta_data <- if (!is.null(fasta)) read_plasmid_fasta(fasta) else NULL

  features <- prepare_plasmid_features(
    annotation = annotation_data,
    fasta = fasta_data,
    genome_length = genome_length,
    feature_types = feature_types,
    label_mode = label_mode,
    category_scheme = category_scheme
  )
  genome_length <- attr(features, "genome_length", exact = TRUE)
  name <- name %||% attr(features, "name", exact = TRUE) %||% "plasmid"
  sequence <- attr(features, "sequence", exact = TRUE) %||% ""
  if (is.null(max_labels)) {
    max_labels <- Inf
  }
  if (is.null(label_min_gap_deg)) {
    label_min_gap_deg <- 0
  }
  if (is.null(inner_radius)) {
    inner_radius <- if (layout == "circular" && category_scheme == "phage") {
      0.14
    } else {
      0.38
    }
  }
  if (is.null(gene_height)) {
    gene_height <- if (layout == "circular" && category_scheme == "phage") {
      0.085
    } else if (layout == "circular") {
      0.10
    } else {
      0.13
    }
  }
  if (is.null(gene_linewidth)) {
    gene_linewidth <- if (layout == "circular" && category_scheme == "phage") {
      0.03
    } else if (layout == "circular") {
      0.22
    } else {
      0.24
    }
  }
  if (!is.null(gene_border_linewidth)) {
    gene_linewidth <- gene_border_linewidth
  }
  if (is.null(gene_arrow_head_bp) && layout == "circular" && category_scheme == "phage") {
    gene_arrow_head_bp <- max(60, genome_length / 700)
  }
  if (is.null(gene_gap_bp)) {
    gene_gap_bp <- if (layout == "circular" && category_scheme == "phage") {
      max(25, genome_length / 5000)
    } else {
      0
    }
  }
  if (is.null(gc_skew_height)) {
    gc_skew_height <- if (layout == "circular") 0.10 else 0.11
  }
  if (is.null(gc_content_height)) {
    gc_content_height <- 0.045
  }
  if (is.null(label_text_size)) {
    label_text_size <- if (layout == "circular") 3.4 else 3.3
  }
  if (is.null(legend_position)) {
    legend_position <- if (layout == "circular") "bottom" else "right"
  }
  legend_columns <- validate_legend_columns(legend_columns, legend_position)
  gc_legend_columns <- validate_gc_legend_columns(gc_legend_columns, legend_position)

  gc_skew <- NULL
  if (isTRUE(show_gc_skew)) {
    gc_skew <- prepare_gc_skew_for_plot(
      skew_table = skew_table,
      fasta = fasta_data,
      sequence = sequence,
      genome_length = genome_length,
      circular = layout == "circular",
      window = gc_window,
      step = gc_step
    )
  }

  p <- if (layout == "circular") {
    plot_circular_plasmid(
      features = features,
      genome_length = genome_length,
      name = name,
      phage_topology = phage_topology,
      category_scheme = category_scheme,
      gc_skew = gc_skew,
      show_labels = show_labels,
      label_unknown = label_unknown,
      label_exclude_categories = label_exclude_categories,
      label_pattern = label_pattern,
      max_labels = max_labels,
      label_min_gap_deg = label_min_gap_deg,
      label_line_length = label_line_length,
      mini_label_line_length = mini_label_line_length,
      label_wrap_width = label_wrap_width,
      label_anchor_radius = label_anchor_radius,
      label_lane_spacing = label_lane_spacing,
      label_lane_step = label_lane_step,
      label_line_colour = label_line_colour,
      label_linewidth = label_linewidth,
      label_line_linetype = label_line_linetype,
      min_feature_bp = min_feature_bp,
      show_ruler = show_ruler,
      ruler_major_bp = ruler_major_bp,
      ruler_minor_bp = ruler_minor_bp,
      palette = palette,
      gene_highlight = gene_highlight,
      gene_radius = gene_radius,
      gene_height = gene_height,
      gene_linewidth = gene_linewidth,
      gene_border_linewidth = gene_border_linewidth,
      gene_arrow_head_bp = gene_arrow_head_bp,
      gene_arrow_head_fraction = gene_arrow_head_fraction,
      gene_gap_bp = gene_gap_bp,
      gc_skew_radius = gc_skew_radius,
      gc_skew_height = gc_skew_height,
      gc_content_radius = gc_content_radius,
      gc_content_height = gc_content_height,
      gc_content_linewidth = gc_content_linewidth,
      gc_legend_linewidth = gc_legend_linewidth,
      ruler_radius = gc_content_radius - 0.08,
      ruler_minor_tick = ruler_minor_tick,
      ruler_major_tick = ruler_major_tick,
      ruler_label_radius = ruler_label_radius,
      ruler_linewidth = ruler_linewidth,
      ruler_minor_linewidth = ruler_minor_linewidth,
      ruler_major_linewidth = ruler_major_linewidth,
      label_text_size = label_text_size,
      label_text_colour = label_text_colour,
      ruler_text_size = ruler_text_size,
      center_text_size = center_text_size,
      legend_position = legend_position,
      legend_text_size = legend_text_size,
      legend_columns = legend_columns,
      gc_legend_columns = gc_legend_columns,
      legend_plot_spacing = legend_plot_spacing,
      inner_radius = inner_radius,
      font_family = font_family
    )
  } else {
    plot_linear_plasmid(
      features = features,
      genome_length = genome_length,
      name = name,
      rows = rows,
      category_scheme = category_scheme,
      gc_skew = gc_skew,
      show_labels = show_labels,
      label_unknown = label_unknown,
      label_exclude_categories = label_exclude_categories,
      label_pattern = label_pattern,
      max_labels = max_labels,
      label_min_gap_deg = label_min_gap_deg,
      min_feature_bp = min_feature_bp,
      palette = palette,
      gene_highlight = gene_highlight,
      gene_height = gene_height,
      gene_linewidth = gene_linewidth,
      gene_border_linewidth = gene_border_linewidth,
      gene_arrow_head_bp = gene_arrow_head_bp,
      gene_arrow_head_fraction = gene_arrow_head_fraction,
      gc_skew_height = gc_skew_height,
      gc_content_height = gc_content_height,
      gc_content_linewidth = gc_content_linewidth,
      gc_legend_linewidth = gc_legend_linewidth,
      sequence_linewidth = sequence_linewidth,
      label_text_size = label_text_size,
      label_text_colour = label_text_colour,
      row_label_text_size = row_label_text_size,
      legend_position = legend_position,
      legend_text_size = legend_text_size,
      legend_columns = legend_columns,
      gc_legend_columns = gc_legend_columns,
      legend_plot_spacing = legend_plot_spacing,
      font_family = font_family
    )
  }

  if (is.null(fasta_data) && nzchar(sequence)) {
    fasta_data <- data.frame(
      name = name,
      description = name,
      sequence = sequence,
      length = nchar(sequence),
      stringsAsFactors = FALSE
    )
  }
  attr(p, "ggplasmid_data") <- list(
    annotation = annotation_data,
    fasta = fasta_data,
    features = features,
    gc_skew = gc_skew,
    genome_length = genome_length,
    name = name
  )

  if (!is.null(output)) {
    save_plasmid_map(p, output, width = width, height = height, dpi = dpi)
    return(invisible(p))
  }
  p
}

#' Draw a plasmid map
#'
#' Alias for [ggplasmid()].
#' @inheritParams ggplasmid
#' @return A ggplot object.
#' @export
plot_plasmid_map <- ggplasmid

#' Retrieve data used by a ggplasmid plot
#'
#' @param plot A plot returned by [ggplasmid()] or [plot_phage_map()].
#' @return A list containing `annotation`, `fasta`, `features`, `gc_skew`,
#'   `genome_length`, and `name`.
#' @export
plasmid_data <- function(plot) {
  data <- attr(plot, "ggplasmid_data", exact = TRUE)
  if (is.null(data)) {
    stop("`plot` does not contain ggplasmid data.", call. = FALSE)
  }
  data
}

#' Draw a phage map
#'
#' A light wrapper around [ggplasmid()] that uses the phage category scheme.
#' `phage_topology = "linear"` only adds a terminal boundary line when
#' `layout = "circular"`.
#'
#' @inheritParams ggplasmid
#' @return A ggplot object.
#' @export
plot_phage_map <- function(annotation = NULL, gbk = NULL, fasta = NULL,
                           skew_table = NULL, output = NULL,
                           layout = c("circular", "linear"),
                           phage_topology = c("circular", "linear"),
                           name = NULL, genome_length = NULL,
                           feature_types = c(
                             "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                             "rep_origin", "mobile_element"
                           ),
                           label_mode = c("auto", "product", "gene"),
                           palette = "default", gene_highlight = NULL,
                           rows = 4, show_gc_skew = TRUE,
                           show_labels = TRUE,
                           label_unknown = FALSE,
                           label_exclude_categories = NULL,
                           label_pattern = NULL, max_labels = NULL,
                           label_min_gap_deg = NULL,
                           label_line_length = 0.02,
                           mini_label_line_length = NULL,
                           min_feature_bp = 1,
                           label_wrap_width = 28L,
                           label_anchor_radius = 1.07,
                           label_lane_spacing = 0.24,
                           label_lane_step = NULL,
                           label_line_colour = "grey70",
                           label_linewidth = 0.30,
                           label_line_linetype = "dashed",
                           show_ruler = TRUE,
                           ruler_major_bp = NULL, ruler_minor_bp = NULL,
                           gene_radius = 1,
                           gene_height = NULL,
                           gene_linewidth = NULL,
                           gene_border_linewidth = NULL,
                           gene_arrow_head_bp = NULL,
                           gene_arrow_head_fraction = 0.55,
                           gene_gap_bp = NULL,
                           gc_skew_radius = 0.78,
                           gc_skew_height = NULL,
                           gc_content_radius = 0.67,
                           gc_content_height = NULL,
                           gc_content_linewidth = 0.28,
                           gc_legend_linewidth = 1.8,
                           ruler_minor_tick = 0.015,
                           ruler_major_tick = 0.035,
                           ruler_label_radius = 0.47,
                           ruler_linewidth = 0.28,
                           ruler_minor_linewidth = 0.18,
                           ruler_major_linewidth = 0.38,
                           sequence_linewidth = 0.55,
                           label_text_size = NULL,
                           label_text_colour = "black",
                           row_label_text_size = 3.4,
                           ruler_text_size = 2.25,
                           center_text_size = 4.3,
                           legend_position = NULL,
                           legend_text_size = 8.5,
                           legend_columns = NULL,
                           gc_legend_columns = NULL,
                           legend_plot_spacing = NULL,
                           inner_radius = NULL,
                           font_family = "sans",
                           gc_window = 100L, gc_step = NULL,
                           width = NULL, height = NULL, dpi = 300) {
  ggplasmid(
    annotation = annotation,
    gbk = gbk,
    fasta = fasta,
    skew_table = skew_table,
    output = output,
    layout = match.arg(layout),
    phage_topology = match.arg(phage_topology),
    name = name,
    genome_length = genome_length,
    feature_types = feature_types,
    label_mode = match.arg(label_mode),
    category_scheme = "phage",
    palette = palette,
    gene_highlight = gene_highlight,
    rows = rows,
    show_gc_skew = show_gc_skew,
    show_labels = show_labels,
    label_unknown = label_unknown,
    label_exclude_categories = label_exclude_categories,
    label_pattern = label_pattern,
    max_labels = max_labels,
    label_min_gap_deg = label_min_gap_deg,
    label_line_length = label_line_length,
    mini_label_line_length = mini_label_line_length,
    label_wrap_width = label_wrap_width,
    label_anchor_radius = label_anchor_radius,
    label_lane_spacing = label_lane_spacing,
    label_lane_step = label_lane_step,
    label_line_colour = label_line_colour,
    label_linewidth = label_linewidth,
    label_line_linetype = label_line_linetype,
    min_feature_bp = min_feature_bp,
    show_ruler = show_ruler,
    ruler_major_bp = ruler_major_bp,
    ruler_minor_bp = ruler_minor_bp,
    gene_radius = gene_radius,
    gene_height = gene_height,
    gene_linewidth = gene_linewidth,
    gene_border_linewidth = gene_border_linewidth,
    gene_arrow_head_bp = gene_arrow_head_bp,
    gene_arrow_head_fraction = gene_arrow_head_fraction,
    gene_gap_bp = gene_gap_bp,
    gc_skew_radius = gc_skew_radius,
    gc_skew_height = gc_skew_height,
    gc_content_radius = gc_content_radius,
    gc_content_height = gc_content_height,
    gc_content_linewidth = gc_content_linewidth,
    gc_legend_linewidth = gc_legend_linewidth,
    ruler_minor_tick = ruler_minor_tick,
    ruler_major_tick = ruler_major_tick,
    ruler_label_radius = ruler_label_radius,
    ruler_linewidth = ruler_linewidth,
    ruler_minor_linewidth = ruler_minor_linewidth,
    ruler_major_linewidth = ruler_major_linewidth,
    sequence_linewidth = sequence_linewidth,
    label_text_size = label_text_size,
    label_text_colour = label_text_colour,
    row_label_text_size = row_label_text_size,
    ruler_text_size = ruler_text_size,
    center_text_size = center_text_size,
    legend_position = legend_position,
    legend_text_size = legend_text_size,
    legend_columns = legend_columns,
    gc_legend_columns = gc_legend_columns,
    legend_plot_spacing = legend_plot_spacing,
    inner_radius = inner_radius,
    font_family = font_family,
    gc_window = gc_window,
    gc_step = gc_step,
    width = width,
    height = height,
    dpi = dpi
  )
}

#' Save a plasmid map
#'
#' @param plot A ggplot object returned by [ggplasmid()].
#' @param filename Output filename.
#' @param width,height Plot size in inches. Defaults depend on layout.
#' @param dpi Output resolution.
#' @param bg Background color.
#' @return The normalized output path, invisibly.
#' @export
save_plasmid_map <- function(plot, filename, width = NULL, height = NULL,
                             dpi = 300, bg = "white") {
  layout <- attr(plot, "ggplasmid_layout", exact = TRUE) %||% "circular"
  if (is.null(width)) {
    width <- if (layout == "linear") 18 else 9
  }
  if (is.null(height)) {
    height <- if (layout == "linear") 12 else 9
  }
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  invisible(normalizePath(filename, mustWork = FALSE))
}
