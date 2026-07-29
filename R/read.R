parse_genbank_location <- function(location, genome_length = NA_real_) {
  loc <- gsub("\\s+", "", clean_text(location))
  strand <- if (grepl("complement\\(", loc, ignore.case = TRUE)) "-" else "+"

  range_matches <- gregexpr("<?([0-9]+)\\.\\.>?([0-9]+)", loc, perl = TRUE)
  ranges_text <- regmatches(loc, range_matches)[[1]]

  if (length(ranges_text) && ranges_text[[1]] != "-1") {
    ranges <- do.call(
      rbind,
      lapply(ranges_text, function(piece) {
        nums <- as_number(regmatches(piece, gregexpr("[0-9]+", piece))[[1]])
        c(nums[[1]], nums[[2]])
      })
    )
  } else {
    nums <- as_number(regmatches(loc, gregexpr("[0-9]+", loc))[[1]])
    ranges <- cbind(nums, nums)
  }

  if (!length(ranges)) {
    return(list(start = NA_real_, end = NA_real_, strand = strand, wraps_origin = FALSE))
  }

  ranges <- matrix(as.numeric(t(ranges)), ncol = 2, byrow = TRUE)
  wraps_origin <- FALSE
  if (nrow(ranges) > 1 && is.finite(genome_length)) {
    wraps_origin <- ranges[1, 1] > ranges[nrow(ranges), 2]
  }

  if (wraps_origin) {
    start <- ranges[1, 1]
    end <- ranges[nrow(ranges), 2]
  } else {
    start <- min(ranges, na.rm = TRUE)
    end <- max(ranges, na.rm = TRUE)
  }

  list(start = start, end = end, strand = strand, wraps_origin = wraps_origin)
}

parse_genbank_qualifier <- function(value) {
  value <- clean_text(value)
  value <- sub('^"', "", value)
  value <- sub('"$', "", value)
  value
}

read_genbank_file <- function(gbk, feature_types) {
  lines <- readLines(gbk, warn = FALSE)
  locus <- grep("^LOCUS\\b", lines, value = TRUE)
  locus <- first_or(locus, "")
  locus_parts <- strsplit(clean_text(locus), "\\s+")[[1]]
  record_name <- if (length(locus_parts) >= 2) {
    locus_parts[[2]]
  } else {
    tools::file_path_sans_ext(basename(gbk))
  }

  length_match <- regmatches(locus, regexpr("[0-9]+(?=\\s+bp\\b)", locus, perl = TRUE))
  genome_length <- suppressWarnings(as.numeric(first_or(length_match, NA_real_)))

  origin_index <- grep("^ORIGIN\\b", lines)
  sequence <- ""
  if (length(origin_index)) {
    end_index <- grep("^//", lines)
    end_candidates <- end_index[end_index > origin_index[[1]]]
    end_index <- first_or(end_candidates, length(lines))
    seq_lines <- lines[(origin_index[[1]] + 1):(end_index - 1)]
    sequence <- toupper(gsub("[^A-Za-z]", "", paste(seq_lines, collapse = "")))
    sequence <- gsub("[^ACGTN]", "", sequence)
    if (nzchar(sequence) && !is.finite(genome_length)) {
      genome_length <- nchar(sequence)
    }
  }

  feature_start <- grep("^FEATURES\\b", lines)
  if (!length(feature_start)) {
    stop("No FEATURES section found in GenBank file: ", gbk, call. = FALSE)
  }
  feature_end_candidates <- grep("^ORIGIN\\b|^BASE COUNT\\b|^//", lines)
  feature_end_candidates <- feature_end_candidates[feature_end_candidates > feature_start[[1]]]
  feature_end <- first_or(feature_end_candidates, length(lines))
  feature_lines <- lines[(feature_start[[1]] + 1):(feature_end - 1)]

  features <- list()
  current <- NULL
  last_key <- NULL

  flush_current <- function() {
    if (is.null(current) || !(current$type %in% feature_types)) {
      return(NULL)
    }
    loc <- parse_genbank_location(current$location, genome_length = genome_length)
    quals <- current$qualifiers
    gene_id <- quals$gene %||% quals$locus_tag %||% quals$protein_id %||% ""
    label <- quals$product %||%
      quals$gene %||%
      quals$mobile_element_type %||%
      quals$note %||%
      current$type
    data.frame(
      gene_id = clean_text(gene_id),
      type = current$type,
      start = loc$start,
      end = loc$end,
      strand = loc$strand,
      label = clean_text(label),
      product = clean_text(quals$product %||% label),
      location = current$location,
      wraps_origin = loc$wraps_origin,
      stringsAsFactors = FALSE
    )
  }

  for (line in feature_lines) {
    feature_match <- regmatches(
      line,
      regexec("^     (\\S+)\\s+(.+)$", line, perl = TRUE)
    )[[1]]

    if (length(feature_match) == 3 && !startsWith(trimws(feature_match[[2]]), "/")) {
      flushed <- flush_current()
      if (!is.null(flushed)) {
        features[[length(features) + 1]] <- flushed
      }
      current <- list(
        type = feature_match[[2]],
        location = clean_text(feature_match[[3]]),
        qualifiers = list()
      )
      last_key <- NULL
      next
    }

    if (is.null(current)) {
      next
    }

    qualifier_match <- regmatches(
      line,
      regexec("^\\s+/([^=]+)=(.*)$", line, perl = TRUE)
    )[[1]]
    if (length(qualifier_match) == 3) {
      last_key <- qualifier_match[[2]]
      current$qualifiers[[last_key]] <- parse_genbank_qualifier(qualifier_match[[3]])
      next
    }

    flag_match <- regmatches(line, regexec("^\\s+/([^=]+)\\s*$", line, perl = TRUE))[[1]]
    if (length(flag_match) == 2) {
      last_key <- NULL
      next
    }

    continuation <- clean_text(line)
    if (nzchar(continuation)) {
      if (!is.null(last_key)) {
        current$qualifiers[[last_key]] <- clean_text(
          paste(current$qualifiers[[last_key]], parse_genbank_qualifier(continuation))
        )
      } else {
        current$location <- clean_text(paste0(current$location, continuation))
      }
    }
  }

  flushed <- flush_current()
  if (!is.null(flushed)) {
    features[[length(features) + 1]] <- flushed
  }

  if (!length(features)) {
    stop("No supported features found in GenBank file: ", gbk, call. = FALSE)
  }

  feature_data <- do.call(rbind, features)
  attr(feature_data, "genome_length") <- genome_length
  attr(feature_data, "sequence") <- sequence
  attr(feature_data, "name") <- record_name
  feature_data
}

#' Format FASTA records
#'
#' @param fasta Path to a FASTA file, or a data frame previously returned by
#'   `read_plasmid_fasta()`.
#' @param width Sequence characters per line.
#' @return A character vector with one FASTA-formatted record per element.
#' @export
format_plasmid_fasta <- function(fasta, width = 70L) {
  fasta_data <- if (is.data.frame(fasta)) fasta else read_plasmid_fasta(fasta)
  width <- as.integer(width)
  if (!is.finite(width) || width < 1L) {
    stop("`width` must be a positive integer.", call. = FALSE)
  }
  vapply(
    seq_len(nrow(fasta_data)),
    function(i) {
      sequence <- fasta_data$sequence[[i]]
      starts <- seq.int(1L, nchar(sequence), by = width)
      sequence_lines <- substring(sequence, starts, pmin(starts + width - 1L, nchar(sequence)))
      paste(c(paste0(">", fasta_data$description[[i]]), sequence_lines), collapse = "\n")
    },
    character(1)
  )
}

#' Read a FASTA file as a data frame
#'
#' @param fasta Path to a FASTA file, or a data frame previously returned by
#'   `read_plasmid_fasta()`.
#' @return A data frame with one row per FASTA record and columns `name`,
#'   `description`, `sequence`, `length`, and `fasta`.
#' @export
read_plasmid_fasta <- function(fasta) {
  if (is.data.frame(fasta)) {
    sequence_col <- find_column(fasta, "sequence", label = "FASTA sequence column")
    out <- as.data.frame(fasta, stringsAsFactors = FALSE)
    out$sequence <- toupper(gsub("[^ACGTN]", "", out[[sequence_col]]))
    if (!"name" %in% names(out)) {
      out$name <- paste0("sequence_", seq_len(nrow(out)))
    }
    if (!"description" %in% names(out)) {
      out$description <- out$name
    }
    out$length <- nchar(out$sequence)
    out$fasta <- format_plasmid_fasta(out[c("name", "description", "sequence", "length")])
    return(out[c("name", "description", "sequence", "length", "fasta")])
  }

  if (!is.character(fasta) || length(fasta) != 1L || !file.exists(fasta)) {
    stop("FASTA file not found: ", first_or(fasta, ""), call. = FALSE)
  }

  lines <- readLines(fasta, warn = FALSE)
  header_index <- which(startsWith(lines, ">"))
  if (!length(header_index)) {
    stop("No FASTA records found in: ", fasta, call. = FALSE)
  }
  record_end <- c(header_index[-1L] - 1L, length(lines))

  records <- lapply(seq_along(header_index), function(i) {
    header <- sub("^>", "", lines[[header_index[[i]]]])
    name <- first_or(strsplit(clean_text(header), "\\s+")[[1]], paste0("sequence_", i))
    sequence_lines <- if (record_end[[i]] > header_index[[i]]) {
      lines[(header_index[[i]] + 1L):record_end[[i]]]
    } else {
      character()
    }
    sequence <- toupper(gsub("[^ACGTN]", "", paste(sequence_lines, collapse = "")))
    data.frame(
      name = name,
      description = clean_text(header),
      sequence = sequence,
      length = nchar(sequence),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, records)
  row.names(out) <- NULL
  if (!all(nzchar(out$sequence))) {
    stop("One or more FASTA records contain no DNA sequence: ", fasta, call. = FALSE)
  }
  out$fasta <- format_plasmid_fasta(out)
  out
}

#' Read a FASTA file
#'
#' Convenience alias for [read_plasmid_fasta()].
#'
#' @inheritParams read_plasmid_fasta
#' @return A data frame with one row per FASTA record and columns `name`,
#'   `description`, `sequence`, `length`, and `fasta`.
#' @export
read_fasta <- function(fasta) {
  read_plasmid_fasta(fasta)
}

#' Read a plasmid annotation table or GenBank file
#'
#' @param annotation A data frame or path to a CSV/TSV annotation table.
#' @param gbk Optional path to a GenBank file.
#' @param feature_types GenBank feature types to keep.
#' @return A data frame of features.
#' @export
read_plasmid_annotation <- function(annotation = NULL, gbk = NULL,
                                    feature_types = c(
                                      "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                                      "rep_origin", "mobile_element"
                                    )) {
  if (!is.null(gbk)) {
    return(read_genbank_file(gbk, feature_types = feature_types))
  }
  if (is.data.frame(annotation)) {
    return(as.data.frame(annotation, stringsAsFactors = FALSE))
  }
  if (is.character(annotation) && length(annotation) == 1) {
    if (!file.exists(annotation)) {
      stop("Annotation table not found: ", annotation, call. = FALSE)
    }
    return(read_delim_auto(annotation))
  }
  stop("Provide either `annotation` as a data frame/path or `gbk` as a GenBank path.", call. = FALSE)
}

#' Read an annotation table
#'
#' Convenience wrapper around [read_plasmid_annotation()] for CSV/TSV tables or
#' already loaded data frames.
#'
#' @param annotation A data frame or path to a CSV/TSV annotation table.
#' @return A data frame of features.
#' @export
read_annotation_table <- function(annotation) {
  read_plasmid_annotation(annotation = annotation)
}

#' Read a GenBank file
#'
#' Convenience wrapper around [read_plasmid_annotation()] for GenBank input.
#'
#' @param gbk Path to a GenBank file.
#' @param feature_types GenBank feature types to keep.
#' @return A data frame of features parsed from the GenBank FEATURES section.
#' @export
read_gbk <- function(gbk, feature_types = c(
                       "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                       "rep_origin", "mobile_element"
                     )) {
  read_plasmid_annotation(gbk = gbk, feature_types = feature_types)
}
