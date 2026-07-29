#' Compute GC skew from a DNA sequence
#'
#' @param fasta Optional FASTA file.
#' @param sequence Optional DNA sequence string. Ignored when `fasta` is supplied.
#' @param genome_length Optional length to use from the start of the sequence.
#' @param window Window size in bp.
#' @param step Step size in bp. Use `NULL` to choose an automatic step that
#'   keeps the plotted GC track compact for long genomes.
#' @param circular Whether windows should wrap around the sequence ends.
#' @return A data frame with `position`, `gc_skew`, `gc_content`, `window`,
#'   and `step`.
#' @export
compute_gc_skew <- function(fasta = NULL, sequence = NULL, genome_length = NULL,
                            window = 100L, step = NULL, circular = TRUE) {
  if (!is.null(fasta)) {
    sequence <- read_fasta_sequence(fasta)
  }
  if (is.null(sequence) || !nzchar(sequence)) {
    stop("Provide `fasta` or `sequence` for GC-skew calculation.", call. = FALSE)
  }

  sequence <- toupper(gsub("[^ACGTN]", "", sequence))
  n <- nchar(sequence)
  if (!is.finite(genome_length %||% NA_real_)) {
    genome_length <- n
  }
  genome_length <- min(as.integer(genome_length), n)
  window <- as.integer(window)
  if (!is.finite(window) || window < 1L) stop("`window` must be a positive integer.", call. = FALSE)
  if (is.null(step)) {
    step <- max(1L, ceiling(genome_length / 3000L))
  }
  step <- as.integer(step)
  if (!is.finite(step) || step < 1L) stop("`step` must be a positive integer.", call. = FALSE)

  bases <- strsplit(substr(sequence, 1, genome_length), "", fixed = TRUE)[[1]]
  left <- (window - 1L) %/% 2L
  right <- window - left - 1L

  if (circular) {
    prefix <- if (left > 0L) bases[((genome_length - left + 1L):genome_length - 1L) %% genome_length + 1L] else character()
    suffix <- if (right > 0L) bases[((seq_len(right) - 1L) %% genome_length) + 1L] else character()
    extended <- c(prefix, bases, suffix)
    start_index <- seq_len(genome_length)
  } else {
    extended <- bases
    start_index <- pmax(seq_len(genome_length) - left, 1L)
  }

  is_g <- as.integer(extended == "G")
  is_c <- as.integer(extended == "C")
  csum_g <- c(0L, cumsum(is_g))
  csum_c <- c(0L, cumsum(is_c))

  positions <- seq.int(1L, genome_length, by = step)
  if (circular) {
    starts <- positions
    ends <- positions + window - 1L
  } else {
    starts <- start_index[positions]
    ends <- pmin(positions + right, genome_length)
  }

  g_count <- csum_g[ends + 1L] - csum_g[starts]
  c_count <- csum_c[ends + 1L] - csum_c[starts]
  denom <- g_count + c_count
  skew <- ifelse(denom == 0, 0, (g_count - c_count) / denom)
  gc_content <- denom / window

  data.frame(
    position = positions,
    gc_skew = as.numeric(skew),
    gc_content = as.numeric(gc_content),
    window = window,
    step = step
  )
}

read_gc_skew_table <- function(skew_table) {
  x <- read_delim_auto(skew_table)
  position_col <- find_column(x, "position", label = "GC-skew position column")
  skew_col <- find_column(x, c("gc_skew", "gcskew", "skew"), label = "GC-skew value column")
  gc_col <- find_column(
    x,
    c("gc_content", "gccontent", "gc", "content"),
    required = FALSE,
    label = "GC-content value column"
  )
  out <- data.frame(
    position = as_number(x[[position_col]]),
    gc_skew = as_number(x[[skew_col]])
  )
  if (!is.na(gc_col)) {
    out$gc_content <- as_number(x[[gc_col]])
  }
  out
}
