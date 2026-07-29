\name{label_adjust}
\alias{label_adjust}
\title{Manually shift selected circular labels}
\usage{
label_adjust(label, hjust = 0, vjust = 0)
}
\arguments{
\item{label}{Label text to match. The match is checked against the plotted
label, raw label, display label, and wrapped label text.}

\item{hjust}{Horizontal shift in plot-panel units. Positive values move the
label to the right.}

\item{vjust}{Vertical shift in plot-panel units. Positive values move the
label upward.}
}
\value{
A data frame usable with the \code{label_adjust} argument in
\code{\link[=ggplasmid]{ggplasmid()}} or \code{\link[=plot_phage_map]{plot_phage_map()}}.
}
\description{
Manually shift selected circular labels
}
