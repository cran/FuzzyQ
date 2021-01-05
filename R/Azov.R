#' Helminth communities of so-iuy mullets from the Sea of Azov
#'
#' Abundance of 25 helminth species from 378 so-iuy mullets collected in the Sea
#'     of Azov and the Black Sea. Fish are grouped in 12 surveys.
#'
#' @usage data(Azov)
#'
#' @format A data frame with 378 rows and 26 columns. The first column (sample)
#'     is a survey identifier. The remaining columns correspond to species
#'     abundances. See source for species abbreviations and survey identifiers.
#'
#' @references Llopis-Belenguer, C., Blasco-Costa, I., Balbuena, J.A.,
#'     Sarabeev, V., Stouffer, D.B. (2020), Native and invasive hosts play
#'     different roles in host-parasite networks. Ecography, 43: 559-568.
#'     \doi{10.1111/ecog.04963}.
#'
#' @source
#' Llopis-Belenguer, C. (2019) Replication data for: Native and invasive hosts
#'     play different roles in host-parasite networks, Harvard Dataverse,
#'     \doi{10.7910/DVN/IWIKOL}.
#'
#' @examples
#' data(Azov)
#' # Apply the FuzzyQ algorithm to each survey:
#' fuzzyq.azov <- by(Azov[, -1], Azov[, "sample"], fuzzyq, rm.absent = FALSE)
#' # Get cluster membership, silhouette widths and commonness indices
#' # per sp. per survey:
#' sppsilw.azov <- lapply(fuzzyq.azov, function(x) x$spp)
#' # Get global silhouette withds, commonness indices and Dunn's normalized
#' # partition coefficient per survey:
#' global.azov <- t(sapply(fuzzyq.azov, function(x) x$global))
#'
"Azov"
