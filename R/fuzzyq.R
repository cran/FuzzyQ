#' Fuzzy Quantification of Common and Rare Species in Ecological Communities
#'
#' Perform fuzzy clustering of each species based on their abundance and
#' occupancy.
#' @param M A matrix or data frame of species abundances (columns). Each row
#'   represents a site.
#' @param diss String. Specify the dissimilarity coefficient to be used. Default
#'   is "gower". See \code{\link[cluster]{daisy}} in package \code{cluster} for
#'   other choices.
#' @param rm.absent Logical. Whether or not absent species are to be removed
#'   from the calculations.
#' @param sorting Logical. If \code{TRUE} (the default) species are sorted in
#'   the output by ascending silhouette widths within each cluster, else species
#'   are arranged in the same order as in the input matrix or data frame.
#' @param keep.Diss Logical. Whether or not the species dissimilarity matrix
#'   shoudl be returned. The default is \code{FALSE}.
#' @param std Logical. Whether or not the measurements of occupancy and
#'   abundance are to be standardized before calculating the dissimilarities.
#'   Measurements are standardized for each variable (column), by subtracting
#'   the variable's mean value and dividing by the variable's mean absolute
#'   deviation. It only takes effect if \code{diss} is different from "gower".
#' @param wgts an optional numeric vector of length 2. To be used if diss =
#'   "gower", specifying weights for occupancy and abundance, respectively.
#'   Default is 1 each as in Gower's original formula.
#' @param ... Arguments to be passed to function \code{fanny} in package
#'   \code{cluster}.
#' @export
#' @return A list of class \code{fuzzyq} containing the following: \describe{
#'   \item{\code{A_O}}{Abundance-occupancy information for each species.}
#'   \item{\code{Diss}}{Object of class dist with pairwise dissimilarities among
#'   species based on A_O. (only if \code{keep.Diss = TRUE)}.}
#'   \item{\code{spp}}{Clustering metrics per species: Cluster membership (where
#'   0 and 1 denote allocation to the rare and common category, respectively),
#'   Silhouette Widths and Commonness Indices).} \item{\code{global}}{Community
#'   level clustering metrics: Average silhouette widths per cluster and
#'   globally, Mean commonness indices per cluster and Normalized Dunn's
#'   coefficient.} }
#' @seealso \code{\link[cluster]{fanny}} and \code{\link[cluster]{daisy}} in
#' package \code{cluster}
#' @examples
#' data(antsA)
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)

fuzzyq <- function(M, diss = "gower", rm.absent = FALSE, sorting = TRUE,
                   keep.Diss = FALSE, std = FALSE, wgts = c(1,1), ...) {
  if (length(dim(M)) != 2 || !(is.data.frame(M) || is.numeric(M)))
    stop("M is not a dataframe or a numeric matrix.")
  if (rm.absent == TRUE &&
      length(which(colSums(M) == 0)) != 0)  M <- M[, -which(colSums(M) == 0)]
  if (length(dim(M)) != 2) stop("Insufficent data: only one species")
  abund <- colMeans(M, na.rm = TRUE)
  names(abund) <- colnames(M)
  M[M > 0] <- 1
  occ  <- colMeans(M, na.rm = TRUE)
  A_O <- cbind(occ, abund)
  colnames(A_O) <- c("frq.occ", "m.abund")
  D <- cluster::daisy(A_O, metric = diss, stand = std, weights = wgts)
  # check that there are at sufficient no. of spp for fanny (k >= n/2-1)
  n <- attr(D, "Size")
  if (n < 6) {
    warning("Insufficient number of spp. for fuzzy clustering. NULLs produced")
    sil.w <- NULL
    global <- NULL
    cr.fan <- list(A_O = A_O, spp = sil.w, global = global)
    if (keep.Diss == TRUE) cr.fan <- append(cr.fan, list(Diss = D), 1)
  } else {
    fanclus <- cluster::fanny(D, k = 2, keep.diss = FALSE,
                              keep.data = FALSE, ...)
    sil.w <- fanclus$silinfo$widths[, c(1, 3)]
    # This part ensures that common and rare are consistently related to
    # clusters. (Rarest sp. is assumed to belong to rare category)
    rar.tag <- names(which(A_O[, 1] == min(A_O[, 1])))[1]
    sil.tag <- which(row.names(sil.w) == rar.tag)
    # 0 = rare; 1 = common
    sil.w[which(sil.w[, 1] == sil.w[sil.tag, 1]), 1] <- 0
    sil.w[which(sil.w[, 1] != sil.w[sil.tag, 1]), 1] <- 1
    # Ensure that "membership" reflects "commonness"
    if (fanclus$membership[which(rownames(fanclus$membership) == rar.tag), 1] <=
        fanclus$membership[which(rownames(fanclus$membership) == rar.tag), 2])
      m <- 1 else m <- 2
    # sort prior to cbind by membership
    sil.w <- sil.w[row.names(fanclus$membership), ]
    sil.w <- cbind(sil.w, fanclus$membership[, m])
    colnames(sil.w)[3] <- "Common.I"
    # sorting options
    if (sorting == TRUE) {
      sil.w <- sil.w[order(sil.w[, 1], sil.w[, 2]), ]
      #order A_O rows as sil.w:
      A_O <- A_O[match(rownames(sil.w), rownames(A_O)), ]
    }
    # compute average sil widths and commoness coeffs
    sil.w <- as.data.frame(sil.w)
    silw <- c(mean(sil.w[sil.w$cluster == 0, 2]),
              mean(sil.w[sil.w$cluster == 1, 2]),
              mean(sil.w[, 2]))
    memb <- c(mean(sil.w[sil.w$cluster == 0, 3]),
              mean(sil.w[sil.w$cluster == 1, 3]))
    global <- c(silw, memb, fanclus$coeff[2])
    names(global) <- c("silw.rar", "silw.com", "silw.all",
                       "commI.rar", "commI.com", "N.Dunn")
    cr.fan <- list(A_O = A_O, spp = sil.w, global = global)
    if (keep.Diss == TRUE) cr.fan <- append(cr.fan, list(Diss = D), 1)
  }
  class(cr.fan) <- append(class(cr.fan), "fuzzyq")
  is.sorted <- sorting
  cr.fan$is.sorted <- is.sorted
  return(cr.fan)
}

