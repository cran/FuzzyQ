#' Sort Species by fuzzyq Clustering
#'
#' Sort species in a matrix or data frame to match the resulting species order
#' of a \code{fuzzyq} object. This is useful prior to plotting Commonness
#' Indices derived from bootstrap replicates.
#' @param M A matrix or data frame with information of species in columns.
#' @param fq A list of class \code{fuzzyq} returned by \code{FuzzyQ::fuzzyq}.
#' @export
#' @return A matrix or data frame with information of species in columns sorted
#'   according to \code{fq$spp}.
#' @examples
#' data(antsA)
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)
#' # Compute species Commonness Indices of species of 1,000 bootstrap
#' # replicates:
#' \donttest{BS.FQAnts <- fuzzyqBoot (antsA, N = 1e3, level='spp')}
#'
#' # Compute 95 % confidence intervals, percentile method, default values:
#' \donttest{BS.sppCI1 <- fuzzyqCI(BS.FQAnts)}
#'
#' # Plot Commonness Indices and their respective confidence intervals:
#' \donttest{BS.sppCI1 <- sortClus(BS.sppCI1, FQAnts)}
#' spp <- FQAnts$spp
#' col.RC <- c("brown2", "turquoise3") # two colors to plot rare and common
#' # species
#' \donttest{plot(spp[, 3], cex.axis = 0.8, xaxt= 'n', ylab = "Commoness index",
#'    ylim = c(0, max(BS.sppCI1)), xlab = "Species", col = col.RC[spp[, 1] + 1],
#'    pch = 16, cex = 0.8, las = 1)
#' ebar.int <- seq_len(nrow(spp))
#' arrows(ebar.int, BS.sppCI1["Lower", ], ebar.int, BS.sppCI1["Upper", ],
#'    length= 0, col = col.RC[spp[, 1] + 1])
#' axis(1, at = ebar.int, labels = rownames(spp), las = 2, cex.axis = 0.6)}
#'
sortClus <- function(M, fq) {
  if (length(dim(M)) != 2 || !(is.data.frame(M) || is.numeric(M)))
    stop("M is not a dataframe or a numeric matrix.")
  if ("fuzzyq" %in% class(fq) == FALSE) stop("fq is not a fuzzyq
                                                  object.")
  if (fq$is.sorted == FALSE) stop("Common-rare species are not sorted in M.
                                      Run fuzzyq with sorting = TRUE")
  M <- M[, match(rownames(fq$spp), colnames(M))]
  return(M)
}

#' Abundance Occupancy Plot
#'
#' Plots the abundance-occupancy relationship of species in a community
#' categorized as common or rare by fuzzyq.
#' @param fq A list of class \code{fuzzyq} returned by \code{FuzzyQ::fuzzyq}.
#' @param col.rc A vector specifying two colors to be used to plot common and
#'   rare species. Accept any valid color specification in R.
#' @param opacity Number within [0,1] specifying the opacity of convex hulls
#'   grouping common and rare species.
#' @param log.x Logical. Whether or not the x axis should be in log10 scale.
#' @param log.y Logical. Whether or not the y axis should be in log10 scale.
#' @param xLab String. Title for the x axis.
#' @param yLab String. Title for the y axis.
#' @param ... Other graphical parameters to be passed to
#'   \code{\link[graphics]{plot}}.
#' @importFrom grDevices adjustcolor chull colors palette
#' @importFrom graphics polygon plot
#' @importFrom stats pnorm qnorm quantile
#' @export
#' @return A scatter plot of occupancy vs. abundance of species. Convex hulls
#'   enclose common and rare species.
#' @examples
#' data(antsA)
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)
#' AOplot(FQAnts) # Plor with default values
#'
#' # Alternative with colors specified in Hex code, logarithmic axes and other
#' # point format
#' AOplot(FQAnts, col.rc = c("#013bad","#bd5f69"),
#'        log.x = TRUE, log.y = TRUE, pch = 4)
AOplot <- function(fq, col.rc = c("red", "blue"), opacity = 0.1,
                  log.x = FALSE, log.y = FALSE,
                  xLab = "Fraction of sites occupied", yLab = "Mean abundance",
                    ...) {
  if ("fuzzyq" %in% class(fq) == FALSE) stop("fq is not a fuzzyq object.")
  if (fq$is.sorted == FALSE) stop("Common-rare species are not sorted in M.
                              Run fuzzyq with sorting= TRUE")
  if (length(col.rc) != 2) stop("Wrong col.rc format. Specify two colors")
  # check color format
   is.color <- function(x) {
    if (is.numeric(x)) return(x > 0 & (x %% 1 == 0))
    if (any(grepl("^[0-9]+$", x))) {
      x[grepl("^[0-9]+$", x)] <- palette()[(as.numeric(x[grepl("^[0-9]+$",
                                            x)]) - 1) %% 8 + 1]
    }
    y <- grepl("^\\#[a-fA-F0-9]{6}$", x) | grepl("^\\#[a-fA-F0-9]{8}$",
                                                 x) | (x %in% colors())
    return(y)
   }
  if (all(is.color(col.rc)) == FALSE) stop("Wrong col.rc specification.
                                           Use a valid color call.")
      x <- fq$A_O[, 1]
      y <- fq$A_O[, 2]
      if (log.x == TRUE) {
         x <- log10(x)
         xLab <- bquote(paste(log[10]~"("~.(xLab)~")"))
         }
      if (log.y == TRUE) {
         y <- log10(y)
         yLab <- bquote(paste(log[10]~"("~.(yLab)~")"))
         }
      plot(x, y, col = col.rc[fq$spp[, 1] + 1], xlab = xLab, ylab = yLab, ...)
      #draw convex hulls around common,rare spp
      xy <- cbind(x, y)
      Hrar <- chull(xy[which(fq$spp[, 1] == 0), ])
      Hcom <- chull(xy[which(fq$spp[, 1] == 1), ])
      polygon(xy[which(fq$spp[, 1] == 0), ][Hrar, ], border = col.rc[1],
              col = adjustcolor(col.rc[1], alpha.f = opacity))
      polygon(xy[which(fq$spp[, 1] == 1), ][Hcom, ], border = col.rc[2],
              col = adjustcolor(col.rc[2], alpha.f = opacity))
}
