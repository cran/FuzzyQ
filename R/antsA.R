#' Ant species abundance from Arnan et el. (2011)
#'
#' Abundance of 46 ant species from 99 sites sampled in the Nothern Territory
#' (Australia). This dataset corresponds to Plot A data in Arnan et al. (2011).
#'
#' @usage data(antsA)
#'
#' @format A data frame with 99 rows (sites) and 46 variables (species
#'   abundances)
#'
#' @references Arnan, X., Gaucherel, C., Andersen, A. N. (2011) Dominance and
#' species co-occurrence in highly diverse ant communities: a test of the
#' interstitial hypothesis and discovery of a three-tiered competition cascade.
#' Oecologia, 166: 783-794. \doi{10.1007/s00442-011-1919-y}.
#'
#' @source Calatayud, J., Andivia, E., Escudero, A., Melian, C. J.,
#' Bernardo-Madrid, R., Stoffel, M., ... , Madrigal-Gonzalez, J.(2019) Positive
#' associations among rare species and their persistence in ecological
#' assemblages. Nature Ecology & Evolution, 4: 40-45.
#'     \doi{10.6084/m9.figshare.9906092}.
#'
#' @examples
#' data(antsA)
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)
#'
"antsA"
