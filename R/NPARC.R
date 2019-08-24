#' \code{NPARC} package
#'
#' Non-parametric analysis of response curves
#'
#' See the preprint on
#' \href{https://www.biorxiv.org/content/10.1101/373845v2}{Childs, Bach, Franken et al. (2019): Non-parametric analysis of thermal proteome profiles reveals novel drug-binding proteins}
#'
#' @docType package
#' @name NPARC
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
#' @import BiocParallel
#' @importFrom magrittr %<>%
#' @importFrom stats coefficients fitted resid D vcov integrate predict nobs median mad nls na.omit nls.control runif
#' @importFrom broom augment
#' @importFrom MASS fitdistr
NULL

