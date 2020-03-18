#' KO-MSP mapping matrix
#'
#' A dataset containing KO mapping to gut MSP
#'
#' @format matrix
#' \describe{
#'   \item{row}{MSPs}
#'   \item{column}{KEGG orthology terms}
#' }
#' @source \url{http://www.diamondse.info/}
"mspKoMat"



#' KEGG orthology description table
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format matrix
#' \describe{
#'   \item{KO}{KO ID}
#'   \item{gene}{Gene ID}
#'   \item{desc}{KO description}
#' }
#' @source \url{http://www.diamondse.info/}
"koDescMap"


#' gut MSP taxonomy table
#'
#' A dataset containing the taxonomy of MSPs
#'
#' @format data frame with taxonomy
#' \describe{
#'   \item{MSP}{MSP id}
#'   \item{size}{the number of genes included}
#'   \item{species}{species-level}
#'   ...
#' }
#' @source \url{http://www.diamondse.info/}
"gutTaxoTab"

#' 7,763 functional clusters
#'
#' A dataset containing the functional clusters
#'
#' @format data frame with taxonomy
#' \describe{
#'   \item{MSP}{MSP id}
#'   \item{size}{the number of genes included}
#'   \item{species}{species-level}
#'   ...
#' }
#' @source \url{http://www.diamondse.info/}
"gutTaxoTab"
