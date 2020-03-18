#' KO-MSP mapping matrix
#'
#' A dataset containing KO mapping to gut MSP
#'
#' @format matrix
#' \describe{
#'   \item{row}{MSPs}
#'   \item{column}{KEGG orthology terms}
#' }
"mspKoMat"



#' KEGG orthology description table
#'
#' A dataset containing KEGG orthology term descriptions
#'
#' @format data frame
#' \describe{
#'   \item{KO}{KO ID}
#'   \item{gene}{Gene ID}
#'   \item{desc}{KO description}
#' }
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
"gutTaxoTab"

#' 7,763 functional clusters
#'
#' A dataset containing the functional clusters
#'
#' @format list of functional clusters
#' \describe{
#'   \item{each item}{list of function/phenotype terms of given cluster}
#'   ...
#' }
"funcModules"

#' longitudinal Swedish cohort metadata
#'
#' A dataset containing sample metadata of longitudinal Swedish cohort
#'
#' @format data frame
#' \describe{
#'   \item{sampleId}{each id corresponds to the column of mgsMat}
#'   \item{metadataId}{each subject id}
#'   \item{time}{time passed after first visit (days)}
#' }
"metaTab"

#' MSP abundance matrix of longitudinal Swedish cohort
#'
#' A dataset containing MSP abundance profile
#'
#' @format matrix
#' \describe{
#'   \item{row}{MSPs}
#'   \item{column}{KEGG orthology terms}
#' }
"mgsMat"

