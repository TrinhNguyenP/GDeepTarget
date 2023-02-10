#' OntargetM. #' an object containing Viability matrix after CRISPR-KO; Viability after Drug Treatment; Drug metadata from Broad, mutation matrix, and expression matrix with common celllines and common drugs.
#' The list are as follows:
#' @name OntargetM
#' @docType data
#' @format the format is the list of 5 \code {matrix} \itemize{
#' \item\$DrugMetadata:\code{matrix}  containing 1448 drugs with 5 columns:  broad_id_trimmed, name, target, drug_category, and moa
#' \item\$secondary_prism: \code{matrix} containing Viability scores after Drug Treatment with 1448 drugs across 319 celllines
#' \item\$avana_CRISPR: \code{matrix} containing Viability scores after CRISPR-KO with 18333 genes across 319 celllines
#' \item\$mutations_mat:\code{matrix} containing Viability scores after CRISPR-KO with 19350 genes across 319 celllines
#' \item\$expression_20Q4:\code{matrix} containing Viability scores after CRISPR-KO with 19182 genes across 319 celllines. }

#' DrugMetadata
#'
#' This dataset was from Corsello_supplemental_tables.xlsx:
#' @source depmap \url{ https://depmap.org/repurposing/#:~:text=Corsello_supplemental_tables.xlsx}
#' @keywords Corsello_supplemental_tables.xlsx

#' Secondary_prism
#'
#' Drug Treatment profile ~ 1,448 drugs against 489 cell lines in an 8-step:
#' @source depmap \url{https://depmap.org/portal/download/all/?releasename=PRISM+Repurposing+19Q4&filename=secondary-screen-dose-response-curve-parameters.csv}
#' @keywords secondary-screen-dose-response-curve-parameters.csv
#'
#' avana_CRISPR
#'
#' CRISPR KO profile ( scores) of ~800 cell lines for 20K genes
#' @source depmap \url{https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q4&filename=CRISPRGeneEffect.csv}

#' mutations_mat
#' Mutation binary profile for ~2K cell lines; 0 is WT; 1 is mutated
#' @source depmap \url{https://depmap.org/portal/download/all/?releasename=DepMap+Public+22Q4&filename=OmicsSomaticMutations.csv}
#'
#' expression_20Q4
#' Expression profile of ~2K cell lines for 20K genes; file name: OmicsExpressionProteinCodingGenesTPMLogp1.csv
#' @source depmap \url{https://depmap.org/portal/download/all/}
NULL
