#' correlation_bet_crispr_drug_r
## Compute a correlation between the every gene crispr KO vs each drug response
### call the core.Test revision ( faster than the original version.)
#' @title Compute a correlation between the every gene crispr KO vs each drug response
#' @description to compute the correlations between the viability of cell lines after CRISPR knockout of each gene and the viability of the same cell lines after drug treatment
#' @param infunc_drugName Drug Name
#' @param infunc_drugResponse scores from the drug
#' @param infunc_CRISPRResponse scores from the KO metnohd
#' @return correlation_df  Genes with corelation values and P val associated with drug across the same celllines
#' @examples
#' data (OntargetM)
#' all_drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
#' sample_drugs=sample(all_drugs, 5)
#' metadata.f <- OntargetM$DrugMetadata[sample_drugs,]
#' drug.prism.f <- OntargetM$secondary_prism[sample_drugs,]
#' drugVScrispr.Corr.Features.List=mclapply(sample_drugs,function(x) correlation_bet_crispr_drug_r(x, infunc_drugResponse=drug.prism.f, infunc_CRISPRResponse=OntargetM$avana_CRISPR),mc.cores = 2)
#' names (drugVScrispr.Corr.Features.List) <- sample_drugs
#'saveRDS(drugVScrispr.Corr.Features.List,"drugVScrispr.Corr.Features.List.RDS")
#' @import parallel
#' @export
correlation_bet_crispr_drug_r <- function(infunc_drugName=drugname, infunc_drugResponse= drugResponse,infunc_CRISPRResponse= CRISPRResponse){
  ## make sure that users already use the common cell lines in both datasets >> the result will be precise.
  common.c=intersect(colnames(infunc_drugResponse),
                     colnames(infunc_CRISPRResponse))
  infunc_drugResp=infunc_drugResponse[infunc_drugName,common.c ]
  infunc_CRISPRResp=infunc_CRISPRResponse[,common.c]
  ## calculate corelation.
  correlation_df=sapply(1:nrow(infunc_CRISPRResp),function(x)
    unlist(cor.test_trimmed_v0.default(infunc_drugResp, infunc_CRISPRResp[x,])))
  correlation_df=t(correlation_df)
  row.names(correlation_df) <- row.names(infunc_CRISPRResponse)
  correlation_df
}


