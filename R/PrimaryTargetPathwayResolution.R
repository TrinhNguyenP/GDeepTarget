#' PrimaryTargetatPathwayRes
#' @title  Extract primary Target at a pathway Level; Provides a score for each pathway to be the MOA of a drug"
#' @description  Predicts a Primary Target at a pathwya Level. It next finds the pathways that are most enriched in the genes with high DKS scores. It does this by performing a pathway enrichment test on the ranked gene list by DKS score. The output is a matrix of pathway-level probabilities for each drug to be the primary MOA
#' @param infunc_drugVScrispr_corr_features_list The list of result from "CorCRISPRDrugR" function.
#' @param infunc_drug_metadata meta data from drug
#' @return Pathway level probability to be a primary MOA; gsea_enrichment_all_drugs
#' @examples

#' data (OntargetM)
#' all_drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
#' sample_drugs=sample(all_drugs, 5)
#' metadata.f <- OntargetM$DrugMetadata[sample_drugs,]
#' drug.prism.f <- OntargetM$secondary_prism[sample_drugs,]
#' drugVScrispr.Corr.Features.List=mclapply(sample_drugs,function(x) CorCRISPRDrugR(x, infunc_drugResponse=drug.prism.f, infunc_CRISPRResponse=OntargetM$avana_CRISPR),mc.cores = 2)
#' names (drugVScrispr.Corr.Features.List) <- sample_drugs
#' gsea_enrichment <- PrimaryTargetatPathwayRes ( infunc_drugVScrispr_corr_features_list=drugVScrispr.Corr.Features.List,infunc_drug_metadata = metadata.f)
#' @export
PrimaryTargetatPathwayRes <- function(infunc_drugVScrispr_corr_features_list=drugVScrispr_corr_features_list,infunc_drug_metadata = metadata.Drug){
  target_by_moa=split(infunc_drug_metadata[,'target'], infunc_drug_metadata[,'moa'])
  target_by_moa_unlisted=sapply(target_by_moa, unlist)
  moa_pathways=lapply(target_by_moa_unlisted, function(y){
    unique(unlist(sapply(as.character(y), function(x)
      unlist(strsplit(x, ', ')) )))
  } )
  moa_pathways=sapply(moa_pathways, na.omit)
  drugVScrispr_corr_Strength=sapply(infunc_drugVScrispr_corr_features_list, function(x) x[,2])
  gsea_enrichment_all_drugs=mclapply(1:ncol(drugVScrispr_corr_Strength), function(x) {
    IF_corr_with_gene=unlist(drugVScrispr_corr_Strength[,x])
    names(IF_corr_with_gene)  = rownames(drugVScrispr_corr_Strength)
    IF_corr_with_gene_ordered=sort(IF_corr_with_gene, decreasing = T)
    gsea_enrichment=fgsea(pathways = moa_pathways,
                          stats = IF_corr_with_gene_ordered,
                          minSize=1,
                          maxSize=100)
    gsea_enrichment
  }, mc.cores = detectCores())

  names(gsea_enrichment_all_drugs)=colnames(drugVScrispr_corr_Strength)
  gsea_enrichment_all_drugs
}
