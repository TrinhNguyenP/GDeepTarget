#' DKSKnownTarget
#' @title provide a similarity score ( corelation values) for Known Target of a drug and the best Primary Target for each one.
#' @description It takes the name of a drug as an input and returns a list containing the DKS score and P-value for the correlation between the drug and each gene.
#' @param infunc_drugVScrispr_corr_features_list The list of result from "CorCRISPRDrugR" function.
#' @param infunc_drug_metadata Meta data from the drug.
#' @return knownTarget_predictions object: The list of targeted genes along with drug information and  "best corelation" and "max corelation values" along with P values
#' @examples
#' data (OntargetM)
#' library (parallel)
#' all_drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
#' sample_drugs=sample(all_drugs, 5)
#' metadata.f <- OntargetM$DrugMetadata[sample_drugs,]
#' drug.prism.f <- OntargetM$secondary_prism[sample_drugs,]
#' drugVScrispr.Corr.Features.List=mclapply(sample_drugs,function(x) CorCRISPRDrugR (x, infunc_drugResponse=drug.prism.f, infunc_CRISPRResponse=OntargetM$avana_CRISPR),mc.cores = 2)
#' names (drugVScrispr.Corr.Features.List) <- sample_drugs
#' KnownTarget.predictions <- DKSKnownTarget( infunc_drugVScrispr_corr_features_list=drugVScrispr.Corr.Features.List,infunc_drug_metadata = metadata.f)
#' saveRDS(KnownTarget.predictions, "KnownTarget_predictions.RDS")
#' @import parallel
#' @export
DKSKnownTarget <- function(infunc_drugVScrispr_corr_features_list=drugVScrispr_corr_features_list,infunc_drug_metadata = metadata.Drug){
  drugSubset_map2_annotation <- match(names(infunc_drugVScrispr_corr_features_list), infunc_drug_metadata[,"broad_id_trimmed"])
  ## based on the name of the drug.
  all_drugNames_common = infunc_drug_metadata[,'name'] [drugSubset_map2_annotation]
  all_drugTargets = infunc_drug_metadata[drugSubset_map2_annotation,c('name','target')]
  all_drugMOA = infunc_drug_metadata[,'moa'][drugSubset_map2_annotation]
  all_drugTargets_seperated = str_split(as.character(all_drugTargets$target), ", ")
  ## this is extract cor values and turn to the matrix where genes are rows and drugs are column.
  drugVScrispr_corr_Strength=sapply(infunc_drugVScrispr_corr_features_list, function(x) x[,2])
  corrMat <- drugVScrispr_corr_Strength
  corrMat_P=sapply(infunc_drugVScrispr_corr_features_list, function(x) x[,1])
  ## map to get the name of drug with corelation values based on the genes.
  all_drugTargets_correlation=sapply(1:length(all_drugTargets_seperated),
                                     function(x) {
                                       ret_cor=drugVScrispr_corr_Strength[
                                         match(all_drugTargets_seperated[[x]],
                                               rownames(drugVScrispr_corr_Strength)), x]
                                       names(ret_cor)=err_handle(all_drugTargets_seperated[[x]])
                                       ret_cor
                                     } )

  names(all_drugTargets_correlation) <- as.character(all_drugNames_common)
  ###
  ### Dataframe of the knowntarget_prediction ( get the annotation from broad.)
  all_drugTargets_MAXcorrelation=sapply(all_drugTargets_correlation, function(x) max(x, na.rm=T))
  ## max corr gene is null. maybe because some drug do multiple genes. the one with multiple genes will have the max corrleation gene.
  all_drugTargets_MAXcorrGene=sapply(all_drugTargets_correlation,
                                     function(x) names(which.max(x)))
  all_drugTargets_MAXcorrGene[sapply(all_drugTargets_MAXcorrGene, length)==0]=NA
  KnownTarget_predictions=data.frame(
    drugName=names(all_drugTargets_MAXcorrelation),
    MaxTargetName=unlist(all_drugTargets_MAXcorrGene),
    Maxcorr=all_drugTargets_MAXcorrelation)
  KnownTarget_predictions=KnownTarget_predictions[order(KnownTarget_predictions$Maxcorr, decreasing = T),]
  KnownTarget_predictions$drugBroadID = infunc_drug_metadata$broad_id_trimmed[match(KnownTarget_predictions$drugName, infunc_drug_metadata$name)]
  ## match the drug name, get the max corelation of that drug, and then calcualtes  the p value.
  BestTargetName=apply(corrMat[,match(KnownTarget_predictions$drugBroadID, colnames(corrMat))],
                       2, function(x) rownames(corrMat)[which.max(x)] )
  BestTargetName[sapply(BestTargetName, length)==0]=NA
  KnownTarget_predictions$BestTargetName=unlist(BestTargetName)
  #best hit score
  BestTargetCorr=apply(corrMat[,match(KnownTarget_predictions$drugBroadID, colnames(corrMat))],
                       2, function(x) max(x, na.rm = T) )
  KnownTarget_predictions$BestTargetCorr=BestTargetCorr

  # Best Hit Significance

  KnownTarget_predictions$BestTargetCorrP = sapply(1:nrow(KnownTarget_predictions), function(x)
    err_handle(corrMat_P[KnownTarget_predictions[x,'BestTargetName'], KnownTarget_predictions[x,'drugBroadID']]) )
  ## ## stop here. error because maxtarget name is NA.
  # Known Target Significance
  KnownTarget_predictions$KnownTargetCorrP = sapply(1:nrow(KnownTarget_predictions), function(x)
    err_handle(corrMat_P[KnownTarget_predictions[x,'MaxTargetName'], KnownTarget_predictions[x,'drugBroadID']]) )
  corrMat_FDR=apply(corrMat_P, 2, function(x) fdrcorr(x))
  # Best Hit Significance - FDR corrected
  KnownTarget_predictions$BestTargetCorrFDR = sapply(1:nrow(KnownTarget_predictions), function(x)
    err_handle(corrMat_FDR[KnownTarget_predictions[x,'BestTargetName'], KnownTarget_predictions[x,'drugBroadID']]) )
  # Known Target Significance - FDR corrected
  KnownTarget_predictions$KnownTargetCorrFDR = sapply(1:nrow(KnownTarget_predictions), function(x)
    err_handle(corrMat_FDR[KnownTarget_predictions[x,'MaxTargetName'], KnownTarget_predictions[x,'drugBroadID']]) )
  KnownTarget_predictions

}
