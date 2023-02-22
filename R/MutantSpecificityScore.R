#' MutantSpecificityScore
#' @title Compute interaction between the drug and KO expression in term of mutant
#' @description compute the Mutant Specificity Score (MS Score) of a drug to its known target
#' @param infunc_KnownTarget  The list of targeted genes along with drug information and corelation values "best corelation" and max corelation values along with P values
#' @param infunc_mutant Mutation matrix ( zero and 1)
#' @param  infunc_drugResponse Drug scores.
#' @param infunc_CRISPRResponse scores from KO method
#' @examples
#' data (OntargetM)
#' KnownTarget.predictions <- readRDS ( "KnownTarget_predictions.RDS")
#' interaction.Features.Mutant<- MutantSpecificityScore (infunc_KnownTarget=KnownTarget.predictions,infunc_mutant=OntargetM$mutations_mat,infunc_drugResponse=OntargetM$secondary_prism, infunc_CRISPRResponse=OntargetM$avana_CRISPR)
#' Target_Mutation_specificity=data.frame(mutation_interaction_strength=sapply(interaction.Features.Mutant, function(x) x[1]), mutation_interaction_P=sapply(interaction.Features.Mutant, function(x) x[2]))
#' Target_Mutation_specificity=cbind(KnownTarget_predictions[,1:2], Target_Mutation_specificity)
#' saveRDS(Target_Mutation_specificity,'Target_Mutation_specificity.RDS');
#' @export
MutantSpecificityScore <- function(infunc_KnownTarget=KnownTarget_predictions,infunc_mutant=dat.Mutation.mached,infunc_drugResponse=drug.Res.PRISM.matched, infunc_CRISPRResponse=CRISPR.KO.Res.matched )
  {

  interaction_Features=lapply(1:nrow(infunc_KnownTarget), function(x)
  {
    infunc_drug_matched=err_handle(infunc_drugResponse[infunc_KnownTarget$drugBroadID[x],])
    infunc_avana_matched=err_handle(infunc_CRISPRResponse[infunc_KnownTarget$MaxTargetName[x],])
    infunc_mutation_matched=err_handle(infunc_mutant[infunc_KnownTarget$MaxTargetName[x],])
    err_handle(summary(lm(infunc_drug_matched ~
                            infunc_avana_matched*infunc_mutation_matched))$coefficients[4,c(1,4)])

    }
)

  names(interaction_Features)=infunc_KnownTarget$drugName
  interaction_Features
  ## MOVE THIS PART TO ANALYIS, MAKE IT EASIER TO let user make dieciosn.
  # This is an interaction where the where cell lines with low expression (<2, Trues in above vector) will have "lower" correlation strength than cell lines with high expression (>2, Falses in above vector)
  ## question: is it use the cut-off Pval for 0.2 for true interaction. MAKE IT FLEXIBLE.
  #infunc_KnownTarget$Whether_interaction=sapply(interaction_Features, function(x) x[1]<0 & x[2]<0.2 )
  #infunc_KnownTarget$interaction_strength=sapply(interaction_Features, function(x) x[1])
  #infunc_KnownTarget$interaction_P=sapply(interaction_Features, function(x) x[2])
  #infunc_KnownTarget
}


