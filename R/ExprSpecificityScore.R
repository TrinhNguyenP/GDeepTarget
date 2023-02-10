#' @title Compute interaction between the drug and KO expression in term of lower vs higher expression
#' @description provide each gene a likelihood score of being a secondary target of drug
#' @param infunc_KnownTarget  The list of targeted genes along with drug informaiton and corelation values "best corelation" and max corelation values along with P values
#' @param infunc_Expression Expression matrix
#' @param  infunc_drugResponse Drug scores.
#' @param infunc_CRISPRResponse scores from KO method
#' @param infunc_low_ex_cut_off  desired cut-off for low expression
#' @examples

#' data (OntargetM)
#' KnownTarget.predictions <- readRDS ( "KnownTarget_predictions.RDS")
#' Low_expression_cellLines = sapply(KnownTarget.predictions$MaxTargetName, function(x) err_handle(names(which(dat.Expr.matched[x,]<2))))
#' interaction.Features.Secondary.LowExpr <- ExprSpecificityScore (infunc_KnownTarget=KnownTarget.predictions,infunc_Expression=OntargetM$expression_20Q4,infunc_drugResponse=OntargetM$secondary_prism, infunc_CRISPRResponse=OntargetM$avana_CRISPR,infunc_low_ex_cut_off  = 2)
#' TargetLoweEpr.specificity=data.frame(LowExpr_interaction_strength=sapply(interaction.Features.Secondary.LowExpr, function(x) x[1]),
#'LowExpr_interaction_P=sapply(interaction.Features.Secondary.LowExpr, function(x) x[2]))
#' identical ( names(interaction.Features.Secondary.LowExpr), KnownTarget.predictions$drugName)
#' TargetLoweEpr.specificity=cbind(KnownTarget.predictions[,1:2], TargetLoweEpr.specificity)
#' saveRDS(TargetLoweEpr.specificity, file='Target_LowedEpr_specificity.RDS')

ExprSpecificityScore <- function(infunc_KnownTarget=KnownTarget_predictions,infunc_Expression=expression_matched,infunc_drugResponse=drug.Res.PRISM.matched, infunc_CRISPRResponse=CRISPR.KO.Res.matched,infunc_low_ex_cut_off  = 3 )
  {
  ## record the # of cellLines having the Low_expression.based.max.target.name to the object KnownTarget.predictions to be used later.
  ## extra code it doesn't use.
  # Targets_wd_lowExp_atLeast_FiveCellLines = infunc_KnownTarget$cellLines_withLOWexp > infunc_no_celllines_low_ex_cut_off
  # <!-- Compute interaction -->
  ## need to change the name.
  interaction_Features=lapply(1:nrow(infunc_KnownTarget), function(x)
  {
    infunc_drug_matched=err_handle(infunc_drugResponse[infunc_KnownTarget$drugBroadID[x],])
    infunc_avana_matched=err_handle(infunc_CRISPRResponse[infunc_KnownTarget$MaxTargetName[x],])
    infunc_expression_matched=err_handle(infunc_Expression[infunc_KnownTarget$MaxTargetName[x],]< infunc_low_ex_cut_off )

    err_handle(summary(lm(infunc_drug_matched ~
                            infunc_avana_matched*infunc_expression_matched))$coefficients[4,c(1,4)])

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


