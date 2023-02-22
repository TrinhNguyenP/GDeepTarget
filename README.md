# GDeepTarget
#V1
## Installation
library(devtools)
install_github("TrinhNguyenP/GDeepTarget")
## Examples
library( GDeepTarget)
data ("OntargetM")
library (parallel)
library (stringr)
all_drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
sample_drugs=sample(all_drugs, 5)
drug.prism.f <- OntargetM$secondary_prism[sample_drugs,]
## Compute a correlation between the every gene crispr KO vs each drug response
drugVScrispr.Corr.Features.List=mclapply(sample_drugs,function(x) CorCRISPRDrugR(x, infunc_drugResponse=drug.prism.f, infunc_CRISPRResponse=OntargetM$avana_CRISPR),mc.cores = 2)
names(drugVScrispr.Corr.Features.List)=sample_drugs
saveRDS(drugVScrispr.Corr.Features.List,
        file = paste('drugVScrispr.Corr.Features.List.RDS', Sys.Date(), '.RDS', sep=''))
## provide a similarity score ( corelation values) for Known Target of a drug and the best Primary Target for each one.
metadata.f <- OntargetM$DrugMetadata[sample_drugs,]
KnownTarget.predictions <- DKSKnownTarget(infunc_drugVScrispr_corr_features_list=drugVScrispr.Corr.Features.List,infunc_drug_metadata = metadata.f)
# Extract primary Target at a pathway Level; Provides a score for each pathway to be the MOA of a drug"
library (fgsea)
gsea_enrichment <- PrimaryTargetatPathwayRes ( infunc_drugVScrispr_corr_features_list=drugVScrispr.Corr.Features.List,infunc_drug_metadata = metadata.f)
## 
saveRDS(gsea_enrichment,
        file = paste('gsea_enrichment.RDS', Sys.Date(), '.RDS', sep=''))
## finding the secondary based on low expression group.
Low_expression_cellLines = sapply(KnownTarget.predictions$MaxTargetName, function(x)
  err_handle(names(which(OntargetM$expression_20Q4[x,]<2))))
identical ( names(Low_expression_cellLines),KnownTarget.predictions$MaxTargetName )

drugVScrispr_corr_features_list_secondary = mclapply(1:nrow(KnownTarget.predictions),
                                                     function(x) err_handle(CorCRISPRDrugR(
                                                       infunc_drugName = unlist(KnownTarget.predictions[x,"drugBroadID"]),
                                                       infunc_drugResponse= OntargetM$secondary_prism[,Low_expression_cellLines[[unlist(KnownTarget.predictions[x,"MaxTargetName"])]]],
                                                       infunc_CRISPRResponse= OntargetM$avana_CRISPR[,Low_expression_cellLines[[unlist(KnownTarget.predictions[x,"MaxTargetName"])]]]
                                                     )) , mc.cores = detectCores())
names (drugVScrispr_corr_features_list_secondary ) <- KnownTarget.predictions[,"drugBroadID"]
saveRDS(drugVScrispr_corr_features_list_secondary,
        paste('drugVScrispr_corr_features_list_secondary', Sys.Date(), '.RDS', sep=''))
        
## record the # of cellLines having the Low_expression.based.max.target.name to the object KnownTarget.predictions to be used later.
cellLines_with_LowTarget = sapply(KnownTarget.predictions[,'MaxTargetName'],
                                  function(x)
                                    err_handle(sum(OntargetM$expression_20Q4[x,] < 2)) )
KnownTarget.predictions$No_cellLines_with_LowExpr_Target <- cellLines_with_LowTarget
## Compute interaction between the drug and KO expression in term of lower vs higher expression
interaction.Features.Secondary.LowExpr <- ExprSpecificityScore (infunc_KnownTarget=KnownTarget.predictions,infunc_Expression=OntargetM$expression_20Q4,
                                                                             infunc_drugResponse=OntargetM$secondary_prism, infunc_CRISPRResponse=OntargetM$avana_CRISPR,infunc_low_ex_cut_off  = 2)
TargetLoweEpr.specificity=data.frame(LowExpr_interaction_strength=sapply(interaction.Features.Secondary.LowExpr, function(x) x[1]),
                                     LowExpr_interaction_P=sapply(interaction.Features.Secondary.LowExpr, function(x) x[2]))
identical ( names(interaction.Features.Secondary.LowExpr), KnownTarget.predictions$drugName)
KnownTarget.predictions <- cbind ( KnownTarget.predictions,TargetLoweEpr.specificity )
## whether interaction is true or false based on cut-off. estimate and p val from lm model
KnownTarget.predictions$Whether_interaction_Ex_based=sapply(interaction.Features.Secondary.LowExpr, function(x) x[1]<0 & x[2]<0.2 )
## Compute interaction between the drug and KO expression in term of mutant
interaction.Features.Mutant<- MutantSpecificityScore (infunc_KnownTarget=KnownTarget.predictions,infunc_mutant=OntargetM$mutations_mat,infunc_drugResponse=OntargetM$secondary_prism, infunc_CRISPRResponse=OntargetM$avana_CRISPR)

Target_Mutation_specificity=data.frame(mutation_interaction_strength=sapply(interaction.Features.Mutant, function(x) x[1]), mutation_interaction_P=sapply(interaction.Features.Mutant, function(x) x[2]))
identical ( names(interaction.Features.Mutant), KnownTarget.predictions$drugName)
KnownTarget.predictions <- cbind ( KnownTarget.predictions,Target_Mutation_specificity )
## mutation interaction with P <0.1
KnownTarget.predictions$predicted_resistance_mutation = KnownTarget.predictions$mutation_interaction_P<0.1
saveRDS(KnownTarget.predictions,
           file = paste('KnownTarget_predictions.RDS', Sys.Date(), '.RDS', sep=''))


