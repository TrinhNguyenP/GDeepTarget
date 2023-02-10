## not sure why we can't save to pdf file has errror.
lm_plot_interaction <- function ( infunc_KnownTarget=KnownTarget_predictions,infunc_Expression=expression_matched,infunc_drugResponse=drug.Res.PRISM.matched, infunc_CRISPRResponse=CRISPR.KO.Res.matched,out.dir= "interaction_plot/",  cellLines_withLOWexp=5, infunc_low_ex_cut_off=2){
  infunc_KnownTarget.f <- infunc_KnownTarget[ which (infunc_KnownTarget[,'cellLines_withLOWexp']>cellLines_withLOWexp),]
  for ( i in 1:nrow(infunc_KnownTarget.f)){
    infunc_drug_matched= infunc_drugResponse[infunc_KnownTarget.f$drugBroadID[i],]
    infunc_avana_matched=infunc_CRISPRResponse[infunc_KnownTarget.f$MaxTargetName[i],]
    infunc_expression_matched=infunc_Expression[infunc_KnownTarget.f$MaxTargetName[i],]< infunc_low_ex_cut_off
    identical ( names(infunc_avana_matched), names (infunc_drug_matched))
    identical ( names(infunc_avana_matched), names (infunc_expression_matched))
    out <- lm(infunc_drug_matched ~ infunc_avana_matched*infunc_expression_matched)
    summ(out)
    pdf (paste0(out.dir,"lm_interaction_plot_" ,infunc_KnownTarget.f$drugName[i],".pdf"))
    interact_plot(out, pred = infunc_avana_matched, modx = infunc_expression_matched, interval = TRUE,
                  int.width = 0.8,x.label="CRISPRResponse_scores", y.label = paste("drug Response scores of ", infunc_KnownTarget.f$drugName[i]) )
    dev.off()
  }
}
