

plot.SBS.signature <- function(signature,
                               CI = NULL,
                               palette =c("#084B9D", "#2473D0", "#A4C1ED", "#F2DAAC", "#F29D52", "#E82C36")) {
  #load("~/CompressiveNMF/data/Cosmic_data.rdata")
  names_sig <- names(signature)
  df_plot <- data.frame("Sig" = as.matrix(signature)) %>%
    mutate(Channel = names_sig, 
           Triplet = apply(str_split(names_sig, "", simplify = TRUE), 1, 
                           function(x) paste0(x[c(1,3,7)], collapse = "")),
           Mutation = apply(str_split(names_sig, "", simplify = TRUE), 1, 
                            function(x) paste0(x[c(3,4,5)], collapse = "")))
  if(!is.null(CI)){
    #sign_inferred <- sign_inferred/sum(sign_inferred) # Renormalize for comparison
    df_plot$lowCI <- CI[, 1]
    df_plot$highCI <- CI[, 2]
  }
  
  p <- ggplot(df_plot, aes(x = Triplet, y = Sig, fill = Mutation))+
    geom_bar(stat = "identity", width = 0.65) +
    facet_wrap(~Mutation, scales = "free_x", nrow = 1)+
    theme_minimal()+
    scale_fill_manual(values = palette)+
    theme(aspect.ratio = 1,
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, color = "gray35",
                                     vjust = .5, size = 7, margin = margin(t = -5)), 
          panel.grid.major.x = element_blank())
  
  if(!is.null(CI)){
    p <- p + 
      geom_errorbar(aes(x = Triplet, ymin =lowCI, ymax = highCI), color = "forestgreen", width = 0.5)
  }
  
  return(p)
  
}
