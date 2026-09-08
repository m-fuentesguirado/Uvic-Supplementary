library(ggplot2)

unified_colours<-c(
  "Baseline"="#008B8B",
  "Post-radiotherapy"="#8B475D",
  "BL"="#008B8B",
  "END"="#8B475D",
  "Raw reads"="#8C8C8C",
  "Duplex consensus"="#008B8B",
  "Observed duplex fraction"="#008B8B",
  "Ideal duplex fraction"="#666666",
  "Healthy controls"="#7A7A7A",
  "Patient samples"="#007C83"
)

theme_unified<-function(base_size=10,base_family="sans"){
  theme_classic(base_size=base_size)+
    theme(
      text=element_text(
        family=base_family,
        colour="black"
      ),
      axis.title=element_text(
        size=11,
        colour="black",
        face="plain"
      ),
      axis.text=element_text(
        size=10,
        colour="black"
      ),
      legend.title=element_text(
        size=9,
        colour="black"
      ),
      legend.text=element_text(
        size=8.5,
        colour="black"
      ),
      plot.tag=element_text(
        size=12,
        face="bold",
        colour="black"
      )
    )
}

standardize_plot<-function(p){
  p+theme_unified()
}
