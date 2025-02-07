#To plot the different Partitions



#Function to run plots for each tissue
plot_tissue<-function(tissue,Berez_model,PK_Sim_model,Poulin_model,RR_model,Schmitt_model){

    #Make dataframe partitions
    #merges the input parameters and experimental partitions with the predicted ones
    tissue_df<-as.data.frame(cbind(exp_partition[,c(1:8)],
                                  exp_partition[,tissue],
                                  Berez_model[,tissue],
                                  PK_Sim_model[,tissue],
                                  Poulin_model[,tissue],
                                  RR_model[,tissue],
                                  Schmitt_model[,tissue]))

    colnames(tissue_df)[c(9,10,11,12,13,14)]<-c("Experimental","Berez","PK_Sim","Poulin","RR","Schmitt")

    #transform partition as numeric values
    for (part_nr in 10:14){
          tissue_df[,part_nr]<-as.numeric(tissue_df[,part_nr])}

    #remove empty  columns
    tissue_df<-tissue_df[!is.na(tissue_df$Experimental), ]

    #Calculate R2 for predicitons and experimental dataset
    partition_R2<-round(c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df))$r.squared,
                  "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df))$r.squared,
                  "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df))$r.squared,
                  "RR"= summary(lm(RR~ Experimental, data = tissue_df))$r.squared,
                  "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df))$r.squared),digits=4)

   #Calculate fold differenc ebetween prediction and exp.
     fold_prediction<-data.frame("Berez"=tissue_df$Berez/tissue_df$Experimental,
                       "PK_Sim"=tissue_df$PK_Sim/tissue_df$Experimental,
                       "Poulin"=tissue_df$Poulin/tissue_df$Experimental,
                       "RR"=tissue_df$RR/tissue_df$Experimental,
                       "Schmitt"=tissue_df$Schmitt/tissue_df$Experimental )

    #merge fold predictions to phys-chem
    fold_pred_by_phys_chem<-cbind(fold_prediction,tissue_df[,c(3,5,6,7,8)])


    #calculate percentage over or underpredicted
    cat_fold<-data.frame("Berez"=c(length(which(fold_prediction$Berez>2))/nrow(fold_prediction),
                                   length(which(fold_prediction$Berez<0.5))/nrow(fold_prediction),
                                   length(which(fold_prediction$Berez>0.5 & fold_prediction$Berez<2))/nrow(fold_prediction)),
                         "PK_Sim"=c(length(which(fold_prediction$PK_Sim>2))/nrow(fold_prediction),
                                   length(which(fold_prediction$PK_Sim<0.5))/nrow(fold_prediction),
                                   length(which(fold_prediction$PK_Sim>0.5 & fold_prediction$Berez<2))/nrow(fold_prediction)),
                         "Poulin"=c(length(which(fold_prediction$Poulin>2))/nrow(fold_prediction),
                                   length(which(fold_prediction$Poulin<0.5))/nrow(fold_prediction),
                                   length(which(fold_prediction$Poulin>0.5 & fold_prediction$Berez<2))/nrow(fold_prediction)),
                         "RR"=c(length(which(fold_prediction$RR>2))/nrow(fold_prediction),
                                   length(which(fold_prediction$RR<0.5))/nrow(fold_prediction),
                                   length(which(fold_prediction$RR>0.5 & fold_prediction$Berez<2))/nrow(fold_prediction)),
                         "Schmitt"=c(length(which(fold_prediction$Schmitt>2))/nrow(fold_prediction),
                                   length(which(fold_prediction$Schmitt<0.5))/nrow(fold_prediction),
                                   length(which(fold_prediction$Schmitt>0.5 & fold_prediction$Berez<2))/nrow(fold_prediction)))


      acid_col<-which(tissue_df["A_B_N"]=="acid")
      basic_col<-which(tissue_df["A_B_N"]=="base")
      neutral_col<-which(tissue_df["A_B_N"]=="neutral")
      high_logP<-which(tissue_df["LogP"]>4.5)
      medium_logP<-which(tissue_df["LogP"]<4.5 & tissue_df["LogP"]>1)

      desc_R2<-data.frame(partition_R2)
      if (length(acid_col)>1){
      desc_R2<-cbind(desc_R2,c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                               "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                               "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                               "RR"= summary(lm(RR~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                               "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[acid_col,]))$r.squared))
      } else {desc_R2<-cbind(desc_R2,rep("NA",5)) }

      if (length(basic_col)>1){
        desc_R2<-cbind(desc_R2,c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                                 "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                                 "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                                 "RR"= summary(lm(RR~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                                 "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[basic_col,]))$r.squared))

      } else {desc_R2<-cbind(desc_R2,rep("NA",5)) }


      if (length(neutral_col)>1){
        desc_R2<-cbind(desc_R2,c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                                 "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                                 "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                                 "RR"= summary(lm(RR~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                                 "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[neutral_col,]))$r.squared))

      } else {desc_R2<-cbind(desc_R2,rep("NA",5)) }


      if (length(high_logP)>1){
        desc_R2<-cbind(desc_R2,c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                                 "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                                 "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                                 "RR"= summary(lm(RR~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                                 "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[high_logP,]))$r.squared))

      } else {desc_R2<-cbind(desc_R2,rep("NA",5)) }


      if (length(medium_logP)>1){
        desc_R2<-cbind(desc_R2,c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                                   "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                                   "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                                   "RR"= summary(lm(RR~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                                   "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[medium_logP,]))$r.squared))
      } else {desc_R2<-cbind(desc_R2,rep("NA",5)) }

      colnames(desc_R2)<-c("all","acid","base","neutral","high logP","middle logP")

      #Make plot different
      vect_partition<-colnames(tissue_df)[c(10,11,12,13,14)]
      plot_tissue<-list()

      #for loop to make plots for the different partitions
      for (t in 1:length(vect_partition)){

        column<-vect_partition[t]
        legend_coord<-min(max(tissue_df$Experimental),max(tissue_df[,10:14]))

        plot_tissue[[t]]<-ggplot(data=tissue_df)+
         geom_point(aes(x=Experimental,y=!!sym(column),col=LogP,shape=A_B_N))+
         labs(shape="Type chemical")+
         scale_shape_manual(values=c(3,18,1))+
          ylab("Predicted")+
         ggtitle(column) +
         theme_bw(base_size = 10)+
          geom_text(x=legend_coord*0.8,
                    y=legend_coord*0.8,
                    label=paste("R2=",partition_R2[t]))+
          scale_y_log10(limits = c(0.1, max(tissue_df[,10:14])))+
          scale_x_log10(limits = c(0.1, max(tissue_df$Experimental)))+
          scale_color_gradient(low="blue", high="red")+
         geom_abline(intercept = 0, slope = 1)

      }

      #combined plot
      #I add supress warening because for every point it does not plot it fives warning
      combined_plot<-ggarrange(suppressWarnings(plot_tissue[[1]]),
                           suppressWarnings(plot_tissue[[2]]),
                           suppressWarnings(plot_tissue[[3]]),
                           suppressWarnings(plot_tissue[[4]]),
                           suppressWarnings(plot_tissue[[5]]),
            common.legend=TRUE)

  combined_plot<-annotate_figure(combined_plot, top = text_grob(paste("Partition to",tissue),
                                        face = "bold", size = 20))

  return(list(fold_pred_by_phys_chem,cat_fold,desc_R2,combined_plot))

}


