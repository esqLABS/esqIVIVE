#To plot the different Partitions

source("partition check/Partition functions.R")


# Make the simulations for the different partition methods



#Function to run plots for each tissue
plot_tissue<-function(tissue){

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

#remove empty brain columns
tissue_df<-tissue_df[!is.na(tissue_df$Experimental), ]

#Calculate R2 for predicitons and experimental dataset
partition_R2<-round(c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df))$r.squared,
            "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df))$r.squared,
            "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df))$r.squared,
            "RR"= summary(lm(RR~ Experimental, data = tissue_df))$r.squared,
            "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df))$r.squared),digits=4)

acid_col<-which(tissue_df["A_B_N"]=="acid")
basic_col<-which(tissue_df["A_B_N"]=="base")
neutral_col<-which(tissue_df["A_B_N"]=="neutral")
high_logP<-which(tissue_df["LogP"]>4.5)
medium_logP<-which(tissue_df["LogP"]<4.5 & tissue_df["LogP"]>1)

desc_R2<-data.frame(partition_R2,
                    c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                    "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                    "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                    "RR"= summary(lm(RR~ Experimental, data = tissue_df[acid_col,]))$r.squared,
                    "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[acid_col,]))$r.squared),
                    c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                      "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                      "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                      "RR"= summary(lm(RR~ Experimental, data = tissue_df[basic_col,]))$r.squared,
                      "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[basic_col,]))$r.squared),
                    c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                      "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                      "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                      "RR"= summary(lm(RR~ Experimental, data = tissue_df[neutral_col,]))$r.squared,
                      "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[neutral_col,]))$r.squared),
                    c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                      "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                      "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                      "RR"= summary(lm(RR~ Experimental, data = tissue_df[high_logP,]))$r.squared,
                      "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[high_logP,]))$r.squared),
                    c("Berez"=summary(lm(Berez ~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                      "PK_Sim"=summary(lm(PK_Sim ~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                      "Poulin"=summary(lm(Poulin ~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                      "RR"= summary(lm(RR~ Experimental, data = tissue_df[medium_logP,]))$r.squared,
                      "Schmitt"=summary(lm(Schmitt~ Experimental, data = tissue_df[medium_logP,]))$r.squared))

colnames(desc_R2)<-c("all","acid","base","neutral","high logP","middle logP")
  
  
#Make plot different 
vect_partition<-colnames(tissue_df)[c(10,11,12,13,14)]
plot_tissue<-list()

#for loop to make plots for the different partitions
for (t in 1:length(vect_partition)){
  
  column<-vect_partition[t]
  
  plot_tissue[[t]]<-ggplot(data=tissue_df)+
   geom_point(aes(x=Experimental,y=!!sym(column),col=LogP,shape=A_B_N))+
   labs(shape="Type chemical")+
   scale_shape_manual(values=c(3,18,1))+
    ylab("Predicted")+
   ggtitle(column) +
  geom_text(x=40,y=40,label=paste("R2=",partition_R2[t]))+
   ylim(0,60)+
   xlim(0,60)+
   theme_bw(base_size = 15)+
   geom_abline(intercept = 0, slope = 1)
  
}
  combined_plot<-ggarrange(plot_tissue[[1]],plot_tissue[[2]],
            plot_tissue[[3]],plot_tissue[[4]],
            plot_tissue[[5]], common.legend=TRUE)
  
  annotate_figure(combined_plot, top = text_grob(paste("Partition to ",tissue), 
                                        face = "bold", size = 20))
  
  return(desc_R2)
}


