
library(readxl)
library(tidyr)
library(ospsuite)
library(dplyr)
library(ggplot2)
library(ggpubr)

#read excel partitions
exp_partition_raw<-read_xlsx("partition check/10928_2017_9548_MOESM1_ESM.xlsx",sheet="Data")


#filter rat and transform columns
exp_partition<-exp_partition_raw %>% 
                filter(Species == "Rat" ) %>%
                filter(A_B_N != "Z" ) %>%
                pivot_wider(names_from = Tissue,  values_from = Exp_PC ,values_fn = list(Tissue = ~!is.na(.))) %>%
                mutate(
                  Type_ionization = case_when(
                  A_B_N == "A" ~ -1,
                  A_B_N== "B" ~ 1,
                  A_B_N == "N" ~ 0 )) %>%
                  mutate( A_B_N = case_when(
                      A_B_N == "A" ~ "acid",
                      A_B_N== "B" ~ "base",
                      A_B_N == "N" ~ "neutral" 
                      ))

exp_partition<-exp_partition %>% group_by(Drug) %>%
                summarise(across(c(CAS,A_B_N,Type_ionization,Effect_pKa), ~first(na.omit(.x))), 
                across(c(LogP,logMA,fu,
                                 Brain,Gut,Kidney,Liver,Muscle,Pancreas 
                                 ,Spleen,Adipose,Blood_Cells ,Bone,Heart,Lung      
                                  ,Skin,Testis,Thymus,Stomach)   
                                 , ~mean(.x, na.rm = TRUE))) %>%
                mutate(logMA=ifelse(is.na(logMA),1.294+0.304*LogP,logMA))

#load rat simulations-----------------------------------------------------------

  Rat_model<-function(partitionQSPR,lipophilicity){
    if (partitionQSPR=="Rodger_Rowland") {
      sim1 <- loadSimulation("partition check/Rat-Rodgers and Rowland.pkml", loadFromCache = FALSE)
      
    } else if  (partitionQSPR=="Schmitt") {
      sim1 <- loadSimulation("partition check/Rat-Schmitt.pkml", loadFromCache = FALSE)
      
    } else if  (partitionQSPR=="PKSim") {
      sim1 <- loadSimulation("partition check/Rat-PK-Sim.pkml", loadFromCache = FALSE)
      
    } else if  (partitionQSPR=="Poulin") {
      sim1 <- loadSimulation("partition check/Rat-Poulin.pkml", loadFromCache = FALSE)
      
    } else if  (partitionQSPR=="Berez") {
      sim1 <- loadSimulation("partition check/Rat-Berez.pkml", loadFromCache = FALSE)
      
    }
    
    outputs<-c(
    "Dose"=getParameter("Applications|Daily ingestion|Dissolved formulation|Application_1|ProtocolSchemaItem|DosePerBodyWeight",sim1),
    "brainKp"=getParameter("Neighborhoods|Brain_int_Brain_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),  
    "adiposeKp"=getParameter("Neighborhoods|Fat_int_Fat_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),  
    "liverKp"=getParameter("Neighborhoods|Periportal_int_Periportal_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1), 
    "kidneyKp"=getParameter("Neighborhoods|Kidney_int_Kidney_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "heartKp"=getParameter("Neighborhoods|Heart_int_Heart_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1), 
    "gutKp"=getParameter("Neighborhoods|Lumen_uje_UpperJejunum_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "muscleKp"=getParameter("Neighborhoods|Muscle_int_Muscle_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "pancreasKp"=getParameter("Neighborhoods|Pancreas_int_Pancreas_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "spleenKp"=getParameter("Neighborhoods|Spleen_int_Spleen_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "boneKp"=getParameter("Neighborhoods|Bone_int_Bone_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),  
    "lungKp"=getParameter("Neighborhoods|Lung_int_Lung_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "skinKp"=getParameter("Neighborhoods|Skin_int_Skin_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),
    "testisKp"=getParameter("Neighborhoods|Gonads_int_Gonads_cell|Test_Chemical|Partition coefficient (intracellular/plasma)",sim1),  
    "pInt"=getParameter("Test_Chemical|Specific intestinal permeability (transcellular)", sim1),
    "Permeability"=getParameter("Test_Chemical|Permeability",sim1), 
    "Fu"=getParameter("Test_Chemical|Fraction unbound (plasma)",sim1),
    "massDrug"=getParameter("Test_Chemical|Total drug mass",sim1))
    
    addOutputs(outputs,simulation = sim1)
    
    #inputs
    parameterPaths <- c("Test_Chemical|Fraction unbound (plasma, reference value)",
                        "Test_Chemical|Lipophilicity")#,
                       # "Test_Chemical|pKa value 0",
                        #"Test_Chemical|Compound type 0")
    
    simBatch <- createSimulationBatch(simulation = sim1, parametersOrPaths = parameterPaths)
    
    
    #exp_partition[is.na(exp_partition)]<-0
    
    ##Check if imported correctly
    
    
    #for having different options for lipophiliticy
    if (lipophilicity=="LogP"){
      lipo_values<-as.double(exp_partition[i,"LogP"])
      
    } else if (lipophilicity=="LogMA"){
      
      lipo_values<-as.double(exp_partition[i,"logMA"])
      
    }
    #The number of parameters to vary for each batch
    #needs to correspond to the vector of parameterPaths and in the same order
    for (i in 1:nrow(exp_partition)){
      
    
      parameterValues = c(as.double(exp_partition[i,"fu"]),
                          as.double(exp_partition[i,"LogP"]))
                        #  as.double(exp_partition[i,"Effect_pKa"]),
                         # as.double(exp_partition[i,"A_B_N"]))
      
      
      simBatch$addRunValues(parameterValues = parameterValues)
      
    }
    
    #Simulations------------------------------------------------------------------
    results <- runSimulationBatches(simBatch)
    
    
    pred_partitions<-data.frame(matrix(ncol = 14, nrow = nrow(exp_partition)))
    colnames(pred_partitions)<-c("Drug","Brain","Adipose","Liver","Kidney",
                                  "Heart","Gut","Muscle","Pancreas",
                                  "Spleen","Bone","Lung","Skin","Testis")

    for (j in 1:nrow(exp_partition)){
      #get table results
      outputValues1<-getOutputValues(results[[1]][[j]])
      
      pred_partitions[j,]<-c(exp_partition$Drug[j],
        "Brain"=as.double(outputValues1$data$`Neighborhoods|Brain_int_Brain_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Adipose"=as.double(outputValues1$data$`Neighborhoods|Fat_int_Fat_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Liver"=as.double(outputValues1$data$`Neighborhoods|Periportal_int_Periportal_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Kidney"=as.double(outputValues1$data$`Neighborhoods|Kidney_int_Kidney_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Heart"=as.double(outputValues1$data$`Neighborhoods|Heart_int_Heart_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Gut"=as.double(outputValues1$data$`Neighborhoods|Lumen_uje_UpperJejunum_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Muscle"=as.double(outputValues1$data$`Neighborhoods|Muscle_int_Muscle_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Pancreasp"=as.double(outputValues1$data$`Neighborhoods|Pancreas_int_Pancreas_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Spleen"=as.double(outputValues1$data$`Neighborhoods|Spleen_int_Spleen_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Bone"=as.double(outputValues1$data$`Neighborhoods|Bone_int_Bone_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Sung"=as.double(outputValues1$data$`Neighborhoods|Lung_int_Lung_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Skin"=as.double(outputValues1$data$`Neighborhoods|Skin_int_Skin_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Testis"=as.double(outputValues1$data$`Neighborhoods|Gonads_int_Gonads_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]))
       
    }
    return("pred_partitions"=pred_partitions)
  }

##MAKE IONIZATION WORK###
