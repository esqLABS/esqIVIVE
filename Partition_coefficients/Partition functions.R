

#load rat simulations-----------------------------------------------------------

  Rat_model<-function(partitionQSPR,lipophilicity,ionization){

    #mind that the path for these pkml files it relatively to the QUARTO file
    if (partitionQSPR=="Rodger_Rowland") {
      sim1 <- loadSimulation("Rat-Rodgers and Rowland.pkml", loadFromCache = FALSE)

    } else if  (partitionQSPR=="Schmitt") {
      sim1 <- loadSimulation("Rat-Schmitt.pkml", loadFromCache = FALSE)

    } else if  (partitionQSPR=="PKSim") {
      sim1 <- loadSimulation("Rat-PK-Sim.pkml", loadFromCache = FALSE)

    } else if  (partitionQSPR=="Poulin") {
      sim1 <- loadSimulation("Rat-Poulin.pkml", loadFromCache = FALSE)

    } else if  (partitionQSPR=="Berez") {
      sim1 <- loadSimulation("Rat-Berez.pkml", loadFromCache = FALSE)

    }

   ## If I need to explore the paths to parameters
   # tree <- getSimulationTree(sim1)

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
    "bloodcells"=getParameter("Test_Chemical|Partition coefficient (blood cells/plasma)", sim1),
    "Permeability"=getParameter("Test_Chemical|Permeability",sim1),
    "Fu"=getParameter("Test_Chemical|Fraction unbound (plasma)",sim1),
    "massDrug"=getParameter("Test_Chemical|Total drug mass",sim1))

    addOutputs(outputs,simulation = sim1)

    #indicate what are the inputs
    parameterPaths <- c("Test_Chemical|Fraction unbound (plasma, reference value)",
                        "Test_Chemical|Lipophilicity",
                        "Test_Chemical|pKa value 0",
                        "Test_Chemical|Compound type 0")

    simBatch <- createSimulationBatch(simulation = sim1, parametersOrPaths = parameterPaths)

    #exp_partition[is.na(exp_partition)]<-0

    #for having different options for lipophiliticy
    if (lipophilicity=="LogP"){
      lipo_values<-as.vector(exp_partition[,"LogP"][[1]])

    } else if (lipophilicity=="LogMA"){

      lipo_values<-as.vector(exp_partition[,"logMA"][[1]])

    } else {warning("error in choice of lipophilicity")}

     #option to consider ionization or not
    if (ionization=="considered"){

      pKa<-as.double(exp_partition[,"Effect_pKa"][[1]])
      ionization_type<-as.double(exp_partition[,"Type_ionization"][[1]])
      pKa[is.na(pKa)] <- 0

    } else if (ionization=="ignored"){

      pKa<-rep(0,nrow(exp_partition))
      ionization_type<-rep(0,nrow(exp_partition))

    }else {warning("error in choice of ionization")}

    #The number of parameters to vary for each batch
    #needs to correspond to the vector of parameterPaths and in the same order
    for (i in 1:nrow(exp_partition)){


      parameterValues = c(as.double(exp_partition[i,"fu"]),
                          lipo_values[i],
                          pKa[i],
                          ionization_type[i])

      simBatch$addRunValues(parameterValues = parameterValues)

    }

    #Simulations------------------------------------------------------------------
    results <- runSimulationBatches(simBatch)


    pred_partitions<-data.frame(matrix(ncol = 15, nrow = nrow(exp_partition)))
    colnames(pred_partitions)<-c("Drug","Brain","Adipose","Liver","Kidney",
                                  "Heart","Gut","Muscle","Pancreas",
                                  "Spleen","Bone","Lung","Skin","Testis","Blood_Cells")

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
        "Pancreas"=as.double(outputValues1$data$`Neighborhoods|Pancreas_int_Pancreas_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Spleen"=as.double(outputValues1$data$`Neighborhoods|Spleen_int_Spleen_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Bone"=as.double(outputValues1$data$`Neighborhoods|Bone_int_Bone_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Sung"=as.double(outputValues1$data$`Neighborhoods|Lung_int_Lung_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Skin"=as.double(outputValues1$data$`Neighborhoods|Skin_int_Skin_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Testis"=as.double(outputValues1$data$`Neighborhoods|Gonads_int_Gonads_cell|Test_Chemical|Partition coefficient (intracellular/plasma)`[1]),
        "Blood_Cells"=as.double(outputValues1$data$`Test_Chemical|Partition coefficient (blood cells/plasma)`[1]))

    }
    return("pred_partitions"=pred_partitions)
  }

##MAKE IONIZATION WORK###
