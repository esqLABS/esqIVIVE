calcPartitions<-function(LogP,hlcAt,MW,pKa,typeIonization,assumptions,fu,BC){

  #Calculate ionization
  if(typeIonization=="acid"){
  fNeutral=1/(1+10**(pH-pKa))
  }else if (typeIonization=="base"){
  fNeutral=1/(1+10**(pKa-pH))
  }else{
  fNeutral=1
  }

  kOW=10**logP

  #Calculate air-water partition coefficient
  #Divide henry law constant in atm/(m3*mol) with the temperature in kelvin and gas constant R (j/k*mol)
  kAir=hlcAt/(0.08206*310)

  #Calculate plastic partitioning
  kPlasticFischer=10**(logP*0.47-4.64)
  kPlasticKramer=10**(logP*0.97-6.94)
  kPlastic=mean(kPlasticFischer,kPlasticKramer)*fNeutral

  if ("assumption"=="Poulin and Theil"){
    #Calculate protein partitioning

    #Calculate lipid partitioning
    kNLip= fNeutral*kOW

  fuInvitro=1/(1+kNLip*cLip+kPlastic*saPlasticVolMedium)

  }else if ("assumption"=="PK-Sim® Standard"){

  }else if ("assumption"=="Rodgers & Rowland"){

      Hema=0.45
      kpuBC=(Hema-1+BP)/(Hema*fu)
      kAPLip_1=kpuBC-
        ((1+10**(pKa-pHBC))/(1+10**(pKa-pH))*fiwBC)-
        (kOW*fnlBC+(0.3*kOW+0.7)*fnpBC)/(1+10**(pKa-pH))

      kAPBC=kAPLip_1*(1+10**(pKa-pH))/(0.5*10**(pKa-pHBC))
      kAPLip= (kAPBC*10**(pKa-pH))*fNeutral
      kNLip= fNeutral*kOW
      if(typeIonization=="base"){
          fuInvitro=1/(1+kNLip*cNLip+
                         (kNLip*0.3+0.7)*cNPLip+
                         kAPLip*cAPLip+ #need to check units of APL ( rodgers are in mg/g)
                         kPla*saPlasticVolMedium)
      }else{
          fuInvitro=1/(1+kNLip*cNLip+
                         (kNLip*0.3+0.7)*cNPLip+
                         kPla*saPlasticVolMedium)
      }

  }else if ("assumption"=="Schmitt"){

  }else if ("assumption"=="Berezhkovskiy"){

  }else {}


  #Warning for volatility
  fuAir=fuInvitro*kAir*volAir_L
  if (fuAir>0.1){
    print("Probable evaporation")
  }else {}

  return(fuInvitro)
}
