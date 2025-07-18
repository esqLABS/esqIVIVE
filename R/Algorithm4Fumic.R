#' Algorithm4Fumic
#'
#' @description
#' The different algorithms from literature to calculate Fu_mic
#'
#' @param logLipo LogP or LogMA of the compound
#' @param BP Blood plasma ratio, this parameter is needed for Rodgers and Rowland and Poulin method for basic chemicals
#' @param fu In Vivo Fraction Unbound in plasma from literature
#' @param ionization Vector of length 2 with ionization class, acid, neutral and base, if not input then it is c(0,0)
#' @param pKa vector of length of 2 with pkA of the compound
#' @param cMicro concentration of microsomes mg/mL
#'
#' @return  fuInvitro 
#' @export
#'
#' @examples
#' Turner("acid",3,100,3,1)
#' Halifax("base",3,100,3,1)
#' Halifax(c("base",0),c(3,0),100,3,1)
#' Austin("base",3,100,3,1)
#' Poulin("base",pKa=3,X=100,Y=120,BP=1,fu=0.2,cCellNL=0.03,logLipo=3,cMicro_mgml=1)
#' Poulin("acid",pKa=3,X=100,cCellNL=0.03,logLipo=3,cMicro_mgml=1)
#' 
Turner<-function(ionization,pKa,X,logLipo,cMicro_mgml){
  
  if (ionization[1]=="base" & pKa[1]>7){
    
    fuInvitro<-1/(1+cMicro_mgml*10^(0.58*logLipo-2.02))
    
  } else if (ionization[1]=="acid" && pKa[1]<7){ 
    
    fuInvitro<-1/(1+cMicro_mgml*10^(0.2*logLipo-1.54))
    
  } else {
    
    fuInvitro<-1/(1+cMicro_mgml*10^(0.46*logLipo-1.51))
    
  }
  return(fuInvitro)
}

Halifax<-function(ionization,pKa,X,logLipo,cMicro_mgml){
  
  LogD<-1/(1+X)*10^logLipo
  
  if (ionization[1]=="base" & pKa[1]>7){
    
    fuInvitro<-1/(1+cMicro_mgml*10^(0.072*LogD^2+0.067*LogD-1.126))
  } else {
    fuInvitro<-1/(1+cMicro_mgml*10^(0.072*logLipo^2+0.067*logLipo-1.126))
  }
  return(fuInvitro)
}

Austin<-function(ionization,pKa,X,logLipo,cMicro_mgml){
  
  if (ionization[1]=="base" & pKa[1]>7){
    LogP<-1/(1+X)*10^logLipo
  } else {
    LogP<-logLipo
  }
  
  fuInvitro<-1/(1+cMicro_mgml*10^(0.56*logLipo-1.41))
  return(fuInvitro)
}

Poulin<-function(ionization,pKa,X,Y,BP,fu,cCellNL,logLipo,cMicro_mgml){
  kNL<-10**logLipo

  if (ionization[1]=="base" & pKa[1]>7){
    Fnl_E<-0.0024
    Fapl_E<-0.00057
    Fw_E<-0.63
    K_EA<-(BP-(1-0.45))/0.45/fu
    kAPL<-(K_EA- ((1+Y)*Fw_E+kNL*Fnl_E)/(1+X)) * ((1+X)/(Y*Fapl_E))
    
    fuInvitro<-as.double(1/(1+((kNL * cCellNL+X*kAPL*cCellAPL)/(1+X))))
    
  } else { 
    
    fuInvitro<-as.double(1/(1+((kNL * cCellNL)/(1+X))))
  }
  return(fuInvitro)
}