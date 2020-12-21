VL_matrix = matrix(nrow = cohortSize, ncol = length(seq(0, maxTime, 1)))
bl_number = 13 # number of baseline characteristics
statesNum = 52 # number of states
#proportion of patients in each fibrosis stages (F0, F1, F2, F3, F4) at simulation start (note: simulation starts at time of HCV infection)
start <- sample(c(1,2,3,4,5), cohortSize, c(0.8590, 0.1509, 0.0184, 0.00239, 0.00015), replace=T)

#data = read.table("baseline_data.txt")


p_HIVMSM = 0.5
p_IDU_start = 0.1
p_IDU_stop = 1 # probability of stoping IDU
p_HIV = 0.01 # probability of getting HIV after HCV
################################Functions#########################################
bl = matrix (nrow = cohortSize, ncol = bl_number)
colnames(bl) <- c("Alcohol","Gender","HIV","HIV_time","IDU","IDU_time_start","IDU_time_stop", "Genotype" , "Age", "BirthYear", "sr.time","Viral_diag", "Diag.time")

V0 = sample(real_vl_density$x,cohortSize,real_vl_density$y, replace = T)
############################################################################################
V0 = sample(real_vl_density$x,cohortSize,real_vl_density$y, replace = T)
epsi = runif(ncol(VL_matrix),0,1)
VL = V0 + runif(length(V0), 0, epsi) * rbinom(length(V0), 1, 0.5)
V0 = VL
############################################
#'Baseline spontaneous recovery time
###################################
#' sr.fun()
sr.fun = function()
{
  logisticfunc <- function(min, max, infl, slo) {min + ((max - min)/(1 + (teta / infl)^ slo))}
  
  sr.t <- rep(0, cohortSize)
  for (i in 1:cohortSize)
  {
    teta = runif( 1, 0, 1)
    
    p = logisticfunc( min = 0, max = 1, infl = 0.25, slo = 2.23)
    
    w <- runif(1, 0, 1)
    ifelse( w <= p, sr.t[i] <- teta, sr.t[i] <- 999)
  }
  sr.t
  
}#End Of Spontaneous Recovery

#####################################################################################
#BACKGROUND MORTALITY
# will be used in the hazard_fun file to calculate the harzard function for mortality.
#To match age category with death rate (cuts is the value of the start of the intervall eg for eg 20-24 => cuts=20)
#########################################################
#' @param Cuts the \code{transition.structure}
#' @param x the \code{transition.structure}
#' @param func.values
#' @param tol the \code{transition.structure}
#'pw.eval.ext()
pw.eval.ext = function(Cuts, x, func.values, tol = 0.0001)
{
  
  lengthCuts = length(Cuts)
  unlist(lapply(x,
                function(xp) {
                  
                  if (xp <  max(Cuts))
                    func.values[ which((Cuts[c(1:(lengthCuts - 1))] - rep(xp, lengthCuts - 1))*
                                         (Cuts[c(2:lengthCuts)] - tol  - rep(xp, lengthCuts-1)) <= 0 , arr.ind =  TRUE)[1] ]
                  
                  else func.values[lengthCuts]
                  
                }
                
                
  ))
  
}

###################################################################################
################################################################################################
###################################################################
bl[,"Alcohol"] = sample(0:2, cohortSize, c(37, 37 ,26), replace=T) # Bl alcohol consumption, 26% heavy consumer, and we distribute the others equally
bl[,"Genotype"] = sample(0:3, cohortSize, c(18.50, 3.70, 25.90, 51.90), replace=T) # Genotype distibution: bl[,2]=(0,1,2,3)=>gen(1,mixed,2,3)
Diagnosis_age = sample(1:3, cohortSize, c(3.54, 78,18.46), replace = T)
Diagnosis_age = sample(1:3, cohortSize, c(11, 78,11), replace = T)
bl[,"Age"] <- Diagnosis_age
bl[,"BirthYear"] <- Diagnosis_age 


for(i in 1:cohortSize)
{
  
  
  age_ID = Diagnosis_age[i]
  age = switch (age_ID,
                runif(1, 14, 30) ,
                runif(1, 30, 59) ,
                runif(1, 59, 73)
  )

  bl[i,"BirthYear"] = 2017 -  age
  
  
  inf_year =  bl[i,"BirthYear"]  + runif(1,10,(age + 4))
  bl[i,"Age"] =   inf_year - bl[i,"BirthYear"]
  
  
}

inf_year = bl[,"BirthYear"] + bl[,"Age"]

bl[,"Gender"] <- rbinom(cohortSize, 1, 0.528)
bl[,"HIV"] <- rbinom(cohortSize, 1, 0.002) # 0.9% among non-IDUs
bl[,"IDU"] <- rbinom(cohortSize, 1, 0)
####################################################################
bl[,"HIV_time"] <- (bl[,"HIV"])
#######################
for (i in 1: cohortSize)
{
  bl[i,"HIV_time"] <- 999
  w = runif(1, 0, 1)
  
  if (bl[i,"HIV"] == 1)
    bl[i,"HIV_time"] = 0 # time having hiv
  
  if (bl[i,"HIV"] == 0)
  {
    ID = ifelse( w < p_HIV, w, -1)
    bl[i,"HIV_time"] = ifelse( ID == -1, 999, runif(1, 0.25, maxTime)) # time having hiv
  }
}
############################################################################################
bl[,"IDU_time_start"] <- bl[,"IDU"]
bl[,"IDU_time_stop"] <- bl[,"IDU"]
for (i in 1: cohortSize)
{
  w = runif(1, 0, 1)
  if (bl[i,"IDU"] == 0)
  {
    IDU_ID = ifelse( w < p_IDU_start, w, 0)
    bl[i,"IDU_time_start"] = ifelse( IDU_ID == 0, 999, maxTime * IDU_ID + runif(1, 0, maxTime) * rbinom(1, 1, 0.5)) # time starting idu
  }
  if (bl[i,"IDU"] == 1)
    bl[i,"IDU_time_start"] = 0 # time having hiv
  teta = runif(1, 0, 1)
  bl[i,"IDU_time_stop"] = ifelse ( bl[i,"IDU_time_start"] < maxTime && teta < p_IDU_stop, bl[i,"IDU_time_start"] + 0.2 * teta * maxTime , 999)
}
bl[,"IDU_time_start"] <- bl[,"IDU"]
bl[,"IDU_time_stop"] <- bl[,"IDU"]
##########################################################################################
bl[,"sr.time"] <- sr.fun()
###################################################################################
### Bl diagnosis time (time of diagnosis for each patient)
diag.time.fun = function(cohortSize, maxTime, period, screen_proba, test_type, V0,bl)
{
  bl.atb<-runif(cohortSize, 0, period)
  
  screen.t<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  diag.times<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  w<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  prob.detect.atb<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  prob.detect.alt<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  prob.detect.rna<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  prob.detect.antg<-matrix(nrow=cohortSize, ncol=floor((maxTime/period)+1))
  
  for(j in 1:ncol(screen.t))
  {
    screen.t[,j]<-bl.atb + ((j-1)*period)
    prob.detect.atb[,j]<-ifelse(screen.t[,j]>1, 0.95, 0.95*(screen.t[,j])^0.273)
    prob.detect.rna[,j]<-1
    prob.detect.alt[,j]<-ifelse(screen.t[,j]<=0.25, 0.88, 0)
    prob.detect.antg[,j] <- ifelse (V0 > 3, 0.934,  0.05)
    
    
    if(test_type=="atb") {prob.detect<-prob.detect.atb}
    else if(test_type=="atb+antg") {prob.detect<-prob.detect.antg} # as the minimum time should be chosen for the next stage
    else if(test_type=="atb+pcr") {prob.detect<-prob.detect.rna}
    else if(test_type=="atb-antg") {prob.detect<- (prob.detect.antg * prob.detect.atb)}
    else if(test_type=="atb-pcr") {prob.detect<-(prob.detect.rna * prob.detect.atb)}
    else if(test_type=="antg") prob.detect<-prob.detect.antg
    else prob.detect<-prob.detect.rna
    
    
    Year_Screen = bl[,"BirthYear"] + bl[,"Age"] + (j - 1) * period
    tmp = runif(cohortSize,0,0.1)
    prob.detect[,j] <- ifelse(2015 - Year_Screen>= 0,0, prob.detect[,j] )
    #if (j > 1)
    screen.t[,j] = ifelse((j > 1) & (2015 - (bl[,"BirthYear"] + bl[,"Age"] + (j - 2) * period) < 0),666,screen.t[,j]) # if it is screened last year, 
    # screen.t[,j] <- 666 # then it is not possible to be tested again
    # screen.t[,j] <- ifelse(2015 - Year_Screen< -1,999, screen.t[,j])
    for(i in 1:cohortSize)
      w[i,j]=runif(1,0,1)
    diag.times[,j]<- ifelse(w[,j]<=prob.detect[,j], screen.t[,j], 666)
    
  }
  
  diag.time<-rep(0, cohortSize)
  for(i in 1:cohortSize)
  {
    diag.time[i]<-min(diag.times[i,])
    w=runif(1,0,1)
    if(w>screen_proba) diag.time[i]<-666 # to adjust proba to be screened!
  }
  
  c(diag.time)
}
diag_results = diag.time.fun(cohortSize, maxTime, period=4, screen_proba=1, test_type,  V0,bl)
bl[,"Diag.time"]<- diag_results[1:cohortSize]
bl[,"Viral_diag"] = V0

###################################################################################
#bm <- read.table("Mort_CH_Swiss.txt", header=T) # Read file containing mortality rates by age groups in Switzerland
bm <- read.table("Mort_CH", header=T) # Read file containing mortality rates by age groups in Switzerland
bm_M <-c(bm[,3],rep(1,100)) # to set background mortality of patients >75 years old as = background mortality of patients 75 years old
bm_F <-c(bm[,4],rep(1,100)) # to set background mortality of patients >75 years old as = background mortality of patients 75 years old

###########################################

###################################################################################
