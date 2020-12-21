calendartime_matrix <- matrix(NA, nrow = nrow(Cohort_comb), ncol = 100)
Cohort_calendar <- Cohort_comb[,(1:52)] + Cohort_comb[,"Age"] + Cohort_comb[,"BirthYear"]
Cohort_comb[,(1:52)] = Cohort_comb[,(1:52)] + Cohort_comb[,"Age"] + Cohort_comb[,"BirthYear"]
indexmatrix <- matrix((1:52),nrow = nrow(Cohort_calendar), ncol = 52, byrow = TRUE)
for (year in 1:100) 
{
  calyear <- year + 1969
  temp <- (Cohort_calendar < calyear + 1) * indexmatrix
  calendartime_matrix[,year] <- apply(temp, 1, max, na.rm = TRUE)
}
calendartime_matrix2 = calendartime_matrix
calendartime_matrix2[,1:45] = 0 # we are only interested to look at 2015 and afterward
calendartime_matrix_main = calendartime_matrix
calendartime_matrix = calendartime_matrix2
 
calendartime_matrix = calendartime_matrix[which((Cohort_comb[,1] < 2018 | Cohort_comb[,2] < 2018 | Cohort_comb[,3] < 2018 | Cohort_comb[,4] < 2018) &
                                                  (Cohort_comb[,49] > 2016 | Cohort_comb[,50] > 2016 | Cohort_comb[,51] > 2016 | Cohort_comb[,52] > 2016)),]

calendartime_matrix_main2 = calendartime_matrix_main[which((Cohort_comb[,1] < 2018 | Cohort_comb[,2] < 2018 |
                                                              Cohort_comb[,3] < 2018 | Cohort_comb[,4] < 2018) &
                                                             (Cohort_comb[,49] > 2016 | Cohort_comb[,50] > 2016 |
                                                                Cohort_comb[,51] > 2016 | Cohort_comb[,52] > 2016)),]

Cohort_comb = Cohort_comb[which((Cohort_comb[,1] < 2018 | Cohort_comb[,2] < 2018 |
                                   Cohort_comb[,3] < 2018 | Cohort_comb[,4] < 2018) &
                                  (Cohort_comb[,49] > 2016 | Cohort_comb[,50] > 2016 | Cohort_comb[,51] > 2016 | Cohort_comb[,52] > 2016)),]

enterSVRtime = apply(Cohort_comb[,c(26:30, 36,42,48:52)],1,min,na.rm=TRUE)
Cohort_comb = Cohort_comb[which(enterSVRtime > 2015),]
calendartime_matrix = calendartime_matrix[which(enterSVRtime > 2015),]
enterSVRtime = enterSVRtime[which(enterSVRtime > 2015)]


#################################################################################################################################
############################# Test what happen if we remove False Negative ###############################################################################
# Cohort_comb2 = Cohort_comb
# Cohort_comb2[(which(Cohort_comb[,"Diag.time"] == 666 & enterSVRtime > 2017& !is.infinite(enterSVRtime))),] = NA
# Cohort_comb = Cohort_comb2
################################################################################################################################


cohortSize = nrow(Cohort_comb)
#################################################################################################################################
## define alpha
################################################################################################################################
{
alpha = matrix(0, nrow = nrow(Cohort_comb), ncol = 13)
colnames(alpha) <- c("tF0","tF1","tF2","tF3","tF4","tDC","tHCC","Death","SVR_F0","SVR_F1","SVR_F2","SVR_F3","SVR_F4")

Liver_T = matrix(0, nrow = nrow(Cohort_comb), ncol = 1)
###########################################################################################################################
# F0 to F1
###########################################################################################################################
for(j in 0:4)
  for(i in 1:nrow(Cohort_comb)) # for each patient
  {
    temp1 = max(Cohort_comb[i,c(1+j,6+j,11+j,16+j,21+j)],na.rm=TRUE) 
    if (!is.infinite(temp1)){
      h = c(1+j,6+j,11+j,16+j,21+j)
      tmp = (Cohort_comb[i, c((h[which.max(Cohort_comb[i,c(1+j,6+j,11+j,16+j,21+j)])]+1):52)])
      temp1 = tmp[min(which(!is.na(tmp)))] 
      temp1 = as.numeric(temp1)
    }
    
    temp0 =  min(Cohort_comb[i,c(1+j,6+j,11+j,16+j,21+j)],na.rm=TRUE)
    
    if (!is.infinite(temp1) && length(temp1) > 0  && length(temp0) > 0 && (temp1 - 2015) >= 0)
    {
      temp0 = max((temp0 - 2015),0)
      temp1 = (temp1 - 2015)
      alpha[i,(j+1)] = (exp(-temp0 * log(1.03))/ log(1.03)) - (exp(-temp1 * log(1.03))/ log(1.03))
    }
  }

############################################################################################################################
# LT
############################################################################################################################

for(i in 1:nrow(Cohort_comb)) # for each patient
{
  temp1 = max(Cohort_comb[i,c(43:48)],na.rm=TRUE) 
  if (!is.infinite(temp1)){
    h = c(43:48)
    tmp = (Cohort_comb[i, c((h[which.max(Cohort_comb[i,c(43:48)])]+1):52)])
    temp1 = tmp[min(which(!is.na(tmp)))] 
    temp1 = as.numeric(temp1)
  }
  
  
  temp0 =  min(Cohort_comb[i,c(41:48)],na.rm=TRUE) 
  if (!is.infinite(temp1) && length(temp1) > 0  && length(temp0) > 0 && (temp1 - 2015) >= 0)
  {
    temp0 = max((temp0 - 2015),0)
    temp1 = (temp1 - 2015)
    Liver_T[i,1] = (exp(-temp0 * log(1.03))/ log(1.03)) - (exp(-temp1 * log(1.03))/ log(1.03))
  }
}
###############################################################################################################################
# DC
###############################################################################################################################
for(i in 1:nrow(Cohort_comb)) # for each patient
{
  temp1 = max(Cohort_comb[i,c(31:36)],na.rm=TRUE) 
  if (!is.infinite(temp1)){
    h = c(31:36)
    tmp = (Cohort_comb[i, c((h[which.max(Cohort_comb[i,c(31:36)])]+1):52)])
    temp1 = tmp[min(which(!is.na(tmp)))] 
    temp1 = as.numeric(temp1)
  }
  
  
  temp0 =  min(Cohort_comb[i,c(31:36)],na.rm=TRUE) 
  if (!is.infinite(temp1) && length(temp1) > 0  && length(temp0) > 0 && (temp1 - 2015) >= 0)
  {
    temp0 = max((temp0 - 2015),0)
    temp1 = (temp1 - 2015)
    alpha[i,6] = (exp(-temp0 * log(1.03))/ log(1.03)) - (exp(-temp1 * log(1.03))/ log(1.03))
  }
  ##############################################################################################################################
  # LT
  # ############################################################################################################################
  # temp1 = max(Cohort_comb[i,c(43:48)],na.rm=TRUE) 
  # if (!is.infinite(temp1)){
  #   h = c(43:48)
  #   tmp = (Cohort_comb[i, c((h[which.max(Cohort_comb[i,c(43:48)])]+1):52)])
  #   temp1 = tmp[min(which(!is.na(tmp)))] 
  #   temp1 = as.numeric(temp1)
  # }
  # 
  # 
  # temp0 =  min(Cohort_comb[i,c(41:48)],na.rm=TRUE) 
  # if (!is.infinite(temp1) && length(temp1) > 0  && length(temp0) > 0 && (temp1 - 2015) >= 0)
  # {
  #   temp0 = max((temp0 - 2015),0)
  #   temp1 = (temp1 - 2015)
  #   Liver_T[i,1] = (exp(-temp0 * log(1.03))/ log(1.03)) - (exp(-temp1 * log(1.03))/ log(1.03))
  # }
  # 
  ###################################################################################################################################
  temp1 = max(Cohort_comb[i,c(37:42)],na.rm=TRUE) 
  if (!is.infinite(temp1)){
    h = c(37:42)
    tmp = (Cohort_comb[i, c((h[which.max(Cohort_comb[i,c(37:42)])]+1):52)])
    temp1 = tmp[min(which(!is.na(tmp)))] 
    temp1 = as.numeric(temp1)
  }
  
  temp0 =  min(Cohort_comb[i,c(37:42)],na.rm=TRUE)  
  if (!is.infinite(temp1) && length(temp1) > 0  && length(temp0) > 0 && (temp1 - 2015) >= 0)
  {
    temp0 = max((temp0 - 2015),0)
    temp1 = (temp1 - 2015)
    alpha[i,7] = (exp(-temp0 * log(1.03))/ log(1.03)) - (exp(-temp1 * log(1.03))/ log(1.03))
  }
  ############################################################################################################################# 
  if(! is.na(Cohort_comb[i,52]))
    alpha[i,8] = 1
  
  for (j in 0:4){
    temp1 = Cohort_comb[i,c(26+j)] 
    
    temp1 =  ifelse(max(which(calendartime_matrix[i,] == (26+j))) <100 && !is.na(temp1) &&
                      !is.infinite(max(which(calendartime_matrix[i,] == (26+j)))),
                    Cohort_comb[i,calendartime_matrix[i,max(which(calendartime_matrix[i,] == (26+j))+1)]],
                    -1)
    
    temp0 =  min(Cohort_comb[i,(26+j)],na.rm=TRUE) 
    if (!is.infinite(temp1) && length(temp1) > 0  && length(temp0) > 0 && (temp1 - 2015) >= 0)
    {
      temp0 = max((temp0 - 2015),0)
      temp1 = (temp1 - 2015)
      alpha[i,(9 + j)] = (exp(-temp0 * log(1.03))/ log(1.03)) - (exp(-temp1 * log(1.03))/ log(1.03))
    }
    
  }
}

##############################################################################################################
}
###############################################################################################################
## Number of diagnosis
################################################################################################################
diagnosistime=apply(Cohort_comb[,(11:15)],1,min,na.rm=TRUE)

new_diag <- vector("numeric", length=100)
for (year in 45:100)
{
  diag = count(diagnosistime < (year + 1970) & diagnosistime >= (year + 1969))
  if(length(which(diag$x == TRUE))>0)
    new_diag[year] <- diag$freq[which(diag$x == TRUE)]
  
}
new_diag2 <- vector("numeric", length=100)

if(length(Cohort_comb[which(!is.na(Cohort_comb[,33])&is.na(Cohort_comb[,32])&is.na(Cohort_comb[,15])),33]) > 0)
{
  tmp = count(floor(Cohort_comb[which(!is.na(Cohort_comb[,33])&is.na(Cohort_comb[,32])&is.na(Cohort_comb[,15])),33]))
  year = tmp$x - 1969
  diag2[year] = tm$freq 
}

if(length(Cohort_comb[which(!is.na(Cohort_comb[,39])&is.na(Cohort_comb[,38])&is.na(Cohort_comb[,15])),39]) > 0)
{
  tmp = count(floor(Cohort_comb[which(!is.na(Cohort_comb[,39])&is.na(Cohort_comb[,38])&is.na(Cohort_comb[,15])),39]))
  year = tmp$x - 1969
  diag2[year] = tm$freq 
}
if(length(Cohort_comb[which(!is.na(Cohort_comb[,45])&is.na(Cohort_comb[,44])&is.na(Cohort_comb[,15])), 45]) > 0)
{
  tmp = count(floor(Cohort_comb[which(!is.na(Cohort_comb[,45])&is.na(Cohort_comb[,44])&is.na(Cohort_comb[,15])),45]))
  year = tmp$x - 1969
  diag2[year] = tm$freq 
}
new_diag3 =  new_diag + new_diag2
####################################################################################

#####################################################################################
cost_stages_matrix = rep.row(cost_stages,nrow(alpha))
utility_stages_matrix = rep.row(utility_stages,nrow(alpha))

time_cost = sum(alpha * cost_stages_matrix)

time_utility = sum(alpha * utility_stages_matrix)

############################################################################################################################
# New diagnosis
##########################################################################################################################  
total_diagnosed_total <- vector("numeric", length=100)
total_infected_total <- vector("numeric", length=100)

year <- 1
total_diagnosed_total[year] <- sum((calendartime_matrix[,year]> 0) *
                                     ((calendartime_matrix[,year]>=11 & calendartime_matrix[,year]<=20) |
                                        calendartime_matrix[,year]==33 | calendartime_matrix[,year]==39))
for(year in 2:100)
  total_diagnosed_total[year] <- sum(((calendartime_matrix[,year]>10 & calendartime_matrix[,year]<26) 
                                      | (calendartime_matrix[,year]>=33 & calendartime_matrix[,year]<=35) |
                                        (calendartime_matrix[,year]>=39 & calendartime_matrix[,year]<=41) | 
                                        (calendartime_matrix[,year]>=45 & calendartime_matrix[,year]<=47)))

for(year in 45:100)
  total_infected_total[year] <- sum(((calendartime_matrix[,year]>0 & calendartime_matrix[,year]<26) |
                                       (calendartime_matrix[,year]>30 & calendartime_matrix[,year]<48 & 
                                          calendartime_matrix[,year]!=36 & calendartime_matrix[,year]!=42)))


######################################################################################################################
# The number of treated people
#######################################################################################################################
new_treated <- vector("numeric", length=100)
year = 1
new_treated[year] <- sum(
  (calendartime_matrix[,year]>=16 & calendartime_matrix[,year]<=20) |  
    (calendartime_matrix[,year]>=21 & calendartime_matrix[,year]<=25))

for (year in 2:100)
{
  new_treated[year] <- sum(
    ( 
      (calendartime_matrix[,year-1] < 16 &
         (calendartime_matrix[,year]>=16 & calendartime_matrix[,year]<=20)) |
        
        (calendartime_matrix[,year-1] < 21 &
           (calendartime_matrix[,year]>=21 & calendartime_matrix[,year]<=25))
    ))
}

new_treated2 =  0 * new_treated
for (year in 2:100)
{
  new_treated2[year] <- sum(  (calendartime_matrix[,year]> 0) *
                                ((calendartime_matrix[,year-1] == 33 & calendartime_matrix[,year]== 34) |
                                   (calendartime_matrix[,year-1] == 39 & calendartime_matrix[,year]== 40) | 
                                   (calendartime_matrix[,year-1] == 45 & calendartime_matrix[,year]== 46) |
                                   (calendartime_matrix[,year-1] == 34 & calendartime_matrix[,year]== 35) |
                                   (calendartime_matrix[,year-1] == 40 & calendartime_matrix[,year]== 41) |
                                   (calendartime_matrix[,year-1] == 46 & calendartime_matrix[,year]== 47) 
                                 
                                ) 
  )
}
new_treated3 =  new_treated2 +  new_treated
################################################################################
# New Infection
##################################################################################

inftime=apply(Cohort_comb[,(1:25)],1,min,na.rm=TRUE)

enterTreatment = apply(Cohort_comb[,c(16:25, 34,40,46,35,41,47)],1,min,na.rm=TRUE)
aa = (which(!is.infinite(enterTreatment)))
bb = enterSVRtime[aa] # treat cases
enterSVRtime[aa] = Inf


total_inf <- vector("numeric", length=100)
for (year in 46:100)
  total_inf[year] = sum(year + 1969 < enterSVRtime & inftime <= year + 1969)


#Total_population =  total_inf[46]/ 0.054 #(5.4 is the prevalence)
Total_population =  total_inf[46]/ 0.0007 #(5.4 is the prevalence)
TN = 0.946 * Total_population
if(ans_IDU == "Yes")
{
  Total_population =  total_inf[46]/ 0.662 #(66.2 is the prevalence)
  #Total_population =  total_inf[46]/ 0.94 #(94% in high prevalence situation:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3072734/)
  TN = 0.338 * Total_population
}
########################################################################################
## Death
######################################################################################
death_2015_2016 = length(which(floor(Cohort_comb[,49]) == 2015 | floor(Cohort_comb[,52]) == 2015 | 
                                 floor(Cohort_comb[,49]) == 2016 | floor(Cohort_comb[,52]) == 2016 |
                                 floor(Cohort_comb[,50]) == 2015 |floor(Cohort_comb[,50]) == 2016|
                                 floor(Cohort_comb[,51]) == 2015 |floor(Cohort_comb[,51]) == 2016))


#####################################################################################################################
#################################### Number of FN, beta  and time ##################################################################

inftime=apply(Cohort_comb[,(1:25)],1,min,na.rm=TRUE)
nFN = length(which(Cohort_comb[,"Diag.time"] == 666 & enterSVRtime > 2017& !is.infinite(enterSVRtime)))

a = (which(Cohort_comb[,"Diag.time"] == 666))
b = (which(enterSVRtime > 2017 & !is.infinite(enterSVRtime)))
c = intersect(a, b)
nFN = length(c)



diagnosistime2=apply(Cohort_comb[,c(11:15, 33,39,45)],1,min,na.rm=TRUE)
t = (diagnosistime2 - 2015)

beta =  1/(1.03 ^ t)
##################################################################################################################
treattime=apply(Cohort_comb[,(16:25)],1,min,na.rm=TRUE)
t1 = (treattime - 2015)
t1[which(is.infinite(t1))] = 0
beta1 = 1/(1.03 ^ t1)
beta1[which(beta1 == 1)] = 0