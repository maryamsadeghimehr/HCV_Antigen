ns = 49
bl_number = 5
number_blocks = 6
fs=5

pp = as.matrix(read.table(sprintf("Total_results_%s.txt",ans)))

cohort@time.to.state = as.data.frame(pp[, 1:49])
cohort@baseline = as.matrix(pp[, 50:61])


# cohort.cirr includes only times for incident cirrhosis and not those who had cirrhosis in an earlier block
cohort.cirr.inc<-cohort
cirr<-cbind(cohort@time.to.state[,5], cohort@time.to.state[,10], cohort@time.to.state[,15],
            cohort@time.to.state[,20], cohort@time.to.state[,25], cohort@time.to.state[,30])

cirr.inc <- cirr

for(j in 2:5)
  for(i in j:5)
    cirr.inc[,i][cirr.inc[,j - 1]!="NA"] <- NA

cohort.cirr.inc@time.to.state[,5] <- cirr.inc[,1]
cohort.cirr.inc@time.to.state[,10] <- cirr.inc[,2]
cohort.cirr.inc@time.to.state[,15]<-cirr.inc[,3]
cohort.cirr.inc@time.to.state[,20]<-cirr.inc[,4]
cohort.cirr.inc@time.to.state[,25]<-cirr.inc[,5]
cohort.cirr.inc@time.to.state[,30]<-cirr.inc[,6]
######################################################
### we use cohort.cirr.inc to avoid counting several times those in F4 who move from one block to the next
cohortize = nrow(cohort.cirr.inc@time.to.state)
perc_states <- 1:ns
for(i in 1:ns){
  perc_states[i]<-(sum(!is.na(cohort.cirr.inc@time.to.state[,i]))*100/cohortize)
  }
#for(i in 1:ns){perc_states[i]<-(sum(!is.na(cohort[,i]))*100/cohortize)}

cinc <- cbind(states=c(1:ns), "cumulative incidence"=perc_states)
print("Percentage of people that went through each states by the end of the simulation:")
print(cinc)
xx <- readline(prompt="Do you want to save the results as a table (y/n)?")
#results$newcolumn<- results$cinc[,2]
if (xx == "y")
{
  write.table(cinc,"results")
}
perc.cirr<-cinc[5,2]
cirr <- cinc[5,2] + cinc[10,2] + cinc[15,2] + cinc[20,2] + cinc[25,2] + cinc[30,2]
perc.DC <- cinc[31,2] + cinc[32,2] + cinc[33,2] + cinc[34,2] + cinc[35,2] + cinc[36,2]
perc.HCC <- cinc[37,2] + cinc[38,2] + cinc[39,2] + cinc[40,2] + cinc[41,2] + cinc[42,2]
perc.LT <- cinc[43,2] + cinc[44,2] + cinc[45,2] + cinc[46,2] + cinc[47,2] + cinc[48,2]

message("cirrhosis")
print(cirr)
message("DC")
print(perc.DC)
message("HCC")
print(perc.HCC)
message("LT")
print(perc.LT)




####################
n_endpoints = 6
ns = 49
bl_number = 5
number_blocks = 6
fs = 5
maxTime = 100
#####################
cohort = new("ArtCohort")
##########################
pp = as.matrix(read.table(ans))
pp = list(pp)
#########################################
est<-matrix(NA,length(pp),n_endpoints)
lower<-matrix(NA,length(pp),n_endpoints)
upper<-matrix(NA,length(pp),n_endpoints)

####################################################
for ( k in 1: length(pp))
{

  cohort@time.to.state = as.data.frame(pp[[k]][,1:ns])

  # cohort.cirr is a sub cohort that includes only times for incident cirrhosis and not those who had cirrhosis in an earlier block
  cohort.cirr.inc<- cohort
  cirr <- cbind(cohort@time.to.state[,5], cohort@time.to.state[,10], cohort@time.to.state[,15],
                cohort@time.to.state[,20], cohort@time.to.state[,25], cohort@time.to.state[,30]) # all F4

  cirr.inc <- cirr

  for(j in 2:5)
    for(i in j:5)
      cirr.inc[,i][cirr.inc[,j-1]!="NA"] <- NA

  cohort.cirr.inc@time.to.state[,5]<-cirr.inc[,1]
  cohort.cirr.inc@time.to.state[,10]<-cirr.inc[,2]
  cohort.cirr.inc@time.to.state[,15]<-cirr.inc[,3]
  cohort.cirr.inc@time.to.state[,20]<-cirr.inc[,4]
  cohort.cirr.inc@time.to.state[,25]<-cirr.inc[,5]
  cohort.cirr.inc@time.to.state[,30]<-cirr.inc[,6]
  #########################################################
  cinc <- cumulativeIncidence(cohort.cirr.inc, seq(0,maxTime,1))

  # plot(cinc, main= "Cumulative incidence", ci=TRUE) #par(mar=c(1,1,1,1))
  ###########################################################################
  est.F0 <- cinc@probabilities[maxTime,1] + cinc@probabilities[maxTime,6]+cinc@probabilities[maxTime,11] +
    cinc@probabilities[maxTime,16] + cinc@probabilities[maxTime,21] + cinc@probabilities[maxTime,26]
  ####################################################################################################
  est.F1 <- cinc@probabilities[maxTime,2] + cinc@probabilities[maxTime,7]+cinc@probabilities[maxTime,12] +
    cinc@probabilities[maxTime,17] + cinc@probabilities[maxTime,22] + cinc@probabilities[maxTime,27]
  ################################################################################################
  est.F2 <- cinc@probabilities[maxTime,3] + cinc@probabilities[maxTime,8]+cinc@probabilities[maxTime,13] +
    cinc@probabilities[maxTime,18] + cinc@probabilities[maxTime,23] + cinc@probabilities[maxTime,28]
  #######################################################################################################
  est.F3 <- cinc@probabilities[maxTime,4] + cinc@probabilities[maxTime,9]+cinc@probabilities[maxTime,14] +
    cinc@probabilities[maxTime,19] + cinc@probabilities[maxTime,24] + cinc@probabilities[maxTime,29]
  #######################################################################################################
  #######################################################################################################
  est.cirr<-cinc@probabilities[maxTime,5] + cinc@probabilities[maxTime,10]+cinc@probabilities[maxTime,15] +
    cinc@probabilities[maxTime,20] + cinc@probabilities[maxTime,25] + cinc@probabilities[maxTime,30]
  lower.cirr<-cinc@lower[maxTime,5]+cinc@lower[maxTime,10]+cinc@lower[maxTime,15]+cinc@lower[maxTime,20]+cinc@lower[maxTime,25]+ cinc@lower[maxTime,30]
  upper.cirr<-cinc@upper[maxTime,5]+cinc@upper[maxTime,10]+cinc@upper[maxTime,15]+cinc@upper[maxTime,20]+cinc@upper[maxTime,25]+ cinc@upper[maxTime,30]
 ############################################################################################################
  est.dc<- cinc@probabilities[maxTime,31] + cinc@probabilities[maxTime,32] + cinc@probabilities[maxTime,33] + cinc@probabilities[maxTime,34] + cinc@probabilities[maxTime,35] + cinc@probabilities[maxTime,36]
  lower.dc<-cinc@lower[maxTime,31] + cinc@lower[maxTime,32] + cinc@lower[maxTime,33] + cinc@lower[maxTime,34] + cinc@lower[maxTime,35] + cinc@lower[maxTime,36]
  upper.dc<-cinc@upper[maxTime,31] + cinc@upper[maxTime,32] + cinc@upper[maxTime,33] + cinc@upper[maxTime,34] + cinc@upper[maxTime,35] + cinc@upper[maxTime,36]
  #############################################################################################################
  est.hcc<-cinc@probabilities[maxTime,37] + cinc@probabilities[maxTime,38] + cinc@probabilities[maxTime,39] + cinc@probabilities[maxTime,40] + cinc@probabilities[maxTime,41] + cinc@probabilities[maxTime,42]
  lower.hcc<-cinc@lower[maxTime,37] + cinc@lower[maxTime,38] + cinc@lower[maxTime,39] + cinc@lower[maxTime,40] + cinc@lower[maxTime,41] + cinc@lower[maxTime,42]
  upper.hcc<-cinc@upper[maxTime,37] + cinc@upper[maxTime,38] + cinc@upper[maxTime,39] + cinc@upper[maxTime,40] + cinc@upper[maxTime,41] + cinc@upper[maxTime,42]
 #############################################################################################################

  est.lt<-cinc@probabilities[maxTime,43] + cinc@probabilities[maxTime,44] + cinc@probabilities[maxTime,45] + cinc@probabilities[maxTime,46] + cinc@probabilities[maxTime,47] + cinc@probabilities[maxTime,48]
  lower.lt<-cinc@lower[maxTime,43] + cinc@lower[maxTime,44] + cinc@lower[maxTime,45] + cinc@lower[maxTime,46] + cinc@lower[maxTime,47] + cinc@lower[maxTime,48]
  upper.lt<-cinc@upper[maxTime,43] + cinc@upper[maxTime,44] + cinc@upper[maxTime,45] + cinc@upper[maxTime,46] + cinc@upper[maxTime,47] + cinc@upper[maxTime,48]
  ###############################################################################################################################
  est.death<-cinc@probabilities[maxTime,49]
  lower.death<-cinc@lower[maxTime,49]
  upper.death<-cinc@upper[maxTime,49]

  #####################################
  ### liver related deaths
  #####################################
  # cohort.ldeath is a subcohort where state 49 is only liver realted death and not total deaths
  cohort.ldeath <- cohort
  cohort.ldeath@time.to.state[is.na(cohort.ldeath@time.to.state)]<- 666
  for (i in 1: nrow(cohort.ldeath@time.to.state))
  {

    cohort.ldeath@time.to.state[i,49][
      cohort.ldeath@time.to.state[i,31] == 666 && cohort.ldeath@time.to.state[i,32] == 666 &&
        cohort.ldeath@time.to.state[i,33] == 666 && cohort.ldeath@time.to.state[i,34] == 666 &&
        cohort.ldeath@time.to.state[i,35] == 666 && cohort.ldeath@time.to.state[i,36] == 666 &&
        cohort.ldeath@time.to.state[i,37] == 666 && cohort.ldeath@time.to.state[i,38] == 666 &&
        cohort.ldeath@time.to.state[i,39] == 666 && cohort.ldeath@time.to.state[i,40] == 666 &&
        cohort.ldeath@time.to.state[i,41] == 666 && cohort.ldeath@time.to.state[i,42] == 666
      ]<- NA

  }

  cincldeath<-cumulativeIncidence(cohort.ldeath, seq(0,maxTime,1))
  est.ldeath<-cincldeath@probabilities[maxTime,49]
  lower.ldeath<-cincldeath@lower[maxTime,49]
  upper.ldeath<-cincldeath@upper[maxTime,49]
  ####################################################################################
  ## liver related deaths happening without replicating HCV (After HCV clearance, in block undetectable HCV)
  ## cohort.ldeath.sub is a cohort where column 49 (death) includes only liver related death occuring in patients without detectable HCV

  cohort.ldeath.sub<-cohort.ldeath

  for (i in 1: nrow(cohort.ldeath.sub@time.to.state))
  {
    cohort.ldeath.sub@time.to.state[i,49][cohort.ldeath.sub@time.to.state[i,20] == 666 ||cohort.ldeath.sub@time.to.state[i,15] == 666  ]<-NA }

  cincldeath.sub<-cumulativeIncidence(cohort.ldeath.sub, seq(0,maxTime,1))
  est.ldeath.sub<-cincldeath.sub@probabilities[maxTime,49]
  lower.ldeath.sub<-cincldeath.sub@lower[maxTime,49]
  upper.ldeath.sub<-cincldeath.sub@upper[maxTime,49]
  #################################################################################################################################

  est[k,]<-c(est.cirr, est.dc, est.hcc, est.death, est.ldeath, est.ldeath.sub)
  colnames(est)<-c("F4", "DC", "HCC", "death", "liv.death.tot", "liv.death.noHCV")
  rownames(est)<-c("treatF3")

  lower[k,]<-c(lower.cirr, lower.dc, lower.hcc, lower.death, lower.ldeath, lower.ldeath.sub)
  colnames(lower)<-c("F4", "DC", "HCC", "death", "liv.death.tot", "liv.death.noHCV")
  rownames(lower)<-c("treatF3")

  upper[k,]<-c(upper.cirr, upper.dc, upper.hcc, upper.death, upper.ldeath, upper.ldeath.sub)
  colnames(upper)<-c("F4", "DC", "HCC", "death", "liv.death.tot", "liv.death.noHCV")
  rownames(upper)<-c("treatF3")
  #############################################################################################################
}
write.table(est, "est_bc")
write.table(lower, "lower_bc")
write.table(upper, "upper_bc")
###########################################################################################################
#matrix of 1 if the patient went to the state, 0 if not:
who = 0*(cohort@time.to.state)
who[cohort@time.to.state >= 0] = 1
who[is.na(cohort@time.to.state)] = 0
##number of patients in each state:
how.many = apply(who,2,sum)/cohortSize
how.many2 = apply(who,2,sum)
advanced = 100*cbind(sum(how.many[c(31: 36)]),sum(how.many[c(37:41)]) , sum(how.many[c(5,10,15,20,25,30)]))
colnames(advanced) = c("DC","HCC","Cirrhosis")
print(advanced)
##############################################################################################
advanced2 = 100*cbind((how.many[c(31: 36)]),(how.many[c(37:42)]) ,(how.many[c(43:48)]), (how.many[c(49:52)]),(how.many[c(5,10,15,20,25,30)]))
colnames(advanced2) = c("DC","HCC","LT","Death","Cirrhosis")
rownames(advanced2) = c("acute","chronic","diag","treated","Second_treatd","SVR + SR")
print("The percentage of population that have the experience of being at each stage:")
print(advanced2)
#############################################################################################
advanced3 = 100*cbind((how.many[1]) + (how.many[6]) + (how.many[11]) + (how.many[16]) + (how.many[21]) + (how.many[26])/6,
                      (how.many[2]) + (how.many[7]) + (how.many[12]) + (how.many[17]) + (how.many[22]) + (how.many[27])/6,
                      (how.many[3]) + (how.many[8]) + (how.many[13]) + (how.many[18]) + (how.many[23]) + (how.many[28])/6,
                      (how.many[4]) + (how.many[9]) + (how.many[14]) + (how.many[19]) + (how.many[24]) + (how.many[29])/6,
                      (how.many[5]) + (how.many[10]) + (how.many[15]) + (how.many[20]) + (how.many[25]) + (how.many[30])/6)

colnames(advanced3) = c("F0","F1","F2","F3","F4")
rownames(advanced3) = c("%Ave")
print("The percentage of population that have the experience of being at each stage:")
print(advanced3)
##############################################################################################
#############################################################################################
advanced4 = 100*cbind((how.many2[1]) + (how.many2[6]) + (how.many2[11]) + (how.many2[16]) + (how.many2[21]) + (how.many2[26])/6,
                      (how.many2[2]) + (how.many2[7]) + (how.many2[12]) + (how.many2[17]) + (how.many2[22]) + (how.many2[27])/6,
                      (how.many2[3]) + (how.many2[8]) + (how.many2[13]) + (how.many2[18]) + (how.many2[23]) + (how.many2[28])/6,
                      (how.many2[4]) + (how.many2[9]) + (how.many2[14]) + (how.many2[19]) + (how.many2[24]) + (how.many2[29])/6,
                      (how.many2[5]) + (how.many2[10]) + (how.many2[15]) + (how.many2[20]) + (how.many2[25]) + (how.many2[30])/6)

colnames(advanced4) = c("F0","F1","F2","F3","F4")
rownames(advanced4) = c("Freq")
print("The number of population that have the experience of being at each stage:")
print(advanced4)



 Primary = 100 * cbind((how.many[c(5,10,15,20,25,30)]),(how.many[c(4,9,14,19,24,29)]), (how.many[c(3,8,13,18,23,28)]), (how.many[c(2,7,12,17,22,27)]), (how.many[c(1,6,11,16,21,26)]))
 colnames(Primary) = c("F4","F3","F2","F1","F0")
 rownames(Primary) = c("acute","chronic","diag","treated","Second_treatd","SVR + SR")
 print("The percentage of population that have the experience of being at each stage:")
 print(Primary)
 ####################################
times =  seq(0, maxTime, length = 100)
posteriorProbabilities = cumulativeIncidence
post <- posteriorProbabilities(cohort.cirr.inc, times)
##cirrhosis !!!WRONG, NEED TO SUBSTRACT THOSE WHO MOVE HORIZONTALLY FROM ONE CIRRHOITC STAGE TO THE NEXT!!!!!!!!!
post.cirr<-post@probabilities[,5]+post@probabilities[,10]+post@probabilities[,15]+post@probabilities[,20]+post@probabilities[,25] + post@probabilities[,30]
#### End Stage Liver Disease
post.esld <- post@probabilities[,31]
for (i in 32:42)
  post.esld <- post@probabilities[,i] + post.esld

post.treat <- post@probabilities[,16]
for (i in 17:20)
  post.treat<- post@probabilities[,i] + post.treat
###To plot cirr
plot(times, post.cirr, type="l", col="blue", main="Cirrhosis", lty="solid")
lines(times, post.cirr, type="l", col="blue")

#legend(0.1,0.4, c(" yearly Atb"), lty=c(1,1), col=c("blue"), cex=0.8)

### To plot treatment
plot(times, post.treat, type="l", col="blue", main="Treated")
lines(times, post.treat, type="l", col="blue")
# ##########################################################################################
# #Spliting by baseline
# #Only alcoholic patiens
# cohort.alcohol.2 = cohort[which(bl[,"Alcohol"] == 2),]
# post.alcohol.2 <- posteriorProbabilities(cohort.alcohol.2, times)
#
# cohort.alcohol.1 = cohort[as.numeric(which(bl[,"Alcohol"] == 1)),]
# post.alcohol.1 <- posteriorProbabilities(cohort.alcohol.1, times)
#
# cohort.alcohol.0 = cohort[which(bl[,"Alcohol"] == 0),]
# post.alcohol.0 <- posteriorProbabilities(cohort.alcohol.0, times)
# #########################################################################################
#
# # Genotype
# cohort.g1 = cohort[which(bl[,"Genotype"] == 0),]
# post.g1 <- posteriorProbabilities(cohort.g1, times)
#
# cohort.g4 = cohort[which(bl[,"Genotype"] == 1),]
# post.g4 <- posteriorProbabilities(cohort.g4, times)
#
# cohort.g23 = cohort[which(bl[,"Genotype"] == 2 | bl[,"Genotype"] == 3 ),]
# post.g23 <- posteriorProbabilities(cohort.g23, times)
#
# #Age (agelow= <=20. agemed: >20 && <60, agehigh: >=60)
# cohort.agelow = cohort[which(bl[,"Age"] <= 40),]
# post.agelow <- posteriorProbabilities(cohort.agelow, times)
# cohort.agemed = cohort[which(bl[,"Age"] > 0 & bl[,"Age"] < 60 ),]
# post.agemed <- posteriorProbabilities(cohort.agemed, times)
#
# cohort.agehigh = cohort[which(bl[,"Age"] > 40),]
# post.agehigh <- posteriorProbabilities(cohort.agehigh, times)
# ##NOT TO FORGET EXAMPLE WITH COMBINATIONS OF BASELINE CHAR
# #If you want to plot some given states out a cohort "cohort.x"
# #subset of states you want to plot: for example a mixture between crue states and grouped
#
# ##################### Endpoints with Bl = Alcohol #################################3
# #HCC
# postHCC.alcohol.2= post.alcohol.2@probabilities[,37] + post.alcohol.2@probabilities[,38] + post.alcohol.2@probabilities[,39] + post.alcohol.2@probabilities[,40] + post.alcohol.2@probabilities[,41]+ post.alcohol.2@probabilities[,42]
# postHCC.alcohol.1= post.alcohol.1@probabilities[,37] + post.alcohol.1@probabilities[,38] + post.alcohol.1@probabilities[,39] + post.alcohol.1@probabilities[,40] + post.alcohol.1@probabilities[,41]+ post.alcohol.1@probabilities[,42]
# postHCC.alcohol.0= post.alcohol.0@probabilities[,37] + post.alcohol.0@probabilities[,38] + post.alcohol.0@probabilities[,39] + post.alcohol.0@probabilities[,40] + post.alcohol.0@probabilities[,41]+ post.alcohol.0@probabilities[,42]
#
# postHCC = post@probabilities[,37] + post@probabilities[,38] + post@probabilities[,39] + post@probabilities[,40] + post@probabilities[,41] + post@probabilities[,42]
# ###########################################################
# #DC
# postDC.alcohol.2= post.alcohol.2@probabilities[,31] + post.alcohol.2@probabilities[,32] + post.alcohol.2@probabilities[,33] + post.alcohol.2@probabilities[,34] + post.alcohol.2@probabilities[,35]+ post.alcohol.2@probabilities[,36]
# postDC.alcohol.1= post.alcohol.1@probabilities[,31] + post.alcohol.1@probabilities[,32] + post.alcohol.1@probabilities[,33] + post.alcohol.1@probabilities[,34] + post.alcohol.1@probabilities[,35]+ post.alcohol.1@probabilities[,36]
# postDC.alcohol.0 = post.alcohol.0@probabilities[,31] + post.alcohol.0@probabilities[,32] + post.alcohol.0@probabilities[,33] + post.alcohol.0@probabilities[,34] + post.alcohol.0@probabilities[,35]+ post.alcohol.0@probabilities[,36]
#
# postDC = post@probabilities[,31] + post@probabilities[,32] + post@probabilities[,33] + post@probabilities[,34] + post@probabilities[,35] + post@probabilities[,36]
# #############################################################################
# # esld
# post.esld.alcohol.2 = postHCC.alcohol.2 + postDC.alcohol.2
# post.esld.alcohol.1 = postHCC.alcohol.1 + postDC.alcohol.1
# post.esld.alcohol.0= postHCC.alcohol.0 + postDC.alcohol.0
# ###########################################################################
# #  CIRRHOSIS
# post.cirr.alcohol.1= apply(post.alcohol.1@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# post.cirr.alcohol.0 = apply(post.alcohol.0@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# postcirr = apply(post@probabilities[,c(5,10,15,20,25)], 1, sum)
# ###########################################################################
# ### Plot with Bl alcohol
# plot(times, post.cirr.alcohol.1, type="l", col="red", main="Cirrhosis")
# lines(times, post.cirr.alcohol.1, type="l", col="red")
# ##################### Endpoints with Bl = Genotype #################################3
# # HCC
# postHCC.g1= post.g1@probabilities[,37] + post.g1@probabilities[,38] + post.g1@probabilities[,39] +
#   post.g1@probabilities[,40] + post.g1@probabilities[,41] + post.g1@probabilities[,42]
# postHCC.g4 = post.g4@probabilities[,37] + post.g4@probabilities[,38] + post.g4@probabilities[,39] +
#   post.g4@probabilities[,40] + post.g4@probabilities[,41] + post.g4@probabilities[,42]
# postHCC.g23 = post.g23@probabilities[,37] + post.g23@probabilities[,38] + post.g23@probabilities[,39] +
#   post.g23@probabilities[,40] + post.g23@probabilities[,41] + post.g23@probabilities[,42]
#
# # #//// DC
# # postDC.g1= post.g1@probabilities[,27]
# # postDC.g4 = post.g4@probabilities[,27]
# # postDC.g4 = post.g23@probabilities[,27]
# # postDC = post@probabilities[,27]
# #
# # CIRRHOSIS
# postcirr.g1= apply(post.g1@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# postcirr.g4 = apply(post.g4@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# postcirr.g23 = apply(post.g23@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# postcirr = apply(post@probabilities[,c(5,10,15,20,25,30)], 1, sum)
#
# # ### Plot with Bl genotype
# # plot(times, postcirr.g4, type="l", col="red", main="Cirrhosis")
# # lines(times, postcirr.g4, type="l", col="red")
# # # add 2nd cohort
# # lines(times, postcirr.g1, type="l", col="blue")
# # lines(times, postcirr.g23, type="l", col="black")
# # legend(0.1,0.4, c("Genotype 1","Genotype 4", "Genotype 2/3"), lty=c(1,1), col=c("blue", "red","black"), cex=0.8)
# #
# #
# # ##################### Endpoints with Bl = Age at infection #################################3
# # #//// HCC
#
# # postHCC.agelow= post.agelow@probabilities[,27]
# # postHCC.agemed = post.agemed@probabilities[,27]
# # postHCC.agehigh = post.gagehigh@probabilities[,27]
# #
# # #//// DC
# # postDC.agelow= post.agelow@probabilities[,27]
# # postDC.agemed = post.agemed@probabilities[,27]
# # postDC.agehigh = post.agehigh@probabilities[,27]
# #
# #
# # CIRRHOSIS
# post.cirr.agelow= apply(post.agelow@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# post.cirr.agemed = apply(post.agemed@probabilities[,c(5,10,15,20,25,30)], 1, sum)
# post.cirr.agehigh = apply(post.agehigh@probabilities[,c(5,10,15,20,25,30)], 1, sum)
#
# ### Plot with Bl Age
# plot(times, post.cirr.agehigh, type="l", col="red", main="Cirrhosis based on the age at HCV infection")
# lines(times, post.cirr.agehigh, type="l", col="red")
# # add 2nd cohort
# lines(times, post.cirr.agelow, type="l", col="blue")
# lines(times, post.cirr.agemed, type="l", col="black")
# legend(0.1,0.6, c("Age at infection >= 40", "Age at infection < 40"), lty=c(1,1), col=c("red","blue"), cex=0.8)
