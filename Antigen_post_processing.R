#######################################################################################################################################
#
# Parameters and defining the dataframes and tabels
#
#
########################################################################################################################################
scenario = c("Baseline","S1","S2","S3","S4")
cat("\014")
# S1: The cost of PCR test was decreased by 20%.
# S2: The cost of antigen test was decreased by 30%.
# S3: Change the cost of liver disease after SVR to zero
# S4: Different set of liver disease costs (based on data from Iran).
# S5: Low prevalence among non-IDUs (looking at only non-IDUs cohort).
# S6: High prevalence among IDUs (looking at only IDUs cohort).

##########################################################################################################################################
current_dir = ("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources/Results/")
wd =  getwd()
for(scenario_ind in 1:5)
{
  cat("\014")
  
{
  if(scenario_ind == 1)
  {
    i = 1 # modify for the cost_stage that you are interested! --> Liver-related costs
    kkk = 1 # modify for COst_tests that you are interested! --> Laboratory costs
    cat("Baseline Scenario")
  }
  if(scenario_ind == 2)
  {
    i = 1 # modify for the cost_stage that you are interested! --> Liver-related costs
    kkk = 2 # modify for COst_tests that you are interested! --> Laboratory costs
    cat("S1: The cost of PCR test was decreased by 20%.")
    
  }
  if(scenario_ind == 3)
  {
    i = 1 # modify for the cost_stage that you are interested! --> Liver-related costs
    kkk = 3 # modify for COst_tests that you are interested! --> Laboratory costs
    cat("S2: The cost of antigen test was decreased by 30%.")
    
  }
  if(scenario_ind == 4)
  {
    i = 2 # modify for the cost_stage that you are interested! --> Liver-related costs
    kkk = 1 # modify for COst_tests that you are interested! --> Laboratory costs
    cat("S3: Change the cost of liver disease after SVR to zero")
    
  }
  if(scenario_ind == 5)
  {
    i = 9 # modify for the cost_stage that you are interested! --> Liver-related costs
    kkk = 1 # modify for COst_tests that you are interested! --> Laboratory costs
    cat("S4: Different set of liver disease costs (based on data from France).")
  }

  period = 4
  x =  0.02
  nFP = 0
 
  
  
  ###############################################################################################################################
  COst_tests <- matrix(0,4,10)
  colnames(COst_tests) <- c("baseline","atb-ag","ag","atbandag","atbandpcr","pcr","pre-treatment","blood_collection","antibody","treatment")
  COst_tests[1,] = c(40.26,22.87,20.87,22.87,40.26,38.26, 115,13,2,3248.33 * 12)
  COst_tests[2,] =  c(32.608,22.87,20.87,22.87,32.608,30.608, 115,13,2,3248.33 * 12) # S1:if I decreased PCR test by 20%: PCR: $30.608
  COst_tests[3,] =  c(40.26,16.609,14.609,16.609,40.26,38.26, 115,13,2,3248.33 * 12) # S2:if I decreased antigen test by 30%: antgen: $14.609

  ans_test_cost <- "No"

  
  PATH1 = setwd("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources")
  PATH2 = setwd("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources/Results/")

  #############################################################################################################################
  Cost_dollor = COst_tests[kkk,1:9]
  TreatCost = COst_tests[kkk,10] #3248.33 * 12# EUro https://onlinelibrary.wiley.com/doi/full/10.1111/liv.13240#support-information-section 
  
  cost_stage  <- matrix(0,11,13) 
  colnames(cost_stage) <- c("CF0","CF1","CF2","CF3","CF4","CDC","C_HCC","C_Death", "CF0_SVR","CF1_SVR","CF2_SVR","CF3_SVR","CF4_SVR")
  cost_stage[1,1:13] = c(447,447,447,447,578,1984,2474,0,447/2,447/2,447/2,447/2,578)
  cost_stage[2,1:13] = c(447,447,447,447,578,1984,2474,0,0,0,0,0,0)    
  
  cost_stage[3,1:13] = c(267, 267, 267, 267,627, 794, 961,0, 267/2, 267/2, 267/2, 267/2,627)
  cost_stage[4,1:13] = c(267, 267, 267, 267,627, 794, 961,0,0,0,0,0,627)
  
  cost_stage[5,1:13] = c(259, 259, 259, 259,1506,5013, 37091, 0, 259/2, 259/2,259/2,259/2, 1506)
  cost_stage[6,1:13] = c(259, 259, 259, 259,1506,5013, 37091, 0,0,0,0,0,0)
  
  cost_stage[7,1:13] = c(221,221,1152,1152,1830,14657,13062,0,221/2,221/2,1152/2,1152/2,1830)
  cost_stage[8,1:13] = c(221,221,1152,1152,1830,14657,13062,0,0,0,0,0,0)
  
  
  cost_stage[9,1:13] = c(391, 391, 391, 456, 1712, 9238, 13311,0, 391/2, 391/2, 391/2, 456/2, 1712)
  cost_stage[10,1:13] = c(391, 391, 391, 456, 1712, 9238, 13311,0,0,0,0,0,0)
  
  LT_cost = 62641.84

  utility_stages  <- matrix(0,1,13) 
  colnames(utility_stages) <- c("CF0","CF1","CF2","CF3","CF4","CDC","C_HCC","C_Death", "CF0_SVR","CF1_SVR","CF2_SVR","CF3_SVR","CF4_SVR")
  utility_stages[1,1:13] = c(0.82,0.82,0.82,0.76,0.76,0.6,0.6,0,0.95,0.95,0.85,0.85,0.85)

  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }

  
}
########################################################################################################################################
# Define the data set of interest by choosing the PATH
# 
# Run the following codes to get the alt_table for analyses

#######################################################################################################################################
all_cases = matrix(NA,1,5)
alt_table_all = matrix(NA,1,5)
for (gg in 1:2)
{
for (path_index in 1:1)
  {
    ans_IDUs = ifelse(gg == 1, "Yes","No")

  if (ans_IDUs == "Yes")
  {
    if (path_index == 1)
      PATH = "../Outputs/Group5_from_July28_IDUS_cluster/"
    
  }
  ########################################################################################################################################## 
  # for non-IDUs
  if (ans_IDUs == "No")
  {

    if (path_index == 1)
      PATH = "../Outputs/Group5_from_code_July29_17_cluster_NOIDUS/"
    
  }
  
  ###########################################################################################################################################
  Table_results = read.table(sprintf("%s/Table_results.csv",PATH),row.names = NULL)
  Table_results = Table_results[,2:ncol(Table_results)]
  Table_results = Table_results[,c(1:4,6)]
  ########################################################################################
  
  alt_table = matrix(NA,5,5)
  colnames(alt_table) <- c("atb->pcr", "atb->antg","antg","atb+antg","pcr")
  rownames(alt_table) <- c("liver disease costs","testing costs for positives (=simulated)","testing costs for negatives","treatment costs","QALYs")
  
  delta_cost_utility = matrix(NA, 20,4)
  colnames(delta_cost_utility) <- c("atb-pcr", "atb-antg","antg","pcr")
  rownames(delta_cost_utility) <- c("Cost","Utility","delta_Cost", "delta_Utility","ICER","Cost","Utility","delta_Cost", "delta_Utility","ICER",
                                    "Cost","Utility","delta_Cost", "delta_Utility","ICER", "Cost","Utility","delta_Cost", "delta_Utility","ICER")
  
  
  
  ###################################################################################
  #for (i = ?)
  {
    cost_stages = cost_stage[i,]
    
    ######################################################################################################################################################
   setwd(dir = PATH)
    # Cohort_comb = read.table("atb_pcr.csv",row.names = NULL)
    # Cohort_comb = Cohort_comb[,c(2:66)]
    # {
    #   calendartime_matrix <- matrix(NA, nrow = nrow(Cohort_comb), ncol = 100)
    #   Cohort_calendar <- Cohort_comb[,(1:52)] + Cohort_comb[,"Age"] + Cohort_comb[,"BirthYear"]
    #   Cohort_comb[,(1:52)] = Cohort_comb[,(1:52)] + Cohort_comb[,"Age"] + Cohort_comb[,"BirthYear"]
    #   indexmatrix <- matrix((1:52),nrow = nrow(Cohort_calendar), ncol = 52, byrow = TRUE)
    #   for (year in 1:100) 
    #   {
    #     calyear <- year + 1969
    #     temp <- (Cohort_calendar < calyear + 1) * indexmatrix
    #     calendartime_matrix[,year] <- apply(temp, 1, max, na.rm = TRUE)
    #   }
    #   calendartime_matrix2 = calendartime_matrix
    #   calendartime_matrix2[,1:45] = 0 # we are only interested to look at 2015 and afterward
    #   calendartime_matrix_main = calendartime_matrix
    #   calendartime_matrix = calendartime_matrix2
    #   
    #   calendartime_matrix = calendartime_matrix[which((Cohort_comb[,1] < 2018 | Cohort_comb[,2] < 2018 | Cohort_comb[,3] < 2018 | Cohort_comb[,4] < 2018) &
    #                                                     (Cohort_comb[,49] > 2016 | Cohort_comb[,50] > 2016 | Cohort_comb[,51] > 2016 | Cohort_comb[,52] > 2016)),]
    #   
    #   calendartime_matrix_main2 = calendartime_matrix_main[which((Cohort_comb[,1] < 2018 | Cohort_comb[,2] < 2018 |
    #                                                                 Cohort_comb[,3] < 2018 | Cohort_comb[,4] < 2018) &
    #                                                                (Cohort_comb[,49] > 2016 | Cohort_comb[,50] > 2016 |
    #                                                                   Cohort_comb[,51] > 2016 | Cohort_comb[,52] > 2016)),]
    #   
    #   Cohort_comb = Cohort_comb[which((Cohort_comb[,1] < 2018 | Cohort_comb[,2] < 2018 |
    #                                      Cohort_comb[,3] < 2018 | Cohort_comb[,4] < 2018) &
    #                                     (Cohort_comb[,49] > 2016 | Cohort_comb[,50] > 2016 | Cohort_comb[,51] > 2016 | Cohort_comb[,52] > 2016)),]
    # }
    setwd(dir = sprintf("%s/Scenario1_atb_pcr/",PATH))
    alpha = read.table("alpha.csv",row.names = NULL)
    alpha = alpha[,2:ncol(alpha)]
    beta = read.table("beta.csv",row.names = NULL)
    beta = beta[,2:ncol(beta)]
    beta1 = read.table("beta1.csv",row.names = NULL)
    beta1 = beta1[,2:ncol(beta1)]
    cost_stages_matrix = rep.row(cost_stages,nrow(alpha))
    utility_stages_matrix = rep.row(utility_stages,nrow(alpha))
    Liver_T = read.table("Liver_T.csv")

    
    # weight <- rep(0, length=nrow(Cohort_comb))
    # diag_age = rep(0, length=nrow(Cohort_comb))
    # 
    #   diag_age = ifelse(Cohort_comb[,"Diag.time"] <200, Cohort_comb[, "Diag.time"],0)
    #   weight = ifelse(diag_age >30 & diag_age<=59, 3, 1)
    # 
    # alpha = alpha * weight
    # beta = beta * weight
    # beta1 = beta1 * weight
    # Liver_T = Liver_T * weight

    n_liver = length(which(Liver_T >0)) * LT_cost
    treat_cost =  sum(beta1 * TreatCost)
    confirm_SVR_cost = (length(beta1) - length(which(beta1 == 0))) * Cost_dollor[6] # number of people who treated
    
    time_cost = sum(alpha * cost_stages_matrix)
    time_utility = sum(alpha * utility_stages_matrix)
    diag_cost = (beta * (Cost_dollor[1] + 2 * Cost_dollor[8])) # add the blood collection cost for antibody and pcr
    test_cost = (Cost_dollor[1] + 2 * Cost_dollor[8])
    
    
    TN = Table_results[4,1]
    nFN = Table_results[1,1]
    cohortSize = Table_results[10,1]
    
    treat_cost =  sum(beta1 * TreatCost)

    alt_table[1,1] = (time_cost +n_liver) /cohortSize # replace the column number 1 with that of the current strategy
    
    alt_table[2,1] = sum(diag_cost)/cohortSize
    
    alt_table[3,1] = (x*TN*test_cost + (1-x)*TN*(Cost_dollor[8]+Cost_dollor[9]))/cohortSize #note for strategies from 3 onwards this should be TN*test_cost/cohortSize
    
    alt_table[4,1] = (sum(treat_cost) + confirm_SVR_cost)/cohortSize
    
    alt_table[5,1] = time_utility/cohortSize
    
    ######################################################################################################################################
    ############## Atb --> AG
    ######################################################################################################################################
    
    setwd(dir = sprintf("%s/Scenario2_atb_ag/",PATH))
    alpha = read.table("alpha.csv",row.names = NULL)
    alpha = alpha[,2:ncol(alpha)]
    beta = read.table("beta.csv",row.names = NULL)
    beta = beta[,2:ncol(beta)]
    beta1 = read.table("beta1.csv",row.names = NULL)
    beta1 = beta1[,2:ncol(beta1)]
    cost_stages_matrix = rep.row(cost_stages,nrow(alpha))
    utility_stages_matrix = rep.row(utility_stages,nrow(alpha))
    Liver_T = read.table("Liver_T.csv")

    n_liver = length(which(Liver_T >0)) * LT_cost
    
    treat_cost =  sum(beta1 * TreatCost)
    confirm_SVR_cost = (length(beta1) - length(which(beta1 == 0))) * Cost_dollor[3] 
    
    time_cost = sum(alpha * cost_stages_matrix)
    time_utility = sum(alpha * utility_stages_matrix)
    
    test_cost = (Cost_dollor[2] + 2 * Cost_dollor[8])
    diag_cost = (beta * (Cost_dollor[2] + 2 * Cost_dollor[8])) # add the blood collection cost for antibody and pcr
    
    nTN_atb_pos = 12427 * cohortSize / 13456 # the number in the cohort
    
    
    TN = Table_results[4,2]
    nFN = Table_results[1,2]
    cohortSize = Table_results[10,2]


    
    alt_table[1,2] = (time_cost +n_liver)/cohortSize # replace the column number 1 with that of the current strategy
    
    alt_table[2,2] = sum(diag_cost)/cohortSize
    
    alt_table[3,2] = (x*TN*test_cost + (1-x)*TN*(Cost_dollor[8]+Cost_dollor[9]))/cohortSize #note for strategies from 3 onwards this should be TN*test_cost/cohortSize
    
    alt_table[4,2] = (sum(treat_cost) + confirm_SVR_cost)/cohortSize
    
    alt_table[5,2] = time_utility/cohortSize
    
    
    
    ###############################################################################################################################################
    ##### Antigen
    ################################################################################################################################################
    setwd(dir = sprintf("%s/Scenario3_ag/",PATH))
    alpha = read.table("alpha.csv",row.names = NULL)
    alpha = alpha[,2:ncol(alpha)]
    beta = read.table("beta.csv",row.names = NULL)
    beta = beta[,2:ncol(beta)]
    beta1 = read.table("beta1.csv",row.names = NULL)
    beta1 = beta1[,2:ncol(beta1)]
    cost_stages_matrix = rep.row(cost_stages,nrow(alpha))
    utility_stages_matrix = rep.row(utility_stages,nrow(alpha))
    
    Liver_T = read.table("Liver_T.csv")

    n_liver = length(which(Liver_T >0)) * LT_cost
    
    treat_cost =  sum(beta1 * TreatCost)
    confirm_SVR_cost = (length(beta1) - length(which(beta1 == 0))) * Cost_dollor[3] 
    
    time_cost = sum(alpha * cost_stages_matrix)
    time_utility = sum(alpha * utility_stages_matrix)
    diag_cost = (beta * (Cost_dollor[3] + 1 * Cost_dollor[8])) # add the blood collection cost for antibody and pcr
    test_cost = (Cost_dollor[3] + 1 * Cost_dollor[8])
    treat_cost =  sum(beta1 * TreatCost)
    nTN = 12414 * cohortSize / 13456 # the number in the cohort
    
    TN = Table_results[4,3]
    nFN = Table_results[1,3]
    cohortSize = Table_results[10,3]
    ###########################################################################################  

    time_utility_test = time_utility  
    
    ###########################################################################################################################################

    
    alt_table[1,3] = (time_cost +n_liver)/cohortSize # replace the column number 1 with that of the current strategy
    
    alt_table[2,3] = sum(diag_cost)/cohortSize
    
    alt_table[3,3] = (TN*test_cost)/cohortSize #note for strategies from 3 onwards this should be TN*test_cost/cohortSize
    
    alt_table[4,3] = (sum(treat_cost) + confirm_SVR_cost)/cohortSize
    
    alt_table[5,3] = time_utility/cohortSize
    
    ###############################################################################################################################################
    ## Antibody and Antigen
    #############################################
 #    setwd(sprintf("%s/Scenario4_antibodyANDantigen/",PATH))
 #    alpha = read.table("alpha.csv",row.names = NULL)
 #    alpha = alpha[,2:ncol(alpha)]
 #    beta = read.table("beta.csv",row.names = NULL)
 #    beta = beta[,2:ncol(beta)]
 #    beta1 = read.table("beta1.csv",row.names = NULL)
 #    beta1 = beta1[,2:ncol(beta1)]
 #    cost_stages_matrix = rep.row(cost_stages,nrow(alpha))
 #    utility_stages_matrix = rep.row(utility_stages,nrow(alpha))
 #    
 #    Liver_T = read.table("Liver_T.csv")
 #    # tmp1 = count(Liver_T)
 #    # n_liver = (sum(tmp1$freq) - tmp1$freq[1]) * LT_cost
 # #   n_liver = 0
 #    n_liver = length(which(Liver_T >0)) * LT_cost
 #    
 #    treat_cost =  sum(beta1 * TreatCost)
 #    confirm_SVR_cost = (length(beta1) - length(which(beta1 == 0))) * Cost_dollor[3] 
 #    
 #    time_cost = sum(alpha * cost_stages_matrix)
 #    time_utility = sum(alpha * utility_stages_matrix)
 #    diag_cost = (beta * (Cost_dollor[4] + 1 * Cost_dollor[8])) # add the blood collection cost for antibody and pcr
 #    test_cost = (Cost_dollor[4] + 1 * Cost_dollor[8])
 #    treat_cost =  sum(beta1 * TreatCost)
 #    
 #    
 #    nTN = 12414 * cohortSize / 13456 # the number in the cohort
 #    
 #    
 #    TN = Table_results[4,4]
 #    nFN = Table_results[1,4]
 #    cohortSize = Table_results[10,4]
 #    
 #    # 
 #    # time_cost_test = time_cost +
 #    #   treat_cost + 
 #    #   sum(diag_cost) + TN * test_cost +
 #    #   nFP * Cost_dollor[7] + nFN * test_cost + confirm_SVR_cost + n_liver
 #    time_utility_test = time_utility   
 # 
 #    
 #    alt_table[1,4] = (time_cost +n_liver)/cohortSize # replace the column number 1 with that of the current strategy
 #    
 #    alt_table[2,4] = sum(diag_cost)/cohortSize
 #    
 #    alt_table[3,4] = (TN*test_cost)/cohortSize #note for strategies from 3 onwards this should be TN*test_cost/cohortSize
 #    
 #    alt_table[4,4] = (sum(treat_cost) + confirm_SVR_cost)/cohortSize
 #    
 #    alt_table[5,4] = time_utility/cohortSize
 #    
    ##############################################################################################################################################
    ## PCR test
    ##############################################################################################################################################
    setwd(dir = sprintf("%s/Scenario5_pcr/",PATH))
    alpha = read.table("alpha.csv",row.names = NULL)
    alpha = alpha[,2:ncol(alpha)]
    beta = read.table("beta.csv",row.names = NULL)
    beta = beta[,2:ncol(beta)]
    beta1 = read.table("beta1.csv",row.names = NULL)
    beta1 = beta1[,2:ncol(beta1)]
    cost_stages_matrix = rep.row(cost_stages,nrow(alpha))
    utility_stages_matrix = rep.row(utility_stages,nrow(alpha))
    
    Liver_T = read.table("Liver_T.csv")
 
    n_liver = length(which(Liver_T >0)) * LT_cost
    
    time_cost = sum(alpha * cost_stages_matrix)
    time_utility = sum(alpha * utility_stages_matrix)
    treat_cost =  sum(beta1 * TreatCost)
    confirm_SVR_cost = (length(beta1) - length(which(beta1 == 0))) * Cost_dollor[6] 
    
    TN = Table_results[4,5]
    nFN = Table_results[1,5]
    cohortSize = Table_results[10,5]
    
    diag_cost = (beta * (Cost_dollor[6] + 1 * Cost_dollor[8])) # add the blood collection cost for antibody and pcr
    
    test_cost = (Cost_dollor[6] + 1 * Cost_dollor[8])
    treat_cost =  sum(beta1 * TreatCost)
    
    
    nTN = 12414 * cohortSize / 13456 # the number in the cohort
    
    cohortSize = 30761
    TN = Table_results[4,5]
    nFN = Table_results[1,5]
    cohortSize = Table_results[10,5]
    ############################################################################################################  

    time_utility_test = time_utility  

    
    alt_table[1,5] = (time_cost +n_liver)/cohortSize # replace the column number 1 with that of the current strategy
    
    alt_table[2,5] = sum(diag_cost)/cohortSize
    
    alt_table[3,5] = (TN*test_cost)/cohortSize #note for strategies from 3 onwards this should be TN*test_cost/cohortSize
    
    alt_table[4,5] = (sum(treat_cost) + confirm_SVR_cost)/cohortSize
    
    alt_table[5,5] = time_utility/cohortSize
    
    
  }
  write.table(alt_table, sprintf("alt_table_%s.csv", path_index))
  if (path_index == 6)
    all_alt_table = alt_table
  
  all_cases = rbind(as.matrix(all_cases),as.matrix(read.table(sprintf("alt_table_%s.csv",path_index))))

  
}
alt_table_sd = alt_table
# if (gg == 1)
# {
#   for (h in 1:5)
#     for (g in 1:5)
#     {
#       alt_table[g,h]  = mean(all_cases[c((g+1),(g+6),(g+11),(g+16),(g+21)),h])
#       alt_table_sd[g,h] = sd(all_cases[c((g+1),(g+6),(g+11), (g+16),(g+21)),h])
#     }
# }
# 
#   if (gg == 2)
#   {
#     for (h in 1:5)
#       for (g in 1:5)
#       {
#         alt_table[g,h]  = mean(all_cases[c((g+1),(g+6)),h])
#         alt_table_sd[g,h] = sd(all_cases[c((g+1),(g+6)),h])
#       }
#  }

alt_table_all = rbind(alt_table_all,alt_table) # the first alt_table will be for IDUs, the second for non-IDUS
}
alt_table_all = alt_table_all[c(2:11),]
w1 = 0.25
w2 = 1 -w1
alt_table = alt_table_all[c(1:5),] * w1 + alt_table_all[c(6:10),] * w2
#ans_results  <- readline("Do you want to work with the whole simulated cohort (1) or the average value of the 5 cohorts (2)? (1/2)")
ans_results = "2"
if (ans_results == "1")
  alt_table = all_alt_table
########################################################################
alt_table = alt_table[,c(1:3,5)]
##############################################################
##### Analyse 1: Testing cost
#############################################################################################################
minpositive = function(x) min(x[x > 0])
{
  Testing = alt_table[c(2,5),] +alt_table[c(3,5),]
  Testing[2,] = alt_table[5,]
  ###########
  #1. sort the costs
  ############
  Testing2 = Testing
  a = sort(Testing[1,])
  Testing2[2,1:4] = Testing[2,names(a)]
  Testing2[1,1:4] = Testing[1,names(a)] # the cost is sorted
  colnames(Testing2)[1:4] = names(a)
  colnames(delta_cost_utility)[1:4] = names(a)
  #######################################################
  best_cost = Testing2[1,1]
  best_utility = Testing2[2,1]
  best_ICER = best_cost/best_utility
  Testing2[2, which(Testing2[1,] == best_cost)] = NA
  Testing2[1, which(Testing2[1,] == best_cost)] = NA
  
  for(hh in 0:3)
  {
    for(kk in 1: 4)
    {
      
      delta_cost_utility[(1 + (5 * hh)),kk] = Testing2[1,kk]
      delta_cost_utility[(2 + (5 * hh)),kk] = Testing2[2,kk]
      delta_cost_utility[(3 + (5 * hh)),kk]  = Testing2[1,kk] - best_cost
      delta_cost_utility[(4 + (5 * hh)),kk]= Testing2[2,kk] - best_utility
      delta_cost_utility[(5 + (5 * hh)),kk] = delta_cost_utility[(3 + (5 * hh)),kk] / delta_cost_utility[(4 + (5 * hh)),kk]
      
    }
    
    
    best_ICER = minpositive(delta_cost_utility[(5 + (5 * hh)),which(!is.na(delta_cost_utility[(5 + (5 * hh)),]))])
    next_best = delta_cost_utility[(1 + (5 * hh)),which(delta_cost_utility[(5 + (5 * hh)),] == best_ICER)]
    best_cost = Testing2[1, which(Testing2[1,] == next_best)]
    best_utility  = Testing2[2, which(Testing2[1,] == next_best)]
    Testing2[2, which(Testing2[1,] == next_best)] = NA
    Testing2[1, which(Testing2[1,] == next_best)] = NA
    
  }
}
#############################################################################################################
##### Analyse 2: Costs of testing and liver disease
#############################################################################################################
{
  delta_cost_utility_liver = matrix(NA, 20,4)
  colnames(delta_cost_utility_liver) <- c("atb-pcr", "atb-antg","antg","pcr")
  rownames(delta_cost_utility_liver) <- c("Cost","Utility","delta_Cost", "delta_Utility","ICER","Cost","Utility","delta_Cost", "delta_Utility","ICER",
                                          "Cost","Utility","delta_Cost", "delta_Utility","ICER", "Cost","Utility","delta_Cost", "delta_Utility","ICER")
  #############################################################################################################
  minpositive = function(x) min(x[x > 0])
  Testing = alt_table[c(2,5),] + alt_table[c(3,5),] + alt_table[c(1,5),]
  Testing[2,] = alt_table[5,]
  ###########
  #1. sort the costs
  ############
  Testing2 = Testing
  a = sort(Testing[1,])
  Testing2[2,1:4] = Testing[2,names(a)]
  Testing2[1,1:4] = Testing[1,names(a)] # the cost is sorted
  colnames(Testing2)[1:4] = names(a)
  colnames(delta_cost_utility_liver)[1:4] = names(a)
  #######################################################
  best_cost = Testing2[1,1]
  best_utility = Testing2[2,1]
  best_ICER = best_cost/best_utility
  Testing2[2, which(Testing2[1,] == best_cost)] = NA
  Testing2[1, which(Testing2[1,] == best_cost)] = NA
  
  for(hh in 0:3)
  {
    for(kk in 1: 4)
    {
      if(length(best_cost))
      {
        delta_cost_utility_liver[(1 + (5 * hh)),kk] = Testing2[1,kk]
        delta_cost_utility_liver[(2 + (5 * hh)),kk] = Testing2[2,kk]
        delta_cost_utility_liver[(3 + (5 * hh)),kk]  = Testing2[1,kk] - best_cost
        delta_cost_utility_liver[(4 + (5 * hh)),kk]= Testing2[2,kk] - best_utility
        delta_cost_utility_liver[(5 + (5 * hh)),kk] = delta_cost_utility_liver[(3 + (5 * hh)),kk] / delta_cost_utility_liver[(4 + (5 * hh)),kk]
      }
    }
    
    
    best_ICER = minpositive(delta_cost_utility_liver[(5 + (5 * hh)),which(!is.na(delta_cost_utility_liver[(5 + (5 * hh)),]))])
    next_best = delta_cost_utility_liver[(1 + (5 * hh)),which(delta_cost_utility_liver[(5 + (5 * hh)),] == best_ICER)]
    best_cost = Testing2[1, which(Testing2[1,] == next_best)]
    best_utility  = Testing2[2, which(Testing2[1,] == next_best)]
    Testing2[2, which(Testing2[1,] == next_best)] = NA
    Testing2[1, which(Testing2[1,] == next_best)] = NA
    
  }
  
}
#############################################################################################################
##### Analyse 3: All costs
#############################################################################################################
{
  delta_cost_utility_all = matrix(NA, 20,4)
  colnames(delta_cost_utility_all) <- c("atb-pcr", "atb-antg","antg","pcr")
  rownames(delta_cost_utility_all) <- c("Cost","Utility","delta_Cost", "delta_Utility","ICER","Cost","Utility","delta_Cost", "delta_Utility","ICER",
                                        "Cost","Utility","delta_Cost", "delta_Utility","ICER", "Cost","Utility","delta_Cost", "delta_Utility","ICER")
  #############################################################################################################
  minpositive = function(x) min(x[x > 0])
  Testing = alt_table[c(2,5),] + alt_table[c(3,5),] + alt_table[c(1,5),] + alt_table[c(4,5),] 
  Testing[2,] = alt_table[5,]
  ###########
  #1. sort the costs
  ############
  Testing2 = Testing
  a = sort(Testing[1,])
  Testing2[2,1:4] = Testing[2,names(a)]
  Testing2[1,1:4] = Testing[1,names(a)] # the cost is sorted
  colnames(Testing2)[1:4] = names(a)
  colnames(delta_cost_utility_all)[1:4] = names(a)
  #######################################################
  best_cost = Testing2[1,1]
  best_utility = Testing2[2,1]
  best_ICER = best_cost/best_utility
  Testing2[2, which(Testing2[1,] == best_cost)] = NA
  Testing2[1, which(Testing2[1,] == best_cost)] = NA
  
  for(hh in 0:3)
  {
    for(kk in 1: 4)
    {
      if(length(best_cost))
      {
        delta_cost_utility_all[(1 + (5 * hh)),kk] = Testing2[1,kk]
        delta_cost_utility_all[(2 + (5 * hh)),kk] = Testing2[2,kk]
        delta_cost_utility_all[(3 + (5 * hh)),kk]  = Testing2[1,kk] - best_cost
        delta_cost_utility_all[(4 + (5 * hh)),kk]= Testing2[2,kk] - best_utility
        delta_cost_utility_all[(5 + (5 * hh)),kk] = delta_cost_utility_all[(3 + (5 * hh)),kk] / delta_cost_utility_all[(4 + (5 * hh)),kk]
      }
    }
    
    
    best_ICER = minpositive(delta_cost_utility_all[(5 + (5 * hh)),which(!is.na(delta_cost_utility_all[(5 + (5 * hh)),]))])
    next_best = delta_cost_utility_all[(1 + (5 * hh)),which(delta_cost_utility_all[(5 + (5 * hh)),] == best_ICER)]
    best_cost = Testing2[1, which(Testing2[1,] == next_best)]
    best_utility  = Testing2[2, which(Testing2[1,] == next_best)]
    Testing2[2, which(Testing2[1,] == next_best)] = NA
    Testing2[1, which(Testing2[1,] == next_best)] = NA
    
  }
}

cat("\014") 

print("Analysis 1:\n")

print(delta_cost_utility)
cat("*********************************************************************************************************************")

print(alt_table[c(2,5),] +alt_table[c(3,5),])
cat("*********************************************************************************************************************")

print("Analysis 2:\n")

print(delta_cost_utility_liver)
print(alt_table[c(2,5),] + alt_table[c(3,5),] + alt_table[c(1,5),])
print("Analysis 3:\n")

print(delta_cost_utility_all)
cat("*********************************************************************************************************************")
print(alt_table[c(2,5),] + alt_table[c(3,5),] + alt_table[c(1,5),] + alt_table[c(4,5),])
setwd(current_dir)
today <- Sys.Date()
dd = format(today, format = "%B%d%Y")
subDir <- sprintf("Outcome_%s_%s",dd,scenario[scenario_ind])
PATH2 = wd
dir.create(file.path(PATH2, subDir), showWarnings = FALSE)
wd = setwd(file.path(PATH2, subDir))
#write.table(delta_cost_utility,"delta_cost_utility.csv")
setwd(file.path(PATH2, subDir))
write.table(delta_cost_utility, sprintf("delta_cost_utility_%s_%s.csv",dd, scenario[scenario_ind]))
write.table(delta_cost_utility_all, sprintf("delta_cost_utility_all_%s_%s.csv",dd, scenario[scenario_ind]))
write.table(delta_cost_utility_liver, sprintf("delta_cost_utility_liver_%s_%s.csv",dd, scenario[scenario_ind]))
write.table(alt_table, sprintf("alt_table_%s_%s.csv",dd, scenario[scenario_ind]))
time <- readline("Do you need time to check (enter the time)")
Sys.sleep(time)


}

############################################################################################################################################
setwd(dir = "S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources/Codes/")
###                                                                   Part 2

### Analysing the simulated data and calculate the proportion of people who live with different liver stages.
### For this part you need to have access to the simulated cohorts for each strategy and you need to modify Cohort_comb
##############################################################################################################################################
#continue <- readline(prompt="Part 2: Do you want to continue part 2 (you should have the simulated cohort data in order to continue this part)? (Yes/No) ")
continue = "NO"
require(MASS)
require(msm)
require(mstate)
require(splines)
require(plyr)
require(pracma)
library(pracma)
require(compare)
library(compare)
require(RColorBrewer)
require(gems)
#rm(list =  ls())
if (continue == "Yes")

{
  all_outcome = matrix(NA,1,11)
  colnames(all_outcome)<-c("F0","F1","F2","F3","F4", "DC", "HCC", "death", "liv.death.tot", "liv.death.noHCV","LT")
  for (gg in 1:2)
  {
    if (gg == 1)
      PATH = "S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources/Outputs/Group5_from_July28_IDUS_cluster"
    if (gg == 2)
      PATH = "S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources/Outputs/Group5_from_code_July29_17_cluster_NOIDUS/"
    
  

  for(test_index in 1:4)
  {
    if(test_index == 1)
      Cohort_comb = read.table(sprintf("%s/antg.csv",PATH), row.names = NULL)
    if(test_index == 2)
      Cohort_comb = read.table(sprintf("%s/atb_antg.csv",PATH), row.names = NULL)
    if(test_index == 3)
      Cohort_comb = read.table(sprintf("%s/pcr.csv",PATH), row.names = NULL)
    if(test_index == 4)
      Cohort_comb = read.table(sprintf("%s/atb_pcr.csv",PATH), row.names = NULL)
    Cohort_comb = Cohort_comb[,2:66]
  
  #########################################################################################################

  #########################################################################################################
  ## Model simulations
  diag_scenario = c("Baseline","Birth","IDU", "ex_IDU")
  diag_ind = 1
  ###############################################################################
  {
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
}
  #Cohort_comb = all_cases
  
  ##########################################################################################################
  ################################################################################################################
  ###### Post processing
  
  
  ################################################################################################################
  wd =  getwd()
  setwd(wd)
  
  
  n_endpoints=6
  ns=52
  bl_number = 16
  number_blocks=6
  fs = 5
  maxTime = 100
  
  cohort = new("ArtCohort")
  
  pp0 = Cohort_comb
  pp = list(pp0)
  est<-matrix(NA,length(pp),n_endpoints)
  lower<-matrix(NA,length(pp),n_endpoints)
  upper<-matrix(NA,length(pp),n_endpoints)
  
  #test = as.matrix(read.table("treatF2_test"))
  
  for ( k in 1: length(pp))
  {
    
    cohort@time.to.state = as.data.frame(pp[[k]][,1:ns])
################################################################################################################################
# Cirrhosis
################################################################################################################################
    
    cirr<-cbind(cohort@time.to.state[,5], cohort@time.to.state[,10], cohort@time.to.state[,15],
                cohort@time.to.state[,20], cohort@time.to.state[,25], cohort@time.to.state[,30])
    
    cirr.inc<-cirr
    
    for(j in 2:6)
      for(i in j:6)
        cirr.inc[,i][cirr.inc[,j-1]!="NA"]<-NA
    
   est.cirr =  sum(!is.na(cirr.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(cirr.inc[,2]))/nrow(Cohort_comb) * 100 +
      sum(!is.na(cirr.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(cirr.inc[,4]))/nrow(Cohort_comb) * 100 + 
      sum(!is.na(cirr.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(cirr.inc[,6]))/nrow(Cohort_comb) * 100 
##################################################################################################################################    
   dc<-cbind(cohort@time.to.state[,31], cohort@time.to.state[,32], cohort@time.to.state[,33],
               cohort@time.to.state[,34], cohort@time.to.state[,35], cohort@time.to.state[,36])
   
   dc.inc<-dc
   
   for(j in 2:6)
     for(i in j:6)
       dc.inc[,i][dc.inc[,j-1]!="NA"]<-NA
   
   est.dc =  sum(!is.na(dc.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(dc.inc[,2]))/nrow(Cohort_comb) * 100 +
     sum(!is.na(dc.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(dc.inc[,4]))/nrow(Cohort_comb) * 100 + 
     sum(!is.na(dc.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(dc.inc[,6]))/nrow(Cohort_comb) * 100 
################################################################################################################################
   ##################################################################################################################################    
   hcc<-cbind(cohort@time.to.state[,37], cohort@time.to.state[,38], cohort@time.to.state[,39],
             cohort@time.to.state[,40], cohort@time.to.state[,41], cohort@time.to.state[,42])
   
   hcc.inc<-hcc
   
   for(j in 2:6)
     for(i in j:6)
       hcc.inc[,i][hcc.inc[,j-1]!="NA"]<-NA
   
   est.hcc =  sum(!is.na(hcc.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(hcc.inc[,2]))/nrow(Cohort_comb) * 100 +
     sum(!is.na(hcc.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(hcc.inc[,4]))/nrow(Cohort_comb) * 100 + 
     sum(!is.na(hcc.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(hcc.inc[,6]))/nrow(Cohort_comb) * 100 

    
   ##################################################################################################################################    
   lt<-cbind(cohort@time.to.state[,43], cohort@time.to.state[,44], cohort@time.to.state[,45],
              cohort@time.to.state[,46], cohort@time.to.state[,47], cohort@time.to.state[,48])
   
   lt.inc<-lt
   
   for(j in 2:6)
     for(i in j:6)
       lt.inc[,i][lt.inc[,j-1]!="NA"]<-NA
   
   est.lt =  sum(!is.na(lt.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(lt.inc[,2]))/nrow(Cohort_comb) * 100 +
     sum(!is.na(lt.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(lt.inc[,4]))/nrow(Cohort_comb) * 100 + 
     sum(!is.na(lt.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(lt.inc[,6]))/nrow(Cohort_comb) * 100 
   
 ###################################################################################################################################### 
# Death   
##################################################################################################################################    
   death<-cbind(cohort@time.to.state[,49], cohort@time.to.state[,50], cohort@time.to.state[,51],
              cohort@time.to.state[,52])
   
   death.inc<-death
   
   est.death =  sum(!is.na(death.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(death.inc[,2]))/nrow(Cohort_comb) * 100 +
     sum(!is.na(death.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(death.inc[,4]))/nrow(Cohort_comb) * 100 
   

   
#################################################################################################################################
    ### liver related deaths
################################################################################################################################
# cohort.ldeath is a subcohort where state 28 is only liver realted death and not total deaths
    
    cohort.ldeath<-cohort
    cohort.ldeath@time.to.state[is.na(cohort.ldeath@time.to.state)]<-666
    for (i in 1: nrow(cohort.ldeath@time.to.state))
    {
      cohort.ldeath@time.to.state[i,52][cohort.ldeath@time.to.state[i,31] == 666 &&
                                          cohort.ldeath@time.to.state[i,32] == 666&&
                                          cohort.ldeath@time.to.state[i,33] == 666&&
                                          cohort.ldeath@time.to.state[i,34] == 666&&
                                          cohort.ldeath@time.to.state[i,35] == 666&&
                                          cohort.ldeath@time.to.state[i,36] == 666&&
                                          cohort.ldeath@time.to.state[i,37] == 666&&
                                          cohort.ldeath@time.to.state[i,38] == 666&&
                                          cohort.ldeath@time.to.state[i,39] == 666&&
                                          cohort.ldeath@time.to.state[i,40] == 666&&
                                          cohort.ldeath@time.to.state[i,41] == 666&&
                                          cohort.ldeath@time.to.state[i,42] == 666]<-NA
    }
    

    est.ldeath <- sum(!is.na(cohort.ldeath@time.to.state[,52]))/nrow(Cohort_comb) * 100

    ## liver related deaths happening without replicating HCV (After HCV clearance, in block undetectable HCV)
    ## cohort.ldeath.sub is a cohort where column 28 (death) includes only liver related death occuring in patients without detectable HCV
    
    cohort.ldeath.sub<-cohort.ldeath
    
    for (i in 1: nrow(cohort.ldeath.sub@time.to.state))
    {
      cohort.ldeath.sub@time.to.state[i,52][cohort.ldeath.sub@time.to.state[i,30] == 666 ||
                                              cohort.ldeath.sub@time.to.state[i,25] == 666 ||
                                              cohort.ldeath.sub@time.to.state[i,20] == 666 ||
                                              cohort.ldeath.sub@time.to.state[i,15] == 666  ]<-NA
    }
    

    est.ldeath.sub<-sum(!is.na(cohort.ldeath.sub@time.to.state[,52]))/nrow(Cohort_comb) * 100

###################################################################################################################################    
    
    est[k,]<-c(est.cirr, est.dc, est.hcc, est.death, est.ldeath, est.ldeath.sub)
    colnames(est)<-c("F4", "DC", "HCC", "death", "liv.death.tot", "liv.death.noHCV")
    rownames(est)<-c("treatFg")

  }
  
  #write.table(est, "est_bc")
#################################################################################################################################
  
  # cohort.cirr is a sub cohort that includes only times for incident cirrhosis and not those who had cirrhosis in an earlier block
  cohort.F3.inc<-cohort
  F3<-cbind(cohort@time.to.state[,4], cohort@time.to.state[,9], cohort@time.to.state[,14],
            cohort@time.to.state[,19], cohort@time.to.state[,24], cohort@time.to.state[,29])
  
  F3.inc<-F3
  
  for(j in 2:6)
    for(i in j:6)
      F3.inc[,i][F3.inc[,j-1]!="NA"]<-NA

  
  est.cirrF3 <- sum(!is.na(F3.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(F3.inc[,2]))/nrow(Cohort_comb) * 100 +
    sum(!is.na(F3.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(F3.inc[,4]))/nrow(Cohort_comb) * 100 + 
    sum(!is.na(F3.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(F3.inc[,6]))/nrow(Cohort_comb) * 100 
  
  
  cohort.F2.inc<-cohort
  F2<-cbind(cohort@time.to.state[,3], cohort@time.to.state[,8], cohort@time.to.state[,13],
            cohort@time.to.state[,18], cohort@time.to.state[,23], cohort@time.to.state[,28])
  
  F2.inc<-F2
  
  for(j in 2:6)
    for(i in j:6)
      F2.inc[,i][F2.inc[,j-1]!="NA"]<-NA
  
  est.cirrF2 <- sum(!is.na(F2.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(F2.inc[,2]))/nrow(Cohort_comb) * 100 +
    sum(!is.na(F2.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(F2.inc[,4]))/nrow(Cohort_comb) * 100 + 
    sum(!is.na(F2.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(F2.inc[,6]))/nrow(Cohort_comb) * 100 
  
  
  cohort.F1.inc<-cohort
  F1<-cbind(cohort@time.to.state[,2], cohort@time.to.state[,7], cohort@time.to.state[,12],
            cohort@time.to.state[,17], cohort@time.to.state[,22], cohort@time.to.state[,27])
  
  F1.inc<-F1
  
  for(j in 2:6)
    for(i in j:6)
      F1.inc[,i][F1.inc[,j-1]!="NA"]<-NA
  
  est.cirrF1 <- sum(!is.na(F1.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(F1.inc[,2]))/nrow(Cohort_comb) * 100 +
    sum(!is.na(F1.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(F1.inc[,4]))/nrow(Cohort_comb) * 100 + 
    sum(!is.na(F1.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(F1.inc[,6]))/nrow(Cohort_comb) * 100 
  
  
  cohort.F0.inc<-cohort
  F0<-cbind(cohort@time.to.state[,1], cohort@time.to.state[,6], cohort@time.to.state[,11],
            cohort@time.to.state[,16], cohort@time.to.state[,11], cohort@time.to.state[,16])
  
  F0.inc<-F0
  
  for(j in 1:6)
    for(i in j:6)
      F0.inc[,i][F0.inc[,j-1]!="NA"]<-NA
  
  est.cirrF0 <- sum(!is.na(F0.inc[,1]))/nrow(Cohort_comb) * 100 + sum(!is.na(F0.inc[,2]))/nrow(Cohort_comb) * 100 +
    sum(!is.na(F0.inc[,3]))/nrow(Cohort_comb) * 100 + sum(!is.na(F0.inc[,4]))/nrow(Cohort_comb) * 100 + 
    sum(!is.na(F0.inc[,5]))/nrow(Cohort_comb) * 100 + sum(!is.na(F0.inc[,6]))/nrow(Cohort_comb) * 100 
  
  outcome = c(est.cirrF0,est.cirrF1, est.cirrF2, est.cirrF3, est,est.LT)
  all_outcome = rbind(all_outcome, outcome)
  
  #write.table(outcome, sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type))
  
}
 
  }
  rownames(all_outcome)<-c("NA","Antigen","Atb-Antigen","PCR","Atb-PCR","Antigen","Atb-Antigen","PCR","Atb-PCR")
  w1 = 0.25
  w2 = 1 -w1
  fibrosis_results = all_outcome[c(2:5),] * w1 + all_outcome[c(6:9),] * w2
  print(fibrosis_results)
}
Table_results1 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources_codes/Outputs/Group5_from_code_July29_17_cluster_NOIDUS/Table_results.csv", row.names = NULL)
Table_results2 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_paper_sources_codes/Outputs/Group5_from_July28_IDUS_cluster/Table_results.csv", row.names = NULL)
(0.25* (Table_results2[7,2:7]/Table_results2[10,2:7]) + 0.75*(Table_results1[7,2:7]/Table_results1[10,2:7]))*100000
