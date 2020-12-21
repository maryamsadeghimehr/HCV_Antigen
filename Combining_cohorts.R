z1 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Data/Antigen_April17_results_general_population/NO_outputDirectoryFg_Georgia_April232019_atb+antg/treatFg_baseline1_April232019",row.names = NULL)
z1 <- z1[(1:nrow(z1)),(2:ncol(z1))]
z2 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Data/Antigen_April17_results_general_population/NO_outputDirectoryFg_Georgia_April232019_atb+antg1/treatFg_baseline1_April232019",row.names = NULL)
z2 <- z2[(1:nrow(z2)),(2:ncol(z2))]
z3=rbind(z1,z2)

write.table(z3, "treatFg_baseline1_April232019_1")


z1 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Data/Antigen_April17_results_general_population/NO_outputDirectoryFg_Georgia_April172019_atb-pcr/treatFg_baseline1_April172019",row.names = NULL)
z1 <- z1[(1:nrow(z1)),(2:ncol(z1))]
z2 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Data/Antigen_April17_results_general_population/NO_outputDirectoryFg_Georgia_April242019_atb-pcr/treatFg_baseline1_April242019",row.names = NULL)
z2 <- z2[(1:nrow(z2)),(2:ncol(z2))]
z3=rbind(z1,z2)
write.table(z3, "treatFg_baseline1_April232019")

#######################################################################################
z1 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Code/Antigen_April11-noretest_general - Copy/Yes_outputDirectoryFg_Georgia_April172019_atb-antg/treatFg_baseline1_April172019",row.names = NULL)
z1 <- z1[(1:nrow(z1)),(2:ncol(z1))]
z2 = atb_antg
z3=rbind(z1,z2)
  
write.table(z3, "treatFg_baseline1_April172019_1")

z1 = read.table("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Code/Antigen_April11-noretest_general - Copy/Yes_outputDirectoryFg_Georgia_April172019_atb-pcr/treatFg_baseline1_April172019",row.names = NULL)
z1 <- z1[(1:nrow(z1)),(2:ncol(z1))]
z2 = atb_antg
z3=rbind(z1,z2)

write.table(z3, "treatFg_baseline1_April172019_1")
######################################################################################


  #setwd("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_antg/")
  test_scenario = c("pcr","antg","atb+antg", "atb+pcr", "atb-antg", "atb-pcr")

 # Cohort_comb = data_atb_pcr
  outcome_cost_matrix <-matrix(NA,7,10)
  colnames(outcome_cost_matrix)<-c("F0","F1","F2","F3","F4", "DC", "HCC", "death", "liv.death.tot", "liv.death.noHCV")
  rownames(outcome_cost_matrix)<-c("pcr","antg","atb+antg", "atb+pcr", "atb-antg", "atb-pcr","")
  case = 1
  z1 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_antg/treatFg_baseline%s_November202018",case),row.names = NULL)
  z1 <- z1[(2:nrow(z1)),(2:ncol(z1))]
  data_antg = z1
  for (case in 2:403)
  {
    print(case)
  z2 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_antg/treatFg_baseline%s_November202018",case),row.names = NULL)
  z2 <- z2[(2:nrow(z2)),(2:ncol(z2))]
  data_antg=rbind(data_antg,z2)
  }
  Cohort_comb = data_antg
  index_test_type = 2
  test_type <- test_scenario[index_test_type]
  source("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/post_processing.R")
  case = 1
  z1 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb-antg/treatFg_baseline%s_November202018",case),row.names = NULL)
  z1 <- z1[(2:nrow(z1)),(2:ncol(z1))]
  data_atb_antg = z1
  for (case in 2:403)
  {
    print(case)
    z2 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb-antg/treatFg_baseline%s_November202018",case),row.names = NULL)
    z2 <- z2[(2:nrow(z2)),(2:ncol(z2))]
    data_atb_antg=rbind(data_atb_antg,z2)
  }
  Cohort_comb = data_antg
  index_test_type = 5
  test_type <- test_scenario[index_test_type]
  source("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/post_processing.R")
  
  case = 1
  z1 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb-pcr/treatFg_baseline%s_November202018",case),row.names = NULL)
  z1 <- z1[(2:nrow(z1)),(2:ncol(z1))]
  data_atb_pcr = z1
  for (case in 2:403)
  {
    print(case)
    z2 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb-pcr/treatFg_baseline%s_November202018",case),row.names = NULL)
    z2 <- z2[(2:nrow(z2)),(2:ncol(z2))]
    data_atb_pcr=rbind(data_atb_pcr,z2)
  }
  Cohort_comb = data_atb_pcr
  index_test_type = 6
  test_type <- test_scenario[index_test_type]
  source("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/post_processing.R")
  case = 1
  z1 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb+antg/treatFg_baseline%s_November202018",case),row.names = NULL)
  z1 <- z1[(2:nrow(z1)),(2:ncol(z1))]
  data_atbantg = z1
  for (case in 2:403)
  {
    print(case)
    z2 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb+antg/treatFg_baseline%s_November202018",case),row.names = NULL)
    z2 <- z2[(2:nrow(z2)),(2:ncol(z2))]
    data_atbantg=rbind(data_atbantg,z2)
  }
  Cohort_comb = data_atbantg
  index_test_type = 3
  test_type <- test_scenario[index_test_type]
  source("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/post_processing.R")
  
  case = 1
  z1 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb+pcr/treatFg_baseline%s_November202018",case),row.names = NULL)
  z1 <- z1[(2:nrow(z1)),(2:ncol(z1))]
  data_atbpcr = z1
  for (case in 2:403)
  {
    print(case)
    z2 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_atb+pcr/treatFg_baseline%s_November202018",case),row.names = NULL)
    z2 <- z2[(2:nrow(z2)),(2:ncol(z2))]
    data_atbpcr=rbind(data_atbpcr,z2)
  } 
  Cohort_comb = data_atbpcr
  index_test_type = 4
  test_type <- test_scenario[index_test_type]
  source("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/post_processing.R")
  
  case = 1
  z1 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_pcr/treatFg_baseline%s_November202018",case),row.names = NULL)
  z1 <- z1[(2:nrow(z1)),(2:ncol(z1))]
  data_pcr = z1
  for (case in 2:403)
  {
    print(case)
    z2 = read.table(sprintf("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/outputDirectoryFg_November202018_pcr/treatFg_baseline%s_November202018",case),row.names = NULL)
    z2 <- z2[(2:nrow(z2)),(2:ncol(z2))]
    data_pcr=rbind(data_pcr,z2)
  }
  Cohort_comb = data_pcr
  index_test_type = 1 
  test_type <- test_scenario[index_test_type]
  source("S:/nas02/Fac/SHARE/GKeiser/HCV_results/Antigen/Antigen_codes/post_processing.R")
 
 