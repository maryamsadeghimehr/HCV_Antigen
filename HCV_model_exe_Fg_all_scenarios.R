#'@name HCV_model_exe
library(MASS)
library(msm)
library(mstate)
library(survival)
library(splines)
library(plyr)
library(gems)
library(readr)
#########################################################
start.time <- Sys.time()
current_dir = getwd()
wd =  getwd()
########################################################
cohortSize <<- 200
blocks_nb <<- 6
blocks_size <<- 5
statesNum <<- 52
maxTime <<- 100 # specify duration of simulation in years
ans <- "Fg"
test_type <- "pcr"
diag_scenario = c("Baseline","Birth","IDU", "ex_IDU")
diag_ind <<- 1
###############################################################################################
### Baseline viral load prepration
### Using the real data for viral load
##############################################################################################
# patient_vlvl_biopsy = read.csv("patient_vlvl_biopsy.csv")
# p_count = count(patient_vlvl_biopsy$pat.id)
# 
# k = 0
# tmp3 = matrix(NA,nrow(p_count),2 )
# for(i in 1:nrow(p_count))
# {
#   tmp = p_count$freq[i]
#   tmp3[i,1] = patient_vlvl_biopsy$pat.id[1 + k]
#   tmp3[i,2] = patient_vlvl_biopsy$hcvrna.level[1 + k]
#   k = k + tmp 
# }
# plot(log10(tmp3[,2]), type = "o", lwd = 3, ylim = c(0,10))
# plot(density(log10(tmp3[,2])), main = "Density of the first viral load")
# real_vl_density = density(log10(tmp3[,2]))
# df <- approxfun(density(log10(tmp3[,2])))
###########################################################################################
# viral load from the whole population
###########################################################################################
patient_biopsy <- read.csv("~/Antigen_git/Antigen/Serology_data.csv", sep="")
pid = count(patient_biopsy$pat.id)
k = 1
vl = matrix(0,nrow(pid),1)
for(i in 1:nrow(pid))
{
  a = which(patient_biopsy[,"pat.id"] == pid$x[i])
  vl[k,1] = patient_biopsy[a[1],"hcvrna.level"]
  k = k + 1
}
plot(density(na.omit(log10(vl[,1]))), main = "Density of the first viral load")
real_vl_density <- density(na.omit(log10(vl[,1])))
#df <- approxfun(density(na.omit(log10(vl[,1]))))
###########################################################################################
bm <<- readr::read_table("Mort_CH_Swiss.txt")
data = read.table("baseline_data.txt")
case = 7
source("baseline4.R")
all_cases_Diagtime = matrix(NA,cohortSize,6)
colnames(all_cases_Diagtime) <- c("antg","pcr","atb+antg", "atb+pcr", "atb-antg", "atb-pcr")
for(i in 1:6)
{
test_type_ind = i
test_scenario = c("antg","pcr","atb+antg", "atb+pcr", "atb-antg", "atb-pcr")
test_type <- test_scenario[test_type_ind]


all_cases = matrix(NA,1,52 + 13)

#for (case in 1:nrow(data))
  for (case in 7:7)
  {
    case <<- case
    #####################Source Codes ###########################
    diag_results = diag.time.fun(cohortSize, maxTime, period=1, screen_proba=1, test_type,  V0,bl)
    bl[,"Diag.time"]<- diag_results[1:cohortSize] 
    all_cases_Diagtime[,test_type_ind] = bl[,"Diag.time"]
    
    source("index_block_fs_stage.R")
    source("hazard_fun.R") # call file specifying the hazard functions
    source("params_cov.R") #choose parameters file
    source("TransitionTime.R")
    source("sim_cohort.R") # call the cohort simulator
    mainDir <- wd
    today <- Sys.Date()
    dd = format(today, format = "%B%d%Y")
    subDir <- sprintf("outputDirectoryFg_%s_%s",dd,test_type)
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    wd = setwd(file.path(mainDir, subDir))
    mainDir <- wd
    setwd(file.path(mainDir, subDir))
    write.table(cohort@time.to.state, sprintf("treatFg_%s_%s",case, dd))
    ########################################################
    
    write.table(cbind(cohort@time.to.state, bl), sprintf("treatFg_baseline%s_%s",case, dd)) 
    all_cases = rbind(as.matrix(all_cases),as.matrix(read.table(sprintf("treatFg_baseline%s_%s",case, dd))))
    setwd(current_dir)
    source("post_processing.R")

  }
}
# setwd("../Figures/")
# hist(bl[,"Viral_diag"], col = "blue", xlab = "Initial viral load", main = "Histogaram of initial viral load")
# hist(bl[which(bl[,"Diag.time"] < 1000),"Viral_diag"], col = "grey", xlab = "viral load", main = "Histogaram of viral load at the time of diagnosis")
# plot(w, lwd = "3", col = "blue", xlab ="Year of birth", main = "Year of birth distribution")
# w = density(bl[,"BirthYear"])
# plot(w, lwd = "3", col = "blue", xlab ="Year of birth", main = "Year of birth distribution")
# setwd("../R")
# 
# tmp = Cohort[which(bl[,"Diag.time"] < 1000),]

antibodypcr = read.table(sprintf("outcome_%s_%s_atb+pcr",case,dd))
ant_pcr = read.table(sprintf("outcome_%s_%s_atb-pcr",case,dd))
antibodyantigen = read.table(sprintf("outcome_%s_%s_atb+antg",case,dd))
ant_antigen = read.table(sprintf("outcome_%s_%s_atb-antg",case,dd))
pcr = read.table(sprintf("outcome_%s_%s_pcr",case,dd))
antigen = read.table(sprintf("outcome_%s_%s_antg",case,dd))

case = 3
test_type = "atb+pcr"
antibodypcr = read.table(sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type))
test_type = "atb-pcr"
ant_pcr = read.table(sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type))
test_type = "atb+antg"
antibodyantigen = read.table(sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type))
test_type = "atb-antg"
ant_antigen = read.table(sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type))
test_type = "pcr"
pcr = read.table(sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type))
test_type = "antg"
antigen = read.table(sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_%s_%s_%s",case,dd, test_type)) 

plot(antibodypcr$x[1:10], col = "red")
lines(antibodyantigen$x[1:10])
lines(ant_antigen$x[1:10], col = "green")
lines(ant_pcr$x[1:10],  col = "yellow")
lines(pcr$x[1:10], col = "grey")
lines(antigen$x[1:10], col = "blue")

outcome_table <- antibodypcr
outcome_table$antibodyantigen <- antibodyantigen$x
outcome_table$ant_antigen <- ant_antigen$x
outcome_table$ant_pcr <- ant_pcr$x
outcome_table$antigen <- antigen$x
outcome_table$pcr <- pcr$x
outcome_table$antibodypcr <- antibodypcr$x
outcome_table$x <- c("F0","F1","F2","F3","F4","DC","HCC","death","liv.death.tot","liv.death.noHCV")
write.csv(outcome_table, sprintf("C:/Users/sadeghim/Documents/Antigen_git/Outcomes/outcome_p%s.csv",case))


plot(density(V0), main = "Density of the first viral load",xlab = "Viral load (log10)", col = "blue", lwd = 3)
polygon(density(V0), col = "blue")

 # Time to detection #########################################################################################
bl[,"Diag.time"]

plot(density(bl[,"Diag.time"]), main = "Density of diagnosis time",xlab = "Time (year)", col = "grey", lwd = 3)
polygon(density(bl[,"Diag.time"]), col = "grey")
