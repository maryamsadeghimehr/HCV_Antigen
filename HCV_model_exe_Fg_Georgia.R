#'@name HCV_model_exe
library(MASS)
library(msm)
library(mstate)
library(survival)
library(splines)
library(plyr)
library(gems)
library(readr)
#require(deSolve)
#########################################################
start.time <- Sys.time()
current_dir = getwd()
wd =  getwd()
########################################################
cohortSize <<- 50000
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
###########################################################################################
# viral load from the whole population
###########################################################################################

patient_biopsy <- read.csv("Serology_data.csv", sep="")
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
###########################################################################################
#bm <<- readr::read_table("Mort_CH_Swiss.txt")
bm <<- readr::read_table2("Mort_CH")
data = read.table("baseline_data.txt")

case = 1
#ans_IDU  <- readline(prompt="Simulating an IDU cohort? (Yes/No) ")

ans_IDU  <- c("Yes", "NO")
for(rr in 1:1)
{
if(ans_IDU[rr] == "Yes")
  source("baseline_Georgia_IDUs.R")
if(ans_IDU[rr] == "NO")
  source("baseline_Georgia_noIDUs2.R")

all_cases_Diagtime = matrix(NA,cohortSize,4)
colnames(all_cases_Diagtime) <- c("antg","pcr", "atb-antg", "atb-pcr")
for(ii in 1:4)
{
 # ii = 1
  test_type_ind = ii
  test_scenario = c("antg","pcr", "atb-antg", "atb-pcr")
  test_type <- test_scenario[test_type_ind]
  
  
  all_cases = matrix(NA,1,52 + 13)
 # for (case in 1:nrow(data))
    #for (case in 1:2)
  {
    #case <<- case
    #####################Source Codes ###########################
    diag_results = diag.time.fun(cohortSize, maxTime, period=4, screen_proba=1, test_type,  V0,bl)
    bl[,"Diag.time"]<- diag_results[1:cohortSize] 
    all_cases_Diagtime[,test_type_ind] = bl[,"Diag.time"]
    setwd(current_dir)
    source("index_block_fs_stage.R")
    source("hazard_fun.R") # call file specifying the hazard functions
    source("params_cov.R") #choose parameters file
    source("TransitionTime.R")
    source("sim_cohort.R") # call the cohort simulator
    mainDir <- wd
    today <- Sys.Date()
    dd = format(today, format = "%B%d%Y")
    subDir <- sprintf("%s_outputDirectoryFg_Georgia_%s_%s",ans_IDU[rr],dd,test_type)
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    wd = setwd(file.path(mainDir, subDir))
    mainDir <- wd
    setwd(file.path(mainDir, subDir))
    write.table(cohort@time.to.state, sprintf("treatFg_%s_%s",case, dd))
    ########################################################
    
    write.table(cbind(cohort@time.to.state, bl), sprintf("treatFg_baseline%s_%s",case, dd)) 
    all_cases = rbind(as.matrix(all_cases),as.matrix(read.table(sprintf("treatFg_baseline%s_%s",case, dd))))
    setwd(current_dir)
    #source("post_processing.R")
    
  }
}
}# end for rr
