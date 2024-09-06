#### Population simulation 

## Simulating a realistic (Canadian) population 
# setwd("G:/My Drive/MLSurvey/Supp")
# source("A01_SimulationCLSA.R") 
setwd("G:/My Drive/MLSurvey")
save(Ppl.Data1, file="Ppl.Data1.RData")




## simple version 








## 
## simulation:  SimSurvey: An R package for comparing the design and analysis of surveys by simulating spatially-correlated population
## too complicated 

# library(SimSurvey)
# 
# set.seed(438) 
# pop<- sim_abundance()%>%sim_distribution()
# 
# a<-pop%>% sim_survey(n_sims=5, set_den = 1/1000, lengths_cap=100,  ages_cap = 5)
# 


#load("G:/My Drive/MLSurvey/simulatedPPl.RData")


## nhanes as an illustration 
#https://cran.r-project.org/web/packages/nhanesA/vignettes/Introducing_nhanesA.html
library(nhanesA)
nhanesTables('EXAM', 2005)








