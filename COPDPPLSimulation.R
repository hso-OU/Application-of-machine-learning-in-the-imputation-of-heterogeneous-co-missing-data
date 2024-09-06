##COPD Population and Sample Simulation 



load(file="CLSAPpl1.Rdata")
rm(list=setdiff(ls(), "Ppl.Data"))
gc()
Ppl.Data2<-  Ppl.Data[order(Ppl.Data$entity_id),] 
rm(Ppl.Data)
gc()

library(dplyr)
Temp<- Ppl.Data2 %>% group_by(Sampling_Strata, WGHTS_PROV_TRM )%>% summarize( sample.size = max(Sample_size))%>% 
  group_by( WGHTS_PROV_TRM ) %>% summarize( Sample.size.Prov = sum(sample.size))
Ppl.Data2<- left_join(Ppl.Data2,  Temp, by = "WGHTS_PROV_TRM")


source("COPDSimulationFunction.R")



## under this situation. MI have less value thatn the  naive estimate most of time
hid.ratio<- function(data= Ppl.Data2) { c(0.75, 0.65)[data$Disease1 + 1] }  
FEV1FVC.func <- function(hid.ratio, FVC ){ hid.ratio* FVC  + rnorm(length(FVC))* 0} 
Ppl.Data2 <-PPL.Disease.Missing.func(Ppl.Data2,hid.ratio, FEV1FVC.func, Disease.prob.func,
                                     Miss.prob.func, miss.lev= 0.10, rand.seed= c( 23456, 1234567 ))


save(Ppl.Data2, file="MedMissingPpl.Rdata")


Strat_sample.func(Ppl.Data2, sample.frac = 0.5,  
                  file.name="Miss_Med_Size_Med_Strat_Sample.RData", iter=10000, rand.seed=1)

Strat_sample.func(Ppl.Data2, sample.frac = 0.1,  
                  file.name="Miss_Med_Size_Sml_Strat_Sample.RData", iter=10000, rand.seed=4)

Strat_sample.func(Ppl.Data2, sample.frac = 1,  
                  file.name="Miss_Med_Size_Lrg_Strat_Sample.RData", iter=10000, rand.seed=7)
###############################################################################################

hid.ratio<- function(data= Ppl.Data2) { c(0.75, 0.65)[data$Disease1 + 1] }  
FEV1FVC.func <- function(hid.ratio, FVC ){ hid.ratio* FVC  + rnorm(length(FVC))* 0} 
Ppl.Data2 <-PPL.Disease.Missing.func(Ppl.Data2,hid.ratio, FEV1FVC.func, Disease.prob.func,
                                     Miss.prob.func, miss.lev= 0.05, rand.seed= c( 23456+1, 1234567+1 ))


save(Ppl.Data2, file="SmlMissingPpl.Rdata")


Strat_sample.func(Ppl.Data2, sample.frac = 0.5,  
                  file.name="Miss_Sml_Size_Med_Strat_Sample.RData", iter=10000, rand.seed=1+1)

Strat_sample.func(Ppl.Data2, sample.frac = 0.1,  
                  file.name="Miss_Sml_Size_Sml_Strat_Sample.RData", iter=10000, rand.seed=4+1)

Strat_sample.func(Ppl.Data2, sample.frac = 1,  
                  file.name="Miss_Sml_Size_Lrg_Strat_Sample.RData", iter=10000, rand.seed=7+1)
###############################################################################################

hid.ratio<- function(data= Ppl.Data2) { c(0.75, 0.65)[data$Disease1 + 1] }  
FEV1FVC.func <- function(hid.ratio, FVC ){ hid.ratio* FVC  + rnorm(length(FVC))* 0} 
Ppl.Data2 <-PPL.Disease.Missing.func(Ppl.Data2,hid.ratio, FEV1FVC.func, Disease.prob.func,
                                     Miss.prob.func, miss.lev= 0.2, rand.seed= c( 23456+2, 1234567+2 ))


save(Ppl.Data2, file="LrgMissingPpl.Rdata")


Strat_sample.func(Ppl.Data2, sample.frac = 0.5,  
                  file.name="Miss_Lrg_Size_Med_Strat_Sample.RData", iter=10000, rand.seed=1+2)

Strat_sample.func(Ppl.Data2, sample.frac = 0.1,  
                  file.name="Miss_Lrg_Size_Sml_Strat_Sample.RData", iter=10000, rand.seed=4+2)

Strat_sample.func(Ppl.Data2, sample.frac = 1,  
                  file.name="Miss_Lrg_Size_Lrg_Strat_Sample.RData", iter=10000, rand.seed=7+2)