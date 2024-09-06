#### Simulatuon  Estimation 




#### Estimation initialization
#setwd("~/COPDSimulation")

#args = as.numeric( commandArgs(trailingOnly=TRUE) )

source("COPDSimulationFunction2a.R")

##### Estimation setting
miss.n <-  1  ## 1 = Sml, 2 = Med, 3 = Lrg 
sample.n<- 3  ## 1 = Sml, 2 = Med, 3 = Lrg
# print( args )
level.vec <-c("Sml","Med","Lrg","XSm")

load(file=paste(level.vec[miss.n], "MissingPpl.Rdata",sep="") )
strat_sample.list<-read.table(file= paste("Miss_", level.vec[miss.n], "_Size_" , level.vec[sample.n], "_Strat_Sample.csv",sep="")
                              ,sep="," ,header=TRUE)#,  skip=3, nrows = 1)



set.seed( miss.n*1000 + sample.n*100 + 3)



i<-441

#i<- 448  #1

while( i <= 1000){
  
  if( i == 442){ i<- 448}
  
  strat_sample <- strat_sample.reconstruct(Ppl.Data2 , strat_sample.list[i,])
  print ( paste( "Now working on i=", i ) )
  
  res <- try(
  Simulation.wrap( strat_sample, rand.seed=NULL, var.char, num.imp.method = c("pmm","sample", "norm.nob","norm.predict"),
                   file.prefix= paste("Miss_", level.vec[miss.n], "_Size_", level.vec[sample.n],sep=""), simindex=i)
  )
  
  ### try it once more
  if(inherits(res, "try-error")){
    print( "Error occured. Try it again.")
    res <- try(
      Simulation.wrap( strat_sample, rand.seed=NULL, var.char, num.imp.method = c("pmm","sample", "norm.nob","norm.predict"),
                       file.prefix= paste("Miss_", level.vec[miss.n], "_Size_", level.vec[sample.n],sep=""), simindex=i)
    )
 
  }
  
  
  if(inherits(res, "try-error")){
    print( "Error occured again. Skip this iteration.") 
       #error handling code, maybe just skip this iteration using
    next; 
  }
  #rest of iteration for case of no error
  i<- i+1 
}



