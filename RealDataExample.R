### Real Data Analysis 
### Real data from  
### https://www.kaggle.com/datasets/rabieelkharoua/asthma-disease-dataset
library(dplyr)
library(readr)
#asthma <- read_csv("asthma_disease_data.csv")

## https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx#core
## 9A. Spirometry (Raw Curves) (June 2001)
## This data release, Series 11 No. 9A contains the NHANES III raw spirometry data file and documentation. This data release does not replace the previous NHANES III data releases
##NH3SPIRO <- read_csv("NH3SPIRO.CSV")

## Source and data dictionary NHANES 2012: 
## Spirometry data 
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/SPX_G.htm
## merged with demographic and the body measurement 
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DEMO_G.htm
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BMX_G.htm
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/SMQ_G.htm
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/SMQFAM_G.htm
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/SMQRTU_G.htm
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/HIQ_G.htm
## https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/HSQ_G.htm

library(haven)
SPX_G<-read_xpt("SPX_G.XPT")   ## Spirometry Data 
DEMO_G<-read_xpt("DEMO_G.XPT") ## Demographic Data 
BMX_G<-read_xpt("BMX_G.XPT")   ## Body Measurement Data 
SMQ_G<-read_xpt("SMQ_G.XPT")       ## Smoking - Cigarette Use Data 
SMQFAM_G<-read_xpt("SMQFAM_G.XPT") ## Smoking - Household Smokers Use Data 
SMQRTU_G<-read_xpt("SMQRTU_G.XPT") ## Smoking - Recent Tobacco Use 
HIQ_G<-read_xpt("HIQ_G.XPT") ## Health Insurance 
HSQ_G<-read_xpt("HSQ_G.XPT") ## Current Health Status 
RDQ_G<-read_xpt("RDQ_G.XPT") ## Respiratory Health 



## Create the dataset (only consider the baseline one for brevity)
Real_data<- SPX_G[,c("SEQN", "SPXNFVC","SPXNFEV1")]
#The proportion of missingness: 
mean(is.na(Real_data[,2]) )
# 0.112475
# Proportion of people having lung function impairment
mean(Real_data$SPXNFEV1/Real_data$SPXNFVC <=0.7,na.rm=TRUE)



### Merging other variables": 
Real_data <- Real_data %>% left_join( 
  DEMO_G[,c(
            "RIAGENDR",  ### gender  
            "RIDAGEYR",  ### age at screen
            "RIDRETH3",  ### Race/Hispanic origin w/ NH Asian
            "SIAINTRP",  ### Interpreter used in SP Interview?
            "INDFMPIR",  ### Ratio of family income to poverty
            "SEQN" )] , by =c("SEQN") )


Real_data <- Real_data %>% left_join( 
  BMX_G[,c( 
          "BMXWT"   ,  ## Weight (kg)
          "BMXHT"   ,  ## Standing Height (cm)
        #"BMXBMI"  ,  ## Body Mass Index (kg/m**2)
          "SEQN" )] , by =c("SEQN") )

## number of yes to the questions
RDQ_G$Num_yes <- apply( 
  RDQ_G[,c("RDQ031", ## Coughing most days - over 3 mo period
           "RDQ050", ## Bring up phlegm most days - 3 mo period
           "RDQ070", ## Wheezing or whistling in chest - past yr
           "RDQ100", ## Chest sound wheezy during exercise
           "RDQ134", ## Doctor prescribe wheezing medication
           "AGQ030"  ## Episode of hay fever in past 12 months?
)]==1,1,sum,na.rm=T)

Real_data <- Real_data %>% left_join( 
  RDQ_G[,c( 
    "Num_yes"   ,  ## number of yes to the above questions 
    "SEQN" )] , by =c("SEQN") )

#### Tabacco exposure   

Real_data <- Real_data %>% left_join( 
  SMQFAM_G[,c( 
    "SMD410", ###  Does anyone smoke inside home?
    "SEQN" )] , by =c("SEQN") )

### % of missing data in each variable 
lapply(Real_data, function(x){mean(is.na(x))} )

### all the variables (except  FVC and FEV1) have missingness less than 10% 
### use MI to fill in their  missingness first. 

### Define the categorical variables: 
Real_data$RIAGENDR <- as.factor( c("Male","Femnale")[Real_data$RIAGENDR])
Real_data$RIDRETH3 <- as.factor( c("Mexican American",
                          "Other Hispanic",	
                          "Non-Hispanic White"  ,	
                          "Non-Hispanic Black"	,"",
                          "Non-Hispanic Asian"	,	
                          "Other Race - Including Multi-Racial")[
                            Real_data$RIDRETH3
                          ])

Real_data$SIAINTRP<- as.factor( c("Yes","No")[Real_data$SIAINTRP])
Real_data$SMD410  <- as.factor( c("Yes","No")[Real_data$SMD410])

## Use the default method to impute
library(mice)
Real_data_temp<-mice::mice(Real_data[,-c(1:3)], m=1, maxit =5
                        , seed= 1234,print=FALSE )

## discard the ID variable as not longer necessary
Real_data_Mod <- cbind( Real_data[,2:3], complete(Real_data_temp,1))

Real_data_Mod$Miss.ind <-  is.na(Real_data_Mod$SPXNFVC)
Real_data_Mod<- Real_data_Mod%>% rename( FEV1 = SPXNFEV1)
Real_data_Mod<- Real_data_Mod%>% rename( FVC = SPXNFVC)
### Source the codes
source("COPDSimulationFunction2b.R")

# var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", "Smoking",
#               "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
#               "FVC", "FEV1" )


Result<-list()
var.char <- c("RIDRETH3", "RIAGENDR", "SIAINTRP", "SMD410",
               "RIDAGEYR", "BMXHT" , "BMXWT", "Num_yes",
              "FVC","FEV1" )

num.imp.methods= c("pmm","sample", "norm.nob","norm.predict")
set.seed(1234)

for( i in 1:4){
num.imp.method<- num.imp.methods[i]
Result[[1+(i-1)*3]] <- Est.func( strat_sample=Real_data_Mod, rand.seed=NULL, var.char=var.char,  num.imp.method = num.imp.method ) ## no clustering
Result[[2+(i-1)*3]] <-Est.func2( strat_sample=Real_data_Mod, rand.seed=NULL, var.char=var.char,  num.imp.method = num.imp.method ) ## with k prototype clustering
Result[[3+(i-1)*3]] <-Est.func3( strat_sample=Real_data_Mod, rand.seed=NULL, var.char=var.char,  num.imp.method = num.imp.method ) ## with DBscan clustering
}
## 
Result[[13]] <-Est.func4( strat_sample=Real_data_Mod, rand.seed=NULL, var.char=var.char,  train.sample.size = 500 ) ## Random forest
## Warnings indicates the some resampled datasets in training are not suitable for nnet (may be due to the number of "yes" in SIAINTRP is too small)
Result[[14]] <-Est.func5( strat_sample=Real_data_Mod, rand.seed=NULL, var.char=var.char,  train.sample.size = 500 )  ##  USE ANN from nnet 
Result[[15]] <-Est.func6( strat_sample=Real_data_Mod, rand.seed=NULL, var.char=var.char,  train.sample.size = 500 )  ## Use rf (Random forest)






