
##Constructing the strata 
Prov.vec   <- c( "AB", "BC", "MB","NB","NL","NS","ON","PE","QC","SK")
Prov.DCS   <- c( "AB", "BC", "MB",     "NL","NS","ON",     "QC") 
Prov.NonDCS<- c( "NB", "PE", "SK")     

Prov.prop  <- c(71000,79000,32000,14000,11000,16000,82000,2800,2e+05,12000)

Age.grp    <- c("45-54","55-64","65-74","75-85")
Sex.grp    <- c("female","male")

## Education   1:11
ED_UDR11.prop<- c( 8000, 17000, 12000, 66000,30000,58000,11000,
                   40000,100000,75000,7000)
## SEX  Female Male
Sex_Ask.prop<- c(F=306000, M=212000)

## non_DCS  FASLE TRUE
DCS.prop<- c( DCS=250000, non_DCS=270000 )

## Age  45:85 as mixture of two normal  
Age.gen.func <- function(n=1){ 
  temp.vec <- rbinom(n,size=1,prob=0.27)+1
  temp.vec <- rnorm(n, mean= c(61,65)[temp.vec],sd=c(7.7,2.7)[temp.vec])
  temp.vec <- round(temp.vec,0)
  out.vec  <- which(! temp.vec%in% c(45:85))
  out.n    <- length(out.vec)
  if (out.n>0) temp.vec[out.vec] <- Age.gen.func( n= out.n)
  temp.vec
}


#Response 
## startlanguage  en fr 
startlanguage.prop<- c(en=320000,fr=200000)

## Experience dry mouth 
ORH_EXP_DRM.Prop<-rbind( No=c(F=230000,M=170000),
                         Yes=c(F= 80000,M= 40000))
## Afraid to walk alone after dark in local area
ENV_AFRDWLK.Prop<-rbind( Strongly_Agree = c(F=  6900,M=  1700),
                         Agree = c(F= 31000,M= 15000), 
                         Disagree = c(F=170000,M=130000), 
                         Strongly_Disagree = c(F= 76000,M= 62000))
#For age group c("45-48","49-54","55-64","65-74","75-85") 
ENV_AFRDWLK.Age.Prop<-rbind( Strongly_Agree = c(  0 ,  0,   0 ,500, 500),
                             Agree = c(  0 ,  0, 500 ,500,   0), 
                             Disagree = c(  0 ,500, 500 ,  0,   0), 
                             Strongly_Disagree = c(500 ,500,   0 ,  0,   0))

## Current marital status 
WEA_MRTL_CURRENT.Prop<-rbind( Single  = c(F=  33000, M= 20000), 
                              Married = c(F= 180000, M=170000),
                              Widowed = c(F=  37000, M=  6000),
                              Divorced =c(F=  38000, M=  8300),
                              Seperated=c(F=   5000, M=  6200))

## Generating function based on SEX
Gen.Sex.func<- function( Prob.Mat, Sex.vec ){
  Tags<- rownames(Prob.Mat)
  No.Tags<- length(Tags)
  ans<-sapply(Sex.vec, function(x)sum(1:No.Tags*rmultinom(1,1,prob=Prob.Mat[,x])))
  Tags[ans]
}

Gen.Sex.Age.func<- function( Prob.Mat, Sex.vec, Prob.Age.Mat, Age.vec ){
  Tags <- rownames(Prob.Mat)
  No.Tags <- length(Tags)
  Len <- length(Sex.vec)
  Age.Grp.Num<- 1 + (Age.vec>=49) + (Age.vec>=55) + (Age.vec>=65)+ (Age.vec>=75)
  #age group=c("45-48","49-54","55-64","65-74","75-85")
  ans<-sapply(1:Len, function(x)sum(1:No.Tags*rmultinom(1,1,prob=Prob.Mat[,Sex.vec[x]]+Prob.Age.Mat[,Age.Grp.Num[x]]) ) )
  Tags[ans]
}

## Height 
HWT_DHT_M_TRM.func<-function( sex.vec ){
  sex.vec1 <- (sex.vec =="M") + 1 ; n<-length(sex.vec)
  temp.vec <- rbinom(n,size=1,prob=c(0.97,0.13)[sex.vec1])+1
  temp.vec <- rnorm(n, mean= c(1.76,1.63,1.73,1.83)[temp.vec+2*sex.vec1-2]
                    ,sd= c(0.10,0.06,0.05,0.05)[temp.vec+2*sex.vec1-2])
  temp.vec 
}
#####HWT_DHT_M_TRM.func(c("F","M","M"))
## Weight 
HWT_DWT_K_TRM.func<- function(height,sex.vec){ 
  rnorm(length(height), mean= -66.66+84.55*height+3.21*(sex.vec=="M"), sd=14.3 )
}

####HWT_DWT_K_TRM.func( HWT_DHT_M_TRM.func(c("F","M","M")), c("F","M","M"))



##Generating a (pseudo) population (520000)
set.seed(1234)
#N<-520000
N <- 14e6  #13655060
Ppl.Data<- data.frame(
  entity_id = 1:N
  ,WGHTS_PROV_TRM   = Prov.vec[ colSums(1:10*rmultinom(N,size=1,prob =Prov.prop)) ]
  ,SEX_ASK_TRM =c("F","M")[rbinom(n=N, size=1, prob= Sex_Ask.prop[2]/sum(Sex_Ask.prop))+1]
  ,AGE_NMBR_TRM= Age.gen.func(n=N) 
  ,DCS.groups  =c("Non_DCS","DCS1","DCS2")[rbinom(n=N, size=1, prob= DCS.prop[1]/sum(DCS.prop))*
                                             (rbinom(n=N, size=1, prob= 0.5)+1 ) +1]  ## Simulation of  DCS1 or DCS2  for BC ON and QC 
  ,ED_UDR11_TRM      = colSums(1:11*rmultinom(N,size=1,prob =ED_UDR11.prop) )
  ,startlanguage_MCQ = c("en","fr")[rbinom(n=N, size=1, prob=startlanguage.prop[2]/sum(startlanguage.prop))+1]
) 

## Change DCS to Non_DCS for Non-DCS province
Ppl.Data$DCS.groups[which(Ppl.Data$WGHTS_PROV_TRM %in% Prov.NonDCS)]<-"Non_DCS"

## Combining  DCS1 and DCS2 in  AB MB NS NL
Ppl.Data$DCS.groups[with(Ppl.Data, which(WGHTS_PROV_TRM %in% c("AB","MB","NS","NL") & DCS.groups !="Non_DCS"))]<- "DCS1"

## Removing  DCS1 and DCS2
Ppl.Data$DCS.vec<- as.character(Ppl.Data$DCS.groups)
Ppl.Data$DCS.vec[which(Ppl.Data$DCS.vec !="Non_DCS")] <- "DCS"


Ppl.Data$ENV_AFRDWLK_MCQ = with(Ppl.Data, Gen.Sex.Age.func(ENV_AFRDWLK.Prop,SEX_ASK_TRM, ENV_AFRDWLK.Age.Prop,AGE_NMBR_TRM))
Ppl.Data$ORH_EXP_DRM_MCQ = Gen.Sex.func(ORH_EXP_DRM.Prop     , Ppl.Data$SEX_ASK_TRM) 
Ppl.Data$WEA_MRTL_CURRENT= Gen.Sex.func(WEA_MRTL_CURRENT.Prop, Ppl.Data$SEX_ASK_TRM)  
Ppl.Data$HWT_DHT_M_TRM   = HWT_DHT_M_TRM.func(Ppl.Data$SEX_ASK_TRM)
Ppl.Data$HWT_DWT_K_TRM   = HWT_DWT_K_TRM.func(Ppl.Data$HWT_DHT_M_TRM,
                                              Ppl.Data$SEX_ASK_TRM)

## Education 
temp.vec <- c("Low Education" , "Medium Education", 
              "Higher Education lower" , "Higher Education upper"  ) 
temp.vec2 <- c(1,1,1,2,2,2,3,3,4,4,2)
Ppl.Data$Education =  temp.vec[temp.vec2[Ppl.Data$ED_UDR11_TRM] ] 

Ppl.Data$Sample_Age_Gpr<-  Age.grp[ with(Ppl.Data, (45<=AGE_NMBR_TRM)+
                                           (55<=AGE_NMBR_TRM)+(65<=AGE_NMBR_TRM)+ (75<=AGE_NMBR_TRM) )]

Ppl.Data$WGHTS_GEOSTRAT_TRM<- with(Ppl.Data, paste(WGHTS_PROV_TRM,DCS.vec,sep="_") )

Ppl.Data$Sampling_Strata<- with(Ppl.Data, 
                                paste(WGHTS_PROV_TRM, DCS.vec, SEX_ASK_TRM
                                      , Sample_Age_Gpr,  sep="_") )

 

####### Basic  Variables Names 
##Constructing the strata 
Prov.vec   <- c( "AB", "BC", "MB","NB","NL","NS","ON","PE","QC","SK")
Prov.DCS   <- c( "AB", "BC", "MB",     "NL","NS","ON",     "QC") 
Prov.NonDCS<- c( "NB", "PE", "SK")     

Prov.prop  <- c(71000,79000,32000,14000,11000,16000,82000,2800,2e+05,12000)

Age.grp    <- c("45-54","55-64","65-74","75-85")
Sex.grp    <- c("female","male")

## Education   1:11
ED_UDR11.prop<- c( 8000, 17000, 12000, 66000,30000,58000,11000,
                   40000,100000,75000,7000)
## SEX  Female Male
Sex_Ask.prop<- c(F=306000, M=212000)

## non_DCS  FASLE TRUE
DCS.prop<- c( DCS=250000, non_DCS=270000 )


########Sample Size##########
#table(Ppl.Data$Sampling_Strata)

Sampling.mat<-expand.grid(Prov=Prov.vec, DCS.vec  =c("DCS","Non_DCS")
                          ,SEX =c("F","M"), Age.grp = Age.grp )  
Sampling.mat<-Sampling.mat[-which(Sampling.mat$Prov%in% Prov.NonDCS 
                                  & Sampling.mat$DCS.vec=="DCS"),]
rownames(Sampling.mat)<-1:136
Prov.Size   <- c(  AB= 90 , BC=100,MB=90,NB=60,NL=60,NS=60,ON=85,PE=55,QC=200,SK=50)
DCS.Size    <- c( DCS= 420, NON.DCS = 280)/700
#Sex.Size <- c(   F= 510,       M = 350)/860
Sex.Size <- c(   F= 500,       M = 490)/990

Sampling.mat$Sample_size<- round( with(Sampling.mat,Prov.Size[Prov]*
                                         Sex.Size[SEX] * DCS.Size[DCS.vec]^(1-Prov%in% Prov.NonDCS)  ) )
Sampling.mat$Sampling_Strata<-  with(Sampling.mat, 
                                     paste(Prov, DCS.vec, SEX, Age.grp,  sep="_") )
temp.vec <- table(Ppl.Data$Sampling_Strata)
Sampling.mat<-merge(Sampling.mat
                    , data.frame(Sampling_Strata=names(temp.vec ), N_h=c(temp.vec) )
                    , by="Sampling_Strata" ) 


library(dplyr)
## using the actual number of CLSA participants
CLSAParticipant <- read.csv("CLSAParticipant.csv")
### Removing Space
CLSAParticipant$Prov    <- as.factor( gsub(" ", "",as.character(CLSAParticipant$Prov) , fixed = TRUE) )
CLSAParticipant$Age.grp <- as.factor( gsub(" ", "",as.character(CLSAParticipant$Age.grp) , fixed = TRUE) )
CLSAParticipant$Sampling_Strata<-  with(CLSAParticipant,
                                        paste(Prov, DCS.vec, SEX, Age.grp,  sep="_") )

Sampling.mat$Sample_size<-NULL
Sampling.mat <- merge(Sampling.mat
                      , CLSAParticipant[,c("Sampling_Strata","Sample_size")]
                      , by="Sampling_Strata" )


## sample sizes 
Sampling.mat$n_h <- Sampling.mat$Sample_size
## 
Sampling.mat$GEOSTRAT_TRM <- with(Sampling.mat, paste( Prov, DCS.vec,sep="_") ) 

## basic design weight, initial weight
Sampling.mat$basic_weight <- Sampling.mat$N_h /Sampling.mat$n_h 
head(Sampling.mat)

### Merge sampling information back to the population data 
Ppl.Data<-merge(x=Ppl.Data ,
                y=Sampling.mat[,c("Sampling_Strata","Sample_size", "N_h", "n_h", "basic_weight","GEOSTRAT_TRM")],
                by= "Sampling_Strata")


#setwd("~/google-drive/MLSurvey")
save(Ppl.Data,file="CLSAPpl.Rdata")

#### Simulation start from now on ...  










