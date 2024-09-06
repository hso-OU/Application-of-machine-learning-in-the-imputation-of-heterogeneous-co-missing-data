##COPD Simulation functions 
library(dplyr)

PPL.Disease.Missing.func<- function( Ppl.Data2,hid.ratio, FEV1FVC.func, Disease.prob.func,
                                     Miss.prob.func, miss.lev= 0.10, rand.seed){
  ### modify the Ppl first before the sampling: 
  ### Creat a scenario to predict something:  
  set.seed(rand.seed[1])
  #Ppl.Data2$Gen1 <- rbinom( , 1, prob= 0.1 + 0.2*Ppl.Data2$WGHTS_PROV_TRM %in% c("BC","ON","QC") )
  Ppl.Data2$Gen1        <-  rpois( dim(Ppl.Data2)[1], lambda=0.5 );
    Ppl.Data2$Gen1[which(Ppl.Data2$Gen1>2)]<-2
  Ppl.Data2$Risk.factor <-  rnorm( dim(Ppl.Data2)[1], mean=c(170,150,190)[Ppl.Data2$Gen1+1], sd=1)  
  
  Ppl.Data2$Env1 <- rbinom( dim(Ppl.Data2)[1], 1, prob= 0.1 + 0.2*Ppl.Data2$WGHTS_PROV_TRM %in% c("BC","ON","QC") )
  ## smoking @ 2017: https://uwaterloo.ca/tobacco-use-canada/adult-tobacco-use/smoking-provinces
  smoking.prov.vec <- c( "BC", "AB", "SK","MB", "ON","QC", "NB", "NS","PE","NL" )
  smoking.prevalence<-c(15.6,18.9,17.8,14.5,12.9,15.7,13.7,18.5,11.8,20.1)/100
  names(smoking.prevalence)<-smoking.prov.vec
  
  ## order it to alphabetical order
  smoking.prevalence<-smoking.prevalence[order(smoking.prov.vec)]
  smoking.prov.vec  <-smoking.prov.vec[order(smoking.prov.vec)]
  
  Ppl.Data2$Smoking <- rbinom( dim(Ppl.Data2)[1], 1, prob=smoking.prevalence[Ppl.Data2$WGHTS_PROV_TRM])
  
  
  Ppl.Data2$Disease1 <- rbinom( dim(Ppl.Data2)[1], 1, prob =  Disease.prob.func(Ppl.Data2))

  
  ### getting Disease or not is hidden but can be measured by the ratio of FEV1 and FVC 
  ### say FEV1  to FVC < 0.7 have a disease 
  Ppl.size<- dim(Ppl.Data2)[1]
  set.seed(rand.seed[2])
  #Ppl.Data2$hid.ratio<- runif(Ppl.size)*c(0.05, 0.05)[Ppl.Data2$Disease1 + 1] + c(0.7, 0.65)[Ppl.Data2$Disease1 + 1]
  Ppl.Data2$hid.ratio<- hid.ratio(data=Ppl.Data2 )
  # Ppl.Data2$hid.ratio<- runif(Ppl.size)*c(0.15, 0.15)[Ppl.Data2$Disease1 + 1] + c(0.7, 0.55)[Ppl.Data2$Disease1 + 1] -
  #                                                                               c(  0, 0.10)[Ppl.Data2$DiseaseA + 1] 
  
  ### mean taken as Caucasian  Age  60 , 170 cm from https://www.cdc.gov/niosh/topics/spirometry/refcalculator.html 
  Ppl.Data2$FVC <-   (Ppl.Data2$SEX_ASK_TRM=="M") +   (3.23-2.51)/2* rnorm(Ppl.size) +3.23/2 
  Ppl.Data2$FEV1<- FEV1FVC.func(Ppl.Data2$hid.ratio, Ppl.Data2$FVC)
  Ppl.Data2$ratio <- Ppl.Data2$FEV1/Ppl.Data2$FVC
  
  
  #### missing  
  
  ( root.temp <-uniroot(function(x){mean(Miss.prob.func(x, data=Ppl.Data2)) -miss.lev},  interval= c(-5,5)) )
  miss.par<- root.temp$root
  
  Ppl.Data2$Miss.prob<- Miss.prob.func(x=miss.par, data=Ppl.Data2 )
  
  Ppl.Data2$Miss.ind <- rbinom( Ppl.size, size=1, prob=Ppl.Data2$Miss.prob) 
  
  
  
  ## prevalence of disease one (our target)
  print( paste( " Ppl. prevalence (hidden): ",mean(Ppl.Data2$Disease1)))
  print( paste( " Ppl. prevalence (realized): ",  mean(Ppl.Data2$ratio[which(Ppl.Data2$Miss.ind==0)]< .7) ))
  
  print( paste( " Ppl. missing probability range: ",  paste( round(range(Ppl.Data2$Miss.prob),4), collapse = "," )))
  print( paste( " Ppl. missing (proportion): ",    mean(Ppl.Data2$Miss.ind) ))
  
  
  print( paste( " Ppl. disease vs missing: "))
  print( tab<- xtabs( ~ Disease1 + Miss.ind, data =Ppl.Data2) ) 
  print( paste( " Ppl. prevalence (realized, complete case only):" , tab[2,1]/sum(tab[,1])  )   )
  
  
  
  Ppl.Data2<- Ppl.Data2 %>%   group_by(WGHTS_PROV_TRM) %>% mutate(Ppl_Size_Prov=n() ) %>% ungroup()
  
  Ppl.Data2<- Ppl.Data2 %>%   group_by(Sampling_Strata) %>% mutate(Ppl_Size_Strata=n() ) %>% ungroup()
  
  
  ### data cleaning 
  Ppl.Data2$ED_UDR11_TRM <- ordered(Ppl.Data2$ED_UDR11_TRM, levels = c(1:6,11,7:10) )
  Ppl.Data2$Education <- ordered( Ppl.Data2$Education, levels= c("Low Education", "Medium Education", 
                                                                 "Higher Education lower", "Higher Education upper") )
  Ppl.Data2$entity_id <- as.factor(Ppl.Data2$entity_id)
  Ppl.Data2$ENV_AFRDWLK_MCQ  <- ordered(Ppl.Data2$ENV_AFRDWLK_MCQ , levels= c("Strongly_Disagree", "Disagree",
                                                                              "Agree", "Strongly_Agree") ) 
  Ppl.Data2$Sample_Age_Gpr <- ordered(Ppl.Data2$Sample_Age_Gpr, levels=c( "45-54", "55-64",  "65-74",  "75-85") )
  
  Ppl.Data2$Smoking <- as.factor( Ppl.Data2$Smoking )
  
  for(  i in 3:(dim(Ppl.Data2)[2])){ if (class( Ppl.Data2[[i]])[1] == "character" ) { Ppl.Data2[[i]]<-as.factor(Ppl.Data2[[i]]) } }
  
  Ppl.Data2 
}


## sampling function 
strat_sample.func<- function( Ppl.Data2, sample.frac = 1, rand.seed=NULL){
  if(!is.null(rand.seed)){set.seed(rand.seed[1])}
  
  strat_sample <- Ppl.Data2 %>%   group_by(WGHTS_PROV_TRM) %>%  
    #slice_sample(n=Sample_size)
    sample_n(round(Sample.size.Prov*sample.frac))
  
  
  ### create missing data
  strat_sample[ which(strat_sample$Miss.ind==1), "FVC" ] <-NA
  strat_sample[ which(strat_sample$Miss.ind==1), "FEV1" ] <-NA
  strat_sample$Gen1 <- NULL 
  strat_sample$Disease1 <- NULL
  
  strat_sample$hidden.ratio <-strat_sample$ratio
  
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  
  strat_sample$Hidden.Miss.prob <- strat_sample$Miss.prob
  
  #mean( strat_sample$ratio<0.7 , na.rm=TRUE)
  
  
  # Change  the sampling vector
  strat_sample$Sampling_Strata <-  strat_sample$WGHTS_PROV_TRM 
  
  #Sample.Data <- Ppl.Data %>% group_by(Sampling_Strata) %>% sample_n(Sample_size)
  
  
  ## inflation weight
  strat_sample$inflation.weight<-  strat_sample$Ppl_Size_Prov / (strat_sample$Sample.size.Prov*sample.frac) #strat_sample$basic_weight
  ## analytic weight
  strat_sample <- strat_sample %>%  left_join(
    strat_sample %>% group_by(Sampling_Strata) %>% summarise(sum.analytic.weight = sum(inflation.weight))
    , by= "Sampling_Strata" )%>%  
    mutate( analytic.weight = inflation.weight/ sum.analytic.weight* Sample_size )  
  
  strat_sample
}

strat_sample.func.slim<- function( Ppl.Data2, sample.frac = 1, rand.seed=NULL){
  if(!is.null(rand.seed)){set.seed(rand.seed[1])}
  
  strat_sample <- Ppl.Data2 %>%   group_by(WGHTS_PROV_TRM) %>%  
    #slice_sample(n=Sample_size)
    sample_n(round(Sample.size.Prov*sample.frac))
  
  as.integer(strat_sample$entity_id)
}

## sample reconstruction function from entity_id
strat_sample.reconstruct<- function( Ppl.Data2, entity_id){
  
  strat_sample <- Ppl.Data2[which(Ppl.Data2$entity_id %in% as.numeric(entity_id)), ]
  
  ### create missing data
  strat_sample[ which(strat_sample$Miss.ind==1), "FVC" ] <-NA
  strat_sample[ which(strat_sample$Miss.ind==1), "FEV1" ] <-NA
  strat_sample$Gen1 <- NULL 
  strat_sample$Disease1 <- NULL
  
  strat_sample$hidden.ratio <-strat_sample$ratio
  
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  
  strat_sample$Hidden.Miss.prob <- strat_sample$Miss.prob
  
  #mean( strat_sample$ratio<0.7 , na.rm=TRUE)
  
  
  # Change  the sampling vector
  strat_sample$Sampling_Strata <-  strat_sample$WGHTS_PROV_TRM 
  
  #Sample.Data <- Ppl.Data %>% group_by(Sampling_Strata) %>% sample_n(Sample_size)
  
  strat_sample<- strat_sample %>%   group_by(Sampling_Strata) %>% mutate(Sample.size.Prov=n()) %>% ungroup()
  
  ## inflation weight
  strat_sample$inflation.weight<-  strat_sample$Ppl_Size_Prov / (strat_sample$Sample.size.Prov) #strat_sample$basic_weight
  
  ## analytic weight
  strat_sample<- strat_sample %>% group_by(Sampling_Strata) %>% mutate(sum.inflation.weight = sum(inflation.weight)) %>% ungroup()
  strat_sample$analytic.weight <-  with(strat_sample, inflation.weight / Sample.size.Prov * Sample.size.Prov)
  
  
  strat_sample
}


## Warpper function with mice 

Prev.MI.pool <-function(imp.est){
  ## between  imputation variance 
  m<- length(imp.est)
  MI.est = mean(imp.est)
  V.B <- var(imp.est)
  MI.se= sqrt(0 + V.B*(1+ 1/m  ))
  MI.df.old = (m -1)  ##  lambda  =1 # https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html#eq:riv  
  MI.df.observed  = 0             ##  lambda  =1  # (9.10)
  
  list( imp.est=imp.est,
        MI.est = MI.est,
        MI.se= MI.se,
        MI.df.old = MI.df.old ,
        MI.CI =  MI.est + qt(c(0.025,0.975), df= MI.df.old)*MI.se
  )
}


Est.func.MI<-function( strat_sample, rand.seed=NULL,var.char,   num.imp.method = "pmm" ){ 
  
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", 
  #               "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
  #               "FVC", "FEV1" )
  
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", "Smoking",
  #               "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
  #               "FVC", "FEV1" )
  
  if(!is.null(rand.seed)){
    mice.result<-mice::mice(strat_sample[,var.char], meth = c( "polyreg", rep("logreg",3), 
                                                               rep(num.imp.method,6))
                            , seed= rand.seed[1],print=FALSE )
  }else{
    mice.result<-mice::mice(strat_sample[,var.char], meth = c( "polyreg",rep("logreg",3), 
                                                               rep(num.imp.method,6))
                            ,print=FALSE )
  }
  
  # imp.vec<-  rep(NA, 5)
  # for( i in 1:5){
  #   temp.data <- mice::complete( mice.result,i)  
  #   imp.vec[i] <- mean(temp.data$FEV1/temp.data$FVC <0.7)
  # } 
  (imp.vec2<-with(mice.result, exp=mean(FEV1/FVC<0.7)))
  
  
  # https://bookdown.org/mwheymans/bookmi/rubins-rules.html  
  #  (imp.vec2<-with(mice.result, exp=  survey::svymean(  (FEV1/FVC<0.7)),   survey::svydesign(  )     )
  #  summary(mice::pool(imp.vec2))
  
  
  ## Using Rubin's rule. no within imputation error
  ## between  imputation veriance 
  imp.est<- unlist( imp.vec2$analyses )
  # MI.est = mean(imp.est)
  # V.B <- var(imp.est)
  # MI.se= sqrt(0 + V.B*(1+ 1/ mice.result$m) )
  # MI.df.old = (mice.result$m -1)  ##  lambda  =1 # https://bookdown.org/mwheymans/bookmi/measures-of-missing-data-information.html#eq:riv  
  # MI.df.observed  = 0             ##  lambda  =1  # (9.10)
  
  #   list( imp.est=imp.est,
  #        MI.est = MI.est,
  #         MI.se= MI.se,
  #         MI.df.old = MI.df.old ,
  #         MI.CI =  MI.est + qt(c(0.025,0.975), df= MI.df.old)*MI.se
  # )
  
  imp.est
}



Est.func<-function( strat_sample, rand.seed=NULL, var.char,  num.imp.method = "pmm" ){ 
  # ## sample prevalence (our target)
  #  print( paste( "sample prevalence (if no missing):",  mean( strat_sample$hidden.ratio<0.7), sep=" ")  )
  # 
  # ## naive estimation is biased 
  # strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  # print(paste(  "sample prevalence (complete case only):", mean( strat_sample$ratio<0.7 , na.rm=TRUE) ) )
  # 
  # ## with "perfect" propensity score (which is fine)
  # print( paste( "sample prevalence (with actual propensity score):", 
  # sum( (strat_sample$ratio <0.7) / (1-strat_sample$Hidden.Miss.prob),na.rm = T)/
  #   sum( (strat_sample$ratio >0 )  / (1-strat_sample$Hidden.Miss.prob),na.rm = T) ) )
  
  
  ## with "estimated" propensity score as Gen1 is hidden (which is okay but may be biased/ subjected to high variance)
  
  # Propensity.glm <- glm( Miss.ind ~ SEX_ASK_TRM  + Smoking+  WGHTS_PROV_TRM+Risk.factor
  #                        + AGE_NMBR_TRM +HWT_DHT_M_TRM + HWT_DWT_K_TRM +startlanguage_MCQ ,
  #                        data = strat_sample, family= binomial)
  
  Propensity.glm <- glm( Miss.ind ~ RIAGENDR + SMD410+RIDRETH3+RIDAGEYR+BMXHT+BMXWT+SIAINTRP,
                         data = strat_sample, family= binomial)
  
  strat_sample$Est.prob  <- predict(Propensity.glm, type="response")
  
  # print( paste( "sample prevalence (with estimated propensity score):", 
  # sum( (strat_sample$ratio <0.7) / (1-strat_sample$Est.prob),na.rm = T)/
  #   sum( (strat_sample$ratio >0 )  / (1-strat_sample$Est.prob),na.rm = T)
  # ))
  ## use MICE
  ### variable number included 
  
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM"  , "AGE_NMBR_TRM","startlanguage_MCQ", 
  #               "WEA_MRTL_CURRENT",  "Risk.factor",
  #               "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM", "Education" , "FVC", "FEV1" )
  
  
  MI.result.pooled<-list()
  
  method.len<- length(num.imp.method)
  
  i <- 1 
  
  while( i <=  method.len){
    imp.est<-Est.func.MI( strat_sample, rand.seed=rand.seed, var.char= var.char,  num.imp.method = num.imp.method[i]) 
    MI.result.pooled[[i]] <- Prev.MI.pool(imp.est)
    i<- i+1
  }
  
  
  names(MI.result.pooled)<- num.imp.method
  
  
  
  
  
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  
  list(  
    #strat_sample = strat_sample, 
    MI.result =   MI.result.pooled  ,
    Naive.result = mean( strat_sample$ratio<0.7 , na.rm=TRUE),
    Sample.result = mean( strat_sample$hidden.ratio<0.7 ),
    Propensity.est =( sum( (strat_sample$ratio <0.7) / (1-strat_sample$Est.prob),na.rm = T)/
                        sum( (strat_sample$ratio >0 )  / (1-strat_sample$Est.prob),na.rm = T)) #, 
    #Prefect.propensity.est = sum( (strat_sample$ratio <0.7) / (1-strat_sample$Hidden.Miss.prob),na.rm = T)/
    #  sum( (strat_sample$ratio >0 )  / (1-strat_sample$Hidden.Miss.prob),na.rm = T)
  )
}

### cluster by k prototype ; the naive result is after clustering 
Est.func2<-function( strat_sample, rand.seed=NULL, var.char,  num.imp.method = "pmm"){ 
  # ## sample prevalence (our target)
  #  print( paste( "sample prevalence (if no missing):",  mean( strat_sample$hidden.ratio<0.7), sep=" ")  )
  # 
  # ## naive estimation is biased 
  # strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  # print(paste(  "sample prevalence (complete case only):", mean( strat_sample$ratio<0.7 , na.rm=TRUE) ) )
  # 
  # ## with "perfect" propensity score (which is fine)
  # print( paste( "sample prevalence (with actual propensity score):", 
  # sum( (strat_sample$ratio <0.7) / (1-strat_sample$Hidden.Miss.prob),na.rm = T)/
  #   sum( (strat_sample$ratio >0 )  / (1-strat_sample$Hidden.Miss.prob),na.rm = T) ) )
  
  
  ## use MICE
  ### variable number included 
  
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM"  , "AGE_NMBR_TRM","startlanguage_MCQ", 
  #               "WEA_MRTL_CURRENT",  "Risk.factor",
  #               "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM", "Education" , "FVC", "FEV1" )
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", "Smoking",
  #               "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
  #               "FVC", "FEV1"
  #               )
  if(!is.null(rand.seed)) set.seed(rand.seed)
  ## too slow  
  #  cluster.results<-clustMixType::validation_kproto(method= "mcclain", data=data.frame(strat_sample[, var.char[1:7] ]) ) 
  
  sample.size <- dim(strat_sample)[1]
  
  ### randomly select 1000 (at most) for clustering ,  stop if wss improvement < 10% or reach 10 cluster 
  cluster.result.list<-list( )
  wss<- 1:10 ;  
  for( k in 2:10){ 
    cluster.result.list[[k]]<-clustMixType::kproto( 
      (strat_sample[ sample( sample.size , min(1000, sample.size) ) , var.char[- length(var.char)+ 0:1] ]) , k , nstart = 3
      , keep.data=FALSE ,      verbose=FALSE)
    wss[k]<-   cluster.result.list[[k]]$tot.withinss
    
    if( k > 3 ) {
      if ( (wss[k]-wss[k-1])/(wss[k-1]-wss[2]) < 0.1 ) { print(paste(" Stop at number of cluster =",k)); break;break; }
    }
  }
  
  ## optimal  k obtained   
  opt.k<- which(  wss [2:k] ==min(wss [2:k]) ) +1
  print(paste(" the selected cluster number is",opt.k));
  rm(cluster.result.list)  
  
  
  imp.able <-TRUE  #  first assume that the data is imputable 
  while( imp.able){ 
  cluster.result<- clustMixType::kproto( strat_sample[ , var.char[- length(var.char)+ 0:1]  ] ,  opt.k , nstart = 3
                                         , keep.data=FALSE ,      verbose=FALSE)
  imp.able <- FALSE  ### to break the for loop 
  
  ## 
  check.vec<-sapply( 1:opt.k,  function(x) length(mice:::find.collinear( strat_sample[ which(cluster.result$cluster==x) ,
                                                                    var.char[- length(var.char)+ 0:1]] )))
  if (max(check.vec) >=  length(var.char)-2 ) {  
    imp.able <- TRUE;  
    opt.k<- max( opt.k-1, 2) ;   ## reduce the number of cluster until imputable, bounded by 2
    print( paste("Constant cluster detected, the new opt.k is",opt.k) ) }
  
  }
  
#  cluster.label <- 1: length(cluster.result$size)
  cluster.label<-names(cluster.result$size)
  
  Singular.cluster<- which(min(cluster.result$size)<=2)   ## usually for "cluster 0" 
  ###  assign the element to other cluster randomly 
  if( length(Singular.cluster) >0) { 
    cluster.result$cluster[ which(cluster.result$cluster %in% cluster.label[Singular.cluster])  ] <-
      as.numeric ( cluster.label[-Singular.cluster][ 
        apply( rmultinom(length(Singular.cluster),1, prob=cluster.result$size[-Singular.cluster]),
               2, function(x) which(x==1))
      ] )
    ### update the labels
    cluster.result$size<- table(cluster.result$cluster)
    cluster.label<-names(cluster.result$size)
  } 
  
  
  
  naive.est<- sum(  sapply(cluster.label, function(x) with(strat_sample[which(cluster.result$cluster == x),] ,
                               mean(FEV1/FVC < 0.7, na.rm=TRUE)* cluster.result$size[x] /nrow(strat_sample) ))  )
  

  
  MI.result.pooled<-list()
  
  method.len<- length(num.imp.method)
  
  i <- 1 
  
  while( i <=  method.len){
    
    MI.result<- sapply(cluster.label, function(x) 
      Est.func.MI( strat_sample[which(cluster.result$cluster == x), var.char ], rand.seed=NULL, 
                   var.char=var.char,  num.imp.method = num.imp.method[i])  
    )
    
    ## combine the clusterd imputed result 
    imp.est <- colSums( t(MI.result) * c( cluster.result$size/sum(cluster.result$size) ) )
    MI.result.pooled[[i]] <- Prev.MI.pool(imp.est)
    i<- i+1
  }
  
  
  names(MI.result.pooled)<- num.imp.method
  
  
  ## naive estimation is biased 
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  
  list(  
    #strat_sample = strat_sample, 
    MI.result =  MI.result.pooled ,
    Cluster.size = cluster.result$size,
    Naive.result = naive.est,   
    Sample.result = mean( strat_sample$hidden.ratio<0.7 ) ### it was the true ratio for checking 
    
  )
}

## clustered by dbscan
Est.func3<-function( strat_sample, rand.seed=NULL, var.char, num.imp.method = "pmm"){ 
  # ## sample prevalence (our target)
  # print( paste( "sample prevalence (if no missing):",  mean( strat_sample$hidden.ratio<0.7) ) )
  # 
  # ## naive estimation is biased 
  # strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  # print( paste( "sample prevalence (complete case only):", mean( strat_sample$ratio<0.7 , na.rm=TRUE)) )
  # 
  # ## with "perfect" propensity score (which is fine)
  # print( paste( "sample prevalence (with actual propensity score):", 
  #        sum( (strat_sample$ratio <0.7) / (1-strat_sample$Hidden.Miss.prob),na.rm = T)/
  #          sum( (strat_sample$ratio >0 )  / (1-strat_sample$Hidden.Miss.prob),na.rm = T)) )
  
  if(!is.null(rand.seed)) set.seed(rand.seed) 
  
  ## use MICE
  ### variable number included 
  
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM"  , "AGE_NMBR_TRM","startlanguage_MCQ", 
  #               "WEA_MRTL_CURRENT",  "Risk.factor",
  #               "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM", "Education" , "FVC", "FEV1" )
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", 
  #               "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
  #               "FVC", "FEV1" 
  # )
  
  
  cat.vec <- which(   sapply(strat_sample[, var.char[- length(var.char)+ 0:1] ], class) != "numeric")
  
  temp.data <-  fpc::cat2bin(  strat_sample[, var.char[- length(var.char)+ 0:1] ],  categorical= cat.vec) 
  temp.data$data  <- apply(  temp.data$data, 2 , as.numeric )
  
  i<-1
  while ( i <=  length(var.char)-2 ){
    if ( i %in%  cat.vec) {
      temp.data$data[,  temp.data$variableinfo[[i]]$varnum] <-  
        temp.data$data[,  temp.data$variableinfo[[i]]$varnum]*
        fpc::distancefactor(temp.data$variableinfo[[i]]$ncat, dim(temp.data$data)[1],type="categorical")
    }else{
      temp.data$data[,  temp.data$variableinfo[[i]]$varnum] <-  
        scale(temp.data$data[,  temp.data$variableinfo[[i]]$varnum])
    }
    i<- i+1
  }
  
  #from https://iopscience.iop.org/article/10.1088/1755-1315/31/1/012012/pdf 
  #choosing k = 300 as we are doing MI later (roughly 20 parameters are needed )
  eps_plot <-  sort( dbscan::kNNdist(temp.data$data, k=400) )
  # slope of 1% diff 
  slope.vec<- eps_plot -dplyr::lag(eps_plot)
  slope.vec<- (slope.vec-dplyr::lag(slope.vec))/dplyr::lag(slope.vec)
  best.eps<- eps_plot[max(which( slope.vec< 0.01))]
  
  imp.able <-TRUE  #  first assume that the data is imputable 
  while( imp.able){ 
  cluster.result<- fpc::dbscan(temp.data$data, eps=  best.eps, MinPts = 400, scale = TRUE, 
                               method = "hybrid", seeds = !is.null(rand.seed), showplot = FALSE, countmode = NULL)
  
  cluster.result$size<- table(cluster.result$cluster)
  cluster.label<-names(cluster.result$size)
  
  imp.able <- FALSE  ### to break the for loop 
  
  ## 
  check.vec<-sapply( 1:length(cluster.result$size)
          ,  function(x) length(mice:::find.collinear( 
            strat_sample[ which(cluster.result$cluster == cluster.label[x])
                          ,  var.char[- length(var.char)+ 0:1] ] )))
 
  if (max(check.vec) >=  length(var.char)-2 ) {  
    imp.able <- TRUE;  
    best.eps<- best.eps*1.1 ;   ## reduce the number of cluster until imputable
    print( paste("Constant cluster detected, the new best.eps is",best.eps) ) }
  
  }
  
  Singular.cluster<- which(min(cluster.result$size)<=2)   ## usually for "cluster 0" 
  ###  assign the element to other cluster randomly 
  if( length(Singular.cluster) >0) { 
    cluster.result$cluster[ which(cluster.result$cluster %in% cluster.label[Singular.cluster])  ] <-
    as.numeric ( cluster.label[-Singular.cluster][ 
      apply( rmultinom(length(Singular.cluster),1, prob=cluster.result$size[-Singular.cluster]),
             2, function(x) which(x==1))
             ] )
    ### update the labels
    cluster.result$size<- table(cluster.result$cluster)
    cluster.label<-names(cluster.result$size)
    } 
  
  
  naive.est<- sum(  sapply(1: length(cluster.result$size), function(x) with(strat_sample[which(cluster.result$cluster == cluster.label[x]),] ,
                                                                            mean(FEV1/FVC < 0.7, na.rm=TRUE)* cluster.result$size[x] /nrow(strat_sample) ))  )
  
  
  
  MI.result.pooled<-list()
  
  method.len<- length(num.imp.method)
  
  i <- 1 
  while( i <=  method.len){
    MI.result<- sapply(cluster.label, function(x) 
      Est.func.MI( strat_sample[which(cluster.result$cluster == x), var.char ], rand.seed=NULL
                   , var.char=var.char,  num.imp.method = num.imp.method[i])  
    )
    
    ## combine the clusterd imputed result 
    imp.est <- colSums( t(MI.result) * c( cluster.result$size/sum(cluster.result$size) ) )
    MI.result.pooled[[i]] <- Prev.MI.pool(imp.est)
    i<- i+1
  }
  
  
  names(MI.result.pooled)<- num.imp.method
  
  
  
  ## naive estimation is biased 
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  
  list(  
    #strat_sample = strat_sample, 
    MI.result = MI.result.pooled ,
    Naive.result = naive.est,
    Sample.result = mean( strat_sample$hidden.ratio<0.7 )
  )
}

## clustered by random forest
library(randomForest)
library(caret)
library(e1071)

## train.sample.size 500 for 20 sec; 1000 for  45 sec 
Est.func4<-function( strat_sample, rand.seed=NULL, var.char,  train.sample.size = 500){ 
  # ## sample prevalence (our target)
  # print( paste( "sample prevalence (if no missing):",  mean( strat_sample$hidden.ratio<0.7) ) )
  # 
  # ## naive estimation is biased 
  # strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  # print( paste( "sample prevalence (complete case only):", mean( strat_sample$ratio<0.7 , na.rm=TRUE)) )
  # 
  # ## with "perfect" propensity score (which is fine)
  # print( paste( "sample prevalence (with actual propensity score):", 
  #               sum( (strat_sample$ratio <0.7) / (1-strat_sample$Hidden.Miss.prob),na.rm = T)/
  #                 sum( (strat_sample$ratio >0 )  / (1-strat_sample$Hidden.Miss.prob),na.rm = T)) )
  
  if(!is.null(rand.seed)) set.seed(rand.seed) 
  
  ## use caret to train rf https://www.guru99.com/r-random-forest-tutorial.html
  ### variable number included 
  
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM"  , "AGE_NMBR_TRM","startlanguage_MCQ", 
  #               "WEA_MRTL_CURRENT",  "Risk.factor",
  #               "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM", "Education" , "FVC", "FEV1" )
  # var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", 
  #               "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
  #               "FVC", "FEV1" 
  # )
  
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  
  trControl <- caret::trainControl(method = "cv",   number = 5)

  ###   sample  1000 for fast training 
  
  missing.vec <- which(is.na(strat_sample$FEV1) )
  
  FVC.vec <- strat_sample$FVC
  FEV1.vec <- strat_sample$FEV1
  

  rf.sample <- sample(1:dim(strat_sample)[1],min(dim(strat_sample)[1],train.sample.size))
  
  
  #https://stats.stackexchange.com/questions/115470/mtry-tuning-given-by-caret-higher-than-the-number-of-predictors
  CC.df<- data.frame( na.omit(strat_sample[rf.sample , var.char[-length(var.char)] ]))

  rf_FVC.train <- caret::train( CC.df[ , - dim(CC.df)[2]]  , CC.df[, dim(CC.df)[2]] , 
                          method = "rf", 
                          metric= "RMSE", trControl = trControl, tuneGrid = NULL)
  
  
  CC.df<- data.frame( na.omit(strat_sample[rf.sample , var.char[-length(var.char)+1] ]))
  
  rf_FEV1.train <- caret::train( CC.df[ , - dim(CC.df)[2]]  , CC.df[, dim(CC.df)[2]] , 
                          method = "rf", 
                          metric= "RMSE", trControl = trControl, tuneGrid = NULL)
  
  
  
  imp.est<- 1:5 *NA 
  
  ## remove the missing value for further imputation 
  FVC.vec[ missing.vec]<- NA
  FEV1.vec[ missing.vec] <-NA
  

  imp.iter<-0
  #https://stackoverflow.com/questions/49186277/caret-method-rf-warning-message-invalid-mtry-reset-to-within-valid-rang
  while( imp.iter < 5 ){
    imp.iter<- imp.iter + 1
  
    rf_FVC<- randomForest::randomForest(FVC~., data=na.omit(strat_sample[rf.sample
                                             , var.char[-length(var.char)] ]), mtry = rf_FVC.train$bestTune$mtry)
    
    rf_FEV1<- randomForest::randomForest(FEV1~., data=na.omit(strat_sample[rf.sample
                                                 , var.char[-length(var.char)+1] ]), mtry = rf_FVC.train$bestTune$mtry) 
    
    FVC.vec[ missing.vec]<-  predict(rf_FVC, newdata = strat_sample[ missing.vec, var.char[-length(var.char)] ])
    FEV1.vec[ missing.vec] <-  predict(rf_FEV1, newdata = strat_sample[missing.vec , var.char[-length(var.char)] ])
    
  
  imp.est[imp.iter]<-  mean( FEV1.vec/FVC.vec < 0.7  )
  
  
  FVC.vec[ missing.vec]<- NA
  FEV1.vec[ missing.vec] <-NA
  }
  # 
  # Class<- (strat_sample$FEV1/strat_sample$FVC <0.7) 
  # temp.data.pred<- as.matrix( temp.data$data[which(is.na(Class)), ] )
  # temp.data.train<-as.matrix(  temp.data$data[which(!is.na(Class)), ] ) 
  # 
  # # temp.data$Class<- (strat_sample$FEV1/strat_sample$FVC <0.7) 
  # # temp.data.pred<- temp.data[which(is.na(temp.data$Class)), ]
  # # temp.data.train<- temp.data[which(!is.na(temp.data$Class)), ]    
  # 
  # 
  # fitControl <- caret::trainControl(## 10-fold CV
  #   method = "repeatedcv",
  #   number = 5,
  #   ## repeated five times
  #   repeats = 5)
  # 
  # #https://stackoverflow.com/questions/70453557/caret-xgbtree-warning-ntree-limit-is-deprecated-use-iteration-range-instea
  # 
  # #https://xgboost.readthedocs.io/en/stable/R-package/xgboostPresentation.html
  # 
  # dtrain <- xgboost::xgb.DMatrix(data = temp.data.train ,  label =  Class[which(!is.na(Class))] )
  # 
  # bstSparse <- xgboost::xgboost(data =   dtrain , max.depth = 2, eta = 1, 
  #                               nthread = 1, nrounds = 2, objective = "binary:logistic")
  # 
  # 
  ### use the sample average as the sample is not balanced 
  # (predict(  bstSparse, temp.data.pred) >mean( Class[which(!is.na(Class))]))
  

  RF.result = Prev.MI.pool(imp.est)
  RF.result$FVC.bestTune <- rf_FVC.train$bestTune
  RF.result$FEV1.bestTune <- rf_FEV1.train$bestTune  
  
  
  ## naive estimation is biased 

  
  list(  
    #strat_sample = strat_sample, 
    RF.result = RF.result,
    Naive.result = mean( strat_sample$ratio<0.7 , na.rm=TRUE),
    Sample.result = mean( strat_sample$hidden.ratio<0.7 )
  )
}


library(nnet)
##  USE ANN from nnet 
Est.func5<-function( strat_sample, rand.seed=NULL,var.char,  train.sample.size = 500){ 

  if(!is.null(rand.seed)) set.seed(rand.seed) 

  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 

  trControl <- caret::trainControl(method = "cv",   number = 5)
  ###   sample  1000 for fast training 
  missing.vec <- which(is.na(strat_sample$FEV1) )
  
  FVC.vec <- strat_sample$FVC
  FEV1.vec <- strat_sample$FEV1
  
  
  ###   sample  500 for fast training 
  sample.size<- dim(strat_sample)[1]
  
  nnet.sample <- sample(1:sample.size,min(sample.size,train.sample.size))

  #https://stats.stackexchange.com/questions/115470/mtry-tuning-given-by-caret-higher-than-the-number-of-predictors
  ##https://cran.r-project.org/web/packages/caret/vignettes/caret.html

  fitControl <- trainControl(method = "repeatedcv", number = 3,    repeats = 3)
  
  nnetGrid <-  expand.grid(size = seq(from = 2, to = 10, by = 2),
                           decay = exp( seq(from = -5, to = -1, by = 1)) )
  
  CC.df<- data.frame( na.omit(strat_sample[nnet.sample , var.char[-length(var.char)] ]))
  nnetFit.FVC <- caret:: train(FVC ~ ., 
                                data =   CC.df,
                                method = "nnet",
                                # metric = "ROC",
                                trControl = fitControl,
                                tuneGrid = nnetGrid,
                                trace = FALSE,
                                verbose = FALSE)

  CC.df<- data.frame( na.omit(strat_sample[nnet.sample , var.char[-length(var.char)+1] ]))
  nnetFit.FEV1 <- caret:: train(FEV1 ~ ., 
                                data =   CC.df,
                                method = "nnet",
                                # metric = "ROC",
                                trControl = fitControl,
                                tuneGrid = nnetGrid,
                                trace = FALSE,
                                verbose = FALSE)
  
  ## bootstrap for MI
  
  FVC.vec <- strat_sample$FVC
  FEV1.vec <- strat_sample$FEV1
  

    
  CC.df<- data.frame( na.omit(strat_sample[    , var.char]))
  nnet.FVC <-nnet::nnet(FVC ~ . ,
                        data =  CC.df[,-length(var.char)], 
                        size = nnetFit.FVC$bestTune$size,  decay = nnetFit.FVC$bestTune$decay, maxit = 1000
                        ,trace = FALSE)
  nnet.FEV1 <-nnet::nnet(FEV1 ~ . ,
                        data =  CC.df[,-length(var.char)+1], 
                        size = nnetFit.FEV1$bestTune$size,  decay = nnetFit.FEV1$bestTune$decay, maxit = 1000
                        ,trace = FALSE)
  
  FVC.vec[missing.vec] <- predict(nnet.FVC, newdata = strat_sample[ missing.vec, var.char[-length(var.char)]  ])
  FEV1.vec[missing.vec]<- predict(nnet.FEV1,newdata = strat_sample[ missing.vec, var.char[-length(var.char)+1]])
  
 

  
  ANN.result =  list (  pt.est = mean( FEV1.vec/FVC.vec < 0.7) 
                      , FEV1.bestTune = nnetFit.FEV1$bestTune
                      , FVC.bestTune = nnetFit.FVC$bestTune 
                      )

  ## naive estimation is biased 
  
  list(  
    #strat_sample = strat_sample, 
    ANN.result = ANN.result,
    Naive.result = mean( strat_sample$ratio<0.7 , na.rm=TRUE),
    Sample.result = mean( strat_sample$hidden.ratio<0.7 )
  )
}


# system.time(
#   Est.func5( strat_sample, rand.seed=NULL, train.sample.size = 500)
# )



## estimating probabilities directly, rf ann
logit<-function(p){ log(p/(1-p))}
expit<-function(x) {1/(1+exp(-x))}
dexpit<-function(x) {exp(-x)/(1+exp(-x))^2}

predict.nnet.1 <- function(nnet.obj, train.data, new.data, wght=nnet.obj$wts ) {
  len<- length(attr(nnet.obj$terms,"variables"))
  x.var <-  as.character(attr(nnet.obj$terms,"variables")[3:len])
  old.x <- as.data.frame(train.data[,x.var])
  new.x <- new.data[,x.var]
  new.len<- nrow(new.x)
  model.x<- model.matrix( as.formula(paste( c( "~",  nnet.obj$terms[[3]]),collapse = "" )  ),
                          data= rbind(new.x,old.x) )[1: new.len, ] 
  
  #n.layer<- length(nnet.obj$n)
  hidden.n<- (prod(nnet.obj$n[1:2])+ nnet.obj$n[2])
  hidden.x <-
    model.x %*% matrix( wght[1:hidden.n], ncol= nnet.obj$n[2])
  hidden.x <- 1/(1+exp(- hidden.x))
  out.x <-  cbind(1,  hidden.x)  %*%  matrix( wght[- (1:hidden.n)], ncol= nnet.obj$n[3])
  
  1/(1+exp(-out.x))
  
}

predict.nnet.grad <- function(nnet.obj, train.data, new.data, wght=nnet.obj$wts ) {
  len<- length(attr(nnet.obj$terms,"variables"))
  x.var <-  as.character(attr(nnet.obj$terms,"variables")[3:len])
  old.x <- as.data.frame(train.data[,x.var])
  new.x <- new.data[,x.var]
  new.len<- nrow(new.x)
  model.x<- matrix( model.matrix( as.formula(paste( c( "~",  nnet.obj$terms[[3]]),collapse = "" )  ),
                                  data= rbind(new.x,old.x) )[1: new.len, ], nrow = new.len )
  
  #n.layer<- length(nnet.obj$n)
  hidden.n<- (prod(nnet.obj$n[1:2])+ nnet.obj$n[2])
  hidden.x1 <-   model.x %*% matrix( wght[1:hidden.n], ncol= nnet.obj$n[2])
  hidden.x <- 1/(1+exp(- hidden.x1))
  out.x <-  cbind(1,  hidden.x)  %*%  matrix( wght[- (1:hidden.n)], ncol= nnet.obj$n[3])
  
  #1/(1+exp(-out.x))
  
  temp.hidden  <- t( t( exp(- hidden.x1)/(1+exp(- hidden.x1))^2) *  wght[- (1:hidden.n)][-1] )
  temp.hidden1 <-  model.x*temp.hidden[,1]
  i<-2
  while( i <=  nnet.obj$n[2]) {
    temp.hidden1 <- cbind(temp.hidden1,
                          model.x*temp.hidden[,i])
    i<-i+1
  } 
  
  list( grad =   c(exp(-out.x)/(1+exp(-out.x))^2) * cbind(  
    temp.hidden1, 1,  hidden.x), 
    raw  =  1/(1+exp(-out.x))
  ) 
  
}

Prev.CI.nnet<-function( nnet.obj, train.data, new.data, resp.var="FEV1", coverage= c(0.025,0.975)){
  
  alpha<- nnet.obj$decay
  temp.J1<- predict.nnet.grad(nnet.obj,train.data=train.data, new.data= train.data)
  
  J1<-   temp.J1$grad
  
  ## (21) 
  #H1<-    J1%*%solve(t(J1)%*%J1+alpha*diag(length(nnet.obj$wts)))%*%t(J1)
  
  temp.JJ <- t(J1)%*%J1
  temp.JJ.1 <- solve(temp.JJ+alpha*diag(length(nnet.obj$wts)))
  
  # RSS <-  sum( (train.data[,c(resp.var)] - predict(nnet.obj, type="raw"))^2)
  RSS <-  sum( (train.data[,c(resp.var)] - temp.J1$raw)^2)
  ## s^2   (23)
  train.n<-dim(train.data)[1]  ## number of obs in training 
  s.2 <- RSS/( train.n - 2*sum( diag( temp.JJ.1%*%temp.JJ ) ) +
                 sum( diag( temp.JJ.1%*%temp.JJ%*%temp.JJ.1%*%temp.JJ ) )  )
  
  ## var(y|x) (20)
  temp.g<- predict.nnet.grad(nnet.obj,train.data=train.data,  new.data= new.data, wght=nnet.obj$wts )
  g0<- temp.g$grad
  
  nn <- dim(new.data)[1]
  ## var(\bar y|x) = sum (var(y|x))/ nn^2 (2) 
  mean.g<-  colMeans(g0[1:nn,])
  var.avg.y<- s.2 * ( 1/nn  + mean.g%*%temp.JJ.1%*%temp.JJ%*%temp.JJ.1%*%(mean.g))
  sd.avg.y<- sqrt(var.avg.y)
  
  # Prev.nnet.noweight <- predict(nnet.obj, newdata = new.data, type="raw")
  
  Prev.nnet.noweight <- temp.g$raw
  
  mean(Prev.nnet.noweight) + qnorm(coverage)* c(sd.avg.y)
  
}

### because of the imbalanced data, direct estimation of COPD may result in biased result 
## train.sample.size 500 for 20 sec; 1000 for  45 sec 
Est.func6<-function( strat_sample, rand.seed=NULL, var.char, train.sample.size = 500){ 
  if(!is.null(rand.seed)) set.seed(rand.seed) 
  
  strat_sample$ratio <- strat_sample$FEV1/strat_sample$FVC 
  strat_sample$COPD  <- as.factor(c("NO","YES")[1+as.numeric(strat_sample$ratio < 0.7)])

  ###   sample  1000 for fast training 
  missing.vec <- which(is.na(strat_sample$COPD) )
  
  sample.size <- dim(strat_sample)[1] 
  train.sample <- sample(1:sample.size, min(sample.size,train.sample.size))
  
  #https://stats.stackexchange.com/questions/115470/mtry-tuning-given-by-caret-higher-than-the-number-of-predictors
  CC.df<- data.frame( na.omit(strat_sample[train.sample , c(var.char[-length(var.char)+0:1],"COPD") ]))
  
  trControl <- caret::trainControl(method = "cv",   number = 5, classProbs = TRUE)

  rf_COPD.train <- caret::train( CC.df[ , - dim(CC.df)[2]]  , CC.df[, dim(CC.df)[2]] , 
                                method = "rf", 
                                metric= "Accuracy", trControl = trControl, tuneGrid = NULL)

  COPD.vec <- strat_sample$COPD  
  ## remove the missing value for further imputation 
  COPD.vec[ missing.vec]<- NA
  imp.est<- 1:5 *NA 
  imp.iter<-0
  #https://stackoverflow.com/questions/49186277/caret-method-rf-warning-message-invalid-mtry-reset-to-within-valid-rang
  while( imp.iter < 5 ){
    imp.iter<- imp.iter + 1
    
    rf_COPD <- randomForest::randomForest(COPD~., data=na.omit(strat_sample[, c(var.char[-length(var.char)+0:1],"COPD")])
                                          , mtry = rf_COPD.train$bestTune$mtry)
    
    COPD.vec[ missing.vec]<-  predict(rf_COPD, newdata = strat_sample[ missing.vec, var.char[-length(var.char)-0:1]])

    imp.est[imp.iter]<-  mean( COPD.vec =="YES")

    COPD.vec[ missing.vec]<- NA
  }
 
  
  RF.Result = Prev.MI.pool(imp.est)
  RF.Result$bestTune <- rf_COPD.train$bestTune
 
  ####################################################################################
  # trControl <- caret::trainControl(method = "cv",   number = 5, classProbs = TRUE)
  # fitControl <- trainControl(method = "repeatedcv", number = 3,    repeats = 3)
  
  nnetGrid <-  expand.grid(size = seq(from = 2, to = 10, by = 2),
                           decay = exp( seq(from = -3, to = 3, by = 1)) )
 
  nnetFit.COPD <- caret:: train(COPD ~ ., 
                               data =   CC.df,
                               method = "nnet",
                               # metric = "ROC",
                               trControl = trControl,
                               tuneGrid = nnetGrid,
                               trace = FALSE,
                               verbose = FALSE)
  
  strat_sample$COPD  <-  as.numeric(strat_sample$ratio < 0.7)
  
  nnet.COPD <-nnet::nnet(COPD ~ . ,
                        data =  na.omit(strat_sample[, c(var.char[-length(var.char)+0:1],"COPD")]), 
                        size = nnetFit.COPD$bestTune$size,  decay = nnetFit.COPD$bestTune$decay, maxit = 1000
                        ,trace = FALSE)
  
  COPD.vec1 <-  as.numeric(strat_sample$ratio < 0.7)
  
  COPD.vec1[missing.vec] <- NA
  
  COPD.vec1[missing.vec] <- predict(nnet.COPD , newdata = strat_sample[ missing.vec, var.char[-length(var.char)-0:1]]) 
  
  
  miss.w <- length(missing.vec)/ sample.size
  
  Prev.CI <- Prev.CI.nnet( nnet.COPD, train.data = na.omit(strat_sample[, c(var.char[-length(var.char)+0:1],"COPD")])
              , new.data=strat_sample[ missing.vec, var.char[-length(var.char)-0:1]]
              , resp.var="COPD", coverage= c(0.025,0.975))
  
  
  Prev.CI2 <- Prev.CI *  miss.w  +  mean( strat_sample$ratio<0.7 , na.rm=TRUE) *(1- miss.w )
  
  NNet.Result<- list( pt.est =  mean( COPD.vec1)
                     ,Prev.CI = Prev.CI2  
                     ,bestTune= nnetFit.COPD$bestTune)
 #####################################################################################################################
  
  strat_sample$COPD  <- as.numeric(strat_sample$ratio < 0.7)
  logit.model <- glm( COPD ~ ., 
                      family=binomial(link = "logit"), data=strat_sample[, c(var.char[-length(var.char)+0:1],"COPD")])
  
  temp<- model.matrix( ~ ., data=strat_sample[missing.vec, c(var.char[-length(var.char)+0:1])])
  
  ## the point estimate by averaging predicted probability
  (est1.pt <- mean( expit( temp%*%coef(logit.model))))
  
  ## SE by delta's method 
  #sqrt(diag(temp%*%vcov(logit.model)%*%t(temp))) ## std error of the linear predictor 
  
  (est1.se <-c(sqrt(colMeans( c(dexpit(temp%*%coef(logit.model)))*temp)%*%vcov(logit.model)  %*%  colMeans( c( dexpit(temp%*%coef(logit.model))) * temp ) )))
  
  #Not this one (est1.se <- as.numeric( sqrt( colSums(temp)%*%vcov(logit.model)%*%colSums(temp))) ) 
  
  ## the Wald CI of the prevalence
  (est1.CI <-  (est1.pt+ qnorm(p=c(0.025,0.975))*est1.se ) )
  
  
  
  Logit.Result<- list( pt.est =  est1.pt* miss.w +  mean( strat_sample$ratio<0.7 , na.rm=TRUE) *(1- miss.w )
                      ,se.est =  est1.se
                      ,Prev.CI = est1.CI* miss.w +  mean( strat_sample$ratio<0.7 , na.rm=TRUE) *(1- miss.w )  
                      ,miss.w=  miss.w)
  #####################################################################################################################
  
  strat_sample$non.Miss.ind <- 1- strat_sample$Miss.ind
  COPD.vec1 <-  as.numeric(strat_sample$ratio < 0.7)
  
  Propensity.model <- glm( non.Miss.ind ~ ., 
                      family=binomial(link = "logit"), 
                      data=strat_sample[, c(var.char[-length(var.char)+0:1],"non.Miss.ind")])
  
  temp<- model.matrix(non.Miss.ind   ~ ., 
                      data=strat_sample[which(COPD.vec1==1), 
                                        c(var.char[-length(var.char)+0:1],"non.Miss.ind")])
  
  ## the point estimate by averaging predicted probability
  Prop.Score <- ( expit( temp%*%coef( Propensity.model )))
  
# switched to an easier method 
#  Prop.est <-   sum(strat_sample[-missing.vec,"COPD"]/Prop.Score)/sum(1/Prop.Score)
  (Prop.est <-   sum(1/Prop.Score)/sample.size)
  
  
  ## SE by delta's method 
  #sqrt(diag(temp%*%vcov(logit.model)%*%t(temp))) ## std error of the linear predictor 
  
  temp.vec <-  -colMeans(temp * c(exp(- temp %*% coef(Propensity.model))))* (dim(temp)[1] /sample.size)
  
  (Prop.se <-c(sqrt( temp.vec %*%  vcov( Propensity.model) %*% temp.vec )))
  
  #Not this one (est1.se <- as.numeric( sqrt( colSums(temp)%*%vcov(logit.model)%*%colSums(temp))) ) 
  
  ## the Wald CI of the prevalence
  (Prop.CI <-  (Prop.est + qnorm(p=c(0.025,0.975))*Prop.se  ) )
  
  Propensity.Result<- list( pt.est = Prop.est
                       ,se.est =  Prop.se 
                       ,Prev.CI = Prop.CI 
                       ,COPD.cases=  dim(temp)[1] )
  
  ##################################################################################################################### 
  
  
  
  ## naive estimation is biased 
  
  list(  
    #strat_sample = strat_sample, 
    RF.COPD.result = RF.Result,
    ANN.COPD.result = NNet.Result,
    Logit.COPD.Result =  Logit.Result,
    Prop.COPD.Result =  Propensity.Result,
    Naive.result = mean( strat_sample$ratio<0.7 , na.rm=TRUE),
    Sample.result = mean( strat_sample$hidden.ratio<0.7 )
  )
}



##################################################################################
#strat_sample  wrap up function 

Strat_sample.func<- function(Ppl.Data2, sample.frac = 0.5,  
                             file.name="Miss_Med_Size_Med_Strat_Sample.csv", iter=2, rand.seed=NULL){
  if(!is.null(rand.seed)) set.seed(rand.seed) 
  strat_sample.list<-list()
  i<-1
  while( i <= iter){ 
    strat_sample.list[[i]]<-strat_sample.func.slim( Ppl.Data2= Ppl.Data2, sample.frac = sample.frac, rand.seed=NULL)#$entity_id
    i<-i+1
  }
  
  strat_sample.df<- do.call("rbind" ,(strat_sample.list))
  #save(strat_sample.list, file=file.name)
  write.csv(strat_sample.df, file=file.name)
}

#simulation wrap up function 
Simulation.wrap <-  function(strat_sample, rand.seed=NULL, var.char, num.imp.method = c( "pmm", "norm.nob" ), 
                             file.prefix= "Result",simindex=1){
  if(!is.null(rand.seed)) set.seed(rand.seed) 
  
  Result<- list()
  
  Result[[1]] <- Est.func( strat_sample=strat_sample, rand.seed=NULL, var.char=var.char,  num.imp.method = num.imp.method )
  Result[[2]] <-Est.func2( strat_sample=strat_sample, rand.seed=NULL, var.char=var.char,  num.imp.method = num.imp.method )
  Result[[3]] <-Est.func3( strat_sample=strat_sample, rand.seed=NULL, var.char=var.char,  num.imp.method = num.imp.method )
  
  
  i<- 1 
  while ( i  <= length(Result)){
    Est.name<- names(Est.vec<- c(simindex, unlist(Result[[i]]) ))
    Est.name[1]<-"Sim.index"
    names(Est.vec)<- Est.name
    write.table( t(Est.vec), file=paste(file.prefix,i,"Est.csv",sep=""), sep="," , append = TRUE, row.name =FALSE,
                 col.names = (Est.vec[1]==1) )
    #write.csv( file=paste(file.prefix,i,"Name.csv",sep=""), Est.name, append = TRUE, row.name =FALSE)
    i<-i+1
  }
  
}
#simulation wrap up function 
Simulation.wrap2 <-  function(strat_sample, rand.seed=NULL, var.char, num.imp.method = c( "pmm", "norm.nob" ), 
                             file.prefix= "Result",simindex=1){
  if(!is.null(rand.seed)) set.seed(rand.seed) 
  
  Result<- list()
  
  Result[[1]] <-Est.func4( strat_sample=strat_sample, rand.seed=NULL, var.char=var.char,  train.sample.size = 500 )
  Result[[2]] <-Est.func5( strat_sample=strat_sample, rand.seed=NULL, var.char=var.char,  train.sample.size = 500 )
  Result[[3]] <-Est.func6( strat_sample=strat_sample, rand.seed=NULL, var.char=var.char,  train.sample.size = 500 )
  
  
  i<- 1 
  while ( i  <= length(Result)){
    Est.name<- names(Est.vec<- c(simindex, unlist(Result[[i]]) ))
    Est.name[1]<-"Sim.index"
    names(Est.vec)<- Est.name
    write.table( t(Est.vec), file=paste(file.prefix,i+3,"Est.csv",sep=""), sep="," , append = TRUE, row.name =FALSE,
                 col.names = (Est.vec[1]==1) )
    #write.csv( file=paste(file.prefix,i,"Name.csv",sep=""), Est.name, append = TRUE, row.name =FALSE)
    i<-i+1
  }
  
}

##########################Other functions for disease missing and relationship ########################################

hid.ratio<- function(data= Ppl.Data2) { c(0.75, 0.65)[data$Disease1 + 1] }  
FEV1FVC.func <- function(hid.ratio, FVC ){ hid.ratio* FVC  + rnorm(length(FVC))* 0} 
Miss.prob.func<-function(x, data=Ppl.Data2 ){  with(data,  (x*( 0.05  + (Env1==1)*0.1  + (Env1==1)*(Gen1==1)*0.15 +
                                                                  +(Gen1==1)*0.1 +  HWT_DHT_M_TRM* 0.01 +#(Gen1==2)*0.15 +
                                                                  (SEX_ASK_TRM=="M")*0.1 + 0.1*Smoking)) )}
Disease.prob.func<-function( data=Ppl.Data2 ){  with(data,  
                                                     ( ( 0.05  + (Env1==1)*0.1  + (Env1==1)*(Gen1==1)*0.15 +
                                                           (Gen1==1)*0.1 +  #(Gen1==2)*0.15 +
                                                           (SEX_ASK_TRM=="M")*0.1 + 0.1*Smoking)) ) }


####################################### Generating Samples details from Population data###########################
##Constructing the strata 
Prov.vec   <- c( "AB", "BC", "MB","NB","NL","NS","ON","PE","QC","SK")
Prov.DCS   <- c( "AB", "BC", "MB",     "NL","NS","ON",     "QC") 
Prov.NonDCS<- c( "NB", "PE", "SK")     
Prov.prop  <- c(71000,79000,32000,14000,11000,16000,82000,2800,2e+05,12000)

Age.grp    <- c("45-54","55-64","65-74","75-85")
Sex.grp    <- c("female","male")


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


# temp.vec <- table(Ppl.Data$Sampling_Strata)
# Sampling.mat<-merge(Sampling.mat
#                     , data.frame(Sampling_Strata=names(temp.vec ), N_h=c(temp.vec) )
#                     , by="Sampling_Strata" ) 



var.char <- c("WGHTS_PROV_TRM", "SEX_ASK_TRM","startlanguage_MCQ", "Smoking",
              "AGE_NMBR_TRM", "HWT_DHT_M_TRM" , "HWT_DWT_K_TRM",  "Risk.factor",
              "FVC", "FEV1" 
)



