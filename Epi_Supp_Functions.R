##############
# Main estimation functions for  "Far from MCAR: 
#   obtaining population-level estimates of HIV viral suppression"
# #
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Primary Statistician for SEARCH
###############

Run.Supp.Xsect <- function(timept, SL.library='glm'){
  
  data.input <- preprocess.cascade()
  
  # define target population for the X-sectional analysis
  restrict<- get.pop.Xsect(data=data.input, timept=timept)
  data.Xsect <- subset(data.input, restrict)
  
  # get the key indicators for all analyses (HIV, pDx, ART, Supp)
  indicators <- preprocess.serial(data = data.Xsect, timept=timept)
  
  # community is the unit of independence
  id <- data.Xsect$id
  
  # unadjusted
  unadjust <- do.unadjusted(data=indicators, id=id)
  
  # npmle (implemented as IPW)
  npmle <- do.dynamic(OC=data.Xsect, data=indicators, id=id)
  
  # TMLE fully stratified on community
  
  clusters <- unique(id)
  nClust <- length(clusters)
  
  tmle <- data.frame(matrix( NA, nrow=nClust, length(unadjust)))
  
  for(j in 1:nClust){
    
    these.units <- id==clusters[j]

    tmle[j,] <- do.with.detQ(OC= data.Xsect[these.units ,], 
                                    data = indicators[these.units, ],
                                   timept=timept, 
                                    SL.library= SL.library)
    print(j)                
  }
  colnames(tmle)  <- c('N', 'nPos', 'nTstVL',
                              'prev.pt', 'prev.CI.lo', 'prev.CI.hi',
                              'supp.pos.pt', 'supp.pos.CI.lo', 'supp.pos.CI.hi')
  
  # average across communities
  tmle <- get.clust.CI.Xsect(data=tmle)

  
  est <- rbind(unadjust,  npmle,  unlist(tmle))
  rownames(est) <- c('unadjust', 'npmle',  'tmle')
  est
  
}


preprocess.cascade <- function(){
  
  load("outputs-withIntOnly.RData")
  print( dim(outputs) )
  data.input <- outputs
  
  # Exclusions
  data.input <- subset(data.input, !(data_flag | dead_0 | move_0) )
  data.input <- subset(data.input, intervention)

  #	add community number
  id<-  rep(NA, nrow(data.input ))
  comm <- unique(data.input$community_name)
  for(j in 1:32){
    these.units <- data.input$community_name==comm[j]
    id[these.units] <- j
  }	
  data.input <- cbind(data.input, id=id)
  
  # transform pairs to be numeric
  data.input$pair <- as.numeric(as.character(data.input$pair))
  
  print('***preprocessing done***')
  
  data.input
}



get.pop.Xsect <- function(data, timept){
  
  n<- nrow(data)
  
  # have to be 15+, cannot have died or outmigrated at t 
  if(timept==0){
    adult <- data$age_0>14
    alive <-  rep(T, n)
    move <- rep(F,n)
  } else if(timept==1){
    adult <- data$age_0>13
    alive <- !data$dead_1
    move <- data$move_1
  } else if (timept==2){
    adult <- data$age_0>12
    alive <- !data$dead_2
    move <- data$move_2
  } else if (timept==3){
    adult <-  rep(T, n) 
    dead <- move <- rep(F, n)
    dead[which(data$dead_3)] <- T
    alive <- !dead
    move[which(data$outmigrate_3)] <-T
  }
  
  # allow for inmigration 
  if(timept==0){
    resident <- data$resident_0
  } else if (timept==1){ # Y1
    resident <- data$resident_recensus_1   
  } else if(timept==2){
    resident <- data$resident_recensus_2
  } else if(timept==3){ 
    resident <- data$resident_3
  }
  
  restrict<-  adult & alive & !move  & resident 
  restrict
}

preprocess.serial <- function(data, timept){
  
  if(timept==0 | timept==3 ) {
    # at baseline or endline
    HIV.variable <- 'hiv'
    pDx.variable <- 'pdx_vl'
    eART.variable <- 'eart_vl'
    }  else{
    # interim data (intervention arm only)
    HIV.variable	<- 'all_hiv'
    pDx.variable <- 'hiv_preCHC'
    eART.variable <- 'i_eart_vl'
  }
  print(c(HIV.variable, pDx.variable, eART.variable))
  
  n <- nrow(data)
  
  pDx <- eART <- Delta  <- HIVpos <- TstVL <- Supp <- rep(0, n)
  
  # no evidence of prior Dx or ART == fail
  pDx[ which(data[, paste(pDx.variable, timept, sep='_')] ==1) ] <-1 
  eART[ which(data[, paste(eART.variable, timept, sep='_')] ==1)] <-1 
  # if on ART then previously dx-ed
  pDx[ eART==1] <- 1
  
  # health fair attendance
  chc <- as.numeric(data[, paste('chc', timept, sep='_')]	)
  # follow-up home-based teting
  tr <- as.numeric(data[, paste('tr', timept, sep='_')])
  
  # did we actually know your HIV status?
  hiv.temp <-  data[, paste(HIV.variable, timept, sep='_')]	
  TstHIV <- as.numeric( !is.na(hiv.temp) )
  # Delta - require that we saw you at t & had a known status
  Delta[ which( (chc==1 | tr==1)  & TstHIV==1 )] <-1
  
  # HIV status as =1 if HIV+ and 0 otherwise 
  HIVpos[ which(hiv.temp) ] <- 1
  
  # Suppression - VL.variable='supp' 
  supp.temp <- data[, paste('supp', timept, sep='_')]
  # Set Supp=NA if not HIVpos
  supp.temp[ which( !is.na(supp.temp) & HIVpos==0) ] <- NA
  TstVL[ which(!is.na(supp.temp) ) ] <- 1
  
  #  suppression as =1 if suppressed and 0 otherwise
  Supp[ which(supp.temp) ] <- 1
  
  data.frame(cbind(pDx, eART, chc, tr, TstHIV, Delta, HIVpos, TstVL, Supp))
}



#*====

do.unadjusted <- function(data, id=NULL){
  
  # unadjusted prevalence estimates among those testing at CHC/tracking
  prev <- call.ltmle(A = data$Delta, Y = data$HIVpos, id=id, do.tmle=T)$est
  # unadj equiv: sum(data$Delta*data$HIVpos)/sum(data$Delta)
  colnames(prev) <- paste('prev', colnames(prev), sep='.')
  
  A <- data$TstVL
  Y <- data$Supp
  supp.pos <-  call.ltmle(A = A, Y = Y, id=id, do.tmle=T)$est
  
  colnames(supp.pos) <- paste('supp.pos', colnames(supp.pos), sep='.')
  
  N <- nrow(data)
  est <- data.frame(cbind(N=N, nPos=round(N*prev$prev.pt), nTstVL = sum(A),
                          prev, supp.pos))
  
  est  
}


do.dynamic <- function(OC, data, id=NULL){
  
  # get baseline strata: age-group, sex, community
  pred.A <- cbind(get.X(data=OC, adj.full=F), comm=droplevels(OC$community_name))
  
  #  stratify on age/sex/community
    
  # prevalence 
  gform <- 'A ~ age.20.29*male*comm + age.30.39*male*comm + age.40.49*male*comm + age.50.59*male*comm + age.60.plus*male*comm'
  Prob.HIVpos <- call.ltmle(pred.A = pred.A,  A = data$Delta, Y = data$HIVpos, 
                              id=id, do.tmle=F, gform=gform)
    
  # supppresion coded as longitudinal dynamic regime (1) all tested; (2) get VL if HIV+
  abar <- matrix(nrow=nrow(data), ncol=2)
  abar[,1] <- 1
  abar[,2] <- data$HIVpos>0
  data.temp <- data.frame(pred.A, data[,c('Delta', 'HIVpos','TstVL','Supp')])
  gform <- c(
      'Delta ~ age.20.29*male*comm + age.30.39*male*comm + age.40.49*male*comm + age.50.59*male*comm + age.60.plus*male*comm',
      'TstVL ~ HIVpos*male*comm*age.20.29 + HIVpos*male*comm*age.30.39 + HIVpos*male*comm*age.40.49 + HIVpos*male*comm*age.50.59 + HIVpos*male*comm*age.60.plus'
    )
    
  est.temp <- ltmle(data=data.temp, Anodes=c('Delta', 'TstVL'), Lnodes='HIVpos', Ynodes='Supp',
                    abar=abar, SL.library=NULL, stratify=T, 
                    estimate.time=F,
                    variance.method='ic',  id=id, iptw.only=T,
                    gform=gform)
  
  Prob.supp <- get.iptw.from.ltmle(est.temp)
  
  compile.results(N=nrow(data), nTstVL = sum(data.temp$TstVL), 
                  Prob.HIVpos=Prob.HIVpos, Prob.supp=Prob.supp)
  
}

#*==
# TMLE
do.with.detQ <- function(OC, data,
                         timept, SL.library = NULL, verbose = F) {
  
  
  baseline.pred <-  get.X(data=OC, adj.full=T) 
  
  
  #**************** Prevalence of HIV at timept t: P(Y*=1)
  Prob.HIVpos <- get.prevalence(OC=OC, 
                                baseline.pred=baseline.pred, 
                                data=data, timept=timept, 
                                SL.library=SL.library,verbose=verbose)
  
  
  #****************Suppression at t:  P(Supp*=1, Y*=1) 
  Prob.supp <- get.supp(OC=OC,
                        baseline.pred=baseline.pred,
                        data=data, timept=timept,
                        SL.library=SL.library,
                        verbose=verbose) 
  
  #************************** Compiling the results
  # via Delta Method
  compile.results(N=nrow(OC), nTstVL=Prob.supp$nTstVL, Prob.HIVpos=Prob.HIVpos, Prob.supp=Prob.supp$est)
  
}



#**************** Prevalence of HIV at timept t: P(Y*=1)
get.prevalence <- function(OC, baseline.pred, data,
                           timept, SL.library, 
                           id=NULL, verbose=F){
  
  # outcome as observed HIV status 
  Y <- data$HIVpos 
  # intervention variable: seen at CHC/tracking with HIV status known
  A <- data$Delta 
  
  # specify adjustment set other than baseline covariates
  adj <- get.adjustment.prev(OC = OC, baseline.pred = baseline.pred,
                             timept = timept)
  
  # using deterministic knowledge about known HIV status at prior timept point 
  est <- call.ltmle(pred.A = adj$adj,  A = A, Y = Y, 
                    SL.library = SL.library, id=id,
                    deterministicQ = adj$detQ, do.tmle=T)
  
  est
  
}

# get.adjustment.prev: function to get the adjustment set for prevalence
get.adjustment.prev <- function(OC, baseline.pred, timept){
  
  if(timept==0){
    adj <- baseline.pred
    detQ<- NULL
  } else { 
    # baseline data
    data_0 <- preprocess.serial(data = OC, timept=0)
    TstHIV_0 <- data_0$TstHIV 
    
    if(timept==1 | timept==3 ) {
      # at Y1 or Y3 (not using interim data to be parallel across arms)
      adj <- data.frame(baseline.pred, TstHIV_0, 
                        detQ.variable=data_0$HIVpos)
      
    } else if (timept==2){
      # if timept=2, adjust for baseline and timept=1
      data_1 <- preprocess.serial(data = OC, timept=1)$HIVpos
      adj <- data.frame(baseline.pred, TstHIV_0, 
                        detQ.variable=data_1 )
    } 
    detQ <- deterministicQ_YES
  }
  list(adj=adj, detQ=detQ)
  
  
}




#****************Suppression at t:  P(Supp*=1, Y*=1) 
get.supp <- function(OC, baseline.pred, data, timept,
                     SL.library, id=NULL, verbose=F){  

  A <- data$TstVL 
  Y <- data$Supp*A
  
  # get adjustment set using deterministic knowledge (primary)
  pred.A <- get.adjustment.supp(OC = OC, data=data,
                                baseline.pred = baseline.pred,
                                timept = timept)
  
  est <- call.ltmle(pred.A = pred.A, A = A, Y = Y, 
                    SL.library = SL.library, id=id, 
                    deterministicQ = deterministicQ_NO, 
                    do.tmle=T)
  
  list(est=est, nTstVL=sum(A))
}  



get.adjustment.supp <- function(OC, data,  baseline.pred, timept){
  
  if(timept==0){
    adj <- baseline.pred
  } else {
    data_0 <- preprocess.serial(data = OC, timept=0)
    TstHIV_0 <- data_0$TstHIV
    
    if( timept==1 | timept==3 ) {
      # at Y1 or Y3 (not using interim data to be parallel across arms)
      adj <- data.frame(baseline.pred, TstHIV_0,
                        prior=data_0$Supp)
      
    }else if (timept==2) {
      # if timept=2, adjust for baseline and timept=1
      data_1 <- preprocess.serial(data = OC, timept=1)
      
      adj <- data.frame(baseline.pred, TstHIV_0,
                        prior=data_1$Supp)
    }
  } 
  
  # all adjust for deterministic information on ART start 
  adj <- data.frame(adj, detQ.variable=data$eART)
  adj
  
}

#*-------

# call.ltmle: function to call the ltmle 
call.ltmle<- function(pred.A=NULL, A,  Y, 
                      SL.library=NULL,
                      deterministicQ=NULL, 
                      observation.weights=NULL, 
                      id=NULL, 
                      verbose=F, do.tmle=T,
                      Qform=NULL, gform=NULL){
  
  
  # create temporary data frame
  if( is.null(pred.A) ){
    data.temp<- data.frame(A, Y)
  } else{
    data.temp<- data.frame(pred.A, A, Y)
  }
  
  est.temp<- ltmle(data=data.temp, Anodes='A', Lnodes=NULL, Ynodes='Y',
                   abar=1,
                   stratify = T, SL.library=SL.library,
                   estimate.time=F,
                   variance.method='ic', 
                   deterministic.Q.function= deterministicQ,
                   observation.weights=observation.weights, id=id,
                   Qform=Qform, gform=gform)
  
  if(do.tmle){
    IC<- est.temp$IC$tmle
    est<- data.frame(pt=est.temp$estimate["tmle"], 
                     CI.lo=summary(est.temp)$treatment$CI[1], 
                     CI.hi=summary(est.temp)$treatment$CI[2] 	)
    OUT <-   list(est=est,IC=IC)
  }else{
    OUT <- get.iptw.from.ltmle(est.temp)
  }
  OUT
}

get.iptw.from.ltmle <- function(est.temp){
  IC<- est.temp$IC$iptw
  est<- data.frame(pt=est.temp$estimate["iptw"], 
                   CI.lo=summary(est.temp,'iptw')$treatment$CI[1], 
                   CI.hi=summary(est.temp,'iptw')$treatment$CI[2] 	)
  list(est=est,IC=IC)
}






