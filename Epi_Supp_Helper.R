##############
# Main helper functions for  "Far from MCAR: 
#   obtaining population-level estimates of HIV viral suppression"
# #
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Primary Statistician for SEARCH
###############

get.count.per <- function(this, sf=0){
  paste(sum(this), ' (', round(mean(this)*100,sf), ')', sep='')
}


# TABLE 1
get.consort <- function(){
  
  data <- preprocess.cascade()
  n<- nrow(data)
  get.measures <- function(data, timept){
    data <- subset(data,  get.pop.Xsect(data=data, timept=timept))
    print(nrow(data))
    indicators <- preprocess.serial(data=data, timept=timept)
    data.frame(
      HIV.fair =sum(indicators$Delta),
      HIV.pos = sum(indicators$Delta & indicators$HIVpos),
      VL.known = sum(indicators$HIVpos & indicators$TstVL & indicators$Delta),
      VL.supp= sum(indicators$Supp & indicators$TstVL & indicators$HIVpos & indicators$Delta)
    )
  }
  
  base <- data.frame(
    N=sum(data$age_0>=15 & data$resident_0),
    get.measures(data, timept=0)
  )
  
  yr1 <- data.frame(
    N=sum(data$age_0>=14 & data$resident_recensus_1),
    get.measures(data, timept=1)
  )
  
  
  yr2 <- data.frame(
    N=sum(data$age_0>=13 & data$resident_recensus_2),
    get.measures(data, timept=2)
  )
  
  yr3 <- data.frame(
    N=sum(data$resident_3),
    get.measures(data, timept=3)
  )
  
  
  coverage <- rbind(base, yr1, yr2, yr3)
  coverage
}



get.table1<- function(data, timept, subset.indicator, demographics){
  
  indicators <- preprocess.serial(data = data, timept=timept)
  
  if(subset.indicator=='All'){
    these <- rep(T, nrow(indicators))
  } else if(subset.indicator=='HIVmeasure'){
    these <- indicators$Delta
  } else if(subset.indicator=='HIVpos'){
    these <-  indicators$HIVpos & indicators$Delta 
  } else if(subset.indicator=='VLmeas'){
    these <- indicators$TstVL & indicators$HIVpos & indicators$Delta
  } else if(subset.indicator=='VLno'){
    these <- !indicators$TstVL & indicators$HIVpos & indicators$Delta
  } else if(subset.indicator=='VLsupp'){
    these <- indicators$Supp & indicators$TstVL & indicators$Delta
  }
  indicators <- subset(indicators, subset=as.logical(these))
  
  data <- subset(data, subset=as.logical(these))
  
  esupp <- rep(F, nrow(data))
  if (timept>0){
    supp0<- preprocess.serial(data = data, timept=0)[,'Supp']
    esupp[which(supp0==1)] <- 1
  }
  if(timept>1){
    supp1<- preprocess.serial(data = data, timept=1)[,'Supp']
    esupp[which(supp1==1)] <- 1
  }
  if(timept>2){
    supp2<- preprocess.serial(data = data, timept=2)[,'Supp']
    esupp[which(supp2==1)] <- 1
  }  

  if(demographics){
    baseline.pred <-  get.X(data=data, adj.full=T)
    table1 <- rbind( 
      N = nrow(data),
      region = get.count.per(is.na(data$region_name)),
      SW.U = get.count.per(data$region_name=='Western Uganda'),
      E.U = get.count.per(data$region_name=='Eastern Uganda'),
      K = get.count.per(data$region_name=='Kenya'),
      Female = get.count.per(!data$sex_0),
      Male = get.count.per(data$sex_0),
      age = get.count.per(is.na(data$age_0)),
      age.15.20 = get.count.per(data$age_0 < 20),
      age.20.29 = get.count.per(data$age_0>19 & data$age_0<30),
      age.30.39 = get.count.per(data$age_0>29 & data$age_0<40),
      age.40.49 = get.count.per(data$age_0>39 & data$age_0<50),
      age.50.59 = get.count.per(data$age_0>49 & data$age_0<60),
      age.60.plus = get.count.per(data$age_0>59),
      #reference is missing
      marital=get.count.per(is.na(data$marital_0)),
      single = get.count.per(baseline.pred$single),
      married = get.count.per(baseline.pred$married),
      #   widowed = get.count.per(baseline.pred$widowed),
      w.divorced.separated = get.count.per(baseline.pred$divorced.separated | baseline.pred$widowed),
      
      # education: reference is less than primary or missing
      education=get.count.per(is.na(data$education_0)),
      less.primary = get.count.per(!is.na(data$education_0) & !data$edu_primary_0 & !data$edu_secondary_plus_0),
      primary = get.count.per(data$edu_primary_0),
      secondary.plus = get.count.per(data$edu_secondary_plus_0),
      
      # occupation: reference NA
      occupation = get.count.per(is.na(data$occupation_0)),
      formal.hi = get.count.per(baseline.pred$formal.hi),
      informal.hi = get.count.per(baseline.pred$informal.hi),
      informal.lo = get.count.per(baseline.pred$informal.lo),
      other = get.count.per(!is.na(data$occupation_0) &
                              !baseline.pred$formal.hi &
                              !baseline.pred$informal.hi &
                              !baseline.pred$informal.lo &
                              !baseline.pred$jobless),
      jobless = get.count.per(baseline.pred$jobless),
      
      # reference wealth is NA missing
      wealth = get.count.per(is.na(data$wealth_0)),
      wealth0 = get.count.per(baseline.pred$wealth0),
      wealth1 = get.count.per(baseline.pred$wealth1),
      wealth2 = get.count.per(baseline.pred$wealth2),
      wealth3 = get.count.per(baseline.pred$wealth3),
      wealth4 = get.count.per(baseline.pred$wealth4)
    )
    
  } else{
    table1<- rbind( 
      N = nrow(data),
      Female = get.count.per(!data$sex_0),
      Male = get.count.per(data$sex_0),
      age = get.count.per(is.na(data$age_0)),
      age.15.24 = get.count.per(data$age_0 < 25),
      age.25.plus = get.count.per(data$age_0>24),
      pdx = get.count.per(indicators$pDx),
      eART = get.count.per(indicators$eART),
      esupp = get.count.per(esupp)
    )
  }

  
table1
  
}

get.age.grp<- function(data, time=0){
  
  if(time==0){
    youth <- data$age_0 < 25
  } else if(time==1){
    youth <- data$age_0 < 24
  } else if(time==2){
    youth <- data$age_0 < 23
  } else{
    youth <- data$age_0 < 22
  }
  youth
}


# get baseline covariates
get.X <- function(data, time=3, adj.full=T){
			
	n <- nrow(data)
	
	# age # reference age group <20
	age.20.29 <- age.30.39 <- age.40.49 <- age.50.59 <- age.60.plus <- rep(0, n)
	age.20.29[ which(data$age_0>19 & data$age_0<30) ] <- 1
	age.30.39[ which(data$age_0>29 & data$age_0<40) ] <- 1
	age.40.49[ which(data$age_0>39 & data$age_0<50) ] <- 1
	age.50.59 [which(data$age_0>49 & data$age_0<60) ] <- 1
	age.60.plus[which(data$age_0>59)  ] <- 1
	age.matrix <- data.frame(cbind(age.20.29, age.30.39, age.40.49, age.50.59, age.60.plus))
	
	# reference is missing
	single <- married <- widowed <- 	divorced.separated <- rep(0, n)
	single[ which(data$marital_0==1)] <- 1
	married[ which(data$marital_0==2) ] <-1
	widowed[ which(data$marital_0 ==3)] <-1
	divorced.separated[ which(data$marital_0==4 | data$marital_0==5)] <-1
	marital <- data.frame(single, married, widowed, divorced.separated)
		
	# education: reference is less than primary or missing
	primary <- as.numeric(data$edu_primary_0)
	secondary.plus <- as.numeric(data$edu_secondary_plus_0)
	education <- data.frame(primary, secondary.plus)

	# occupation: reference NA
	formal.hi <- as.numeric(data$formal_hi_occup_0)
	informal.hi <- as.numeric(data$informal_hi_occup_0)
	informal.lo <- as.numeric(data$informal_low_occup_0)
	jobless <- as.numeric(data$jobless_0)
	student <- as.numeric(data$student_0)
	fisherman <- as.numeric(data$fisherman_0)
	occupation<- data.frame(formal.hi, informal.hi, informal.lo, jobless, student, fisherman)	
	
	# alcohol use: ref is NA
	alcohol.yes <- alcohol.no <- rep(0, n)
	alcohol.yes[which(data$alcohol_0) ] <- 1
	alcohol.no[which(!data$alcohol_0) ] <- 1

	# reference wealth is NA missing
	wealth0 <- wealth1<- wealth2 <- wealth3 <- wealth4 <-  rep(0, n) 
	wealth0[ which(data$wealth_0==0)] <- 1
	wealth1[ which(data$wealth_0==1)] <- 1
	wealth2[ which(data$wealth_0==2)] <- 1
	wealth3[ which(data$wealth_0==3)] <- 1
	wealth4[ which(data$wealth_0==4)] <- 1
	wealth <- data.frame(cbind(wealth0, wealth1, wealth2, wealth3, wealth4))

	#mobility indicators
	mobile <- as.numeric(data$mobile_0)
	
	# shifted main residence
	shift.no <- shift.yes <- rep(0,n)
	shift.no[which(!data$shifted_0)] <-1
	shift.yes[which(data$shifted_0)] <-1
	
	# nights home
	nights <- as.numeric(as.character(data$nightsHome_0))
	nights0 <- nights1.2 <- nights3.4 <- nights5 <- rep(0,n)
	nights0[which(nights==0)] <-1
	nights1.2[which(nights==1 | nights==2)] <-1
	nights3.4[which(nights==3 | nights==4)] <-1
	nights5[which(nights==5)] <- 1
	mobility <- data.frame(mobile, shift.no, shift.yes, nights0, nights1.2, nights3.4, nights5)
	
	# health-seeking
	chc.BL <- as.numeric(data$chc_0)
	self.hivtest.yes <- self.hivtest.no <- rep(0,n)
	self.hivtest.yes[which(data$self_hivtest_0)]<-1
	self.hivtest.no[which(!data$self_hivtest_0)] <-1
	health<- data.frame(chc.BL, self.hivtest.yes, self.hivtest.no)

	male <- rep(0,n)
	male[which(data$sex_0)] <- 1

	
	X <- cbind(
		age.matrix, marital,
		education, occupation,
		alcohol.yes, alcohol.no, wealth, mobility, male)
		
  if( !adj.full){
	    X <- cbind(age.matrix, male)
  }
	
	X
	
}





#*-------

compile.results <- function(N, nTstVL, Prob.HIVpos, Prob.supp){
  
  prev <- Prob.HIVpos$est
  colnames(prev) <- paste('prev', colnames(prev), sep='.')
  
  supp.pos <- data.frame(get.var.delta(mu1=Prob.supp$est$pt,
                                       IC1=Prob.supp$IC, 
                                       mu0=Prob.HIVpos$est$pt,
                                       IC0=Prob.HIVpos$IC)$est)
  colnames(supp.pos) <- paste('supp.pos', colnames(supp.pos), sep='.')
  
  est <- data.frame(cbind(N=N, 
                          nPos=round(N*Prob.HIVpos$est$pt),
                          nTstVL=nTstVL,
                          prev, supp.pos))
  
  est
}


# get.var.delta - function to get inference via the delta method
# 	assumes inputed estimators are asymptotically linear
#		i.e. written in first order as an empircal mean of an influence curve (IC)
#	input:  point estimates (mu1, mu0), corresponding influence curves (IC1, IC0)
#		significance level
#	output: point estimate, var, wald-type CI 

get.var.delta <- function(mu1, mu0=NULL, IC1, IC0=NULL, alpha=0.05){
  
  mu1<- unlist(mu1)
  
  if(is.null(mu0)){ 
    # if single TMLE 
    psi <- mu1
    IC <- IC1
    log <- F
    
  } else { 
    # if ratio of TMLEs (i.e. target = psi/psi0)
    mu0 <- unlist(mu0)
    # get inference via the delta method on log scale
    psi <- log(mu1/mu0)
    IC <- 1/mu1*(IC1) - 1/mu0*IC0
    log <- T
  }
  
  # variance of asy lin est is var(IC)/n
  var<- var(IC)/length(IC)
  # testing and CI	
  cutoff <- qnorm(alpha/2, lower.tail=F)
  se <- sqrt(var)
  CI.lo <- psi - cutoff*se
  CI.hi <- psi + cutoff*se
  
  if(log){
    est <- data.frame(pt=exp(psi), CI.lo=exp(CI.lo), CI.hi=exp(CI.hi) ) 
  }else{
    est <- data.frame(pt=psi, CI.lo=CI.lo, CI.hi=CI.hi) 
  }
  
  list(est=est, IC=IC)
}

# get.clust.CI.Xsect
# take obtain inference for cluster-specific estimates of prevalence & suppression

get.clust.CI.Xsect <- function(data){	
  
  A <- rep(1, nrow(data))
  # weight wrt number of indv in population
  alpha.prev <- get.weights(data, weighting='indv', this.col='N')$alpha
  
  prev <- call.ltmle(A=A, Y=data$prev.pt, observation.weights=alpha.prev)$est
  # weight following wrt number of estimated HIV+
  alpha <- get.weights(data, weighting='indv', this.col='nPos')$alpha
  supp <- call.ltmle(A=A, Y=data$supp.pos.pt, observation.weights=alpha)$est
  
  colnames(prev) <- paste('prev', colnames(prev), sep='.')
  colnames(supp) <- paste('supp', colnames(supp), sep='.')
  
  c(N=sum(data$N), nPos=sum(data$nPos), nTstVL=sum(data$nTstVL), prev, supp)
  
}



get.weights<- function(data, weighting, this.col){
  

  nTot <- sum(data[,this.col])
  J <- nrow(data) 

  # aggregated to cluster-level
  if(weighting=='clust'){
     data$alpha  <- 1
  } else{
    data$alpha  <- data[,this.col]*J/nTot
  }
  
  data
  
}



#*-------
#SCREENING ALGORITHMS FOR SUPERLEARNER
# See SuperLearner help file for more info: ?SuperLearner

screen.corRank10 <- function(Y, X, family, ...) screen.corRank(Y, X, family, rank = 10, ...)


# FUNCTIONS TO ENCODE OUR DETERMINISTIC KNOWLEDGE
# See ltmle help file for more info: ?ltmle

# deterministicQ_YES
# if detQ.variable==1, then outcome==1 with probability 1
deterministicQ_YES<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "detQ.variable")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]==1
  return(list(is.deterministic=is.deterministic, Q.value=1))
}

# deterministicQ_NO
# if detQ.variable==0,  then outcome==0 with probability 1
deterministicQ_NO<- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "detQ.variable")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data[,L2.index]== 0
  return(list(is.deterministic=is.deterministic, Q.value=0))
}


#*-------

get.file.name <- function(timept, date=NULL){
  
  if(is.null(date)){
    date <- 	format(Sys.timept(), "%d%b%Y")
  }
  file.name <- paste('Supp',
                     paste('Year', timept, sep=''),
                     paste('v', date,  sep=''), sep="_")
  file.name
}

