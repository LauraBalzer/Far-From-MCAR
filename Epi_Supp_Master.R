##############
# R code for "Far from MCAR: 
#   obtaining population-level estimates of HIV viral suppression"
# 
# Additional sensitivity analyses available by request
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Primary Statistician for SEARCH
#
# This code is largely adapted from the primary & secondary analyses
# for Havlir et al. NEJM 2019: github.com/LauraBalzer/SEARCH_Analysis_Adults
###############

rm(list=ls())


set.seed(1)
library('SuperLearner')
library('ltmle')

source('Epi_Supp_Functions.R') 
source('Epi_Supp_Helper.R') 


DATE <- 'final'

# specify the follow-up year
# options: {0,1,2,3}
timept <- 3


# SL.Library 
SL.library <- list( c('SL.glm', 'screen.corRank10'), 
                    c('SL.glm', 'screen.corP'),
                    c('SL.gam', 'screen.corRank10'), 
                    c('SL.gam', 'screen.corP'),
                    c('SL.mean', 'All'))

#########################################################


file.name <-  get.file.name(timept=timept, date=DATE)

out<- Run.Supp.Xsect(timept=timept, SL.library=SL.library)

save( out, timept, SL.library, file=paste(file.name, 'Rdata', sep='.'))
round(out, 2)

# TABLES 
# table1
table1 <- get.consort()

# table2 
data.input <- preprocess.cascade()
table2 <- NULL; demographics=F

for(timept in 0:3){
  # define target population for the X-sectional analysis
  data <- subset(data.input, get.pop.Xsect(data=data.input, timept=timept))
  table2 <- cbind(table2, 
                  get.table1(data, timept=timept, subset.indicator='VLmeas', demographics=demographics),
                  get.table1(data, timept=timept, subset.indicator='VLno', demographics=demographics)
                  )
}
colnames(table2) <- paste(c('VL', 'noVL'), c(0,0,1,1,2,2,3,3), sep='')
write.csv(t(table2), file=paste('T.Table2', demographics, DATE, 'csv', sep='.'))


# etable1 - baseline
etable1 <- NULL; demographics=T 
data <- subset(data.input, get.pop.Xsect(data=data.input, timept=0))
  
etable1 <- cbind(
  get.table1(data, timept=0, subset.indicator='All', demographics=demographics),
  get.table1(data,timept=0, subset.indicator='HIVmeasure', demographics=demographics),
  get.table1(data,timept=0, subset.indicator='HIVpos', demographics=demographics),
  get.table1(data, timept=0, subset.indicator='VLmeas', demographics=demographics)
)
colnames(etable1) <- c('enumerated', 'HIVmeasre', 'HIVpos', 'VLmeas')
  
library(xtable)
xtable(etable1 ) 
 


