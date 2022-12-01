library(mgcv)
library(mgcViz)
library(itsadug)
library(plyr)
library(tools)

source('tw_getpred.R')

#paths need to be set:
if(.Platform$OS.type=="unix"){
	gamdata_dir <- '/home/tilsen/Data/timewarping_opendata/GAM/'
} else {
	gamdata_dir <- 'M:/Data/timewarping_opendata/GAM/'
}

if (! (file.exists(gamdata_dir))) {
	print('data directory not found')
	stop()
}

#defaults:
fitmethod <- "fREML"
discreteopt <- TRUE
nthreads <- 16
numknots <- 50
numknots_re <- 25
basis <- "cr"
mval <- 2
mvalre <- 1
gamma <- 1

#datatransform <- "orig"
datatransform <- "log"

#fitmethod <- "REML"
#extraspec <- "REML"

extraspec <- ""

modelspec <- paste0(datatransform,extraspec)

files <- c('gamdata_CVCarl_deltas_100_P05.csv',
           'gamdata_CVCacl_deltas_100_P05.csv')

for (i in (1:length(files))){
  
  fname <- file_path_sans_ext(files[i])
  datacsv <- paste0(gamdata_dir,files[i])
  outfname <- paste0(sub(pattern="gamdata_",replacement="",fname),"_",modelspec)
  kcheckcsv <- paste0(gamdata_dir,"gamkcheck_",outfname,".csv")
  predcsv <- paste0(gamdata_dir,"gampred_",outfname,".csv")
  modelfile <- paste0(gamdata_dir,"gammodel_",outfname,".rds")
  
  print(paste("estimating gam for:",outfname))
  
  D <- read.csv(datacsv)
  
  if (datatransform=="log") {
	D$lrrgsd <- log(D$lrrgsd)
  }
  
  D$subj <- as.factor(D$subj)
  D$cond <- as.factor(D$cond)
  
  #sort data for autocorrelation correction:
  D <- D[order(D$subj, D$cond, D$trial, D$time), ]
  D$start.event <- D$time==0
    
  #first fit model without AR coefficient
  m0 <- bam(lrrgsd ~ cond + 
              s(time,by=cond,bs=basis,k=numknots,m=mval) + 
              s(time,subj,by=cond,bs="fs",m=mvalre,k=numknots_re),
            data=D,method=fitmethod,discrete=discreteopt,nthreads=nthreads)
    
  #now refit with autocorrelation
  m0acf <- acf_resid(m0,plot=FALSE)
  rhoval <- m0acf[2]
  m1 <- bam(lrrgsd ~ cond +
              s(time,by=cond,bs=basis,k=numknots,m=mval) + 
              s(time,subj,by=cond,bs="fs",m=mvalre,k=numknots_re),
            rho = rhoval, AR.start = D$start.event,
            data=D,method=fitmethod,discrete=discreteopt,nthreads=nthreads)
  
  #gam.check(m1)
  kcheck <- k.check(m1,subsample=5000,n.rep=400)
  write.csv(kcheck,kcheckcsv)
  
  times <- unique(D$time)
  conds <- seq(1,length(unique(D$cond)))
  subjs <- unique(D$subj)
  
  pred <- get_pred(m1,times,conds,subjs)
  write.csv(pred,predcsv)

  saveRDS(m1,modelfile)
  
}









