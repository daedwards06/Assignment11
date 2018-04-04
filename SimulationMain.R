######################################################
#
# AffyMetrix Gene Expression Microarray Simulator
# Version: 0.9.4
# Author: Andrew Hardin
#
######################################################

source("C:/git/Capstone/drews/SimConfig.R")
#General Steps

# Step #1: Read and Determine Simulation Parameters
# Step #2: Construct Experiment
# Step #3: Construct Arrays
# Step #4: Construct Probes
# Step #5: Apply Fold Change Effects

simulate.array <- function(exp.file,config.file) {

  # Check for experiment file
  # 
  # The experiment file has the details specific to the simulated experiment
  #

  if(missing(exp.file)) {
    stop("Experiment File Not Specified")
  }
  if(!file.exists(exp.file)) {
    stop(paste("Experiment File",exp.file,"Not Found"))
  }

  # Check for configuration file
  #
  # The configuration file has the details specific to the arrays used to
  # build the simulation. In particular they specify the structure of the
  # arrays used to generate the bootstrap samples.
  #

  if(missing(config.file)) {
    stop("Configuration File Not Specified")
  }
  if(!file.exists(config.file)) {
    stop(paste("Configuration File",config.file,"Not Found"))
  }

  # Retrieve Configuration

  exp.list <- read.exp(exp.file)
  param.list <- read.config(config.file)

  verify.config(param.list)
  exp.list <- verify.exp(exp.list)

  #Retrieve High Level Experiment Parameters

  exp <- construct.experiment(param.list,exp.list)
  
  # Simulate the low level probe details
  
  probes <- construct.probes(param.list,exp.list,exp)
  return(probes)
  
  #probes <- apply.factor.changes(probes,param.list,exp.list)

  #
  # This is now deprecated.
  # The probes contain the complete information that is needed.
  
  #create.cel.file(paste(exp.list$General$OutputDir,"/",exp.list$General$ExperimentSuffix,sep = ""),
  #  round(2^probes,0),param.list$General$Header,exp.list$GroupCoding)

}

construct.experiment <- function(param.list,exp.list) {

  # Determine General Experimental Conditions

  # Case 1: 
  
  # Experiment Mean is drawn from a normal distribution with specified mean and variance.
  # Experiment Variance is drawn from a linear relationship with the experiment mean through a linear relation


  # Experiment Mean is drawn from a normal distribution with specified
  # mean and variance

  if(!file.exists(param.list$Bootstrap$BootMetrics)) {
    stop(paste("Invalid File",param.list$Bootstrap$BootMetrics))
  }

  exp.vals      <- read.csv(param.list$Bootstrap$BootMetrics)[,-1]
  exp.mean.mean <- tapply(exp.vals[,2],exp.vals[,1],mean)
  exp.mean.sd   <- tapply(exp.vals[,2],exp.vals[,1],sd)
  exp.sd.mean   <- tapply(exp.vals[,3],exp.vals[,1],mean)
  exp.sd.sd     <- tapply(exp.vals[,3],exp.vals[,1],sd)
   
  exp.mean.res  <- unlist(tapply(exp.vals[,2],exp.vals[,1],function(i) {
     return((i - mean(i))/sd(i)) }))
  exp.sd.res <- unlist(tapply(exp.vals[,3],exp.vals[,1],function(i) {
     return((i - mean(i))/sd(i)) }))

  sample.index <- sample(1:length(exp.mean.mean),1)
  exp.loc.mean     <- exp.mean.mean[sample.index]
  exp.loc.sd       <- exp.mean.sd[sample.index]
  exp.scale.mean   <- exp.sd.mean[sample.index]
  exp.scale.sd     <- exp.sd.sd[sample.index]

  sample.index.mean <- sample(1:length(exp.mean.res),exp.list$NumOfArrays,replace = TRUE)
  sample.index.sd <- sample(1:length(exp.sd.res),exp.list$NumOfArrays,replace = TRUE)
  return(list(loc = (exp.mean.res[sample.index.mean]*exp.loc.sd) + exp.loc.mean,scale = (exp.sd.res[sample.index.sd]*exp.scale.sd) + exp.scale.mean))

}


#
# Step #4: Construct Probes
#

construct.probes <- function(param.list,exp.list,exp.param) {

  # BootConfig contains the location of a file describing
  # the structure of the experiment. 
  #
  # 

  if(!file.exists(param.list$General$PMIndex)) {
    stop(paste("Invalid File",param.list$General$MMIndex))
  }
  
  if(!file.exists(param.list$General$MMIndex)) {
    stop(paste("Invalid File",param.list$General$MMIndex))
  }
  
  nrow = as.integer(param.list$General$ArraySize)
  pm.index = readBin(param.list$General$PMIndex,what="integer",n=10000000)
  mm.index = readBin(param.list$General$MMIndex,what="integer",n=10000000)

  exp.listing = read.csv(param.list$Bootstrap$BootDesign)
  # Select an experimental sample from the samples database and use the data to bootstrap the probe
  # level variability.

  group <- unique(exp.listing[,3])
  exp <- sample(group,1)
  exp.mask <- exp.listing[,3] == exp
  exp.spec <- exp.listing[exp.mask,]
  group.det <- unique(exp.spec[,4])
  exp.det <- sample(group.det,1)
  exp.mask.det <- exp.listing[,4] == exp.det
  exp.mask.full <- exp.mask & exp.mask.det

  #Determine which rows are pm/mm which must be done in pairs, and which are single

  non.index <- 1:nrow
  non.index <- non.index[-c(pm.index,mm.index)]

  #Retrieve the cel intensity and adjust the bootstrap samples to be on the same scale
 
  boot.sam <- sapply(exp.listing[exp.mask.full,2],function(i) {
      return(readBin(paste(param.list$Bootstrap$BootDir,"/",i,".std",sep = ""),what = "numeric", 2000000))
  })
    
  new.boot.sam <- matrix(rep(0,nrow * exp.list$NumOfArrays), ncol = exp.list$NumOfArrays)

  #Build bootstrap by looping through number of arrays

  boot.seq <- 1:ncol(boot.sam)
  l.non <- length(non.index)
  l.oth <- length(pm.index)
  l.non.seq <- 1:l.non
  l.oth.seq <- 1:l.oth

  for(i in 1:exp.list$NumOfArrays) {
    boot.sel <- sample(boot.seq,l.non,replace = T)
    for(j in l.non.seq) {
      ind <- non.index[j]
      new.boot.sam[ind,i] <- (boot.sam[ind,boot.sel[j]] * exp.param$scale[i]) + exp.param$loc[i]
    }
    boot.sel <- sample(boot.seq,l.oth,replace = T)
    for(j in l.oth.seq) {
      ind.pm <- pm.index[j]
      ind.mm <- mm.index[j]
      new.boot.sam[ind.pm,i] <- (boot.sam[ind.pm,boot.sel[j]] * exp.param$scale[i]) + exp.param$loc[i]
      new.boot.sam[ind.mm,i] <- (boot.sam[ind.mm,boot.sel[j]] * exp.param$scale[i]) + exp.param$loc[i]
    }
  }
  return(new.boot.sam)
}

#
# Step #5: Apply Fold Change Effects
#

apply.factor.changes <- function(probes,param.list,exp.list) {

  if(!file.exists(exp.list$General$ExperimentFactor)) {
    stop(paste("Invalid File",exp.list$General$ExperimentFactor))
  }

  factor.list <- read.csv(exp.list$General$ExperimentFactor)
  gene.list <- as.character(factor.list[,1])

  probe.group <- rep(1:length(exp.list$GroupCoding),each=exp.list$GroupCoding)

  # Predicted intensity
  # Based upon spherical curve

  logistic <- function(i,min.x,min.y,max.x,max.y) { 
     s.y <- (max.y - min.y)
     s.x <- (max.x - min.x)/6
     mid.x <- (min.x + max.x)/2
     return( (s.y/(1 + exp(-((i - mid.x)/s.x)))) + min.y)
  }

  logistic.inv <- function(i,min.x,min.y,max.x,max.y) {
     i[i <= min.y] <- min.y+0.001
     i[i >= max.y] <- max.y-0.001
     s.y <- (max.y - min.y)
     s.x <- (max.x - min.x)/6
     mid.x <- (min.x + max.x)/2  
     retval <- mid.x - s.x*log((s.y/(i - min.y)) - 1)
     retval[retval < min.x] <- min.x
     retval[retval > max.x] <- max.x
     return(retval)
  }

  data <- apply(probes,2,function(i) {
    quant <- quantile(i,c(0.413,0.9985))
    #The high range of 512 was chosen arbitrarily from the data, seems to represent 
    #the highest concentration of meaning on the HGU95 and HGU133 platforms, converted to log2 that is a value
    #of 9
    return(c(-4,quant[1],9,quant[2]))
  })

  source(param.list$General$IndexProbes)
  fold.pattern <-  t(apply(factor.list[,2:ncol(factor.list)],1, 
    function(i) { rep(i,each = exp.list$GroupCoding) }))

  quant.max <- runif(1,14,15)

  for(j in 1:ncol(probes)) {

    quant.min <- quantile(probes[,j],c(0.0001))

    vals <- sapply(pm.index,function(i) {
      return(median(probes[i,1]))
    })
    
    quant.med <- quantile(vals,c(0.43,0.99))
    quant.low <- quant.med[1]
    quant.high <- quant.med[2]

    for(i in 1:nrow(factor.list)) {
        if(fold.pattern[i,j] == 0) { next }

        probe.index <- cel.probes[[gene.list[i]]]

        #Data too low or high

        orig.median.int <- median(probes[probe.index,j])
        pred.conc <- logistic.inv(orig.median.int,data[1,j],quant.low,data[3,j],quant.high) + fold.pattern[i,j]
        pred.conc[pred.conc < data[1,j]] <- data[1,j]
        pred.conc[pred.conc > data[3,j]] <- data[3,j]
        pred.median.int <- logistic(pred.conc,data[1,j],quant.low,data[3,j],quant.high)
        mult.ratio <- pred.median.int - orig.median.int 
        probes[probe.index,j] <- probes[probe.index,j] + mult.ratio
    }

    probes[,j]  <- probes[,j] + rnorm(length(probes[,j]),0,0.1)
    probes[probes[,j] > quant.max,j] <- quant.max + (16-quant.max)*(pgamma(probes[probes[,j] > quant.max,j]-quant.max, shape = 3,scale = 1))
    probes[probes[,j] < quant.min,j] <- quant.min - (quant.min-4) *(pgamma(quant.min - probes[probes[,j] < quant.min,j], shape = 3,scale = 1))
  }
  return(probes)
}


