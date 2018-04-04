# Specify the Working Directory

library(affy)

#
# This code retrieves all the .CEL files in the directory (assuming each directory contains an experiment
# with appropriate .CEL files).
#
# They are then standardized by normalizing each individual .CEL file. The results are then
# saved into directories containing the original intensities and the log2 intensities.
#
# This is done to simplify the simulations so it doesn't have to load the .CEL files with all the details
# then standardize afterwards. 

######

# List the directorys containing the experiments
#
# Improved Code Would Allow for Specification of the File Externally
#
standardize_cel_files = function(raw_dirs, ori_loc, std_loc, details_file, metrics_file, bootstrap_design) {
  raw_files_dir = "./Exp_Raw/"
  ori_files_dir = "./Exp_Ori/"
  std_files_dir = "./Exp_Std/"
  #setwd('c:/git/capstone/drews')
  
  #raw_files_dir = raw_dirs
  #ori_files_dir = ori_loc
  #  std_files_dir = std_loc
  
    
  
  dir.names <- list.files(raw_files_dir)
  
  dir.names
  
  # Return list for each directory of the .CEL files
  
  dir.files <- lapply(dir.names, function(i) {
    return(list.files(paste(raw_files_dir,i,sep=""))) })
  
  # Create a Dataframe for each experiment and directory
  
  exp.num <- sapply(dir.files,length)
  exp.list <- data.frame(rep(1:length(dir.files),exp.num),unlist(dir.files),rep(dir.names,exp.num))
  exp.list <- exp.list[order(exp.list[,2]),]
  # Given the Filelist Retrieve and Save the Files
  
  apply(exp.list,1,retrieve.save)
  
  exp.list <- cbind(exp.list,apply(exp.list,1,retrieve.ori.name))
  exp.list <- cbind(exp.list,apply(exp.list,1,retrieve.std.name))
  
  colnames(exp.list) <- c("Exp","Filename","Dir","Orig","Std")
  write.csv(exp.list,details_file)
  
  bootstrap.list <- exp.list[,c("Dir","Filename","Exp")]
  bootstrap.list <- data.frame(lapply(bootstrap.list, function(x) {
                      gsub(".CEL", "", x)
                  }))
  bootstrap.list$Group=1
  colnames(bootstrap.list) <- c("Directory","Name","Experiment","Group")
  
  write.csv(bootstrap.list,bootstrap_design, row.names=FALSE)
  
  # Calculate Means and SD's of Each Array
  
  metrics <- apply(exp.list,1,retrieve.metrics) 
  exp.metrics <- data.frame(cbind(exp.list[,1],t(metrics)))
  
  # Calculate the Experiment Specific Traits
  
  mean.mean <- rep(tapply(exp.metrics[,2],exp.metrics[,1],mean),table(exp.metrics[,1]))
  mean.sd <- rep(tapply(exp.metrics[,2],exp.metrics[,1],sd),table(exp.metrics[,1]))
  
  sd.mean <- rep(tapply(exp.metrics[,3],exp.metrics[,1],mean),table(exp.metrics[,1]))
  sd.sd <-   rep(tapply(exp.metrics[,3],exp.metrics[,1],sd),table(exp.metrics[,1]))
  
  std.mean <- (exp.metrics[,2] - mean.mean)/mean.sd
  std.sd   <- (exp.metrics[,3] - sd.mean)/sd.sd
  
  # Create The Updates Experiment Metrics
  
  exp.metrics <- data.frame(cbind(exp.metrics,std.mean,std.sd))
  
  names(exp.metrics) <- c("Exp","Mean","SD","S.Mean","S.SD")
  
  write.csv(exp.metrics,metrics_file)
}

# Function to save the log2 and standardized log2 intensities

retrieve.save <- function(i) {
  
  # Open the Files
  cel <- ReadAffy(filenames = i[2],celfile.path = paste(raw_files_dir,i[3],sep=""))
  int <- intensity(cel)[,1]
  
  int[int < 1] <- 1
  log2.int <- log2(int)
  split.list <- unlist(strsplit(i[2],"\\."))
  std.log2.int <- (log2.int - mean(log2.int))/sd(log2.int)
  writeBin(log2.int,paste(ori_files_dir,split.list[1],".ori",sep=""))
  writeBin(std.log2.int,paste(std_files_dir,split.list[1],".std",sep=""))
}

# Function to Retrieve the Log2 Intensity Files

retrieve.ori.name <- function(i) {
  split.list <- unlist(strsplit(i[2],"\\."))
  return(paste(split.list[1],".ori",sep=""))
}

# Function to Retrieve the Log2 Standardized Files

retrieve.std.name <- function(i) {
  split.list <- unlist(strsplit(i[2],"\\."))
  return(paste(split.list[1],".std",sep=""))
}

# Create the Metrics For The Experiments

retrieve.metrics <- function(i) {
  data <- readBin(paste(ori_files_dir,i[4],sep = ""),what = "numeric", 2000000)
  return(c(mean(data),sd(data)))
}


#standardize_cel_files('./Exp_Raw/','./Exp_Ori/','./Exp_Std/','ExperimentDetails.csv','ExperimentMetrics.csv', 'bootstrap_design.csv')
