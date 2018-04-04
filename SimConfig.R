read.config <- function(config.file) {
  config.lines <- readLines(config.file)
  param.list <- list(General = list(),Bootstrap = list())
  cursec <- NULL

  for(i in 1:length(config.lines)) {
    if(config.lines[i] == "") { next }
    else if(config.lines[i] == "[General]") {
       cursec <- 1
    }
    else if(config.lines[i] == "[Bootstrap]") {
       cursec <- 2
    }
    else {
      if(is.null(cursec)) {
        stop("Invalid Configuration File")
      }
      split.vals <- strsplit(config.lines[i],"=")
      if(length(split.vals[[1]]) != 2) {
        stop(paste("Invalid Value At Line",i,":",config.lines[i]))
      }
      config.name <- split.vals[[1]][1]
      config.setting <- split.vals[[1]][2]
      param.list[[cursec]] <- append(param.list[[cursec]],config.setting)
      names(param.list[[cursec]])[length(param.list[[cursec]])] <- config.name
    }
  }
  return(param.list)
}

read.exp <- function(exp.file) {
  config.lines <- readLines(exp.file)
  param.list <- list(General = list(),ExpGroups = list())
  cursec <- NULL

  for(i in 1:length(config.lines)) {
    if(config.lines[i] == "") { next }
    else if(config.lines[i] == "[General]") {
       cursec <- 1
    }
    else if(config.lines[i] == "[ExpGroups]") {
       cursec <- 2
    }
    else {
      if(is.null(cursec)) {
        stop("Invalid Configuration File")
      }
      split.vals <- strsplit(config.lines[i],"=")
      if(length(split.vals[[1]]) != 2) {
        stop(paste("Invalid Value At Line",i,":",config.lines[i]))
      }
      config.name <- split.vals[[1]][1]
      config.setting <- split.vals[[1]][2]
      param.list[[cursec]] <- append(param.list[[cursec]],config.setting)
      names(param.list[[cursec]])[length(param.list[[cursec]])] <- config.name
    }
  }
  return(param.list)
}

verify.config <- function(param.list) {
  if(is.null(param.list$Bootstrap)) {
    stop("Invalid Configuration File, Should Contain [Bootstrap]")
  }
  if(is.null(param.list$Bootstrap$BootDir)) {
    stop("Invalid Configuration File, Should Contain BootDir setting under [Bootstrap]")
  }
#  if(is.null())
#  if(is.null(param.list$Bootstrap$BootMetrics)) {
#    stop("Invalid Configuration File, Should Contain BootMetrics setting under [Bootstrap]")
#  }
#  if(is.null(param.list$Bootstrap$BootConfig)) {
#    stop("Invalid Configuration File, Should Contain BootConfig setting under [Bootstrap]")
#  }
#  if(is.null(param.list$General$Header)) {
#     stop("Invalid Configuration File, Should Contain Header setting under [General]")
#  }
#  if(is.null(param.list$General$PMMM)) {
#     stop("Invalid Configuration File, Should Contain PMMM setting under [General]")
#  }
#  if(is.null(param.list$General$IndexProbes)) {
#     stop("Invalid Configuration File, Should Contain IndexProbes setting under [General]")
#  }
}

verify.exp <- function(exp.list) {
  if(is.null(exp.list$General)) {
    stop("Invalid Configuration File, Should Contain [General]")
  }
  if(is.null(exp.list$General$ExperimentSuffix)) {
    stop("Invalid Configuration File, Should Contain ExperimentSuffix setting under [General]")
  }
  if(is.null(exp.list$General$OutputDir)) {
    stop("Invalid Configuration File, Should Contain OutputDir setting under [General]")
  }
  if(is.null(exp.list$General$ExperimentFactor)) {
    stop("Invalid Configuration File, Should Contain ExperimentFactor setting under [General]")
  }

  if(length(exp.list$ExpGroups) == 0) {
    stop("Invalid Configuration File, at least one entry must be specified under ExpGroups")
  }
  numArrays <- 0
  grpCoding <- vector()
  for(i in 1:length(exp.list$ExpGroups)) {
    exp.list$ExpGroups[[i]] <- as.integer(exp.list$ExpGroups[[i]])
    numArrays <- numArrays + exp.list$ExpGroups[[i]]
    if(is.na(exp.list$ExpGroups[[i]])) {
      stop(paste("Entry",i,"in ExpGroups must be an integer"))
    }
    grpCoding <- c(grpCoding,exp.list$ExpGroups[[i]])
  }
  exp.list <- append(exp.list,list(NumOfArrays = numArrays,GroupCoding = grpCoding))
  return(exp.list)
}