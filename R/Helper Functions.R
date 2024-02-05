indicator<-function(condition) ifelse(condition,1,0)

gam.variables <- function(data){
  sapply(1:ncol(data), function(i){
    if (is.numeric(data[,i]) | is.integer(data[, i])){ paste("s", "(",names(data[i]), ")", sep="")}
    else{names(data[i]) }
  })
}

cox.variables <- function(data){
  sapply(1:ncol(data), function(i){
    if (is.numeric(data[,i]) | is.integer(data[, i])){ paste("ns", "(",names(data[i]), ")", sep="")}
    else{names(data[i]) }
  })
}
