plot.space <- function(name){
  tt <- read.table(name)
  x <- 1:(length(as.list(tt$V1)))
  y <- 1:(length(names(tt))-1)
  z <- 1:(length(x)*length(y))
  dim(z) <- c(length(x), length(y))
  for (i in 1:(length(x))){ for (j in 1:(length(y))) { z[i,j] <- tt[i,j]}}
  image(x,y,z) 
}
