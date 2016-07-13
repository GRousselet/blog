# make data frame for single vector with default names "data"
# and "gr"
# GAR, University of Glasgow, 2016-07-07
mkdf1 <- function(x,name=c("data","gr")){
gr <- c(rep('Group 1',length(x)))
data <- data.frame(x,gr)
names(data) <- name
data
}

# make data frame for two vectors with default names "data"
# and "gr"
# GAR, University of Glasgow, 2016-07-09
mkdf2 <- function(x,y,name=c("data","gr")){
data <- c(x,y)
gr <- c(rep('Group 1',length(x)),rep('Group 2',length(y)))
df <- data.frame(data,gr)
names(df) <- name
df
}
