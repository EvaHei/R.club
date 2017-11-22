# works also with categorical data

#kNN erklärt: https://www.youtube.com/watch?v=jU5j7Jd6hy0

##example

library(agridat)
summary(ilri.sheep)
str(ilri.sheep)
sheep2 <-ilri.sheep[,c(1,5,6,7,8,10)]
summary(sheep2) ##NAs due to death

##impute

library(VIM)

?kNN

# kNN(data, variable = colnames(data)# default all, k=5)
#k should be between 5 and 10

sheepimp <- kNN(ilri.sheep, variable = c("weanwt","weanage"), k=6)

