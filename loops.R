library(rbenchmark)
library(dplyr)
library(reshape2)

###############################################################################################################
## Daten einlesen
TDL.pop <- read.table("TDL-jan.csv", sep=",", header=TRUE)

##BENCHMARK####################################################################################################
benchmark(
###############################################################################################################

###############################################################################################################
## No-Loop v1.0
##BENCHMARK####################################################################################################
"noloop" = {
###############################################################################################################  
early.TDL2 <- TDL.pop %>% arrange(DTF30) %>% group_by(E, DL, T) %>% slice(1:20) %>% mutate(no = row_number())
early.TDL2 <- subset(early.TDL2, select=c(9,3:5,2))
early.TDL2 <- dcast(early.TDL2, no ~ E + T + DL, value.var = "Genotyp")[,2:17]

late.TDL2 <- TDL.pop %>% arrange(desc(DTF30)) %>% group_by(E, DL, T) %>% slice(1:20) %>% mutate(no = row_number())
late.TDL2 <- subset(late.TDL2, select=c(9,3:5,2))
late.TDL2 <- dcast(late.TDL2, no ~ E + T + DL, value.var = "Genotyp")[,2:17]
###BENCHMARK###################################################################################################
},
###############################################################################################################
## Dreifach-Loop v2.0
##BENCHMARK####################################################################################################
"janloop" = {
###############################################################################################################  
nu <- 20
early.TDL1 <- data.frame(row.names=1:nu)
late.TDL1 <- data.frame(row.names=1:nu)

for (t in 1:2) {
  for (dl in 1:2) {
    for (e in 1:4) {
      
      TDL.subset <- TDL.pop[TDL.pop$E == e & TDL.pop$DL == dl & TDL.pop$T == t,]
      subset.order <- TDL.subset[with(TDL.subset, order(DTF30)),]$Genotyp
      subset.length <- length(subset.order)
      
      early.TDL1[,paste("DL",dl,"E",e,"T",t)] <- subset.order[1:nu]
      late.TDL1[,paste("DL",dl,"E",e,"T",t)] <- subset.order[(subset.length-(nu-1)):subset.length]
   
      }}}
###BENCHMARK###################################################################################################
  },
###############################################################################################################
## Das Original
##BENCHMARK####################################################################################################
"evaloop" = {
###############################################################################################################
nu <- 20 #number of GT in early or late
early.TDL <- data.frame(nrow=nu)
late.TDL <- data.frame(nrow=nu)
early.pos <- c()
late.pos <- c()

for (t in 1:2) {
  for (dl in 1:2) {
    for (e in 1:4) {
      coln <- paste0(e,"DL",dl,"T",t)
      sub.TDL <- subset(TDL.pop, TDL.pop$E == e & TDL.pop$DL == dl & TDL.pop$T==t)
      GT.ord <- order(sub.TDL$DTF30
                      ,decreasing = F, na.last = NA)
      TDL.early <- (head(sub.TDL[GT.ord,]$Genotyp, n=nu))
      TDL.late <- (tail(sub.TDL[GT.ord,]$Genotyp, n=nu))
      early.TDL <- cbind(early.TDL,TDL.early)
      late.TDL <- cbind(late.TDL,TDL.late)
      colnames(early.TDL)[length(early.TDL)] <- coln
      colnames(late.TDL)[length(late.TDL)] <- coln
      
      early.pos <- c(early.pos,TDL.early)
      late.pos <- c(late.pos,TDL.late)
    }
  }
}
###BENCHMARK###################################################################################################
}, replications = 100, columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self"))
###############################################################################################################