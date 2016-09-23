###data generation for simulation

library(pcIRT)

repnumb <- 1000
load("trueitpar3-15.RData")

for (i in 1:repnumb){
  res <- simMPRM(truep, 1000)
  nam <- paste0("datTAMMP", i)
  save(res, file=nam)
}




