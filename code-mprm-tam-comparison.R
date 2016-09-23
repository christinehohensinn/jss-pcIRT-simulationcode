
#######main function comparison MPRM TAM

require(TAM)
require(pcIRT)

load("trueitpar3-15.RData")
load("desA.RData")
load("desB.RData")

repnumb <- 1000
ncat <- nrow(truep)
nitems <- ncol(truep)

#create objects for saving results
parTAM <- matrix(NA, length(as.vector(truep[-ncat,])), repnumb)
parPC <- matrix(NA, length(as.vector(truep[-ncat,])), repnumb)

seTAM <- matrix(NA, length(as.vector(truep[-ncat,])), repnumb)
sePC <- matrix(NA, length(as.vector(truep[-ncat,])), repnumb)

iterTAM <- vector(mode="numeric", length=repnumb)
convPC <- vector(mode="numeric", length=repnumb)

mseTAM <- vector(mode="numeric", length=repnumb)
msePC <- vector(mode="numeric", length=repnumb)

timeTAM <- vector(mode="list", length=repnumb)
timePC <- vector(mode="list", length=repnumb)


for(i in 1:repnumb){

naml <- paste0("datTAMMP", i)
load(naml)

#colnames for  TAM
daten <- res$datmat
colnames(daten) <- paste0("I", 1:ncol(daten))

tp <- system.time(mod.pc <- MPRM(daten))
            
tt <- system.time(mod.tam <- tam.mml(resp=daten , A=A2 , B=B2 , control= list( maxiter=1000, progress=FALSE) ))

#normalization to 0 for comparison to pcIRT
itpTAM <- rbind(t(mod.tam$AXsi)[1,]-mean(t(mod.tam$AXsi)[1,]),t(mod.tam$AXsi)[2,]-mean(t(mod.tam$AXsi)[2,]))

mse.tam <- sum((as.vector(itpTAM)-(as.vector(truep[-ncat,]))*(-1))^2)
mse.pc <- sum((as.vector(mod.pc$itempar[-ncat,])-(as.vector(truep[-ncat,])))^2)

##saving results

#item parameters
parTAM[,i] <- as.vector(itpTAM)
parPC[,i] <- as.vector(mod.pc$itempar[-ncat,])

seTAM[,i] <- as.vector(t(mod.tam$se.AXsi[,-ncat]))
sePC[,i] <- as.vector(mod.pc$itempar_se[-ncat,])

iterTAM[i] <- mod.tam$iter
convPC[i] <- mod.pc$convergence

#MSE
mseTAM[i] <- mse.tam
msePC[i] <- mse.pc

tr1 <- paste0("paraTAM", i)
tr2 <- paste0("modPC", i)
save(itpTAM, file=tr1)
save(mod.pc, file=tr2)

#times
timeTAM[[i]] <- tt
timePC[[i]] <- tp
 
outnam <- paste0("res", i)
save(itpTAM, mod.pc, tt, tp, file=outnam)
rm(mod.tam, mod.pc, itpTAM, tt, tp, daten, naml)
}
save.image("results.RData")


