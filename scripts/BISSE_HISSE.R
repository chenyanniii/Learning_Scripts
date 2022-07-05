## BISSE & HISSE
## The following script and data is borrowed from ilhabela workshop in 2015.
# https://lukejharmon.github.io/ilhabela/2015/07/05/BiSSE-and-HiSSE/

# other learning material:   
# https://github.com/lukejharmon/ilhabela/blob/master/workingFiles/diversification/diversification.md

# https://lukejharmon.github.io/ilhabela/lectures/alfarinho/state_dependent_diversification_ilhabela_2015.pdf


###############################################
#################### BiSSE ####################
###############################################

library(ape)
library(TreeSim)

library(diversitree)

simTree1 <- sim.bd.age(age = 10, numbsim = 1, lambda = 0.4, mu = 0)[[1]]

# first fit a Yule model
pbModel <- make.yule(simTree1)
pbMLFit <- find.mle(pbModel, 0.1)

# next fit a Birth-death model
bdModel <- make.bd(simTree1)
bdMLFit <- find.mle(bdModel, c(0.1, 0.05), method = "optim", lower = 0)

# compare models
anova(bdMLFit, pure.birth = pbMLFit)

bdSamples <- mcmc(bdModel, bdMLFit$par, nsteps = 1e+05, lower = c(0, 0), upper = c(Inf,  Inf), w = c(0.1, 0.1), fail.value = -Inf, print.every = 10000)

postSamples <- bdSamples[c("lambda", "mu")]
profiles.plot(postSamples, col.line = c("red", "blue"), las = 1, legend = "topright")

postSamples$r <- with(bdSamples, lambda - mu)
postSamples$eps <- with(bdSamples, mu/lambda)

profiles.plot(postSamples[, c("r", "eps")], col.line = c("red", "blue"), las = 1,  legend = "topright")




simPars <- c(0.4, 0.2, 0.05, 0.05, 0.05, 0.05)
set.seed(3)
simBisseData <- tree.bisse(simPars, max.t = 14, x0 = 0)

hist <- history.from.sim.discrete(simBisseData, 0:1)
plot(hist, simBisseData)

nbModel <- make.bisse(simBisseData, simBisseData$tip.state)
p <- starting.point.bisse(simBisseData)
nbModelMLFit <- find.mle(nbModel, p)

rbind(real = simPars, estimated = round(coef(nbModelMLFit), 2))


cnbModel <- constrain(nbModel, lambda1 ~ lambda0)
cnbModel <- constrain(cnbModel, mu1 ~ mu0)

cnbModelMLFit <- find.mle(cnbModel, p[c(-1, -3)])


anova(nbModelMLFit, constrained = cnbModelMLFit)

prior <- make.prior.exponential(1/(2 * 0.4))

mcmcRun <- mcmc(nbModel, nbModelMLFit$par, nsteps = 1000, prior = prior, w = 0.1, print.every = 100)


col <- c("blue", "red")
profiles.plot(mcmcRun[, c("lambda0", "lambda1")], col.line = col, las = 1, legend = "topright")


profiles.plot(mcmcRun[, c("mu0", "mu1")], col.line = col, las = 1, legend = "topright")


# looks like speciation rate differs, but not extinction. can we confirm?
sum(mcmcRun$lambda0 > mcmcRun$lambda1)/length(mcmcRun$lambda1)


sum(mcmcRun$mu0 > mcmcRun$mu1)/length(mcmcRun$mu1)

###############################################
#################### HiSSE ####################
###############################################
library(hisse)
library(geiger)

#setwd("~/Desktop/Research Projects/Macroevolution/Demo_Scripts/data")
grunt.tree <- read.tree("./data/grunts.phy")

grunt.data <- read.csv("./data/grunts.csv", row.names = 1)
head(grunt.data)

name.check(grunt.tree, grunt.data)

hd<-cbind(rownames(grunt.data), grunt.data[,"habitat"])

#####????
trans.rates.hisse <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
trans.rates.hisse[!is.na(trans.rates.hisse) & !trans.rates.hisse == 0] = 1
null.2.hisse <- hisse(grunt.tree, hd, hidden.states=TRUE, turnover=c(1,1,2,2), eps=c(1,1,2,2), trans.rate=trans.rates.hisse)

null.logL <- null.2.hisse$loglik
null.AIC <- null.2.hisse$AIC


trans.rates.bisse <- TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse

pp.bisse.no.hidden <- hisse(grunt.tree, hd, hidden.states=FALSE, turnover=c(1,2), eps=c(1,2), trans.rate=trans.rates.bisse)

pp.bisse.no.hidden

alt.logL <- pp.bisse.no.hidden$loglik
alt.AIC <- pp.bisse.no.hidden$AIC

logL <- c(null.logL, alt.logL)
AIC <- c(null.AIC, alt.AIC)

res <- as.data.frame(logL, row.names = c("null", "BiSSE"))
res$AIC <- AIC
res

trans.rates.bisse <- TransMatMaker.old(hidden.states=FALSE)
bisse.null <- hisse(grunt.tree, hd, hidden.states=FALSE, turnover=c(1,1), eps=c(1,1), trans.rate=trans.rates.bisse)

bisse.null

bisse.null.logL <- bisse.null$loglik
bisse.null.AIC <- bisse.null$AIC


res <- rbind(data.frame(logL = bisse.null.logL, AIC = bisse.null.AIC), res)

row.names(res) <- c("bisse null", "hisse null", "habitat dependent")
res


trans.rates.hisse <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(2,3,5,7,8,9,10,12))


hisse.grunts <- hisse(grunt.tree, hd, hidden.states=TRUE, turnover=c(1,2,0,3), eps=c(1,2,0,3), trans.rate=trans.rates.hisse)


hisse.logL <- hisse.grunts$loglik
hisse.AIC <- hisse.grunts$AIC
res <- rbind(res, data.frame(logL = hisse.logL, AIC = hisse.AIC, row.names = "habitat with hidden state"))

res
