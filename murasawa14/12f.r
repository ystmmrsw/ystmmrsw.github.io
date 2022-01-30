# Preparation

setwd("D:/Documents/Dropbox/r/12f")
#library(Matrix)
library(MCMCpack)
library(mvtnorm)
#options(digits = 10)
set.seed(1)
start <- proc.time()

# Functions

ApplyOLS <- function(mY, mX) {
	ct  <- nrow(mY)
	mxx <- crossprod(mX)
	mxy <- crossprod(mX, mY)
	mb  <- t(solve(mxx, mxy))
	me  <- mY - mX %*% t(mb)
	mee <- crossprod(me)
	ms  <- mee / ct
	return(list(mB = mb, mXX = mxx, mEE = mee, mS = ms))
}

SpecifyPrior <- function(mY, iP, dLambda, dN_0) {
	cn   <- ncol(mY)
	mm_0 <- matrix(0, cn, iP * cn)
	vs2  <- NULL
	for (i in 1:cn) {
		dy.ar <- ar(mY[, i], FALSE, iP, method = "ols", demean = FALSE)
		vs2   <- c(vs2, dy.ar$var.pred)
	}
	md_0 <- (1 / dLambda ^ 2) * diag((1:iP) ^ 2) %x% diag(vs2)
	ms_0 <- (dN_0 - cn - 1) * diag(vs2)
	return(list(mM_0 = mm_0, mD_0 = md_0, dN_0 = dN_0, mS_0 = ms_0))
}

SpecifyPosterior <- function(mM_0, mD_0, dN_0, mS_0, mY, mX) {
	ct       <- nrow(mY)
	lols     <- ApplyOLS(mY, mX)
	mphi_OLS <- lols$mB
	mxx      <- lols$mXX
	mee      <- lols$mEE
	mxxInv   <- solve(mxx)
	md_1     <- mxx + mD_0
	md_1Inv  <- solve(md_1)
	mm_1     <- (mphi_OLS %*% mxx + mM_0 %*% mD_0) %*% md_1Inv
	dn_1     <- ct + dN_0
	md_0Inv  <- solve(mD_0)
	ms_1     <- (mphi_OLS - mM_0) %*% solve(mxxInv + md_0Inv) %*% t(mphi_OLS - mM_0) + mee + mS_0
	return(list(mM_1 = mm_1, mD_1 = md_1, dN_1 = dn_1, mS_1 = ms_1))
}

LogMultiGamma <- function(cN, dAlpha) {
	(cN * (cN - 1) / 4) * log(pi) + sum(lgamma(dAlpha - (1: cN - 1) / 2))
}

ComputeLogMarginalLikelihood <- function(mD_0, dN_0, mS_0, mD_1, dN_1, mS_1) {
	cn          <- nrow(mS_0)
	dlprior     <- -(dN_0 / 2) * log(det(mS_0 / 2)) - (cn / 2) * log(det(mD_0)) + LogMultiGamma(cn, dN_0 / 2)
	dlposterior <- -(dN_1 / 2) * log(det(mS_1 / 2)) - (cn / 2) * log(det(mD_1)) + LogMultiGamma(cn, dN_1 / 2)
	dlml        <- dlposterior - dlprior
	return(dlml)
}

ComputeSDRatio <- function(mD_0Inv, dN_0, mS_0, mM_1, mD_1Inv, dN_1, mS_1) {
	cn    <- nrow(mM_1)
	ip    <- ncol(mM_1) / cn
	mR    <- diag(c(t(matrix(1, 1, ip) %x% (matrix(1, cn, cn) - diag(cn)))))
	mR    <- unique(mR)[-1,]
	vmu_1 <- c(t(mM_1))
	vx_1  <- c(mR %*% vmu_1)
	vx_0  <- 0 * vx_1
	dp_0  <- 0
	dp_1  <- 0
	cR    <- 1000
	for (r in 1:cR) {
		msigma_0 <- riwish(dN_0, mS_0)
		msigma_1 <- riwish(dN_1, mS_1)
		dp_0     <- dp_0 + dmvnorm(vx_0,, mR %*% (msigma_0 %x% mD_0Inv) %*% t(mR))
		dp_1     <- dp_1 + dmvnorm(vx_1,, mR %*% (msigma_1 %x% mD_1Inv) %*% t(mR))
	}
	return(dp_1 / dp_0)
}

GetCompanionMatrix <- function(mPhi) {
	cn <- nrow(mPhi)
	ip <- ncol(mPhi) / cn
	ma <- rbind(mPhi, cbind(diag((ip - 1) * cn), matrix(0, (ip - 1) * cn, cn)))
	return(ma)
}

ComputeInitStateVariance <- function(mPhi, mSigma) {
	ma     <- GetCompanionMatrix(mPhi)
	cm     <- nrow(ma)
	cn     <- nrow(mPhi)
	mb     <- rbind(t(chol(mSigma)), matrix(0, cm - cn, cn))
	mc     <- diag(cm ^ 2) - ma %x% ma
#	mc     <- Matrix(mc)
	mgamma <- matrix(solve(mc) %*% c(mb %*% t(mb)), cm, cm)
	return(mgamma)
}

DrawPhi <- function(mSigma, mM_1, mD_1Inv) {
	cn          <- nrow(mM_1)
	ip          <- ncol(mM_1) / cn
	ck          <- ip * cn ^ 2
	vmean       <- c(t(mM_1))
	mvariance   <- mSigma %x% mD_1Inv
	dlambda_max <- 1
	while (dlambda_max >= 1) {
		vphi        <- rmvnorm(1, vmean, mvariance)
		mphi        <- t(matrix(vphi, ip * cn, cn))
		ma          <- GetCompanionMatrix(mphi)
		ma.eigen    <- eigen(ma)
		vlambda     <- ma.eigen$values
		dlambda_max <- max(abs(vlambda))
	}
	return(mphi)
}

ComputePmove <- function(mGamma, mGamma_star, vS_0) {
	dlnf0      <- dmvnorm(vS_0,, mGamma,      log = TRUE)
	dlnf0_star <- dmvnorm(vS_0,, mGamma_star, log = TRUE)
	dalpha     <- min(exp(dlnf0_star - dlnf0), 1)
	return(dalpha)
}

# Data

iStart <- 1948
iEnd   <- 2011.75
iFreq  <- 4
mY     <- read.csv("usmacroq.csv")
mY     <- ts(mY, start = iStart, frequency = iFreq)
mdY    <- diff(mY)
#plot(mY)
#plot(mdY)
summary(mdY)
mdY <- scale(mdY, scale = F)
cN  <- ncol(mdY)
cT  <- nrow(mdY)
iP  <- 12

# OLS

ar(mdY, FALSE, iP, method = "ols", demean = FALSE)
mx <- NULL
for (i in 1:iP) {
	mx <- ts.union(mx, lag(mdY, -i))
}
ms_p       <- window(mdY, iStart + .25, iStart + iP / iFreq)
vs_p       <- c(t(ms_p))
my         <- window(mdY, iStart + .25 + iP / iFreq, iEnd)
my         <- matrix(my, cT - iP, cN)
mx         <- window(mx, iStart + .25 + iP / iFreq, iEnd)
mx         <- matrix(mx, cT - iP, iP * cN)
lOLS       <- ApplyOLS(my, mx)
mPhi_OLS   <- lOLS$mB
mSigma_OLS <- lOLS$mS

# Empirical Bayes

lnML <- function(x) {
	dlambda    <- x[1]
	dn_0       <- x[2]
#	dlambda    <- x
#	dn_0       <- cN + 2
	lprior     <- SpecifyPrior(mdY, iP, dlambda, dn_0)
	mm_0       <- lprior$mM_0
	md_0       <- lprior$mD_0
	dn_0       <- lprior$dN_0
	ms_0       <- lprior$mS_0
	lposterior <- SpecifyPosterior(mm_0, md_0, dn_0, ms_0, my, mx)
	md_1       <- lposterior$mD_1
	dn_1       <- lposterior$dN_1
	ms_1       <- lposterior$mS_1
	dlml       <- ComputeLogMarginalLikelihood(md_0, dn_0, ms_0, md_1, dn_1, ms_1)
	return(-dlml)
}
lhyper  <- optim(c(1, cN + 2), lnML, method = "L-BFGS-B", lower = c(.0001, cN - 1))
dLambda <- lhyper$par[1]
dN_0    <- lhyper$par[2]
#lhyper  <- optimize(lnML, c(0, 100))
#dLambda <- lhyper$minimum
#dN_0    <- cN + 2

# Prior

lPrior  <- SpecifyPrior(mY, iP, dLambda, dN_0)
mM_0    <- lPrior$mM_0
mD_0    <- lPrior$mD_0
dN_0    <- lPrior$dN_0
mS_0    <- lPrior$mS_0
mD_0Inv <- solve(mD_0)

# Posterior

lPosterior <- SpecifyPosterior(mM_0, mD_0, dN_0, mS_0, my, mx)
mM_1       <- lPosterior$mM_1
mD_1       <- lPosterior$mD_1
dN_1       <- lPosterior$dN_1
mS_1       <- lPosterior$mS_1
mD_1Inv    <- solve(mD_1)

# Log marginal likelihood

dlnML <- ComputeLogMarginalLikelihood(mD_0, dN_0, mS_0, mD_1, dN_1, mS_1)

# Bayes factor of AR vs VAR (S-D density ratio)

dbf <- ComputeSDRatio(mD_0Inv, dN_0, mS_0, mM_1, mD_1Inv, dN_1, mS_1)

# Initial value

msigma <- riwish(dN_1, mS_1)
mphi   <- DrawPhi(msigma, mM_1, mD_1Inv)
#mgamma <- ComputeInitStateVariance(mphi, msigma)

# M-H algorithm

cBurn    <- 0
cR       <- 10000
mphi_s   <- mcmc(matrix(, cR, iP * cN ^ 2))
msigma_s <- mcmc(matrix(, cR, cN * (cN + 1) / 2))
for (r in (1 - cBurn):cR) {

	# Draw parameters

	msigma_star <- riwish(dN_1, mS_1)
	mphi_star   <- DrawPhi(msigma_star, mM_1, mD_1Inv)
#	mgamma_star <- ComputeInitStateVariance(mphi_star, msigma_star)
#	dalpha      <- ComputePmove(mgamma, mgamma_star, vs_p)
#	if (runif(1) <= dalpha) {
		mphi   <- mphi_star
		msigma <- msigma_star
#		mgamma <- mgamma_star
#	}

	# Save draws

	if (r >= 1) {
		mphi_s[r,]   <- c(t(mphi))		
		msigma_s[r,] <- vech(msigma)
	}
	print(r)
}

# Diagnostics

#plot(mphi_s)
#plot(msigma_s)
summary(mphi_s)
summary(msigma_s)
#autocorr.plot(mphi_s)
#autocorr.plot(msigma_s)
geweke.diag(mphi_s)
geweke.diag(msigma_s)
rejectionRate(mphi_s)
rejectionRate(msigma_s)

# B-N Decomposition

amphi   <- array(dim = c(cN, iP * cN, cR))
amsigma <- array(dim = c(cN, cN, cR))
for (r in 1:cR) {
	amphi[,, r]   <- t(matrix(mphi_s[r,], iP * cN, cN))
	amsigma[,, r] <- xpnd(msigma_s[r,], cN)
}
apply(amphi, 1:2, mean)
apply(amsigma, 1:2, mean)
ms <- NULL
for (i in 0:(iP - 1)) {
	ms <- ts.union(ms, lag(mdY, -i))
}
ms      <- window(ms, iStart + iP / iFreq, iEnd)
ms      <- matrix(ms, cT - iP + 1, iP * cN)
amgap   <- array(dim = c(cN, cT - iP + 1, cR))
amdgap  <- array(dim = c(cN, cT - iP, cR))
amd2gap <- array(dim = c(cN, cT - iP - 1, cR))
amcorr  <- array(dim = c(cN, cN, cR))
mC      <- cbind(diag(cN), matrix(0, cN, (iP - 1) * cN))
for (r in 1:cR) {
	ma            <- GetCompanionMatrix(amphi[,, r])
	mw            <- -mC %*% solve(diag(iP * cN) - ma) %*% ma
	amgap[,, r]   <- mw %*% t(ms)
	amdgap[,, r]  <- t(diff(t(amgap[,, r])))
	amd2gap[,, r] <- t(diff(t(amgap[,, r]), 2))
	amcorr[,, r]  <- cor(t(amgap[,, r]))
}
mgap_med   <- t(apply(amgap, 1:2, median))
mgap_lower <- t(apply(amgap, 1:2, quantile, prob = .025))
mgap_upper <- t(apply(amgap, 1:2, quantile, prob = .975))
mgap_med   <- ts(mgap_med, start = iStart + iP / iFreq, frequency = iFreq)
mgap_lower <- ts(mgap_lower, start = iStart + iP / iFreq, frequency = iFreq)
mgap_upper <- ts(mgap_upper, start = iStart + iP / iFreq, frequency = iFreq)
mY         <- window(mY, iStart + iP / iFreq, iEnd)
mnr_med    <- mY - mgap_med
mnr_lower  <- mY - mgap_upper
mnr_upper  <- mY - mgap_lower
amgapDI    <- (amgap > 0)
mgapDI     <- t(apply(amgapDI, 1:2, mean))
mgapDI     <- ts(mgapDI, start = iStart + iP / iFreq, frequency = iFreq)
amets      <- (amd2gap[, 1:(cT - iP - 3),] > 0) * (amdgap[, 2:(cT - iP - 2),] > 0) * (amgap[, 3:(cT - iP - 1),] > 0) * (amdgap[, 3:(cT - iP - 1),] < 0) * (amd2gap[, 3:(cT - iP - 1),] < 0)
amcts      <- (amd2gap[, 1:(cT - iP - 3),] < 0) * (amdgap[, 2:(cT - iP - 2),] < 0) * (amgap[, 3:(cT - iP - 1),] < 0) * (amdgap[, 3:(cT - iP - 1),] > 0) * (amd2gap[, 3:(cT - iP - 1),] > 0)
mets       <- t(apply(amets, 1:2, mean))
mcts       <- t(apply(amcts, 1:2, mean))
mets       <- ts(mets, start = iStart + iP / iFreq + .5, frequency = iFreq)
mcts       <- ts(mcts, start = iStart + iP / iFreq + .5, frequency = iFreq)

# Stats

mY[, 4]         <- exp(mY[, 4])
mnr_med[, 4]    <- exp(mnr_med[, 4])
mnr_lower[, 4]  <- exp(mnr_lower[, 4])
mnr_upper[, 4]  <- exp(mnr_upper[, 4])
mgap_med[, 4]   <- mY[, 4] - mnr_med[, 4];
mgap_lower[, 4] <- mY[, 4] - mnr_upper[, 4];
mgap_upper[, 4] <- mY[, 4] - mnr_lower[, 4];
summary(mgap_med)

# Plots

mTP <- read.csv("usrecess.csv")
ShadeRecession <- function(l, u, mTP) {
	rect(mTP$Peak, l, mTP$Trough, u, col = "gray", border = "NA")
}
DrawNaturalRate <- function(vy_ar, vy_nr, stitle) {
	ts.plot(vy_ar, vy_nr, main = stitle, lty = c(1, 2))
	ShadeRecession(-1, 10, mTP)
	par(new = T)
	ts.plot(vy_ar, vy_nr, main = stitle, lty = c(1, 2))
}
DrawGap <- function(vy500, vy025, vy975, stitle) {
	ts.plot(vy500, vy025, vy975, 0, main = stitle, lty = c(1, 2, 2, 1))
	ShadeRecession(-1, 1, mTP)
	par(new = T)
	ts.plot(vy500, vy025, vy975, 0, main = stitle, lty = c(1, 2, 2, 1))
}
DrawScatter <- function(vx, vy, sxlab, sylab) {
	plot(vx, vy, xlim = c(-.03, .03), ylim = c(-.03, .03), xlab = sxlab, ylab = sylab)
	grid(2, 2, lty = 1)
}
DrawPDF <- function(vx, stitle) {
	plot(density(vx), xlim = c(-1, 1), ylim = c(0, 5), main = stitle)
	grid(2, 1, lty = 1)
}
DrawProbPositiveGap <- function(vy, stitle) {
	ts.plot(vy, .5, ylab = "probability of positive gap",
		main = stitle, lty = c(1, 1))
	ShadeRecession(-1, 2, mTP)
	par(new = T)
	ts.plot(vy, .5, ylab = "probability of positive gap",
		main = stitle, lty = c(1, 1))
}
DrawProbRecession <- function(vy, stitle) {
	ts.plot(vy, .5, ylim = c(0, 1), ylab = "probability of recession",
	main = stitle, lty = c(1, 1))
	ShadeRecession(-1, 2, mTP)
	par(new = T)
	ts.plot(vy, .5, ylim = c(0, 1), ylab = "probability of recession",
	main = stitle, lty = c(1, 1))
}
DrawProbRevival <- function(vy, stitle) {
	ts.plot(vy, .5, ylim = c(0, 1), ylab = "probability of revival",
		main = stitle, lty = c(1, 1))
	ShadeRecession(-1, 2, mTP)
	par(new = T)
	ts.plot(vy, .5, ylim = c(0, 1), ylab = "probability of revival",
		main = stitle, lty = c(1, 1))
}

par(mfrow = c(4, 1))
DrawNaturalRate(mY[, 1], mnr_med[, 1], "output")
DrawNaturalRate(mY[, 2], mnr_med[, 2], "inflation rate")
DrawNaturalRate(mY[, 3], mnr_med[, 3], "real interest rate")
DrawNaturalRate(mY[, 4], mnr_med[, 4], "unemployment rate")
dev.copy2eps(file = "12f_fig1.eps", width = 9, height = 13.5)

DrawGap(mgap_med[, 1], mgap_lower[, 1], mgap_upper[, 1], "output")
DrawGap(mgap_med[, 2], mgap_lower[, 2], mgap_upper[, 2], "inflation rate")
DrawGap(mgap_med[, 3], mgap_lower[, 3], mgap_upper[, 3], "real interest rate")
DrawGap(mgap_med[, 4], mgap_lower[, 4], mgap_upper[, 4], "unemployment rate")
dev.copy2eps(file = "12f_fig2.eps", width = 9, height = 13.5)

par(mfrow = c(3, 2))
DrawScatter(mgap_med[, 1], mgap_med[, 2], "output gap", "inflation rate gap")
DrawPDF(amcorr[2, 1,], "correlation coefficient")
DrawScatter(mgap_med[, 4], mgap_med[, 2], "unemployment rate gap", "inflation rate gap")
DrawPDF(amcorr[2, 4,], "correlation coefficient")
DrawScatter(mgap_med[, 4], mgap_med[, 1], "unemployment rate gap", "output gap")
DrawPDF(amcorr[1, 4,], "correlation coefficient")
dev.copy2eps(file = "12f_fig3.eps", width = 9, height = 13.5)

par(mfrow = c(4, 1))
DrawProbPositiveGap(mgapDI[, 1], "output")
DrawProbPositiveGap(mgapDI[, 2], "inflation rate")
DrawProbPositiveGap(mgapDI[, 3], "real interest rate")
DrawProbPositiveGap(mgapDI[, 4], "unemployment rate")
dev.copy2eps(file = "12f_fig4.eps", width = 9, height = 13.5)

DrawProbRecession(mets[, 1], "output")
DrawProbRecession(mets[, 2], "inflation rate")
DrawProbRecession(mcts[, 3], "real interest rate")
DrawProbRecession(mcts[, 4], "unemployment rate")
dev.copy2eps(file = "12f_fig5.eps", width = 9, height = 13.5)

DrawProbRevival(mcts[, 1], "output")
DrawProbRevival(mcts[, 2], "inflation rate")
DrawProbRevival(mets[, 3], "real interest rate")
DrawProbRevival(mets[, 4], "unemployment rate")
dev.copy2eps(file = "12f_fig6.eps", width = 9, height = 13.5)

write.csv(mgap_med[, 1], "temp.csv")

(proc.time() - start) / 60
