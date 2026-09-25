# Preparation

library(ggplot2)
library(reshape2)
setwd("c:/Users/murasawa/Dropbox/r/15e")

# Functions

X2Y <- function(mX, cN_2) {
  cn   <- ncol(mX)
  cn_1 <- cn - cN_2
  mx_1 <- mX[, 1:cn_1]
  mx_2 <- mX[, (cn_1 + 1):cn]
  cbind(mx_1, diff(mx_2))
}

Y2dY <- function(mY, cN_2) {
  mdy  <- diff(mY)
  mdy0 <- scale(mdy, scale = F)
  cn   <- ncol(mY)
  cn_1 <- cn - cN_2
  cbind(mdy0[, 1:cn_1], mdy[, (cn_1 + 1):cn])
}

ApplyOLS <- function(mY, mX) {
  ct  <- nrow(mY)
  mxx <- crossprod(mX)
  mxy <- crossprod(mX, mY)
  mb  <- t(solve(mxx, mxy))
  me  <- mY - mX %*% t(mb)
  mee <- crossprod(me)
  ms  <- mee / ct
  return(list(mB = mb, mS = ms))
}

GetCompanionMatrix <- function(mPhi) {
  cn <- nrow(mPhi)
  ip <- ncol(mPhi) / cn
  ma <- rbind(mPhi, cbind(diag((ip - 1) * cn), matrix(0, (ip - 1) * cn, cn)))
  return(ma)
}

BN1 <- function(mPhi, mS) {
  cn <- nrow(mPhi)
  ip <- ncol(mPhi) / cn
  mc <- cbind(diag(cn), matrix(0, cn, (ip - 1) * cn))
  ma <- GetCompanionMatrix(mPhi)
  mw <- -mc %*% solve(diag(ip * cn) - ma, ma)
  t(mw %*% t(mS))
}

BN12 <- function(mPhi, mS, cN_2) {
  cn   <- nrow(mPhi)
  ip   <- ncol(mPhi) / cn
  mc   <- cbind(diag(cn), matrix(0, cn, (ip - 1) * cn))
  cn_1 <- cn - cN_2
  mc_1 <- mc[1:cn_1,]
  mc_2 <- mc[(cn_1 + 1):cn,]
  ma   <- GetCompanionMatrix(mPhi)
  miaa <- solve(diag(ip * cn) - ma, ma)
  mw_1 <- -mc_1 %*% miaa
  mw_2 <-  mc_2 %*% miaa %*% miaa
  mw   <- rbind(mw_1, mw_2)
  t(mw %*% t(mS))
}

# Data

iStart <- 1980.25
iEnd   <- 2013.25
iFreq  <- 4
mX     <- read.csv("jpmacror.csv")
#iEnd   <- 2013.5                     # if nominal interest rate
#mX     <- read.csv("jpmacron.csv")   # if nominal interest rate
mX     <- ts(mX, start = iStart - .25, frequency = iFreq)
cN     <- ncol(mX)
cN_2   <- 1
mY     <- X2Y(mX, cN_2)
#mY     <- mX   # if I(1) only
mY     <- window(mY, start = iStart, frequency = iFreq)
mdY    <- Y2dY(mY, cN_2)
#mdY    <- scale(diff(mY), scale = F)   # if I(1) only
plot(mX)
plot(mdY)
summary(mdY)
cT <- nrow(mdY)
iP <- 8

# OLS

ar(mdY, FALSE, iP, method = "ols", demean = FALSE)
mx <- NULL
for (i in 1:iP) {
	mx <- ts.union(mx, lag(mdY, -i))
}
ms_p   <- window(mdY, start = iStart + .25, end = iStart + iP / iFreq)
vs_p   <- c(t(ms_p))
my     <- window(mdY, start = iStart + .25 + iP / iFreq)
my     <- matrix(my, cT - iP, cN)
mx     <- window(mx, start = iStart + .25 + iP / iFreq, end = iEnd)
mx     <- matrix(mx, cT - iP, iP * cN)
lOLS   <- ApplyOLS(my, mx)
mPhi   <- lOLS$mB
mSigma <- lOLS$mS

ma       <- GetCompanionMatrix(mPhi)
ma.eigen <- eigen(ma)
vlambda  <- ma.eigen$values
abs(vlambda)

# B-N Decomposition

ms <- NULL
for (i in 0:(iP - 1)) {
  ms <- ts.union(ms, lag(mdY, -i))
}
ms   <- window(ms, iStart + iP / iFreq, iEnd)
ms   <- matrix(ms, cT - iP + 1, iP * cN)
mgap <- BN12(mPhi, ms, cN_2)
#mgap <- BN1(mPhi, ms)   # if I(1) only
mgap <- ts(mgap, start = iStart + iP / iFreq, frequency = iFreq)
mX   <- window(mX, iStart + iP / iFreq)
mnr  <- mX - mgap

# Stats

colnames(mnr)  <- colnames(mX)
colnames(mgap) <- colnames(mX)
summary(mgap)
cor(mgap)
pairs(data.frame(mgap))

# Plots

mTP <- read.csv("jprecess.csv")

mar <- data.frame(date = c(time(mX)), mX)
mar <- melt(mar, id.vars = "date")
levels(mar$variable)[levels(mar$variable) == "Pi"]  <- "inflation rate"
levels(mar$variable)[levels(mar$variable) == "r"]   <- "real interest rate"
levels(mar$variable)[levels(mar$variable) == "U"]   <- "unemployment rate"
levels(mar$variable)[levels(mar$variable) == "lnY"] <- "log output"
ggplot(mar, aes(x = date, y = value)) +
  annotate("rect", xmin = mTP[, 1], xmax = mTP[, 2], ymin = -Inf, ymax = Inf, fill = "grey") +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  geom_line() +
  theme(axis.title = element_blank())
ggsave(file = "15e_fig0.eps", width = 6, height = 6)

mX    <- data.frame(rate = "actual" , date = c(time(mX)) , mX)
mnr   <- data.frame(rate = "natural", date = c(time(mnr)), mnr)
mrate <- rbind(mX, mnr)
mrate <- melt(mrate, id.vars = c("rate", "date"))
ggplot(mrate, aes(x = date, y = value, linetype = rate)) +
  annotate("rect", xmin = mTP[, 1], xmax = mTP[, 2], ymin = -Inf, ymax = Inf, fill = "grey") +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  geom_line() +
  theme(axis.title = element_blank(), legend.position = "top")

mgap <- data.frame(date = time(mgap), mgap)
mgap <- melt(mgap, id.vars = "date")
levels(mgap$variable)[levels(mgap$variable) == "Pi"]  <- "inflation rate gap"
levels(mgap$variable)[levels(mgap$variable) == "r"]   <- "interest rate gap"
levels(mgap$variable)[levels(mgap$variable) == "U"]   <- "unemployment rate gap"
levels(mgap$variable)[levels(mgap$variable) == "lnY"] <- "output gap"
ggplot(mgap, aes(x = date, y = value)) +
  annotate("rect", xmin = mTP[, 1], xmax = mTP[, 2], ymin = -Inf, ymax = Inf, fill = "grey") +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  geom_hline(yintercept = 0) +
  geom_line() +
  theme(axis.title = element_blank())
#ggsave(file = "15e_fig1.eps", width = 6, height = 6)
ggsave(file = "15e_fig2.eps", width = 6, height = 6)
