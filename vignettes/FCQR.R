####################################
# Simulation1.
####################################
library(FCQR);
library(xtable); library(survival); library(fda);
library(quantreg); library(ggplot2); library(MASS);

n <- 200; v <- 1; sdx <- 0.2; sdz <- 1; lam0 <- 5;

t <- seq(0, 1, 0.01);  m <- length(t); deltat <- (diff(t)[1]);
nk <- 50; c0 <- 1;
b0 <- c(1, 2); lb0 = length(b0);

grd <- seq(0.1, 0.9, 0.1/5);
qv <- c(which(eleq(0.3, grd)), which(eleq(0.4, grd)), which(eleq(0.5, grd)), 
        which(eleq(0.6, grd)), which(eleq(0.7, grd)));

# ---------------------------------------------------------#
# creat B-spline basis.
nbasis <- 5; norder <- 4;
basis <- create.bspline.basis(range(t),nbasis = nbasis, norder = norder)
nbasis <- basis$nbasis
valueBasis <-eval.basis(t, basis)
# ---------------------------------------------------------#

nsim = 10;
Base0.m <- Sigb0.m <- array(0, dim = c(nsim, length(qv), lb0));
Bfun.m <- Sigfun.m <- array(0, dim = c(nsim, m, length(qv)))
cdata <- matrix(0, ncol = nsim, nrow = n);
itr <- 0; tt <- 1;
repeat{
  itr <- itr +1
  # ---------------------------------------------------------#
  # Generate Sim function.
  temp <- simfun(t, v, n, sdx, sdz, b0,  lam0, c0);
  Zb <- temp$Zb # baseline covariate.
  time <- temp$Timedata # survival time and censoring time.
  Zs <- temp$mZ # functional covariate.
  cdata[, itr] <- 1*(time[, 1] <= time[, 2])
  # ---------------------------------------------------------#
  res.itr <- try(Fcqr(time, Zs, Zb, t, grd, qv, valueBasis));
  if(class(res.itr) != 'try-error'){
    Base0.m[itr, , 1:lb0] <- res.itr$base0;
    Sigb0.m[itr, ,1:lb0] <- res.itr$sigb0;
    Bfun.m[itr, , ] <- t(res.itr$Bfun);
    Sigfun.m[itr, , ] <- t(res.itr$sigfun);
    tt <- tt + 1;
  }
  if(tt%%100==0) {print(tt);print(itr);print(Sys.time())}
  if(tt == nsim + 1){
    break
  }
}

cenRate = 1 - mean(colMeans(cdata));
#**********************************#
# The result of coefficient function.
#**********************************#
alfun <- simBX(t, v, nk)$B0
plot.m <- alfun.be(alfun, Sigfun.m, Bfun.m, grd, qv, sdx);
tau.id <- 3;
resplot6 <- alphaplot.be(t, grd, qv, tau.id, plot.m, -10, 1);
resplot6;
#**********************************#
# The result of Beta0
#**********************************#
lgrd <- length(grd);
btautrue <- cbind(rep(b0[1], lgrd), b0[1]*qnorm(grd, 0, sdx))
b0res <- rbind(b0res.be(Base0.m[,,1], Sigb0.m[,,1], btautrue[,1], qv),
               b0res.be(Base0.m[,,2], Sigb0.m[,,2], btautrue[,2], qv));
colnames(b0res) = c("BIAS", "SD", "SE", "CP");
rownames(b0res) = paste0(rep(c("$\\hat\\beta_1(","$\\hat\\beta_2("), 5),
                         rep(grd[qv], each  = 2),rep(")$", 10));
b0res;
####################################
# Simulation2: Mimic the real data.
####################################
rm(list = ls())
library(FCQR);
library(xtable); library(survival); library(fda);
library(quantreg); library(ggplot2); library(MASS);

t <- time; zsm <- MSdata[, 3:98]; 
c0 <- 10;  sdx <- 0.2; sdz <- 1; beta1 <- 0.1; norder <- 4;
# ------------------------------------------------------------
grd <- seq(0.1, 0.9, 0.005/5);
qv <- c(which(eleq(0.3, grd)), which(eleq(0.4, grd)), which(eleq(0.5, grd)), 
        which(eleq(0.6, grd)), which(eleq(0.7, grd)));

bsnb <- c(4, 5, 6, 7, 8);
nsim <- 10; lam0 <- 1 # for 20% censoring rate.
resop = optbasis(nsim, bsnb, zsm, lam0, c0, beta1, sdz, sdx,
                 grd, t, norder, qv);
resoptb = resop$inopt;
rownames(resoptb) = c("IGACV", "IBIC")
colnames(resoptb) = bsnb;
resoptb;

# -------------------------------------------------------------
nbasis <- 6;
btautrue <- beta1*qnorm(grd, 0, sdx);
# -------------------------------------------------------------
# report tau equals 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8;
grd <- seq(0.1, 0.9, 0.1/5);
qv <- c(which(eleq(0.2, grd)), 
        which(eleq(0.3, grd)), which(eleq(0.4, grd)),
        which(eleq(0.5, grd)), which(eleq(0.6, grd)), 
        which(eleq(0.7, grd)), which(eleq(0.8, grd)));
# -------------------------------------------------------------
SIMRES <- Mimicfun(zsm, lam0, c0, beta1, sdz, sdx, nsim,
                   grd, t, nbasis, norder, qv);

Base0.m <- SIMRES$Base0.m;
Sigb0.m <- SIMRES$Sigb0.m;
Sigfun.m <- SIMRES$Sigfun.m;
Bfun.m <- SIMRES$Bfun.m;

b0res <- b0res.be(Base0.m, Sigb0.m, btautrue, qv)
colnames(b0res) = c("BIAS", "SD", "SE", "CP");
rownames(b0res) = seq(0.2, 0.8, 0.1);

n <- nrow(MSdata);
alfun <- SimMimicfun(1, zsm, t, n, sdz, sdx, c0, beta1, lam0)$alp;
plotres <- alfun.be(alfun, Sigfun.m, Bfun.m, grd, qv, sdx);
tau.id <- 4;
mimicplot <- alphaplot.be(t, grd, qv, tau.id, plotres, -16, 16);

b0res;
mimicplot;
####################################
# Real data
####################################
rm(list=ls());
library(FCQR);
# --------------------------------------------------
library(imputeTS); library(survival); library(survminer); library(fda);
library(quantreg); library(statmod); library(xtable); library(Cairo);
# --------------------------------------------------

a0 <- 0; b0 <- 1; gn <- 20;
norder <- 4;
grd <- seq(0.1, 0.9, 0.005/5); lgrd <- length(grd);
qv <- c(which(eleq(0.2, grd)), which(eleq(0.3, grd)), which(eleq(0.4, grd)), which(eleq(0.5, grd)), 
        which(eleq(0.6, grd)), which(eleq(0.7, grd)), which(eleq(0.8, grd)));
# -----------------------------------------------------#
# using GACV method to determine optimal nbasis
CBA <- c(1, 2, 3, 4, 5)
BICR <- GACVR <- matrix(0, ncol=length(CBA), nrow=lgrd)
NBS <- rep(0,length(CBA))

n <- dim(MSdata)[1]
cn <- (log(n))^(1/2)
for(i in 1:length(CBA)){
  cba <- CBA[i];
  res <- predfun(MSdata, a0, b0, gn, time, norder, cba);
  TM <- res$TM; DEL <- res$DEL;
  nbasis <- res$nbasis; ns <- nbasis + 1;
  pre.y <- res$pre.y;
  bicres <- bicfun(n, ns, cn, TM, pre.y, DEL, grd);
  GACVR[, i] <- bicres$gacv;
  BICR[, i] <- bicres$bic;
  NBS[i] <- nbasis;
}

bicres <- optfun(grd, BICR);
gacvres <- optfun(grd, GACVR);
resoptReal <- rbind(gacvres[2, ], bicres[2, ]);
rownames(resoptReal) <- c("IGACV", "IBIC");
colnames(resoptReal) <- 4:8;

resoptReal

# -----------------------------------------------------#
cba <- 4;
covres <- intzfun0(MSdata, a0, b0, gn, time, norder, cba)
covdata <- covres$covdata
Z <- covres$Z
valueBasisT <- covres$valueBasisT;
weight0 <- rep(1, nrow(MSdata))
coefres <- coeffun0(valueBasisT, covdata,  MSdata,
                    grd, qv, weight0, Z, time)
# -----------------------------------------------------
# result for baseline covariate
base0 <- coefres$base0
sigb0 <- coefres$sigb0

c2 <- base0 + 1.96 * sqrt(sigb0)
c1 <- base0 - 1.96 * sqrt(sigb0)
rbind(c1,c2)

ci0 <- paste("(",round(c1,5)," ,",round(c2,5),")",sep="")
resbase0 <- rbind(round(base0, 5), ci0)
rownames(resbase0) <- c("Estimators", "95% CI ")
resbase = resbase0[, 2:5];
colnames(resbase) = seq(0.3, 0.6, 0.1);

resbase

# result for time dependent covariate
Bfun <- coefres$Bfun
sigfun <- coefres$sigfun
DA00 <- DA01 <- DA02 <- matrix(0, ncol = length(time), nrow = length(grd[qv]))
for(j in 1:length(grd[qv])){
  DA00 <- Bfun;
  bt <- DA00[j, ];
  vc <- sigfun[j, ];
  DA01[j, ] <- bt + 1.96 * sqrt(vc);
  DA02[j, ] <- bt - 1.96 * sqrt(vc);
}

tau.id <- 3
 
bt <- Bfun[tau.id, ]
vc <- sigfun[tau.id, ]
da00 <- bt;
da01 <- bt + 1.96 * sqrt(vc);
da02 <- bt - 1.96 * sqrt(vc)

yid <- 0.5
max.y <- yid;
min.y <- -yid
t <- seq(19.25, 43, length.out = 96);
realplot  = bandplot(da00, da01, da02, tau.id, t, min.y, max.y);
realplot

