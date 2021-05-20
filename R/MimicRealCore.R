SimMimicfun <- function(itr, zsm, t, n, sdz, sdx, c0, beta1, lam0){

  set.seed(itr);
  z1m <- apply(zsm, 1, min);
  m <- length(t)
  deltat <- diff(t)[1];

  alp <- 5*sin(0.5*pi + 3.25*pi*t)
  Alp <- matrix(rep(alp, n), nrow = n, byrow = TRUE)

  zrand <- t(replicate(n, runif(m, -c0, c0)));
  Zs0 <- zsm + zrand;
  z1rand <- runif(n, -c0, c0)
  z10 <- z1m + z1rand;

  Zs <- Zs0/(sd(Zs0));
  z1 <- z10/(sd(z10));

  Zsa <- Zs*Alp;
  izs <- apply(Zs, 1, sum)*deltat;
  izsa <- apply(Zsa, 1, sum)*deltat;

  xi <- rnorm(n, 0, sdx);

  logT <- (beta1*z1 + izs)*xi + izsa

  # if(lam0==0){C <- exp(logT)}else{C <- rexp(n, lam0);}
  if(lam0==0){C <- exp(logT)}else{C <- runif(n, 0, lam0);}

  # msTime[i, ] <- ind
  Timedata <- cbind(exp(logT), C);

  return(list(Timedata = Timedata, alp = alp, z1 = z1, Zs = Zs))
}


FcqrMimicfun <- function(time, Zs, Zb, t, grd, qv, valueBasis,
                         int.id = 0, sdid = 0){
  # ----------------------------------------------------------- #
  # lb0: the number of non-zero of beta1
  lb0 = ncol(Zb);
  # deltat: the increment of t.
  deltat <- (diff(t)[1]);
  m <- length(t);
  lgrd <- length(grd);
  #---------------------------------------------------------#
  # event time and censor indicator denoted as dat.yc
  surtime <- cbind(pmin(time[, 1], time[, 2]), (time[, 1] <= time[, 2]))

  # estimate the function using Bspline and integrate
  #int_0 ^T Z(s)Br(s) ds denote ad covMat
  covMat <- Zs%*%valueBasis*deltat

  # ix <- which(apply(abs(covMat), 2, mean)> 1e-5)
  # dat <- cbind(dsurtime,  covMat[,ix])
  # -----------------------------------------------------------
  nbasis <- dim(valueBasis)[2]
  covdata <- cbind(Zb, covMat);
  T <- log(surtime[, 1])

  # whether a the covariates are standarded.
  if(sdid==1){
    std<- matrix(apply(covdata, 2,sd), nrow(covdata), ncol(covdata), byrow = T)
    covdata <- covdata/std
    T <- T/sd(T)
  }

  # estimate the coefficients
  # -----------------------------------------------------------
  # if(lb0==0){dat.cov <- dat[, 3:ncol(dat) ]}else{
  #   dat.cov <- cbind(Zb, dat[, 3:ncol(dat) ])}

  if(int.id==0){
    fit <- crq(Surv(T, surtime[,2])~.-1, data = data.frame(covdata),  taus=grd, method= "PengHuang")
  }else{
    fit <- crq(Surv(T, surtime[,2])~., data = data.frame(covdata),  taus=grd, method= "PengHuang")
  }

  # ------------------------------------------------------------
  # estimate the coefficients
  coefest <- t(coef(fit,grd))
  rescoef <- t(coefest[, 1:(nbasis + 1)])

  coef2 <- rescoef[2:dim(rescoef)[1], ]
  res <- valueBasis%*%coef2;
  beta0 <- rescoef[1, ];
  #-----------------------------------------------#
  # predict value
  pz1 <- (matrix(rep( covdata[, 1], lgrd), ncol = lgrd))*(matrix(rep(beta0, length(Zb)), ncol = lgrd, byrow = TRUE))
  pzs <- covdata[, -1]%*%coef2
  pre.y <- exp(pz1 +  pzs)
  #-----------------------------------------------#
  # for baseline covariate
  sigb0 <- base0 <- rep(NA, length(grd[qv]));
  # for covariate function
  sigfun <- Bfun <- matrix(NA, ncol = m, nrow = length(grd[qv]))

  # m <- length(t)
  if(sum(is.na(coefest)) < (length(coefest)/2)){
    #95% confidence intercal only report tau =0.2-0.8
    fc <- summary(fit, taus = grd[qv], alpha = .05, se="boot", R = 200, covariance=TRUE)
    for(tt in 1:length(fc)){
      sigfun[tt, ] <- diag(valueBasis %*% fc[[tt]]$cov[-1, -1] %*% t(valueBasis))
      sigb0[tt] <- fc[[tt]]$cov[1, 1];
      cf <- fc[[tt]]$coefficients[, 1]
      Bfun[tt, ] <- valueBasis %*% cf[-1];
      base0[tt] <- cf[1]
    }
  }
  # ------------------------------------------------------------
  return(list(pre.y = pre.y, covdata = covdata, T = T,
              base0 = base0, Bfun = Bfun, sigb0 = sigb0, sigfun = sigfun));
}

Mimicfun <- function(zsm, lam0, c0, beta1, sdz, sdx, nsim,
                     grd, t, nbasis, norder, qv){


  # --------------------------
  lgrd <- length(grd);
  m <- length(t);
  n <- dim(zsm)[1];
  cn <- (log(n))^(1/2);
  # --------------------------

  basis <- create.bspline.basis(range(t),nbasis = nbasis, norder = norder)
  nbasis <- basis$nbasis
  valueBasis <-eval.basis(t, basis)

  gacvr <- bicr <- aicr <- matrix(0, nrow = nsim, ncol = lgrd);
  cdata <- matrix(0, ncol = nsim, nrow = n);
  Base0.m <- Sigb0.m <- matrix(0, nrow = nsim, ncol = length(qv));
  Bfun.m <- Sigfun.m <- array(0, dim = c(nsim, m, length(qv)));

  for(itr in 1:nsim){

    resMat <-  matrix(NA, n, m)
    covMat <- matrix(NA, n, nbasis)

    temp <- SimMimicfun(itr, zsm, t, n, sdz, sdx, c0, beta1, lam0);
    time <- temp$Timedata;
    Zs <- temp$Zs;
    Zb <- temp$z1;

    cdata[, itr] <- 1*(time[, 1] <= time[, 2])
    res.itr <- try(FcqrMimicfun(time, Zs, Zb, t, grd, qv, valueBasis));
    # ---------------------------------------------------------#

    if(class(res.itr) != 'try-error'){

      Base0.m[itr, ] <- res.itr$base0;
      Sigb0.m[itr, ] <- res.itr$sigb0;
      Bfun.m[itr, , ] <- t(res.itr$Bfun);
      Sigfun.m[itr, , ] <- t(res.itr$sigfun);

      # simcoef1[itr, 1:lb0, ] <- res.itr$res1;
      # simcoeffun[itr, , ] <- res.itr$res2

      T <- res.itr$T
      TM <- matrix(rep(T, lgrd), ncol = lgrd);
      pre.y <- res.itr$pre.y

      del <- 1*(time[, 1] <= time[, 2]);
      DEL <- matrix(rep(del, lgrd), ncol = lgrd);

      ns <- nbasis + 1
      res <- bicfun(n, ns, cn, TM, pre.y, DEL, grd)

      gacvr[itr, ] <- res$gacv
      bicr[itr, ] <- res$bic
      aicr[itr, ] <- res$aic

    }
    if(itr%%500==0) {print(itr);print(Sys.time())}
  }

  return(list(bicr=bicr, gacvr = gacvr, aicr = aicr, cdata = cdata, Base0.m=Base0.m, Sigb0.m = Sigb0.m,
              Bfun.m=Bfun.m, Sigfun.m = Sigfun.m));

}

optbasis <- function(nsim, bsnb, zsm, lam0, c0, beta1, sdz, sdx,
                     grd, t,  norder, qv){
  n <- dim(zsm)[1];
  lgrd <- length(grd);

  cenrate <- rep(0, length(bsnb))
  GACVR <- BICR <- AICR <- array(0, dim = c(nsim, length(bsnb), lgrd))
  GACVR.IN <- BICR.IN <- AICR.IN <- matrix(0, nrow = nsim, ncol = length(bsnb))
  for(nb in 1:length(bsnb)){

    nbasis <- bsnb[nb];
    SIMRES <- Mimicfun(zsm, lam0, c0, beta1, sdz, sdx, nsim,
                       grd, t, nbasis, norder, qv)

    cenrate[nb] <- 1 - mean(colMeans(SIMRES$cdata));

    aa1 <- SIMRES$gacvr;
    aa2 <- SIMRES$bicr;
    aa3 <- SIMRES$aicr;

    aa1[which(is.infinite(aa1), arr.ind = TRUE)] <- NA;
    aa2[which(is.infinite(aa2), arr.ind = TRUE)] <- NA;
    aa3[which(is.infinite(aa3), arr.ind = TRUE)] <- NA;

    GACVR[, nb, ] <- aa1;
    BICR[, nb, ] <- aa2;
    AICR[, nb, ] <- aa3;

    # integrate for tau
    del.tau <- diff(grd)[1]
    GACVR.IN[, nb] <- apply(aa1*del.tau, 1, sum, na.rm = TRUE);
    BICR.IN[, nb] <- apply(aa2*del.tau, 1, sum, na.rm = TRUE);
    AICR.IN[, nb] <- apply(aa3*del.tau, 1, sum, na.rm = TRUE)
    # print(nb);
  }

  # inopt <- rbind(countint(bsnb, GACVR.IN),
  #                countint(bsnb, BICR.IN),
  #                countint(bsnb, AICR.IN))
  #
  # tauopt <- rbind(counttau(GACVR, lgrd, bsnb),
  #                 counttau(BICR, lgrd, bsnb),
  #                 counttau(AICR, lgrd, bsnb))

  inopt <- rbind(countint(bsnb, GACVR.IN),
                 countint(bsnb, BICR.IN))

  tauopt <- rbind(counttau(GACVR, lgrd, bsnb),
                  counttau(BICR, lgrd, bsnb))

  return(list(inopt = inopt, tauopt = tauopt, cenrate = cenrate))

}


