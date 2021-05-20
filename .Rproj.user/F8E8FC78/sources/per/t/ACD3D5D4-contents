#' Title
#'
#' @param t
#' @param v
#' @param nk
#'
#' @return
#' @export
#'
#' @examples
simBX <- function(t, v, nk){
  kk <- pi*(seq(1, nk, 1) -1);
  phik <- sqrt(2)*cos(t%*%t(kk))
  phik[, 1] <- rep(1, length(t));
  uk <- runif(nk, -sqrt(3), sqrt(3));
  gk <- ((-1)^(seq(2, (nk+1), 1)))*(seq(1, nk, 1)^(-v/2));
  XK <- gk*uk*t(phik);
  X <- abs(colSums(XK));

  bk <- 4*((-1)^(seq(1, nk, 1)))*(seq(1, nk, 1)^(-2));
  B0 <- colSums(bk*t(phik));
  alp <- X*B0
  return(list(X = X, alp = alp, B0 = B0));
}

simfun <- function(t, v, n, sdx, sdz, b0,  lam0, c0){
  # set.seed(itr);

  Z <- alpha <- matrix(NA, n, m);
  Timedata <- matrix(NA, n, 2);

  Z1 <- rnorm(n, 0, sdz);
  Z2 <- runif(n, 0, c0);
  xi <- rnorm(n, 0, sdx)

  for(i in 1:n){
    # set.seed(itr * 1000 + i);
    resbx <- simBX(t, v, nk);
    Z[i, ] <- resbx$X;
    alpha[i, ] <- resbx$alp
  }
  iz <- apply(Z, 1, sum)*deltat;
  iza <- apply(alpha, 1, sum)*deltat;

  mT0 <- b0[1]*Z1 + ( b0[2]*Z2  + iz)*xi + iza

  if(lam0==0){C <- exp(mT0)}else{C <- 0.05*rexp(n, lam0);}

  Timedata <- cbind(exp(mT0), C);
  mZ <- Z;

  Z0 <- cbind(Z1,Z2);

  return(list(Timedata = Timedata, mZ = mZ, Zb = Z0))

}

#-------------------------------------------------#
# eleq function find the position of specified elements.(which function doesn't  work)
eleq <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

Hfun <- function(u){
  hu <- -log(1-u)
  return(hu)
}

# alfun <- function(n, TM, pre.y){
#   tv <- c(0, grd)
#   HH <- diff(Hfun(tv))
#   pre.ynew <- cbind(rep(0, n), pre.y[, -dim(pre.y)[2]])
#   ind.m <- 1*(exp(TM) >=pre.ynew);
#   AA <- ind.m*(matrix(rep(HH, dim(TM)[1]), nrow = dim(TM)[1], byrow = TRUE))
#   ATau <- t(apply(AA, 1, cumsum));
#   # ATau is n by lgrd matrix
#   return(ATau)
# }
#
# bicfun <- function(n, ns, cn, TM, pre.y, DEL){
#   ATau <- alfun(n, TM, pre.y);
#   rtb <- TM - log(pre.y);
#   Rbtau <- apply(rtb*(ATau - 1*(rtb<0)*DEL), 2, sum)/(n - ns)
#
#   gacv <- Rbtau
#   temp <- apply(rtb*(ATau - 1*(rtb<0)*DEL), 2, mean)
#   if(length(which(temp<0))!=0){
#     na.id <- which(temp < 0);
#     temp[na.id] <- NA
#     aic <- log(temp) + ns/n;
#     bic <- log(temp) + (ns*log(n))/n;
#   }else{
#     aic <- log(temp) + ns/n
#     bic <- log(temp) + (ns*log(n))/n
#   }
#
#   return(list(gacv = gacv, aic = aic, bic = bic))
# }

countint <- function(bsnb, BICR.IN){
  idna <- apply(BICR.IN, 1, function(x){
    sum(is.na(x))
  })
  if(sum(idna==length(bsnb))!=0){
    AA <- BICR.IN[-which(idna==length(bsnb)), ]
  }else{
    AA <- BICR.IN
  }

  aa <- apply(AA, 1, function(x){
    x[which(is.na(x))] <- max(x, na.rm = TRUE)
    return(which.min(x))})

  count <- rep(0, length(bsnb))
  for(i in 1:length(bsnb)){
    count[i] <- sum(aa==i)
  }
  return(count)
}
counttau <- function(BICR, lgrd, bsnb){
  count <- rep(0, length(bsnb))
  bb <- rep(0, lgrd)
  for(i in 1:lgrd){
    aa <- BICR[,,i]
    bb[i] <- bsnb[which.max(countint(bsnb, aa))]
  }
  for(i in 1:length(bsnb)){
    count[i] <- sum(bb==bsnb[i])
  }
  return(count)
}

Fcqr <- function(time, Zs, Zb, t, grd, qv,
                 valueBasis, sdid = 0, int.id = 0){
  # ----------------------------------------------------------- #
  # lb0: the number of non-zero of beta1
  lb0 = ncol(Zb);
  # deltat: the increment of t.
  deltat <- (diff(t)[1]);
  # ----------------------------------------------------------- #
  # event time and censor indicator denoted as dat.yc
  surtime <- cbind(pmin(time[, 1], time[, 2]), (time[, 1] <= time[, 2]))

  # estimate the function using Bspline and integrate
  #int_0 ^T Z(s)Br(s) ds denote ad covMat
  covMat <- Zs%*%valueBasis*deltat

  # ix <- which(apply(abs(covMat), 2, mean)> 1e-5)
  # dat <- cbind(dsurtime,  covMat[,ix])
  # ----------------------------------------------------------- #
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
  #   dat.cov <- cbind(Z1, dat[, 3:ncol(dat) ])}

  if(int.id==0){
    fit <- crq(Surv(T, surtime[,2])~.-1, data = data.frame(covdata),  taus=grd, method= "PengHuang")
  }else{
    fit <- crq(Surv(T, surtime[,2])~., data = data.frame(covdata),  taus=grd, method= "PengHuang")
  }

  # ------------------------------------------------------------
  # estimate the coefficients
  coefest <- t(coef(fit,grd))
  rescoef <- t(coefest[, 1:(nbasis + lb0)])

  coef2 <- rescoef[(lb0+1):dim(rescoef)[1], ]
  res <- valueBasis%*%coef2;
  beta0 <- rescoef[1:lb0, ];
  #-----------------------------------------------#
  # predict value
  if(lb0==1){
    pz1 <- t(outer(beta0, covdata[,1]))
  }else{
    res1 <- t(outer(beta0[1,], covdata[,1]));
    res2 <- t(outer(beta0[2,], covdata[,2]));
    pz1 <- res1 + res2
  }
  pzs <- covdata[, -(1:lb0)]%*%coef2
  pre.y <- exp(pz1 +  pzs)
  #-----------------------------------------------#
  #95% confidence intercal only report tau =0.3,0.4,0.5,0.6,0.7
  m <- length(t)

  # for baseline covariate
  if(lb0==1){
    sigb0 <- base0 <- rep(NA, length(grd[qv]));
  }else{
    sigb0 <- base0 <- matrix(NA, ncol = lb0, nrow = length(grd[qv]))
  }

  # for covariate function
  sigfun <- Bfun <- matrix(NA, ncol = m, nrow = length(grd[qv]))

  if(sum(is.na(coefest)) < (length(coefest)/2)){
    fc <- summary(fit, taus = grd[qv], alpha = .05, se="boot", R = 200, covariance=TRUE)

    for(tt in 1:length(fc)){

      covbfun <- fc[[tt]]$cov[-(1:lb0), -(1:lb0)]
      covb0 <- fc[[tt]]$cov[1:lb0, 1:lb0]
      sigfun[tt, ] <- diag(valueBasis %*% covbfun %*% t(valueBasis))

      cf <- fc[[tt]]$coefficients[, 1]
      Bfun[tt, ] <- valueBasis %*% cf[-(1:lb0)];
      if(lb0==1){
        sigb0[tt] <- covb0;
        base0[tt] <- cf[1:lb0];
      }else{
        sigb0[tt, ] <- diag(covb0);
        base0[tt, ] <- cf[1:lb0];
      }
    }
  }
  # ------------------------------------------------------------
  return(list(pre.y = pre.y, covdata = covdata, T = T,
              base0 = base0, Bfun = Bfun, sigb0 = sigb0, sigfun = sigfun));
}


#------simulation function to get the coefficient-------#
simulation.fun <- function(n, ns, v, b0, lam0, c0, sdz, sdx, nsim, BN,
                           lgrd, t, nbasis, norder, qv, int.id, sdid){

  basis <- create.bspline.basis(range(t),nbasis = nbasis, norder = norder)
  nbasis <- basis$nbasis
  valueBasis <-eval.basis(t, basis)

  gacvr <- bicr <- aicr <- matrix(0, nrow = nsim, ncol = lgrd);
  cdata <- matrix(0, ncol = nsim, nrow = n);

  m <- length(t)
  Base0.m <- Sigb0.m <- array(0, dim = c(nsim, length(qv), lb0));
  Bfun.m <- Sigfun.m <- array(0, dim = c(nsim, m, length(qv)))

  itr <- 0;
  tt <- 1;
  repeat{
    itr <- itr +1

    resMat <-  matrix(NA, n, m);
    covMat <- matrix(NA, n, nbasis);
    temp <- simfun(t, v, n, sdx, sdz, b0, lam0, c0);

    Z1 <- temp$Z1 # baseline covariate.
    time <- temp$Timedata # survival time and censoring time.
    Zs <- temp$mZ # functional covariate.

    cdata[, itr] <- 1*(time[, 1] <= time[, 2])
    res.itr <- try(Fcqr(time, Zs, Z1, lb0, grd, qv, ind.id, valueBasis, sdid));
    if(class(res.itr) != 'try-error'){

      Base0.m[itr, , 1:lb0] <- res.itr$base0;
      Sigb0.m[itr, ,1:lb0] <- res.itr$sigb0;
      Bfun.m[itr, , ] <- t(res.itr$Bfun);
      Sigfun.m[itr, , ] <- t(res.itr$sigfun);

      # simcoef1[itr, 1:lb0, ] <- res.itr$res1;
      # simcoeffun[itr, , ] <- res.itr$res2

      T <- res.itr$T
      TM <- matrix(rep(T, lgrd), ncol = lgrd);
      pre.y <- res.itr$pre.y

      del <- 1*(time[, 1] <= time[, 2]);
      DEL <- matrix(rep(del, lgrd), ncol = lgrd);

      ns <- nbasis + lb0
      res <- bicfun(n, ns, cn, TM, pre.y, DEL)

      gacvr[itr, ] <- res$gacv
      bicr[itr, ] <- res$bic
      aicr[itr, ] <- res$aic
      tt <- tt + 1;
    }
    if(tt%%100==0) {print(tt);print(itr);print(Sys.time())}
    if(tt == nsim + 1){
      break
    }
  }
  return(list(bicr=bicr, cdata = cdata, Base0.m=Base0.m, Sigb0.m = Sigb0.m,
              Bfun.m=Bfun.m, Sigfun.m = Sigfun.m));

}
#-------------------------------------------------------#

alfun.be <- function(alfun, Sigfun.m, Bfun.m, grd, qv, sdx){
  # alfun <- simBX(t, v, nk)$B0

  plot.m <- array(0, dim=c(dim(Sigfun.m)[2], 6, dim(Sigfun.m)[3]));

  for(i in 1:dim(Sigfun.m)[3]){

    # da00:true curve, da02: mean of simulation curve
    # da01 & da02 simulation 95% ci. da04&da05 bootstrap 95%ci

    sigma <- apply(Sigfun.m[,,i], 2, mean, na.rm = TRUE)

    da01 <- apply(Bfun.m[,,i], 2, quantile,0.025, na.rm = TRUE);
    da02 <- apply(Bfun.m[,,i], 2, mean, na.rm = TRUE)
    da03 <- apply(Bfun.m[,,i], 2, quantile,0.975, na.rm = TRUE);

    da04 <- da02 + 1.96 * sqrt(sigma);
    da05 <- da02 - 1.96 * sqrt(sigma);

    da00 <-   (alfun)   +  qnorm(grd[qv[i]], 0, sdx)
    plot.m[, , i] <- cbind(da00, da01, da02, da03, da04, da05)
  }
  return(plot.m)
}

b0res.be <- function(Base0.m, Sigb0.m, btautrue, qv){
  # bias:sampling biases
  # sse:sampling stardard error estimator
  # see:mean of standard error estimator
  # cp: 95% cp
  b0true <- btautrue[qv];
  sd.m <- sqrt(Sigb0.m);

  b0hat <- apply(Base0.m, 2, mean, na.rm=TRUE)
  bias <- abs(b0true - b0hat)
  sse <- apply(Base0.m, 2, sd, na.rm=TRUE)
  see <- apply(sd.m, 2, mean, na.rm=TRUE)

  cp <- rep(0, length(qv))
  for(i in 1:length(qv)){
    temp <- Base0.m[, i];
    upper <- temp + 1.96*sd.m[, i]
    lower <- temp - 1.96*sd.m[, i]
    cp[i] <- mean((lower<= b0true[i])*(b0true[i]<=upper),na.rm = TRUE)
  }

  res <- cbind(bias, sse, see, cp)
  return(res)
}

alphaplot.be <- function(t, grd, qv,  tau.id, plot.m, min.y, max.y){

  tauvalue <- grd[qv[tau.id]]
  data <- plot.m[,,tau.id];


  # min.y <- floor(min(data)*10)/10;
  # max.y <- ceiling(max(data)*10)/10;
  # Round to nearest "nid"
  # nid <- 1;
  # min.y <- round(floor(min(data)/nid))*nid;
  # max.y <- round(ceiling(max(data)/nid))*nid;
  #
  # min.y <- -10;
  # max.y <- 1;
  #
  min.x <- min(t);
  max.x <- max(t);

  data <- data.frame(data)

  names(data)[1:6] <- c('y0', 'y1', 'y2', 'y3', 'y4', 'y5');
  da00 <- data.frame( x = t, y = data$y0, type = paste('line1') )
  da01 <- data.frame( x = t, y = data$y1, type = paste('line2') )
  da02 <- data.frame( x = t, y = data$y2, type = paste('line3') )
  da03 <- data.frame( x = t, y = data$y3, type = paste('line4') )
  da04 <- data.frame( x = t, y = data$y4, type = paste('line5') )
  da05 <- data.frame( x = t, y = data$y5, type = paste('line6') )


  name.label <- c("true", "sim025", "simnk", "sim975", "boot025",  "boot975")

  g <- ggplot()

  g <- g + geom_line ( data = da00, aes( x, y, color = type), cex=1.5, linetype = 1 )
  g <- g + geom_line ( data = da01, aes( x, y, color = type), cex=1.5, linetype = 3 )
  g <- g + geom_line ( data = da02, aes( x, y, color = type), cex=1.5, linetype = 5 )
  g <- g + geom_line ( data = da03, aes( x, y, color = type), cex=1.5, linetype = 3 )
  g <- g + geom_line ( data = da04, aes( x, y, color = type), cex=1.5, linetype = 4 )
  g <- g + geom_line ( data = da05, aes( x, y, color = type), cex=1.5, linetype = 4 )

  # g <- g + geom_point( data = da00, aes( x, y, color = type ), shape = 0 )
  # g <- g + geom_point( data = da01, aes( x, y, color = type ), shape = 1 )
  # g <- g + geom_point( data = da02, aes( x, y, color = type ),  shape = 0 )
  # g <- g + geom_point( data = da03, aes( x, y, color = type ),  shape = 0 )
  # g <- g + geom_point( data = da04, aes( x, y, color = type ), shape = 0 )
  # g <- g + geom_point( data = da05, aes( x, y, color = type ),  shape = 0 )

  g <- g + guides( colour = guide_legend( override.aes =
                                            list(linetype = c(1, 3, 5, 3, 4, 4),
                                                 shape =  c(1, 3, 5, 3, 4, 4) )))

  # g <- g + labs(title=bquote(tau == .(tauvalue)))

  g <- g + scale_color_manual( name = NULL, labels = name.label,
                               values = c( "red", "seagreen4", "seagreen4", "seagreen4", "coral4", "coral4"))


  g <- g + 		scale_y_continuous(limits=c(min.y, max.y),breaks=seq(min.y, max.y, 2))+
    scale_x_continuous(limits=c(min.x, max.x),breaks= seq(0, 1,0.2))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face=0.8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
          #axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          #panel.background = element_blank())+ theme(legend.position = c(0.15,0.75))
          panel.background = element_blank())

  # remove the legend of the chart
  g <- g + labs(title=bquote(tau == .(tauvalue)))
  g <- g + theme(legend.position="none")

  # print(g)
  return(g)
}
