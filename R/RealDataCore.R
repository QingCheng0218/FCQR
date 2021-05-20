# integrate function for int_a^b Br(s) Z(s) ds, where Z denoted as msbp
intzfun0 <- function(msdata, a0, b0, gn, time, norder, cba){
  zid <- 3:98;
  msbp0 <- msdata[, zid];

  Z <- apply(msbp0, 1, min);
  Z.M <- matrix(rep(Z, dim(msbp0)[2]), nrow=dim(msbp0)[1]);
  msbp <- msbp0 - Z.M

  if(cba==1){
    basistime <- create.bspline.basis(c(0, 1), breaks=c(min(time),  max(time)), norder = norder);

    nbasis <- basistime$nbasis
    valueBasisT <-eval.basis(time, basistime);
  }else if(cba==2){
    basistime <- create.bspline.basis(c(0, 1),nbasis = cba + 3,  norder = norder);#breaks=c(min(time), 0.5, max(time)),

    nbasis <- basistime$nbasis
    valueBasisT <-eval.basis(time, basistime);
  }else if(cba==3){
    basistime <- create.bspline.basis(c(0, 1), nbasis = cba + 3,  norder = norder); #breaks=c(min(time), 0.25, 0.75, max(time)),

    nbasis <- basistime$nbasis
    valueBasisT <-eval.basis(time, basistime);

  }else if(cba==4){
    basistime <- create.bspline.basis(c(0, 1),nbasis = cba + 3,   norder = norder);#breaks=c(min(time), 0.25, 0.5, 0.75, max(time)),

    nbasis <- basistime$nbasis
    valueBasisT <-eval.basis(time, basistime);
  }else if(cba==5){
    basistime <- create.bspline.basis(c(0, 1), nbasis = cba + 3,  norder = norder);#breaks=c(min(time), 0.25, 0.4, 0.6, 0.75, max(time)),

    nbasis <- basistime$nbasis
    valueBasisT <-eval.basis(time, basistime);
  }
  inte <- gauss.quad(n = gn)
  a1 <- (b0 - a0)/2;
  a2 <- (a0 + b0)/2;
  w <- inte$weights;
  nodes <- inte$nodes;
  neww <- a1*w
  newnodes <- nodes*a1+a2;

  basis <- basistime #create.bspline.basis(range(newnodes), nbasis = nbasis, norder = norder);
  valueBasis <-eval.basis(newnodes, basis);

  ycom <- t(apply(msbp, 1, function(x){
    sp <- smooth.spline(time, x)
    y <- predict(sp, newnodes)$y;
    return(y)
  }))
  # -------------------------------------
  New.W <- matrix(rep(neww, dim(msdata)[1]), nrow = dim(msdata)[1], byrow = TRUE)
  Yint <- New.W*ycom;
  covdata <- Yint%*%valueBasis*a1


  return(list(msdata = msdata,  covdata = covdata, Z = Z, nbasis = nbasis, valueBasisT =valueBasisT))
}

# AFT model to obtain the optimal weight
weifun <- function(msdata, covdata){
  surtime <- msdata[, 1:2]
  su_wei <- Surv(surtime[, 1], surtime[, 2]);
  fit.fun <- survreg(su_wei ~ covdata-1, dist = 'exponential')
  #standardized residuals
  resi <- residuals(fit.fun)/fit.fun$scale
  su_resi <- Surv(resi[surtime[, 2]==1]^2, surtime[surtime[, 2]==1, 2]);
  fit.res<- survreg(su_resi ~ covdata[surtime[, 2]==1, ] - 1, dist = 'exponential')
  sigres <- abs(covdata%*%coef(fit.res))
  weight <- 1/sigres
  return(weight)
}

coeffun0 <- function(valueBasisT, covdata,  msdata,
                     grd, qv, weight0, Z, time){

  surtime <- msdata[, 1:2]
  nbasis <- dim(valueBasisT)[2]
  covdata <- cbind(Z, covdata);
  std0 <- apply(covdata, 2,sd)
  std<- matrix(std0, nrow(covdata), ncol(covdata), byrow = T)
  covdata <- covdata/std
  T <- log(surtime[,1]/365)/sd(log(surtime[,1]/365))

  fit <- crq(Surv(T, surtime[,2])~.-1, data = data.frame(covdata),  taus=grd,method= "PengHuang")
  # estimate the coefficients
  coefest <- t(coef(fit,grd))
  rescoef <- t(coefest[, 1:(nbasis + 1)])

  coef2 <- rescoef[2:dim(rescoef)[1], ]
  res <- valueBasisT%*%coef2;
  beta0 <- rescoef[1, ];
  #-----------------------------------------------#
  # predict value
  pz1 <- (matrix(rep( covdata[, 1], lgrd), ncol = lgrd))*(matrix(rep(beta0, length(Z)), ncol = lgrd, byrow = TRUE))
  pzs <- covdata[, -1]%*%coef2
  pre.y <- exp(pz1 +  pzs)
  #-----------------------------------------------#

  #95% confidence intercal only report tau =0.3,0.4,0.5,0.6,0.7
  m <- length(time)
  fc <- summary(fit, taus = grd[qv], alpha = .05, se="boot", R = 200, covariance=TRUE)
  # for baseline covariate
  sigb0 <- base0 <- rep(0, length(grd[qv]));
  # for covariate function
  sigfun <- Bfun <- matrix(0, ncol = m, nrow = length(grd[qv]))

  for(tt in 1:length(grd[qv])){
    sig0 <- fc[[tt]]$cov[1, 1];
    sigf0 <- fc[[tt]]$cov[-1, -1];
    SD0 <- t(t(std0[-1]))%*%std0[-1]

    sig1 <- sigf0/SD0;

    sigfun[tt, ] <- diag(valueBasisT %*% sig1 %*% t(valueBasisT))
    sigb0[tt] <- sig0/(std0[1]^2);

    cf <- fc[[tt]]$coefficients[, 1]
    Bfun[tt, ] <- valueBasisT %*% (cf[-1]/std0[-1]);
    base0[tt] <- cf[1]/std0[1]
  }


  return(list(rescoef = rescoef, pre.y = pre.y, covdata = covdata, T = T,
              base0 = base0, Bfun = Bfun, sigb0 = sigb0, sigfun = sigfun));
}

bandplot <- function(da00, da01, da02, tau.id, t, min.y, max.y){

  min.x <- min(t);
  max.x <- max(t);

  da01[which(da01>max.y)] <- max.y;
  da02[which(da02<min.y)] <- min.y;
  # x.break <- seq(20, 42, 2)
  # time.label=c( "20:00","22:00", "24:00", "02:00", "04:00","06:00",
  #               "08:00", "10:00", "12:00", "14:00", "16:00", "18:00");

  x.break <- seq(20, 42, 4)
  time.label=c( "20:00", "24:00", "04:00",
                "8:00", "12:00","16:00");


  data <- cbind(t, da00, da01, da02);
  data <- data.frame(data);
  tauvalue <- grd[qv[tau.id]];


  p1 <- ggplot(data, aes(t, da00)) + geom_line(data=data)+
    geom_ribbon(data=data,aes(ymin=da01,ymax=da02),alpha=0.2)

  p1 <- p1 + scale_y_continuous(limits=c(min.y, max.y),
                                breaks=seq(min.y, max.y, 0.2))+
    scale_x_continuous(limits=c(min.x, max.x),
                       breaks= x.break, labels=time.label)

  p1 <- p1 + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face=0.8),
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          axis.text.x = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          #panel.background = element_blank())+ theme(legend.position = c(0.15,0.75))
          panel.background = element_blank())

  # p1 <- p1 + geom_vline(xintercept = 30, linetype="dotted");
  p1 <- p1 + geom_hline(yintercept = 0, linetype="dotted");

  p1 <- p1 + labs(title=bquote(tau == .(tauvalue)))
  p1 <- p1 + theme(legend.position="bottom");

  p1 <- p1 + labs(x = "Time, hours", y =expression(hat(alpha)(s, tau)))
  return(p1);
}

impuplot <- function(y, min.y, max.y){
  i <- 1
  repeat{
    id.out <- which(y< min.y | y>max.y);
    missnu <- length(id.out)
    y[id.out] <- NA;
    y <- na.interpolation(y, option ="spline")
    i <- i+1
    if(sum(y< min.y | y>max.y)==0 | i>5) break
  }
  if(sum(y<min.y)!=0) y[y<min.y] <- min.y
  if(sum(y>max.y)!=0) y[y>max.y] <- max.y
  return(list(missnu = missnu, y = y))
}

Hfun <- function(u){
  hu <- -log(1-u)
  return(hu)
}

alfun <- function(TM, pre.y, grd, n){
  tv <- c(0, grd)
  HH <- diff(Hfun(tv))
  pre.ynew <- cbind(rep(0, n), pre.y[, -dim(pre.y)[2]])
  ind.m <- 1*(exp(TM) >=pre.ynew);
  AA <- ind.m*(matrix(rep(HH, dim(TM)[1]), nrow = dim(TM)[1], byrow = TRUE))
  ATau <- t(apply(AA, 1, cumsum));
  # ATau is n by lgrd matrix
  return(ATau)
}


bicfun <- function(n, ns, cn, TM, pre.y, DEL, grd){
  ATau <- alfun(TM, pre.y, grd, n);
  rtb <- TM - log(pre.y);
  Rbtau <- apply(rtb*(ATau - 1*(rtb<0)*DEL), 2, sum)/(n - cn*ns)

  gacv <- Rbtau
  aic <- log(apply(rtb*(ATau - 1*(rtb<0)*DEL), 2, mean)) + ns/n
  bic <- log(apply(rtb*(ATau - 1*(rtb<0)*DEL), 2, mean)) + (ns*log(n))/n

  return(list(gacv = gacv, aic = aic, bic = bic))
}


predfun <- function(msdata, a0, b0, gn, time, norder, cba){
  # predict value: exp(b^t x)
  del <- msdata[, 2]
  DEL <- matrix(rep(del, lgrd), ncol = lgrd);

  covres <- intzfun0(msdata, a0, b0, gn, time, norder, cba)
  nbasis <- covres$nbasis
  covdata <- covres$covdata;
  Z <- covres$Z;
  valueBasisT <- covres$valueBasisT
  # weight0 <- weifun(msdata, covdata);
  weight0 <- rep(1, nrow(msdata))

  coefres <- coeffun0(valueBasisT, covdata,  msdata, grd, qv, weight0, Z, time)

  T <- coefres$T
  TM <- matrix(rep(T, lgrd), ncol = lgrd);

  pre.y <- coefres$pre.y

  return(list(pre.y=pre.y, TM = TM, DEL = DEL, nbasis = nbasis));
}

# integration result
optfun <- function(grd, BICR){
  de.tau <- diff(grd)[1]
  bic.in <- apply(BICR*de.tau, 2, sum);
  bic.in;
  rank(bic.in)
  # interes <- paste(round(bic.in,5),"(",rank(bic.in),")", seq="")
  interes <- round(bic.in,5);
  # BIC for each tau
  aa <- BICR
  # aa <- BICR[,2:5]
  bb <- apply(aa, 1, function(x){return(which.min(x))})
  tauopt <- rep(0, length(CBA))
  for(i in 1:length(CBA)){
    tauopt[i] <- sum(bb==CBA[i])
  }


  res <- rbind(NBS, interes, tauopt)
  return(res)
}
