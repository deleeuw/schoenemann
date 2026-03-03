library(RSpectra)

schoenemann <- function(dmat, p = 2, fsel = 0) {
  cmat <- -doubleCenterRect(dmat) / 2
  svdd <- svd(cmat, nu = p, nv = p)
  lvdd <- diag(sqrt(svdd$d[1:p]))
  gmat <- svdd$u %*% lvdd
  hmat <- svdd$v %*% lvdd
  fmat <- dmat + 2 * tcrossprod(gmat, hmat)
  ftil <- columnCenter(fmat)
  ggmt <- NULL
  for (s in 1:p) {
    for (t in 1:s) {
      fa <- ifelse(s == t, 1, 2)
      ggmt <- cbind(ggmt, fa * gmat[, s] * gmat[, t])
    }
  }
  kmat <- cbind(columnCenter(ggmt), 2 * columnCenter(gmat))
  if (fsel == 0) {
    fs <- apply(ftil, 1, mean)
  } else {
    fs <- ftil[, fsel]
  }
  alst <- lm.fit(kmat, fs)
  print(alst$residuals)
  if (alst$rank < ncol(kmat)) {
    print("warning: matrix K is singular, powering on nevertheless")
  }
  avec <- alst$coefficients
  mmat <- matrix(0, p, p)
  k <- 1
  for (s in 1:p) {
    for (t in 1:s) {
      mmat[s, t] <- mmat[t, s] <- avec[k]
      k <- k + 1
    }
  }
  emat <- eigen(mmat)
  eval <- emat$values
  evec <- emat$vectors
  if (min(eval) < 1e-15) {
    print("warning: M matrix not psd, nevertheless powering on")
  }
  tmat <- evec %*% diag(sqrt(abs(eval)))
  smat <- solve(t(tmat))
  munu <- solve(tmat, rev(rev(avec)[1:p]))
  xhat <- gmat %*% tmat
  xhat <- xhat + outer(rep(1, dim(xhat)[1]), munu)
  yhat <- hmat %*% smat
  dhat <- outer(diag(tcrossprod(xhat)), diag(tcrossprod(yhat)), "+")
  dhat <- dhat - 2 * tcrossprod(xhat, yhat)
  return(list(
    d = dmat,
    dhat = dhat,
    xhat = xhat,
    yhat = yhat,
    f = ftil,
    k = kmat,
    m = mmat
  ))
}

doubleCenterRect <- function(x) {
  r <- apply(x, 1, mean)
  s <- apply(x, 2, mean)
  m <- mean(x)
  return(x - outer(r, s, "+") + m)
}

columnCenter <- function(x) {
  apply(x, 2, function(z)
    z - mean(z))
}

plotUnfold <- function(h, rowLabs = as.character(1:nrow(h$xhat)), colLabs = as.character(1:nrow(h$yhat))) {
  plot(rbind(h$xhat, h$yhat), type = "n")
  text(h$xhat,
       rowLabs,
       col = "RED",
       cex = 1)
  text(h$yhat,
       colLabs,
       col = "BLUE",
       cex = 1.5)
}
