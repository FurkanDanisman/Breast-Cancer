# Last Updated: 6 November 2006
# Last updated: 17 march 2010 (return: vt)
# Transformed to R on March 9 2024

# Required packages/libraries
library(multipol)
library(pracma)
library(matconv)
library(matlab)
library(Matrix)

# EH Analysis 
EH_analysis_packed <- function(yy, zadd, status, xx, initial, method) {
  
  zadd_original <- zadd
  
  # ANALYSIS
  nn <- nrow(xx)  # nn: no. of subjects
  mxx <- ncol(xx)  # mxx: no. of covariates
  
  nsubtype <- nrow(zadd)
  nprm2 <- ncol(zadd)
  
  zadd <- cbind(zadd[, 1], zadd[, 2], zadd[, 3], zadd[, 4])
  
  # Defining the grand design matrix
  # together (if different design matrices are chosen for different exposure/covariate
  # this part needs to be changed)
  zz <- zadd
  
  for (i in 1:(mxx-1)) {
    if (mxx==1) {
      break
    }
    zz <- blkdiag(zz, zadd) 
  }
  
  # Permuting the rows of the design matrix
  nz <- nrow(zadd)
  mz <- ncol(zadd)
  perm <- as.matrix(seq(1, 1+(mxx - 1)*nz, by = nz)) 
  for (i in 2:nz) {
    perm <- rbind(perm, t(t(seq(i, (i + (mxx - 1) * nz), by = nz))))  
  }
  zz <- zz[perm, ]
  
  # Fitting the two-stage model
  # Defining the response vector yyc, takes values in 0:nsubtype indicating whether
  # the subject is control (0) or the numeric level for the disease subtype
  yyc <- zeros(nn,1)
  for (j in 1:nz) {
    yyc[yy[, j] == 1] <- j
  }
  
  # Defining initial value
  if (initial == 0) {
    theta0 <- zeros(length(zz[1, ]),1)
  }else if (initial == 1) {
    lrp_result <- logistic_regression_poly(yy, xx) 
    theta00 <- lrp_result$beta  
    covmat <- lrp_result$covmat 
    theta00 <- theta00[2:(mxx+1)] 
    theta01 <- zeros(nprm2-1,mxx)
    theta01 <- rbind(t(theta00),theta01) 
    for (i in 1:mxx) {  
      if (i==1) {
        theta0 <- theta01[,i]
      }else{theta0 <- rbind(theta0,theta01[,i])}
    }
  }

  # Calling the program for fitting the data using the method for multiple disease characters (output estimates are grouped by covariates)
  if (method == 1) { 
    result_cmvpolysp <- cmvpolysp(theta0, yy, yyc, xx, zz)
    theta <- result_cmvpolysp$delta
    vt <- result_cmvpolysp$v
    found <- result_cmvpolysp$found
    theta_estimates_with_modified_zadd <- theta
  } else if (method == 0) {
    result_mvpoly <- mvpoly(theta0, yy, xx, zz)
    theta <- result_mvpoly$delta
    vt <- result_mvpoly$v
  }
  
  # Second part
  zadd <- zadd_original
  zz <- zadd
  
  d <- zeros(4,1) 
  d[1] <- theta_estimates_with_modified_zadd[1, 1]
  d[2] <- theta_estimates_with_modified_zadd[2, 1]
  d[3] <- theta_estimates_with_modified_zadd[3, 1]
  d[4] <- theta_estimates_with_modified_zadd[4, 1]
  
  result_cpolyscoreTest <- cpolyscoreTest(d, yy, yyc, xx, zz)
  s <- result_cpolyscoreTest$sc  
  i <- result_cpolyscoreTest$info   
  score_test_value <- t(s) %*% solve(-i) %*% s
  score_test_all <- 0
  g <- 0
  for (j in 2:nrow(i)) {
    s <- s[j]
    seq <- 1:nrow(i)
    seq <- seq[-j]
    cv <- i[j,seq,drop=FALSE]
    ii <- i[-j,-j]
    i_1 <- -solve(i[j,j] - cv%*%solve(ii) %*% t(cv))
    score_test_all[j-1] <- t(s) %*% solve(i_1) %*% s
    g[j-1] <- paste("Theta_1",j,sep = "")
    s <- result_cpolyscoreTest$sc 
  }
  
  score_test_all <- as.matrix(score_test_all)
  score_df <- as.data.frame(rbind(score_test_all,score_test_value))
  colnames(score_df) <- "Score Test Values"
  rownames(score_df) <- c(g,"Theta0") 
  
  
  # OUTPUT
  # Second-stage parameters, st. dev., pvalues and confidence intervals:
  est <- theta 
  sd <- sqrt(diag(vt))  
  pvalue <- 2 * (1 - pnorm(abs(est / sd)))
  orest <- exp(est)
  orl95 <- exp(est - 1.9599 * sd)
  oru95 <- exp(est + 1.9599 * sd)
  logorl95 <- est - 1.9599 * sd
  logoru95 <- est + 1.9599 * sd
  
  # below (i,j) th. entry is for the (second stage)parameter of jth. contrast
  # for the ith. covariate
  
  est_matrix <- matrix(est, nrow = mz, ncol = mxx,  byrow = TRUE)  
  sd_matrix <- matrix(sd, nrow = mz, ncol = mxx,  byrow = TRUE) 
  pvalue_matrix <- matrix(pvalue, nrow = mz, ncol = mxx,  byrow = TRUE) 
  orest_matrix <- matrix(orest, nrow = mz, ncol = mxx,  byrow = TRUE) 
  orl95_matrix <- matrix(orl95, nrow = mz, ncol = mxx,  byrow = TRUE) 
  oru95_matrix <- matrix(oru95, nrow = mz, ncol = mxx,  byrow = TRUE)
  df_output <- data.frame("log(OR) estimates"=est_matrix,"OR estimates"=orest_matrix,"p-values"=pvalue_matrix,
              "Standard Deviation"=sd_matrix,"Lower Bound of 95% CI" = orl95_matrix,
              "Upper Bound of 95% CI" = oru95_matrix)
  cat("TABLE: \n")
  print(df_output)
  
  return(score_df)
}

# Functions

# 1. CMVPOLYSP
cmvpolysp <- function(delta0, y, yc, x, z) {
  iter <- 0
  rerror <- 1.0
  delta <- delta0
  maxiter <- 50
  tolx <- 1e-07
  found <- 0
  
  while ((iter < maxiter) & all(rerror > 1e-05)) {
    result_cpolyscoreTest <- cpolyscoreTest(delta0, y, yc, x, z)
    l <- result_cpolyscoreTest$score
    h <- result_cpolyscoreTest$info
    
    p <- -solve(h) %*% l
    delta <- flnsrch(cpolyscoreTest, delta0, 0.5 * sum(l * l), t(h) %*% l, p, tolx, 1.0, 0.5, y, yc, x, z)
    delta <- delta$gamma
    rerror <- pmax(abs(delta - delta0) / pmax(abs(delta0), 0.1)) 
    delta0 <- delta
    iter <- iter + 1
  }
  
  if (all(rerror < 1e-05)) {
    found <- 1
  }
  # [delta,l,found]=fsolve('cpolyscore',delta0,OPTIMSET('Jacobian','on'),y,yc,x,z);
  
  result_cpolyscoreTest2 <- cpolyscoreTest(delta, y, yc, x, z)
  sc <- result_cpolyscoreTest2$score
  info <- result_cpolyscoreTest2$info
  # y = sparse(y)
  ycase <- matrix(rowSums(y) == 1,byrow = T) 
  x0 <- x[!ycase, , drop = FALSE]
  x1 <- x[ycase, , drop = FALSE]
  y1 <- y[ycase, , drop = FALSE]
  yc1 <- yc[ycase, , drop = FALSE]
  
  n0 <- nrow(x0)
  p <- ncol(x0) 
  n1 <- nrow(x1)
  p <- ncol(x1)
  k <- ncol(y1)
  n1 <- nrow(y1)
  q <- length(delta0)
  
  ###
  
  nn1 <- as.matrix(colSums(y1))
  
  b <- z %*% delta
  beta <- matrix(b, nrow = p, ncol = k)
  xb0 <- x0 %*% beta
  exb0 <- exp(xb0)
  sexb0 <- as.matrix(colSums(exb0)) 
  
  vx0 <- Reshape(x0, n0 * p, 1) 
  rvx0 <- repmat(vx0,1,k) 
  rexb0 <- repmat(exb0,p,1) 
  rxexb0 <- rvx0 * rexb0 
  seyep <- diag(p) 
  onen <- kronecker(seyep, matrix(1,1,n0))
  sxexb0 <- t(onen %*% rxexb0) 
  
  xb1 <- x1 %*% beta
  exb1 <- exp(xb1)
  exb1y1 <- exb1 * y1
  cexb1 <- as.matrix(rowSums(exb1y1)) 
  rcexb1 <- repmat(cexb1,1,p) 
  xcexb1 <- x1 * rcexb1
  
  d <- rcexb1 + repmat(sexb0[yc1,,drop=FALSE],1,p) 
  n <- xcexb1 + sxexb0[yc1, ]
  
  s <- n / d
  sc <- reshape(t(t(y1) %*% (x1 - s)),k * p,1) 

  # Computing meat part of the sandwich variance formula
  hri <- Reshape(t(repmat(as.matrix(1:k),1,p)),k*p,1) 
  hx0 <- repmat(x0,1,k)
  hexb0 <- exb0[, hri,drop=FALSE] 
  rnn1 <- repmat(nn1[hri],n0,1) 
  hs0 <- repmat((sexb0[hri,,drop=FALSE]),n0,1) 
  hs1 <- reshape(repmat(reshape(t(sxexb0),1,k*p),n0,1),n0,k*p) 
  w <- rnn1 * (hexb0 / hs0) * (hx0 - hs1 / hs0)
  C <- t(z) %*% t(w) %*% w %*% z
  
  iinfo <- solve(-info)
  v <- iinfo %*% (-info + C) %*% iinfo
  return(list(delta=delta,v=v,info=info)) 
}


# 2. CPOLYSCORETEST
cpolyscoreTest <- function(delta, y, yc, x, z) {
  ycase <- matrix(rowSums(y) == 1, byrow = T)
  x0 <- x[!ycase, , drop = FALSE]
  x1 <- x[ycase, , drop = FALSE]
  y1 <- y[ycase, , drop = FALSE]
  yc1 <- yc[ycase, , drop = FALSE]
  
  n0 <- nrow(x0)
  p <- ncol(x0)
  n1 <- nrow(x1)
  k <- ncol(y1)
  nn1 <- as.matrix(rowSums(y1))
  
  b <- z %*% delta
  beta <- reshape(b, p,k)
  xb0 <- x0 %*% beta
  exb0 <- exp(xb0)
  sexb0 <- (as.matrix(colSums(exb0)))
  
  vx0 <- reshape(x0, n0 * p, 1)
  rvx0 <- repmat(vx0, 1 , k)
  rexb0 <- repmat(exb0,p,1)
  rxexb0 <- rvx0 * rexb0
  seyep <- diag(p)
  onen <- kronecker(seyep, matrix(rep(1, n0), ncol = n0))
  sxexb0 <- t(onen %*% rxexb0)
  
  xb1 <- x1 %*% beta
  exb1 <- exp(xb1)
  exb1y1 <- exb1 * y1
  cexb1 <- as.matrix(rowSums(exb1y1))
  rcexb1 <- repmat(cexb1,1,p)
  xcexb1 <- x1 * rcexb1
  
  d <- rcexb1 + repmat(sexb0[yc1, ,drop=FALSE],1,p)
  n <- xcexb1 + sxexb0[yc1, , drop = FALSE]
  
  s <- n / d
  sc <- reshape(t(t(y1) %*% (x1 - s)), k * p, 1)
  sc <- t(z) %*% sc
  
  # Computing observed information
  cri <-  repmat(as.matrix((1:p)),p,1) 
  hri <-  reshape((repmat(as.matrix((1:p)),1,p)),p*p,1) 
  m <- 0
  sxxexb0 <- zeros(k,p*p) 
  
  for (i in 1:p) {
    for (j in 1:p) {
      m <- m + 1
      sxxexb0[, m] <- t(repmat(as.matrix(x0[, i] * x0[, j]), 1, k) * exb0) %*% matrix(1,n0,1) 
    }
  }
  
  xx1 <- x1[, cri] * x1[, hri]
  rcexb1 <- repmat(cexb1,1,p*p) 
  xxcexb1 <- xx1 * rcexb1
  
  d1 <- rcexb1 + repmat(sexb0[yc1, ,drop=FALSE],1, p * p)
  n1 <- xxcexb1 + sxxexb0[yc1, , drop = FALSE]
  i1 <- n1 / d1
  
  m <- 0
  y1i2 <- zeros(k,p*p)
  for (i in 1:p) {
    for (j in 1:p) {
      m <- m + 1
      y1i2[, m] <- t(y1) %*% (s[, i] * s[, j])
    }
  }
  
  ikp2 <- y1i2 - t(y1) %*% i1
  ippk <- reshape(ikp2,p, p, k)
  ick <- matrix(list(), k,1)

  for (i in 1:k) {
    ick[[i]] <- t(ippk[,,i])
    if (i==1) {
      i_addition <- ick[[i]]
    }else{
      i_addition <- blkdiag(i_addition,ick[[i]])
    }
  }
  i <- i_addition
  info <- t(z) %*% i %*% z
  
  return(list(score = sc, info = info))
}

# 3. LOGISTIC REGRESSION POLY
logistic_regression_poly <- function(yy, zz) {
  # Design matrix
  design <- cbind(matrix(1,nrow(yy),1),zz)
  
  # Number of columns in design matrix
  mm <- ncol(design)
  # Initialize variables
  beta_old <- zeros(mm,1)
  jj <- 0
  
  while (jj < 999) {
    jj <- jj + 1
    
    # Calculate predicted values and related quantities
    gg <- design %*% beta_old
    pp <- logistic(gg)
    rr <- pp * (1 - pp)
    rr <- max(0.0001, min(0.9999, rr))
    hh <- design * c((rr %*% matrix(1,1,mm)))
    aa <- t(design) %*% hh
    bb <- t(design) %*% (1-pp)

    # Update beta
    beta_new <- beta_old +(solve(aa) %*% bb)
    
    # Calculate the difference between old and new beta
    epss <- sum(abs(beta_old - beta_new))
    
    beta_old <- beta_new
    
    if (epss < 0.00001) {
      jj <- 99999999
    }
    
  }
  
  return(list(beta = beta_new, covmat = solve(aa)))
}


# 4. FLNSRCH
flnsrch <- function(FUN, gamma0, f0, g0, p, tolx, lmax, lmin, y,yc,x,z) {
  stpmax <- 1.0 * length(gamma0) # Length may change based on what gamma0 is.
  found <- 0.0
  sp2 <- sqrt(sum(p^2))
  if (sp2 > stpmax) {
    p <- p * (stpmax / sp2) # scale if attempted step is too large
  }
  slope <- sum(g0 * p)
  if (slope >= 0.0) {
    stop('Roundoff problem in lnsrch')
  }
  rerror <- max(abs(p) / max(abs(gamma0), 0.1))
  alamin <- tolx / rerror
  alam <- 1.0 # try the full step first
  while (found == 0.0) {
    gamma <- gamma0 + alam * p
    result <- FUN(gamma,y,yc,x,z) # Currently it is assumed this form of function will be called.
    l <- result[[1]] # Might require adjustment. Currently it is assumed this form of function will be called.
    h <- result[[2]] # Might require adjustment. Currently it is assumed this form of function will be called.
    f <- 0.5 * sum(l * l)  
    if (alam < alamin) {
      gamma <- gamma0 + alam * p # convergence on x
      found <- 1
      return(list(gamma = gamma, found = found))
    } else if (f <= f0 + 0.0001 * alam * slope) { # sufficient function decrease
      found <- 1
      return(list(gamma = gamma, found = found))
    } else { # backtrack
      if (alam == 1.0) { # first time
        scale <- -slope / (2.0 * (f - f0 - slope))
      } else { # subsequent
        rhs1 <- f - f0 - alam * slope
        rhs2 <- f2 - f0 - alam2 * slope
        a <- ((rhs1 / (alam^2)) - (rhs2 / (alam2^2))) / (alam - alam2)
        b <- ((-alam2 * rhs1 / (alam^2)) + (alam * rhs2 / (alam2^2))) / (alam - alam2)
        if (a == 0.0) {
          scale <- -slope / (2.0 * b)
        } else {
          disc <- b^2 - 3 * a * slope
          if (disc < 0.0) {
            scale <- lmax * alam
          } else if (b <= 0.0) {
            scale <- (-b + sqrt(disc)) / (3.0 * a)
          } else {
            scale <- -slope / (b + sqrt(disc))
          }
        }
        if (scale > lmax * alam) {
          scale <- lmax * alam
        }
      }
    }
    alam2 <- alam
    f2 <- f
    alam <- max(scale, lmin * alam)
  }
}

# 5. LOGISTIC
logistic <- function(gg) {
  return(exp(gg) / (1 + exp(gg)))
}

# 6. MVPOLY
mvpoly <- function(delta0, y, x, z) {
  n <- nrow(x)
  p <- ncol(x)
  k <- ncol(y)
  
  ik <- diag(k)
  i_n <- diag(n)
  uk <- matrix(1,k, k)
  yy <- matrix(y,nrow = n*k,ncol = 1,byrow = F) 
  xx <- kron(ik,x) 
  onek <- matrix(1,k, 1)
  inonek <- kron(onek,i_n) 
  xxz <- xx %*% z
  
  niter <- 1
  rerror <- 1.0
  delta <- as.matrix(delta0)
  
  while ((niter < 100) & (rerror > 0.000001)) {
    beta <- z %*% delta
    lxx <- xx %*% beta
    pxx <- exp(lxx)
    pxx <- matrix(pxx, n, k)
    spxx <- rowSums(pxx)
    pxx <- pxx / (1 + repmat(spxx,1,k))
    pxx <- matrix(pxx, n * k, 1)
    dd <- Matrix(0,n*k,n*k)
    diag(dd) <- pxx
    aa <- Matrix(dd %*% inonek)
    aa2 <- aa %*% t(aa)
    ww <- dd - aa2
    wwyyw <- yy - pxx + ww %*% lxx
    xxzwxxz <- t(xxz) %*% ww %*% xxz  # Information matrix
    ixxzwxxz <- solve(xxzwxxz)
    delta <- ixxzwxxz %*% (t(xxz) %*% wwyyw)
    rerror <- max(abs(delta - delta0) / max(abs(delta0), 0.1))
    niter <- niter + 1
    delta0 <- delta
  }
  
  v <- ixxzwxxz
  return(list(delta = delta, v = v))
}

# Running the EH ANALYSIS
yy <- read.csv("yy.csv",header = F)
yy <- as.matrix(yy)
xx <- read.csv("x.csv",header=F)
xx <- as.matrix(xx)
zz <- read.csv("zz.csv",header = F)
zz <- as.matrix(zz)

status <- 1
initial <- 1
method <- 1

EH_analysis_packed(yy,zz,status,xx,initial,method)

# In R: 1.070748e-27
# In Matlab: 1.1257e-27 
