# Required libraries

library(questionr)
library(dplyr)
library(spida2)
library(nlme)
library(tibble)
library(multipol)
library(pracma)
library(matconv)
library(matlab)
library(Matrix)
library(AlgDesign)

# Required Functions 

# 1. LOGISTIC REGRESSION POLY

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

# 2. MVPOLY

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
    dd <- matrix(0,n*k,n*k)
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

# 3. LOGISTIC
logistic <- function(gg) {
  return(exp(gg) / (1 + exp(gg)))
}


# Integrated Function 

cgml <- function(data, ydep, wdep, zvar, dvar, id, design_matrix,print, ssest = NULL, csest = NULL,
                 status,initial,response_level) {
  
  ydep_names <- ydep
  ydep <- data.frame(data[,ydep,drop=FALSE])
  wdep <- data.frame(data[,wdep,drop=FALSE])
  zvar <- data.frame(data[,zvar,drop=FALSE])
  dvar <- data.frame(data[,dvar,drop=FALSE])
  id <- data.frame(data[,id,drop=FALSE])
  
  stop <- 0
  
  if (length(data) == 0) {
    stop <- 1
    stop("ERROR: DATA is not provided.")
  }
  
  if (length(ydep) == 0) {
    stop <- 1
    stop("ERROR: Dependent variable YDEP is not provided.")
  }
  
  if (length(wdep) == 0) {
    stop <- 1
    stop("ERROR: Dependent variable WDEP is not provided.")
  }
  
  if (length(zvar) == 0) {
    stop <- 1
    stop("ERROR: Independent variable ZVAR is not provided.")
  }
  
  if (length(dvar) == 0) {
    stop <- 1
    stop("ERROR: Independent variable DVAR is not provided.")
  }
  
  if (length(id) == 0) {
    stop <- 1
    stop("ERROR: ID is not provided.")
  }
  
  if (stop > 0) {
    cat("NOTE: Macro stopped because of errors!")
    return("")
  }
  
  if (length(print)==0) {
    print <- "Y" 
  }
  
  
  # NAIVE ESTIMATES
  p <- ncol(zvar)
  q <- ncol(dvar)
  towork <- data
  
  for (i in 1:length(names(towork))) {
    col <- names(towork)[i]
    if (all(towork[,col] == dvar)) {
      names(towork)[i] <- "dvar"
    }else if(towork[,col,drop=FALSE]%in%ydep){
      for (j in 1:ncol(ydep)) {
        if(all(towork[,col] == ydep[,j])) {
          names(towork)[i] <- paste0("ydep",j)
          break
        } 
      }
    }else if(all(towork[,col] == wdep)) {
      names(towork)[i] <- "wdep"
    }else if(all(towork[,col] == zvar)) {
      names(towork)[i] <- "zvar"
    }else if(all(towork[,col] == id)) {
      names(towork)[i] <- "id"
    }
  }
  
  towork <- subset(towork, !is.na(ydep))
  towork <- subset(towork, !is.na(wdep))
  towork <- towork %>% arrange(id) 
  ydep_column_names <- names(towork)[grep("ydep", names(towork))]
  ydep_column_names <- sort(ydep_column_names)
  
  mx <- lme(wdep~ 0 + dvar,random = ~ 1 | id,data = towork,correlation = corAR1(form = ~1 | id))  
  mx <- coef(mx)
  mx <- rownames_to_column(mx,"id")
  zy <- towork %>%
    group_by(id) %>%
    filter(row_number() == 1) %>%
    dplyr::select(id, zvar, ydep_column_names) %>% 
    ungroup()
  
  
  zymx <- merge(zy, mx, by = "id")
  
  # MVPOLY Integration #
  yy <- as.matrix(ydep)
  xx <- as.matrix(dvar)
  zz <- as.matrix(zvar)
  xx <- cbind(xx,zz)
  xx <- as.matrix(xx)
  
  a <- data.frame("yy_factor"=as.factor(yy),"b"=1:nrow(yy))
  w <- which(sapply(a, class) == 'factor');
  df <- as.data.frame.matrix(model.matrix(
    as.formula(sprintf("b ~  0 + %s", paste0(names(w), collapse = "+"))),
    data = a))
  df <- as.matrix(df)
  df <- df[,-1,drop=FALSE]
  
  zz <- design_matrix
  zadd <- zz
  zadd_original <- zadd
  
  # ANALYSIS
  nn <- nrow(xx)  # nn: no. of subjects
  mxx <- ncol(xx)  # mxx: no. of covariates
  
  nsubtype <- nrow(zadd)
  nprm2 <- ncol(zadd)
  
  # zadd <- cbind(zadd[, 1], zadd[, 2], zadd[, 3], zadd[, 4]) Prior knowledge is necessary to include this.  
  
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
    if(nz==1){
      break
    }
    perm <- rbind(perm, t(t(seq(i, (i + (mxx - 1) * nz), by = nz))))  
  }
  
  zz <- zz[perm, ,drop=FALSE]
  
  # Fitting the two-stage model
  # Defining the response vector yyc, takes values in 0:nsubtype indicating whether
  # the subject is control (0) or the numeric level for the disease subtype
  
  yyc <- zeros(nn,1)
  for (j in 1:nz) {
    yyc[df[, j] == 1] <- j
  }
  
  # Defining initial value
  if (initial == 0) {
    theta0 <- zeros(length(zz[1, ]),1)
  }else if (initial == 1) {
    lrp_result <- logistic_regression_poly(df, xx) 
    theta00 <- lrp_result$beta  
    covmat <- lrp_result$covmat 
    theta00 <- theta00[2:(mxx+1),,drop=FALSE]
    if (nprm2==1) {
      theta0 <- as.matrix(theta00)
    }else{
      theta01 <- zeros(nprm2-1,mxx) 
      theta01 <- rbind(t(theta00),theta01) 
      for (i in 1:mxx) {  
        if (i==1) {
          theta0 <- theta01[,i,drop=FALSE]
        }else{theta0 <- rbind(theta0,theta01[,i,drop=FALSE])}
      }
    }
  }
  
  result_mvpoly <- mvpoly(theta0, df, xx, zz)
  theta <- result_mvpoly$delta
  vt <- result_mvpoly$v
  
  level <- response_level - 1
  n_level <- NULL
  if (length(level)==1) {
    if (level==0) {
      n_level <- 0
    } else{
      n_level <- 1:level
    }
  } else{
    for (k in 1:length(level)) {
      n_level <- c(n_level,1:level[k])
    }
  }
  
  iter_ss <- NULL
  iter_cs <- NULL
  
  theta_codes <- "0"
  for (k in 1:length(level)) {
    for (j in 2:nprm2) {
      if (n_level[j-1]==level[k]) {
        theta_codes <- c(theta_codes,paste0(n_level[j-1],"(",k,")"))
        break
      }else{
        theta_codes <- c(theta_codes,paste0(n_level[j-1],"(",k,")"))
      }
    }
  }
  
  nc_flag1=rep(0,nprm2)
  nc_flag2=rep(0,nprm2)
  nc_flag3=rep(0,nprm2)
  nc_flag4=rep(0,nprm2)
  
  for (j in 1:nprm2) {
    
    beta_n <- theta[c(j,j+nprm2),,drop=FALSE]
    
    dataold <- as.matrix(towork[, c("id", "zvar","ydep1", "dvar", "wdep"),drop=FALSE])
    
    rawdata <- as.matrix(towork[, c("zvar", "ydep1", "dvar", "wdep"),drop=FALSE])
    
    ydep <- as.matrix(towork[,"ydep1",drop=FALSE])
    
    id <- as.matrix(towork[,"id",drop=FALSE])
    
    id <- rbind(id,0)
    # Initialize variables
    totobs <- nrow(ydep)
    nsub <- 0
    
    i <- 1
    while(i <= totobs){
      start <-  i
      while ((id[i,] == id[start,]) & (i <= totobs)){
        i <-  i + 1
      }
      if(start == 1){
        ystart=start
      }else{
        ystart <-  rbind(ystart, start)
      }
      nsub = nsub + 1
    }
    
    if(nsub==1){
      yend <- totobs
    }else{
      yend <- rbind(ystart[2:nsub,,drop=FALSE]-1,totobs)
    }
    rr <- 0
    mq <- 0
    # Loop over each subject group
    for (i in 1:nsub) {
      l <- ystart[i,]
      u <- yend[i,]
      
      
      DATAi <- rawdata[l:u, ,drop=FALSE]
      mi <- nrow(DATAi)
      
      
      temp <- matrix(0, ncol = (p+q+2)) 
      
      
      for (k in 1:mi) {
        if (!is.na(DATAi[k, 1])) {
          temp <- rbind(temp, DATAi[k, ,drop=FALSE])
        }
      }
      ntemp <- nrow(temp)
      if (ntemp == 1) {next}
      DATAi <- temp[2:ntemp,,drop=FALSE]
      di <- DATAi[, (p+2):(p+q+1),drop=FALSE]
      wi <- DATAi[, (p+q+2),drop=FALSE]
      
      
      dd <- t(di) %*% di
      detdd <- det(dd)
      if (detdd != 0) {
        ri <- wi - di %*% solve(t(di)%*%di) %*% t(di) %*% wi
        rr <- rr + (t(ri)%*%ri)
        mq <- mq + (mi-q)
      }
    }
    
    sgsqu_n <- rr / mq
    beta_n <-  as.matrix(beta_n)
    theta_n <-  rbind(beta_n,sgsqu_n)
    
    maxiter = 10000
    
    # Sufficiency Score # 
    theta_ss <- theta_n
    theta_max <- theta_n
    theta_ss_old <- theta_ss 
    iter <- 1
    # Start iterative process
    for (iter in 1:maxiter){
      while(theta_ss[(p+q+1),1]<= 0) {
        theta_ss <- (theta_ss + theta_ss_old) / 2
      }
      if (iter==2) {
        dh_max <- dethess
        theta_max <- theta_ss_old
      }
      
      
      theta_ss_old <- theta_ss
      b0 <- theta_ss[1:p,1,drop=FALSE]
      b1 <- theta_ss[(p+1):(p+q),1,drop=FALSE]
      sgu = theta_ss[p+q+1, 1,drop=FALSE]
      ntol  = 0
      score = matrix(0, (p+q+1), 1)
      hess  = matrix(0, (p+q+1), (p+q+1))
      scsc  = matrix(0, (p+q+1), (p+q+1))
      for(i in 1:nsub){
        l = ystart[i,]  
        u = yend[i,]
        DATAi = rawdata[l:u, ,drop=FALSE]
        mi = nrow(DATAi)
        temp = matrix(0, 1, (p+q+2))
        for(k in 1:mi){
          if(!any(is.na(DATAi[k, 1]))) { temp = rbind(temp,DATAi[k, ,drop=FALSE])}
        }
        ntemp = nrow(temp)
        if(ntemp ==1){
          next}
        DATAi = temp[2:ntemp, ,drop=FALSE]
        zi = DATAi[1, 1:p,drop=FALSE]   
        zi = t(zi)
        yi = DATAi[1, (p+1),drop=FALSE]
        di = DATAi[ , (p+2):(p+q+1),drop=FALSE]
        wi = DATAi[ , (p+q+2),drop=FALSE]
        q  = ncol(di)
        dd = t(di)%*%di
        detdd = det(dd)
        if(detdd == 0){next}
        dd     = solve(dd)
        si     = t(di) %*% wi + yi %*% sgu %*% b1
        mu     = t(b0) %*% zi + t(si - sgu %*% b1 / 2) %*% dd %*% b1
        if(mu >= 0){ mu = 1 / (1 + exp(-mu))}
        if(mu < 0){ mu = exp(mu) / (1 + exp(mu))}
        ymu    = yi - mu
        mumu1  = mu * (mu-1)
        ddsub1 = dd %*% (si - sgu %*% b1)
        b1ddb1 = t(b1) %*% dd %*% b1
        udd    = sgu %*% dd
        ddb1   = dd %*% b1
        sgu2   = sgu**2
        ri     = wi - di %*% dd %*% t(di) %*% wi
        rr     = t(ri) %*% ri
        sc1 = ymu %*% zi
        sc2 = ymu %*% ddsub1
        sc3 = - (mi - q) / 2 / sgu + rr / 2 / sgu2 - ymu %*% b1ddb1 / 2 
        sci   = rbind(sc1, sc2, sc3)
        score = score + sci
        hess11 = mumu1 %*% zi %*% t(zi)
        hess12 = mumu1 %*% zi %*% t(ddsub1)
        hess13 = - mumu1 %*% b1ddb1 %*% zi / 2;
        hess22 = mumu1 %*% ddsub1 %*% t(ddsub1) - ymu %*% udd
        hess23 = - mumu1 %*% b1ddb1 %*% ddsub1 / 2 - ymu %*% ddb1
        hess33 = (mi-q)/2/sgu2 - rr/sgu/sgu2 + mumu1*(b1ddb1**2)/4
        hess1 =  cbind(hess11   ,  hess12   ,hess13)
        hess2 = cbind(t(hess12) ,  hess22   ,hess23)
        hess3 = cbind(t(hess13) , t(hess23) ,hess33)
        hess = hess + rbind(hess1,hess2,hess3)
        ntol <- ntol +1
      }
      dethess = det(hess);
      if (iter > 1){
        if (dethess > dh_max){
          dh_max = dethess
          theta_max = theta_ss
        }
      }
      if (abs(dethess) < 1.0e-04){
        theta_ss = theta_max
        break
      }
      theta_ss = (theta_ss_old) - ( solve(hess) %*% score )
      diffmax = 0
      for(k in 1:(p+q+1)){
        diff = abs(theta_ss[k,1] - theta_ss_old[k,1])
        if(diff > diffmax){diffmax = diff}
      }
      if(diffmax < 10**(-6)){break}
      if(iter == maxiter){
        nc_flag1[j]=1
        theta_ss = theta_max
        se = matrix(NA, (p+q+1), 1)
      }
    }
    iter_ss[j] <- iter
    
    
    
    se <- matrix(NA, nrow = p + q + 1, ncol = 1)
    
    if (!identical(theta_ss, matrix(NA, nrow = p + q + 1, ncol = 1))) {
      b0 <- theta_ss[1:p,1,drop=FALSE]
      b1 <- theta_ss[(p + 1):(p + q),1,drop=FALSE]
      sgu <- theta_ss[p + q + 1,1,drop=FALSE]
      
      score <- matrix(0, nrow = p + q + 1, ncol = 1)
      hess <- matrix(0, nrow = p + q + 1, ncol = p + q + 1)
      scsc <- matrix(0, nrow = p + q + 1, ncol = p + q + 1)
      
      for (i in 1:nsub) {
        l <- ystart[i,]
        u <- yend[i,]
        DATAi <- rawdata[l:u, , drop = FALSE]
        
        temp <- matrix(0, nrow = 1, ncol = p + q + 2)
        mi <- nrow(DATAi)
        for (k in 1:mi) {
          if (!is.na(DATAi[k, 1])) {
            temp <- rbind(temp, DATAi[k, ,drop=FALSE])
          }
        }
        ntemp <- nrow(temp)
        if (ntemp == 1) {
          next
        }
        DATAi <- temp[2:ntemp,,drop=FALSE]
        zi <- DATAi[1,1:p,drop=FALSE]
        zi <- t(zi)
        yi <- DATAi[1,p+1,drop=FALSE]
        di <- DATAi[,(p+2):(p+q+1),drop=FALSE]
        wi <- DATAi[,(p+q+2),drop=FALSE]
        q <- ncol(di)
        dd <- t(di)%*%di
        detdd <- det(dd)
        
        if (det(dd) == 0) {
          next
        }
        dd <- solve(dd)  
        si <- t(di) %*% wi + yi * sgu * b1
        mu <- t(b0) %*% zi + (t(si - sgu * b1 / 2) %*% dd %*% b1)
        mu <- ifelse(mu >= 0, 1 / (1 + exp(-mu)), exp(mu) / (1 + exp(mu)))
        ymu <- yi - mu
        y1 <- yi - 1
        y12 <- yi - 0.5
        mumu1 <- mu * (mu - 1)
        ddsub1 <- dd %*% (si - sgu %*% b1)
        b1ddb1 <- t(b1) %*% dd %*% b1
        dwy1ub1 <- t(di) %*% wi + (yi - 1) %*% sgu %*% b1
        dw2y1ub1 <- t(di) %*% wi + (2 * yi - 1) %*% sgu %*% b1
        udd <- sgu %*% dd
        ddb1 <- dd %*% b1
        sgu2 <- sgu^2
        q <- ncol(di)
        ri <- wi - di %*% dd %*% t(di) %*% wi
        rr <- t(ri) %*% ri
        sc1 <- ymu %*% zi
        sc2 <- ymu %*% ddsub1
        sc3 <- - (mi - q) / 2 / sgu + rr / 2 / sgu2 - ymu %*% b1ddb1 / 2
        
        
        sci <- rbind(sc1, sc2, sc3)  
        scsc <- scsc + sci %*% t(sci)
        
        
        hess11 <- mumu1 * zi %*% t(zi)
        hess12 <- mumu1 * zi %*% t(ddsub1)
        hess13 <- mumu1 * y12 %*% b1ddb1 %*% zi
        hess21 <- mumu1 * dd %*% dwy1ub1 %*% t(zi)
        hess22 <- mumu1 * dd %*% dwy1ub1 %*% t(dw2y1ub1) * dd + ymu %*% y1 %*% udd
        hess23 <- mumu1 * y12 * b1ddb1 %*% dd %*% dwy1ub1 + ymu %*% y1 %*% ddb1
        hess31 <- -mumu1 * b1ddb1 %*% t(zi) / 2
        hess32 <- -mumu1 * b1ddb1 %*% t(dw2y1ub1) %*% dd / 2 - ymu %*% t(ddb1)
        hess33 <- (mi - q) / 2 / sgu2 - rr / sgu / sgu2 - mumu1 %*% y12 %*% (b1ddb1^2) / 2
        
        
        hess <- hess + rbind(cbind(hess11, hess12, hess13),
                             cbind(hess21, hess22, hess23),
                             cbind(hess31, hess32, hess33))
      }
      
      sdwh <- -hess
      detsdwh <- det(sdwh)
      adj <- ifelse(abs(detsdwh) > 1e20, 1e-10, 1e10)
      sdwh <- sdwh * adj
      detsdwh <- det(sdwh)
      
      if (detsdwh == 0) {
        nc_flag2[j] <- 1
      } else {
        sdwh <- solve(sdwh)
        sdwh <- sdwh %*% scsc %*% t(sdwh) * adj * adj
        var <- as.matrix(diag(sdwh))
        se <- sqrt(var)
      }
    }
    
    
    p_value <- matrix(NA,nrow = p + q + 1,ncol = 1)
    for (k in 1:(p + q + 1)) {
      p_value[k,] <- 2 * (1 - pnorm(abs(theta_ss[k,] / se[k,])))
    }
    
    
    ssest <- cbind(theta_ss, se, p_value)
    
    
    colnames(ssest) <- c("Estimate", "se", "p_value")
    varnames <- c(sapply(1:p, function(i) toupper(colnames(zvar)[i])),
                  sapply(1:q, function(i) toupper(colnames(dvar)[i])),
                  "sigma_u^2")
    rownames(ssest) <- varnames
    ss_est <- ssest
    
    
    if (j==1) {
      ssest_responses <- NULL
    }
    rownames(ssest) <- paste0(varnames,"_",theta_codes[j])
    
    ssest_responses <- rbind(ssest_responses,ssest) 
    
    
    # Conditional Score
    theta_cs <- theta_n
    theta_max <- theta_n
    iter <- 1
    
    # Start iterative process
    for (iter in 1:maxiter) {
      while (theta_cs[(p + q + 1), 1] <= 0) {
        theta_cs <- (theta_cs + theta_cs_old) / 2
      }
      if (iter == 2) {
        dh_max <- dethess
        theta_max <- theta_cs_old
      }
      
      theta_cs_old <- theta_cs
      b0 <- theta_cs[1:p, 1,drop=FALSE]
      b1 <- theta_cs[(p + 1):(p + q), 1,drop=FALSE]
      sgu <- theta_cs[(p + q + 1), 1,drop=FALSE]
      ntol <- 0
      score <- matrix(0, nrow = p + q + 1, ncol = 1)
      hess <- matrix(0, nrow = p + q + 1, ncol = p + q + 1)
      scsc <- matrix(0, nrow = p + q + 1, ncol = p + q + 1)
      
      for (i in 1:nsub) {
        l <- ystart[i,]
        u <- yend[i,]
        DATAi <- rawdata[l:u, , drop = FALSE]
        mi <- nrow(DATAi)
        
        temp <- matrix(0, nrow = 1, ncol = p + q + 2)
        for (k in 1:mi) {
          if (!any(is.na(DATAi[k, 1]))) {
            temp <- rbind(temp, DATAi[k, ,drop=FALSE])
          }
        }
        
        ntemp <- nrow(temp)
        if (ntemp == 1){next}
        
        DATAi <- temp[2:ntemp, , drop = FALSE]
        zi <- DATAi[1, 1:p,drop=FALSE]
        zi <- t(zi)
        yi <- DATAi[1, (p + 1),drop=FALSE]
        di <- DATAi[, (p + 2):(p + q + 1),drop=FALSE]
        wi <- DATAi[, (p + q + 2),drop=FALSE]
        q <- ncol(di)
        dd <- t(di) %*% di
        if (det(dd) == 0){next}
        
        dd <- solve(dd)
        si <- t(di) %*% wi + yi %*% sgu %*% b1
        mu <- t(b0) %*% zi + t(si - sgu %*% b1 / 2) %*% dd %*% b1
        if(mu >= 0){ mu = 1 / (1 + exp(-mu))}
        if(mu < 0){ mu = exp(mu) / (1 + exp(mu))}
        ymu <- yi - mu
        mumu1 <- mu * (mu-1)
        ddsub1 <- dd %*% (si - sgu %*% b1)
        b1ddb1 <- t(b1) %*% dd %*% b1
        udd <- sgu %*% dd
        ddb1 <- dd %*% b1
        muub1 <- mu %*% sgu %*% b1
        yub1 <- yi %*% sgu %*% b1
        sgu2 <- sgu**2
        ri <- wi - di %*% dd %*% t(di) %*% wi
        rr <- t(ri) %*% ri 
        sc1 <- ymu %*% zi
        sc2 <- ymu %*% dd %*% (si-muub1) 
        sc3 <- -(mi - q)/2/sgu +  rr/2/sgu2 + ymu %*% (1/2 - mu) %*% b1ddb1
        sci <- rbind(sc1, sc2, sc3) 
        score <- score + sci
        
        hess11 <- mumu1 %*% zi %*% t(zi)
        hess12 <- mumu1 %*% zi %*% t(ddsub1)
        hess13 <- - mumu1 %*% b1ddb1 %*% zi / 2
        hess21 <- mumu1 %*% dd %*% (si + yub1 - 2 %*% muub1) %*% t(zi)
        hess22 = mumu1 %*% dd %*% (si+yub1-2%*%muub1)%*%t(ddsub1) - mu%*%ymu%*%udd
        hess23 = - mumu1 %*% b1ddb1%*%dd%*%(si+yub1-2%*%muub1)/2-mu%*%ymu%*%ddb1
        hess31 = mumu1 %*% (yi - 2*mu + 1/2) %*% b1ddb1 %*% t(zi)
        hess32 = mumu1 %*% (yi-2*mu+1/2)%*%b1ddb1%*%t(ddsub1) + ymu%*%(1 - 2*mu)%*%t(ddb1)
        hess33 = (mi-q)/2/sgu2 - rr/sgu/sgu2 - mumu1%*%(yi-2*mu+1/2)%*%(b1ddb1**2)/2
        hess1 = cbind(hess11,hess12, hess13)
        hess2 = cbind(hess21, hess22,hess23)
        hess3 = cbind(hess31, hess32,hess33)
        hess = hess + rbind(hess1,hess2,hess3)
        ntol = ntol + 1
      }
      
      dethess <- det(hess)
      if (iter > 1){
        if(dethess > dh_max){
          dh_max <- dethess
          theta_max <- theta_cs
        }
      }
      
      if (abs(dethess) < 1e-6) {
        theta_cs <- theta_max
        break 
      }
      
      
      theta_cs <- theta_cs_old - (solve(hess) %*% score)
      
      
      diffmax <- 0
      
      for (k in 1:(p+q+1)) {
        diff <- abs(theta_cs[k,1] - theta_cs_old[k,1])
        if (diff>diffmax) {
          diffmax <- diff
        }
      }
      
      if (diffmax < 10**-6) {
        break
      }  
      if (iter == maxiter) {
        nc_flag3[j] <- 1 
        theta_cs <- theta_max
        se <- matrix(NA,(p+q+1),1)
      }
    }
    
    iter_cs[j] <- iter
    
    se <- matrix(NA, nrow = p + q + 1, ncol = 1)
    
    if (!identical(theta_cs, matrix(NA, nrow = p + q + 1, ncol = 1))) {
      b0 <- theta_cs[1:p,1,drop=FALSE]
      b1 <- theta_cs[(p + 1):(p + q),1,drop=FALSE]
      sgu <- theta_cs[p + q + 1,1,drop=FALSE]
      
      score <- matrix(0, nrow = p + q + 1, ncol = 1)
      hess <- matrix(0, nrow = p + q + 1, ncol = p + q + 1)
      scsc <- matrix(0, nrow = p + q + 1, ncol = p + q + 1)
      
      for (i in 1:nsub) {
        l <- ystart[i,]
        u <- yend[i,]
        DATAi <- rawdata[l:u, , drop = FALSE]
        
        temp <- matrix(0, nrow = 1, ncol = p + q + 2)
        mi <- nrow(DATAi)
        for (k in 1:mi) {
          if (!is.na(DATAi[k, 1])) {
            temp <- rbind(temp, DATAi[k, ,drop=FALSE])
          }
        }
        ntemp <- nrow(temp)
        if (ntemp == 1) {
          next
        }
        DATAi <- temp[2:ntemp,,drop=FALSE]
        zi <- DATAi[1,1:p,drop=FALSE]
        zi <- t(zi)
        yi <- DATAi[1,(p+1),drop=FALSE]
        di <- DATAi[,(p+2):(p+q+1),drop=FALSE]
        wi <- DATAi[,(p+q+2),drop=FALSE]
        q <- ncol(di)
        dd <- t(di)%*%di
        detdd <- det(dd)
        
        
        if (det(dd) == 0) {
          next
        }
        dd <- solve(dd)  
        si <- t(di) %*% wi + yi %*% sgu %*% b1
        mu <- t(b0) %*% zi + (t(si - sgu %*% b1 / 2) %*% dd %*% b1)
        mu <- ifelse(mu >= 0, 1 / (1 + exp(-mu)), exp(mu) / (1 + exp(mu)))
        ymu <- yi - mu
        ymu2 <- ymu**2
        y1 <- yi - 1
        y12 <- yi - 0.5
        mumu1 <- mu * (mu - 1)
        muub1 <- mu %*% sgu %*% b1
        b1ddb1 <- t(b1) %*% dd %*% b1
        dwy1ub1 <- t(di) %*% wi + (yi - 1) %*% sgu %*% b1
        dw2y1ub1 <- t(di) %*% wi + (2 * yi - 1) %*% sgu %*% b1
        dw2ymuub1 = t(di) %*% wi + 2 * (yi - mu) %*% sgu %*% b1
        y2mu12 <- yi - 2 * mu + 0.5
        udd <- sgu %*% dd
        ddb1 <- dd %*% b1
        sgu2 <- sgu^2
        ri <- wi - di %*% dd %*% t(di) %*% wi
        rr <- t(ri) %*% ri
        sc1 <- ymu %*% zi
        sc2 <- ymu %*% dd %*% (si - muub1)
        sc3 <- - (mi - q) / 2 / sgu + rr / 2 / sgu2 - ymu %*% (0.5-mu) %*% b1ddb1
        
        sci <- rbind(sc1, sc2, sc3) 
        scsc <- scsc + sci %*% t(sci)
        
        
        hess11 <- mumu1 * zi %*% t(zi)
        hess12 <- mumu1 * zi %*% t(dw2y1ub1) %*% dd
        hess13 <- mumu1 * y12 %*% b1ddb1 %*% zi
        hess21 <- mumu1 * dd %*% dw2ymuub1 %*% t(zi)
        hess22 <- mumu1 * dd %*% dw2ymuub1 %*% t(dw2y1ub1) * dd + ymu2 %*% udd
        hess23 <- mumu1 * y12 %*% b1ddb1 %*% dd %*% dw2ymuub1 + ymu2 %*% ddb1
        hess31 <- mumu1 * y2mu12 %*% b1ddb1 %*% t(zi)
        hess32 <- mumu1 * y2mu12 %*% b1ddb1 %*% t(dw2y1ub1) %*% dd + ymu %*% (1-2*mu) %*% t(ddb1)
        hess33 <- (mi - q) / 2 / sgu2 - rr / sgu / sgu2 + mumu1 %*% y12 %*% y2mu12 %*% (b1ddb1^2)
        
        
        hes1 <- cbind(hess11,hess12,hess13)
        hes2 <- cbind(hess21,hess22,hess23)
        hes3 <- cbind(hess31,hess32,hess33)
        hess <- hess + rbind(hes1,hes2,hes3) 
      }
      
      sdwh <- -hess
      detsdwh <- det(sdwh)
      adj <- ifelse(abs(detsdwh) > 1e20, 1e-10, 1e10)
      sdwh <- sdwh * adj
      detsdwh <- det(sdwh)
      
      if (detsdwh == 0) {
        nc_flag4[j] <- 1
      } else {
        sdwh <- solve(sdwh)
        sdwh <- sdwh %*% scsc %*% t(sdwh) * adj * adj
        var1 <- as.matrix(diag(sdwh))
        se <- sqrt(var1)
      }
    }
    
    p_value <- matrix(NA,(p+q+1),1)
    for (k in 1:(p+q+1)) {
      p_value[k,] <- 2 * (1 - pnorm(abs(theta_cs[k,] / se[k,]))) 
    }
    csest <- cbind(theta_cs, se, p_value)
    
    colnames(csest) <- c("Estimate", "se", "p_value")
    varnames <- c(sapply(1:p, function(i) toupper(colnames(zvar)[i])), 
                  sapply(1:q, function(i) toupper(colnames(dvar)[i])),
                  "sigma_u^2")
    rownames(csest) <- varnames
    cs_est <- csest
    
    if (j==1) {
      csest_responses <- NULL
    }
    
    rownames(csest) <- paste0(varnames,"_",theta_codes[j])
    
    csest_responses <- rbind(csest_responses,csest) 
    
  }
  
  if(print != "N") {
    description <- c("Number of Subjects" = nsub, "Number of Observations" = totobs)
    cat("Model Fitting Information\n")
    print(description)
    
    cat("\nResponse Variable Name\n")
    print(ydep_names)
    start_value <- theta
    cat("\nStarting Values\n")
    initial_thetas <- matrix(theta,ncol = nprm2,nrow = ncol(xx))
    
    rownames(initial_thetas) <- c(varnames[1:ncol(xx)])
    colnames(initial_thetas) <- rep("1",nprm2)
    for (j in 1:nprm2) {
      colnames(initial_thetas)[j] <- paste0("Theta_",theta_codes[j])
    }
    print(initial_thetas)
    
    iterations <- rbind(iter_ss,iter_cs)
    rownames(iterations) <- c("Sufficiency Score","Conditional Score")
    colnames(iterations) <- colnames(initial_thetas)
    
    cat("\nNumber of Iterations for each method\n")
    print(iterations)
    
    for (j in 1:nprm2) {
      text_nlevel <- paste0("For Theta_",theta_codes[j])
      if(nc_flag1[j] == 1) cat("\n",text_nlevel,"Sufficiency score method did not converge\n")
      if(nc_flag2[j] == 1) cat("\n",text_nlevel,"Sufficiency score method stopped because the empirical sandwich matrix is singular -- Standard error is not provided\n")
      if(nc_flag3[j] == 1) cat("\n",text_nlevel,"Conditional score method did not converge\n")
      if(nc_flag4[j] == 1) cat("\n",text_nlevel,"Conditional score method stopped because the empirical sandwich matrix is singular -- Standard error is not provided\n")
      
      
    }
    
  }
  if (print!="N") {
    cat("\nSUFFICIENCY SCORE METHOD\nEstimates of the Model Parameters:\n")
    print(ssest_responses)
    
    cat("\nCONDITIONAL SCORE METHOD\nEstimates of the Model Parameters:\n")
    print(csest_responses)
    
  }
  
}

# Importing the dataset
yy <- read.csv("seker_02122013.csv",header = T)

# Transforming the dataset
testdata <- transform(yy,
                      akstrsqr = 1/sqrt(aks),
                      logaks = log(aks),
                      logpp2 = log(pp2),
                      dogumsekli = dogumsekli-1,                      dogumkilosu01 = as.integer(dogumkilosu > 4000),
                      dogumhaftasi01 = as.integer(dogumhaftasi > 37))


# Dogumsekli (1x1)

design_matrix <- matrix(1,1,1)

cgml(data=testdata, ydep="dogumsekli", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = 1)


# Dc (4x1)

zz <- matrix(c(1,1,1,1,
               1,0,0,0,
               0,1,0,0,
               0,0,1,0),ncol = 4,nrow = 4,byrow = F)

design_matrix <- zz

cgml(data=testdata, ydep="dc", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = 4)


# Dc (2x2)
zz <- matrix(c(1,1,1,1,
               1,1,0,0,
               1,0,1,0),ncol = 3,nrow = 4,byrow = F)

design_matrix <- zz

cgml(data=testdata, ydep="dc", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = c(2,2))


# Yasayan (1x5)
zz <- matrix(c(1,1,1,1,1,
               1,0,0,0,0,
               0,1,0,0,0,
               0,0,1,0,0,
               0,0,0,1,0),ncol = 5,nrow = 5,byrow = F)

design_matrix <- zz
testdata[testdata$yasayan==7,"yasayan"] <- 5

cgml(data=testdata, ydep="yasayan", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = 5)


# Parite (1x6)

zz <- matrix(c(1,1,1,1,1,1,
               1,0,0,0,0,0,
               0,1,0,0,0,0,
               0,0,1,0,0,0,
               0,0,0,1,0,0,
               0,0,0,0,1,0),ncol = 6,nrow = 6,byrow = F)

design_matrix <- zz
testdata[testdata$parite==7,"parite"] <- 6

cgml(data=testdata, ydep="parite", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = 6)



# Parite (3x2)

zz <- matrix(c(1,1,1,1,1,1,
               1,1,0,0,0,0,
               0,0,1,1,0,0,
               0,1,0,1,0,1),nrow = 6,ncol = 4,byrow = F)

design_matrix <- zz

cgml(data=testdata, ydep="parite", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = c(3,2))


# Gebelik (1x9)

zz <- matrix(c(1,1,1,1,1,1,1,1,1,
               0,1,0,0,0,0,0,0,0,
               0,0,1,0,0,0,0,0,0,
               0,0,0,1,0,0,0,0,0,
               0,0,0,0,1,0,0,0,0,
               0,0,0,0,0,1,0,0,0,
               0,0,0,0,0,0,1,0,0,
               0,0,0,0,0,0,0,1,0,
               0,0,0,0,0,0,0,0,1),ncol = 9,nrow = 9,byrow = F)

design_matrix <- zz

cgml(data=testdata, ydep="gebelik", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = c(9))


# Gebelik (3x3)

zz <- matrix(c(1,1,1,1,1,1,1,1,1,
               1,1,1,0,0,0,0,0,0,
               0,0,0,1,1,1,0,0,0,
               1,0,0,1,0,0,1,0,0,
               0,1,0,0,1,0,0,1,0),ncol = 5,nrow = 9,byrow = F)

design_matrix <- zz

cgml(data=testdata, ydep="gebelik", wdep="aks", zvar="yas", dvar="gebelikhaftasi", id="arsivno",
     design_matrix = design_matrix,print="Y", ssest="est1", csest="est2",status = 1,initial = 0,response_level = c(3,3))