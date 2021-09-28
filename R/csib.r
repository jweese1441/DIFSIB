#' Crossing SIBTEST
#'
#' This function assesses non-uniform DIF with one crossing point using the Crossing-SIBTEST procedure (Li & Stout, 1996).
#' The function outputs the non-uniform and uniform DIF magnitudes as well as the significance testing.
#'
#' @param data_ref  The dataset that contains the reference group participants. Dichotomous data must only contain 0's and 1's.
#' @param data_foc  The dataset that contains the focal group participants. Dichotomous data must only contain 0's and 1's.
#' @param minc  The minimum number of individuals in a give total test score category. Initialized to 2 per Shealy and Stout (1993) recommendation
#' @param cusr  The guess correction factor. This is initialized to 0.2, but can be changed based on the items used. The cusr value is the average of guess correction across the items in the dataset.
#' @param idw  The weighting factor. This is initialized to use the reference and focal group weighting recommended by Shealy and Stout (1993)
#' @param suspect_items  the item(s) to be assess for DIF. Listing one item number will assess for DIF on that item, multiple items (e.g. c(2,3,8)) will assess DIF for the bundle of items.
#' @param matching_items  the list of items to match on.
#' @param listwise  initialized to 1, this inidicates that the procedure will not delete any rows with missing data, changing this to 2 will listwise delete.
#' @aliases csib
#' @export csib
#' @references Chalmers, R. P. (2018). Improving the Crossing-SIBTEST statistic for detecting non-uniform DIF. Psychometrika, 83, 2, 376-386.
#' @references DIF-Pack (2021) Measured Progress. https://psychometrics.onlinehelp.measuredprogress.org/tools/dif/
#' @references Li, H.-H. & Stout, W. (1996). A new procedure for detection of crossing DIF. Psychometrika, 61, 647-677.
#' @references Shealy, R. & Stout, W. (1993). A model-based standardization approach that separates true bias/DIF from group ability differences and detect test bias/DTF as well as item bias/DIF. Psychometrika, 58, 159-194.
#' @examples
#' 
#' \dontrun{
#'
#' #perform Crossing-SIBTEST with one suspect item 
#' 
#' csib(data_ref = data_refC, data_foc = data_focC ,minc = 2, 
#' cusr = .2,  suspect_items = c(20), matching_items = c(1:19))
#'
#' #perform Crossing-SIBTEST with one suspect item and no guessing (cusr = 0) 
#'
#' csib(data_ref = data_refC, data_foc = data_focC, minc = 2, 
#' cusr = 0, suspect_items = c(20), matching_items = c(1:19))
#'
#' #perform Crossing-SIBTEST with a bundle of suspect items 
#'
#' csib(data_ref = data_refC, data_foc = data_focC, minc = 2, 
#' cusr = .2, suspect_items = c(16,17,18,19,20), matching_items = c(1:15))
#'
#' }
csib <- function(data_ref, data_foc,minc=2,cusr = .2,idw = 0,suspect_items, matching_items,listwise = 1){
  ############################################################################################################1
  ############################################################################################################1
  # NOTE: THIS CODE IS ADAPTED FROM THE ORIGINAL FORTRAN CODE THAT POWERED THE SIBTEST V1.7 (2005 SOFTWARE
  # RELEASE DATE). THAT CODE HAS BEEN RELEASED BY MEASURED PROGRESS AND IS NOW CONSIDERED OPEN SOURCED.
  # YOU CAN FIND THAT CODE AT THE  FOLLOWING WEBSITE
  # https://psychometrics.onlinehelp.measuredprogress.org/tools/dif/.
  #
  # I take no credit for creating the code, only the translation and it's implementation in R
  # This code will be updated to become vectorized so it is more efficient.
  #
  ############################################################################################################1
  ############################################################################################################1
  if(listwise == 2){
    data_ref <- na.omit(data_ref)
    data_foc <- na.omit(data_foc)
  }
  if (any(is.na(data_ref)) || any(is.na(data_foc))) {
    return("SIBTEST requires that both reference and focal group datasets do not have missing values.
            Please select TRUE to remove missing values to listwise remove missing values")}
  if(ncol(data_ref) != ncol(data_foc)){
    return("Verify that your data is correct. The number of columns in the reference group data does not equal the number of columns in the focal group data
		this will cause an issue when selecting matching and suspect items.")}
  if(any(suspect_items %in% matching_items)){
    return("SIBTEST cannot run when when suspect items and matching subtest items contain the same the same items
            Please check your suspect and matching subtest item selections")}
  sibuni <- function(data_ref,data_foc,xmaj,ymaj,j0maj,xmin,ymin,j0min,itdifj,
                     itdifn,ivalid,n0,m0,ni,minc,cusr,iflg,idw,
                     buni,r1,se,wgt,adjmj,adjmn,xnk,ybarmj,ybarmn,jknmaj,jknmin
                     ,kr20pn,kr20pj,ehr,ehf){
    estaujn <- estaunn <- estauj <-estaun <-NULL
    adjmj <- adjmn <- obprmj<-obprmn<-scores<-iclind <- NULL
    itdfwj = itdifj
    itdfwn = itdifn
    ybarmj <- ybarmn <-ysum2j <- ysum2n <-
      spbmaj <- spbmin <-spbmj2 <-spbmn2 <-rep(0,(ni+1))
    j0maj <- nrow(data_ref)
    j0min <- nrow(data_foc)
    for(j in 1:j0maj){
      spbmaj[xmaj[j]] = spbmaj[xmaj[j]]+ymaj[j]
      spbmj2[xmaj[j]] = spbmj2[xmaj[j]]+(ymaj[j]*ymaj[j])
    }
    for(j in 1:j0min){
      spbmin[xmin[j]] = spbmin[xmin[j]]+ymin[j]
      spbmn2[xmin[j]] = spbmn2[xmin[j]]+(ymin[j]*ymin[j])
    }
    ybarmj[c(2:(m0+1))] <- spbmaj[c(1:(m0))]/n0
    ybarmn[c(2:(m0+1))] <- spbmin[c(1:(m0))]/n0
    ysum2j[c(2:(m0+1))] <- spbmj2[c(1:(m0))]/(n0*n0)
    ysum2n[c(2:(m0+1))] <- spbmn2[c(1:(m0))]/(n0*n0)
    ysdmaj <- ifelse(jknmaj > 1,sqrt(((jknmaj*ysum2j)-(ybarmj^2))/(jknmaj*(jknmaj-1))),0)
    ysdmin <- ifelse(jknmin > 1,sqrt(((jknmin*ysum2n)-(ybarmn^2))/(jknmin*(jknmin-1))),0)
    ybarmj <-ifelse(jknmaj >0, ybarmj/jknmaj, ybarmj)
    ybarmn <-ifelse(jknmin >0,  ybarmn/jknmin, ybarmn)
    adjmj <- ybarmj
    adjmn <- ybarmn
    rj0maj = j0maj
    obprmj <- jknmaj/rj0maj
    rj0min = j0min
    obprmn = jknmin/(rj0min)
    obprmj[is.na(obprmj)] <-0
    obprmn[is.na(obprmn)] <-0
    iclind <- rep(1,(m0+1))
    iclind[1] <- 0
    iclind[m0+1] <- 0
    iclind <- ifelse(jknmaj < (minc),0, iclind)
    iclind <- ifelse(jknmin < (minc),0, iclind)
    lowcts = round(cusr*(m0),0)
    iclind[c(1:(lowcts+1))] <- 0
    ipelimf = 0
    ipelimr = 0
    ntsc = 0
    ntsc <- sum(iclind,na.rm = T)
    ipelimr <- sum(jknmaj[!iclind],na.rm = T)
    ipelimf <- sum(jknmin[!iclind],na.rm = T)
    pcel = (ntsc)/(m0-1)
    pelimr = (ipelimr)/(j0maj)
    pelimf = (ipelimf)/(j0min)
    if(ntsc == 0 | pelimr == j0maj | pelimf == j0min){
      return("No usable data in the datasets")
    }
    for(k in 1:(m0+2)){
      scores[k] = (k-1)
    }
    sss <- sctat(obprmj,scores,m0+1)
    vsavmj <- sss[1]
    vsvrmj <- sss[2]
    vsvrmj <- vsvrmj^2
    sss <- sctat(obprmn,scores,m0+1)
    vsavmn <- sss[1]
    vsvrmn <- sss[2]
    vsvrmn <- vsvrmn^2
    estauj <- estaun <- c(0:(m0+1))/m0
    for(i in 1:m0){
      k0=ivalid[i]
      itdfwj[k0] = (itdfwj[k0]-cusr)/(1-cusr)
      if (itdfwj[k0] < 0){
        itdfwj[k0] = 0
      }
      itdfwn[k0] = (itdfwn[k0]-cusr)/(1-cusr)
      if (itdfwn[k0] < 0){
        itdfwn[k0] = 0
      }
    }
    kk <- c(1:(m0+1))
    estaun <- estauj <- (kk-1)/m0
    kr20pj <- (m0/(m0-1))*(1-((sum(itdfwj[ivalid]*(1-itdfwj[ivalid]),na.rm = T))/vsvrmj))
    kr20pn <- (m0/(m0-1))*(1-((sum(itdfwn[ivalid]*(1-itdfwn[ivalid]),na.rm = T))/vsvrmn))
    kk <- c(1:(ni))
    estauj <- (vsavmj+kr20pj*((kk-1)-vsavmj))/m0
    estaun <- (vsavmn+kr20pn*((kk-1)-vsavmn))/m0
    estauj[is.na(estauj)] <-0
    estaun[is.na(estaun)] <-0
    adjmj  <- shadjm(m=1,m0,adjmj,ybarmj,iclind,estauj,estaun,iflg)
    if(length(ybarmj) == 1){
      stop("There was an error computing the adjusted mean for the reference group")
    }
    adjmn <- shadjm(m=2,m0,adjmn,ybarmn,iclind,estauj,estaun,iflg)
    if(length(ybarmn) == 1){
      stop("There was an error computing the adjusted mean for the reference group")
    }
    xnk <- adjmj - adjmn
    xnk[!iclind] <- 0
    dhatnk <- ((1/jknmaj)*ysdmaj^2)+((1/jknmin)*(ysdmin)^2)
    dhatnk[!iclind] <- 0
    yyy <- bver2(xnk,jknmaj,jknmin,dhatnk,iclind,m0,idw,wgt,iflg)
    se <- yyy[[1]][2]
    buni <- yyy[[1]][1]
    r1 <- yyy[[1]][3]
    wgt <- yyy[[2]]
    se = se*n0
    buni = buni*n0
    if(r1 < 0){
      p <- (pnorm(r1,0,1))*2
    }else{
      p <- (1-pnorm(r1,0,1))*2
    }
    sendback <- list()
    sendback[[1]] <-c(buni,se,r1,p)
    sendback[[2]] <- dhatnk
    sendback[[3]] <- wgt
    sendback[[4]] <- xnk
    return(sendback)
  }
  bver2 <- function(xnk,jknmaj,jknmin,dhatnk,iclind,n,idw,wgt,iflag){
    wgt <- sd_denom <- rep(0,n+1)
    rn = 0.0
    rd = 0.0
    cntmaj=0.0
    cntmin=0.0
    for(k in 1:n){
      if (iclind[k] == 1){
        cntmaj = cntmaj+jknmaj[k]
        cntmin = cntmin+jknmin[k]
      }
    }
    for(k in 1:n){
      if(iclind[k] == 1){
        if(idw == 0){
          wgt[k] = (jknmaj[k] + jknmin[k])/(cntmaj+cntmin)
        }
        if(idw == 1){
          wgt[k] = jknmin[k]/cntmin
        }
      }
      rn = rn + xnk[k]*wgt[k]
      rd = rd + (wgt[k]^2.0)*dhatnk[k]
    }
    buni = rn
    if(rd == 0){
      iflag=6
      se = 0
      r1 = -99
      return(c(buni,se,r1))
    }else{
      se = rd^0.5
      r1 = rn /se
    }
    xx123 <- list()
    xx123[[1]] <-c(buni,se,r1)
    xx123[[2]] <- wgt
    return(xx123)
  }
  shadjm <- function(m,n0,ybardt,sercor,iclind,estauj,
                     estaun,iflg){
    ivec <- NULL
    taumid <- c(rep(0,n0+1))
    kount = 0
    for(k in 1:(n0+1)){
      if(iclind[k] == 1){
        kount = kount + 1
        ivec[kount] = k
      }
    }
    for(i in 1:n0+1){
      if(iclind[i] == 1){
        l = i
        break
      }

    }
    for(i in rev(1:(n0+1))){
      if(iclind[i] == 1){
        r = i
        break
      }
    }
    if(r <= l){
      iflg <- -99
      return(iflg)
    }
    taumid[1] = max(estauj[l],estaun[l])
    taumid[kount] = min(estauj[r],estaun[r])
    for(kc in (2):(kount-1)){
      taumid[kc] = (estauj[ivec[[kc]]]+estaun[ivec[[kc]]])/2
    }
    if(m == 1){etprm1 <- rep(0,n0+1)}
    if(m == 2){etprm2 <- rep(0,n0+1)}
    if (m == 1){
      for(kc in (2):(kount-1)){
        etprm1[kc] = (sercor[ivec[kc+1]] - sercor[ivec[kc-1]]) /(estauj[ivec[kc+1]] - estauj[ivec[kc-1]] )
      }
    }
    if(m == 2){
      for(kc in (2):(kount-1)){
        etprm2[kc] = (sercor[ivec[kc+1]] - sercor[ivec[kc-1]]) /(estaun[ivec[kc+1]] - estaun[ivec[kc-1]] )
      }
    }
    estfix1 <- NULL
    estfix <- 0
    if (m == 1){
      for(kc in (2):(kount-1)){
        estfix = (estauj[ivec[kc]] - taumid[kc])*etprm1[kc]
        ybardt[ivec[kc]] = ybardt[ivec[kc]] - estfix
        estfix1 <- rbind(estfix1,estfix)
        if(ybardt[ivec[kc]] < 0){ ybardt[ivec[kc]]=0 }
        if(ybardt[ivec[kc]] > 1){ ybardt[ivec[kc]]=1 }
      }
    }else{
      if(m == 2){
        for(kc in  (2):(kount-1)){
          estfix = (estaun[ivec[kc]] - taumid[kc])*etprm2[kc]
          ybardt[ivec[kc]] = ybardt[ivec[kc]] - estfix
          if(ybardt[ivec[kc]] < 0){ ybardt[ivec[kc]]=0 }
          if(ybardt[ivec[kc]] > 1){ ybardt[ivec[kc]]=1 }
        }
      }
    }
    if(m == 1){
      ybardt[l] = ethat1(taumid[1],estauj,sercor,n0,l,r,type = 1)
      ybardt[r] = ethat1(taumid[kount],estauj,sercor,n0,l,r,type =2)
    }else{
      if(m == 2){
        ybardt[l] = ethat1(taumid[1],estaun,sercor,n0,l,r,type = 1)
        ybardt[r] = ethat1(taumid[kount],estaun,sercor,n0,l,r,type =2)
      }
    }
    return(ybardt)
  }
  ethat1 <- function( tau,estauj,ybarmj,n0,l,r,type ){
    if (tau <= estauj[l]){
      ethat1 = ybarmj[l]
      return(ethat1)
    }else{
      if (estauj[r] <= tau){
        ethat1 = ybarmj[r]
        return(ethat1)
      }
    }
    if(type == 1){
      kfound = -1
      for(k in (l):(r-1)){
        if (estauj[k] <= tau & tau < estauj[k+1]){
          kfound = k
          alpha = (tau - estauj[kfound]) /( estauj[kfound+1] - estauj[kfound] )
          ethat1 = (alpha)*ybarmj[kfound+1] + (1-alpha)*ybarmj[kfound]
          return(ethat1)
        }
      }
    }
    else{
      if(type == 2){
        kfound = -1
        for(k in (l+1):r){
          if (estauj[k-1] < tau & tau <= estauj[k]){
            kfound = k
            alpha = (tau - estauj[kfound-1]) /(estauj[kfound] - estauj[kfound-1] )
            ethat1 = (1-alpha)*ybarmj[kfound-1] + (alpha)*ybarmj[kfound]
            return(ethat1)
          }
        }
      }
    }
  }
  sctat <- function(rf,sc,nsin,rm,sd,skew,rkurt){
    epsrf = 0.1e-8
    rm = 0.0e0
    ssq = 0.0e0
    scube = 0.0e0
    s4th = 0.0e0
    for(i in 1:nsin){
      if(epsrf <= rf[i]){
        rm = rm + sc[i]*rf[i]
      }
    }
    for(i in 1:nsin){
      if(epsrf <= rf[i]){
        t4 = sc[i] - rm
        t1 = t4*rf[i]
        t2 = t4*t1
        t3 = t4*t2
        ssq = ssq + t2
        scube = scube + t3
        s4th = s4th + t4*t3
      }
    }
    var = ssq
    sd = sqrt(var)
    skew = scube
    rkurt = s4th
    t1 = var*sd
    if(t1 == 0){
      skew = 0
      rkurt = 0
    }else{
      skew = skew/t1
      t1 = t1*sd
      rkurt = rkurt/t1
      rkurt = rkurt - 3.0
    }
    return(c(rm,sd))
  }
  itdifj <- itdifn <- NULL
  nitem = ncol(data_ref)
  xnitem <- nitem
  ni = nitem
  n0 <- netsum <- length(suspect_items)
  isusp <- suspect_items
  nthsum <- length(matching_items)
  m0 = nthsum
  ivalid <- matching_items
  iden <-  idifj <- idifn <-ixpsum <- sad<-  NULL
  iflag = 0
  for(k in 1:nitem){
    idifj[k] = 0
    idifn[k] = 0
    ixpsum[k] = 0
    iden[k] <- 2
  }
  iden[isusp] <- 0

  iden[ivalid] <-1
  jknmaj <- rep(0,nitem+1)
  jknmin <- rep(0,nitem+1)
  ixsum = 0
  ix2sum = 0
  iyjsum = 0
  ixjsum = 0
  ix2jsum = 0
  idifj <- idifn <-
    xmaj <- ymaj <- xmin <- ymin <-  u1 <- rep(0,xnitem)
  for(j in 1:nrow(data_ref)){
    xmaj[j] <- 0
    ymaj[j] <- 0
    xr = 0
    u <- data_ref[j,]
    for(ii in 1:nitem){
      u1[ii] <- u[[ii]]
    }
    for(k in 1:nitem){
      if(u1[k] != 1){
        u1[k] <- 0
      }

      if(iden[k] == 1){
        xmaj[j] = xmaj[j] + u1[k]

      }
      if(iden[k] == 0){
        ymaj[j] = ymaj[j] + u1[k]
      }
      if((iden[k] != 0 & iden[k] != 1)){
        xr = xr + u1[k]
      }

      idifj[k] = idifj[k] + u1[k]
    }
    xx <- xmaj[[j]]
    jknmaj[xx]=jknmaj[xx]+1
    iyjsum = iyjsum + ymaj[j]
    ixjsum = ixjsum + xmaj[j]
    ix2jsum = ix2jsum + xmaj[j]*xmaj[j]
    ix = xmaj[j] + ymaj[j] + xr
    ixsum = ixsum + ix
    ix2sum = ix2sum + ix*ix
    for(k in 1:nitem){
      if(u1[k] == 1){
        ixpsum[k] = ixpsum[k] + ix
      }
    }
  }
  jr = nrow(data_ref)
  xjr = jr
  xbr = ixsum/(xjr)
  sxr = sqrt(ix2sum/xjr - xbr*xbr)
  ix2sumt = ix2sum
  ixsum <- ix2sum <- iynsum <- ixnsum <- ix2nsum <- 0
  for(j in 1:nrow(data_foc)){
    xmin[j] <- 0
    ymin[j] <- 0
    xr = 0
    u <- data_foc[j,]
    for(ii in 1:nitem){
      u1[ii] <- u[[ii]]
    }
    for(k in 1:nitem){
      if(iden[k] == 1){
        xmin[j] = xmin[j] + u1[k]
      }
      if(iden[k] == 0){
        ymin[j] = ymin[j] + u1[k]
      }
      if(!(iden[k] == 0 | iden[k] == 1)){
        xr = xr + u1[k]
      }

      idifn[k] = idifn[k] + u1[k]
    }
    xx <- xmin[[j]]
    jknmin[xx]=jknmin[xx]+1
    iynsum = iynsum + ymin[j]
    ixnsum = ixnsum + xmin[j]
    ix2nsum = ix2nsum + xmin[j]*xmin[j]
    ix = xmin[j] + ymin[j] + xr
    ixsum = ixsum + ix
    ix2sum = ix2sum + ix*ix
    for(k in 1:nitem){
      if(u1[k] == 1){
        ixpsum[k] = ixpsum[k] + ix
      }
    }
  }
  jf = nrow(data_foc)
  xjf = jf
  xbf = (ixsum)/(xjf)
  sxf = sqrt((ix2sum/xjf) - (xbf*xbf))
  ix2sumt = ix2sumt + ix2sum
  pbsflag = 0
  pbs <- rep(0,nitem)
  fnex = xjr + xjf
  ad = xbr - xbf
  sdavg = sqrt((xjr*sxr*sxr + xjf*sxf*sxf)/(fnex))
  if(sdavg > 0){
    sad = ad/sdavg
  }
  sad1 <- sad
  term2 <- 0
  xbar = ((xjr)*xbr + (xjf)*xbf)/(fnex)
  term1 = ((xjr)*xbr + (xjf)*xbf)*((xjr)*xbr + (xjf)*xbf)
  sdx = sqrt((ix2sumt - (term1/(fnex)))/(fnex))
  for(k in 1:nitem){
    iterm4 <- idifn[k] +idifj[k]
    iterm3 <- jr + jf +2 - idifj[k] - idifn[k]
    if(iterm3 > 0 & iterm4 > 0){
      iterm2 = iterm4/iterm3
      if(sdx < 0){
        pbsflag = 1;
        pbs[k] = -99
      }
      if(sdx >= 0){
        xbp = ixpsum[k]/iterm4
        pbs[k] = (xbp-xbar)*(sqrt(iterm2))/sdx
      }
    }
    if(iterm3 <= 0 | iterm4 <= 0){
      pbsflag = 1;
      pbs[k] = -99
    }
  }
  pp <- (idifj+idifn)/fnex
  xbr = (ixjsum)/(xjr)
  varr = (ix2jsum)/(xjr) - xbr*xbr
  sdr = sqrt(varr)
  xbf = (ixnsum)/(xjf)
  varf = (ix2nsum)/(xjf) - xbf*xbf
  sdf = sqrt(varf)
  ad = xbr - xbf
  sdavg = sqrt(((xjr)*varr + (xjf)*varf)/(fnex))
  if(sdavg > 0) {sad = ad/sdavg }
  sad2 <- sad
  for(k in 1:nitem){
    itdifj[k] = (idifj[k])/(xjr)
    itdifn[k] = (idifn[k])/(xjf)
  }
  jknmaj11 <- table(xmaj)
  jknmin11 <-table(xmin)
  jknmaj1 <- jknmin1 <-rep(0,ni+1)
  jknmaj1[as.integer(names(jknmaj11))+1] <- jknmaj11
  jknmin1[as.integer(names(jknmin11))+1] <- jknmin11

  output <- sibuni(data_ref=data_ref,data_foc=data_foc,xmaj,ymaj,jr,xmin,ymin,
                   jf,itdifj,itdifn,ivalid,
                   netsum,nthsum,nitem,minc,cusr,iflag,idw,
                   buni,r1,se,wgt,adjmj,adjmn,xnk,ybarmj,
                   ybarmn,jknmaj1,jknmin1,kr20n,kr20j,ehr,ehf)

  dhatnk <- unlist(output[[2]], recursive = TRUE, use.names = TRUE)
  wgt <- unlist(output[[3]], recursive = TRUE, use.names = TRUE)
  xnk <- unlist(output[[4]], recursive = TRUE, use.names = TRUE)
  dp <- rep(0,(nthsum+1))
  vardp <- rep(0,(nthsum+1))
  weight <-rep(0,(nthsum+1))
  for(k in 1:(nthsum+1)){
    dp[k] = xnk[k]*wgt[k]
    vardp[k] = (wgt[k]*wgt[k])*dhatnk[k]
  }
  vardp[is.na(vardp)] <- 0
  x <- rep(0,1000)
  y <- rep(0,1000)
  for(i in 1:(nthsum+1)){
    if(jknmaj1[i] >  jknmin1[i]){
      weight[i] <- jknmaj1[i]
    }
    if(jknmaj1[i] <=  jknmin1[i]){
      weight[i] = jknmin1[i]
    }
    x[i] = i
    y[i] = xnk[i]
  }
  suma = 0.0
  sumb = 0.0
  sumc = 0.0
  sumd = 0.0
  sume = 0.0
  killr = ifelse((0.01*jr > 0), ifelse((0.01*jr - floor(0.01*jr)) >= .5,
                                     ceiling(0.01*jr),floor(0.01*jr)),
                 ifelse(0.01*jr == 0, 0,
                        ifelse((0.01*jr - floor(0.01*jr)) >= .5,floor(0.01*jr),ceiling(0.01*jr)
                 )))
  killf = ifelse((0.01*jf > 0), ifelse((0.01*jf - floor(0.01*jf)) >= .5,ceiling(0.01*jf),floor(0.01*jf)),
                 ifelse(0.01*jf == 0, 0,
                        ifelse((0.01*jf - floor(0.01*jf)) >= .5,floor(0.01*jf),ceiling(0.01*jf)
                 )))
  for(i in 1:(nthsum+1)){
    if(jknmaj1[i] > killr){
      if(jknmin1[i] > killf){
        suma = suma + weight[i]
        sumb = sumb +  x[i] * weight[i]
        sumc = sumc + x[i]*x[i] * weight[i]
        sumd = sumd + y[i] * weight[i]
        sume = sume + x[i] * y[i] * weight[i]
      }
    }
  }
  beta0=((sumc*sumd)-(sumb*sume))/((suma*sumc)-(sumb^2))
  beta1=((suma*sume)-(sumb*sumd))/((suma*sumc)-(sumb^2))
  if(-beta0/beta1 > 0){
    if(-beta0/beta1 - floor(-beta0/beta1) >= 0.5){
      kcross = ceiling(-beta0/beta1)
    }
    if(-beta0/beta1-floor(-beta0/beta1) < 0.5){
      kcross = floor(-beta0/beta1)
    }
  }
  if(-beta0/beta1 == 0){
    kcross = 0
  }
  if(-beta0/beta1 < 0){
    if(-beta0/beta1 - floor(-beta0/beta1) >= 0.5){
      kcross = floor(-beta0/beta1)
    }
    if(-beta0/beta1-floor(-beta0/beta1) < 0.5){
      kcross = ceiling(-beta0/beta1)
    }
  }
  if (kcross <= 0){ kc=1}
  if (kcross >= (nthsum+1)){ kc=(nthsum+1)}
  if ((kcross > 0) & (kcross < (nthsum+1))){ kc=kcross}
  bb1=0
  bb2=0
  kc3 = kc
  bc32 = 0
  smvbc32 = 0
  smvbc32u = 0
  smvbc32b = 0
  bc32b <- bc32u <-0
  if (kc3 == 1){
    for(k in 1:(nthsum+1)){
      bc32 = bc32+dp[k]
      smvbc32 = smvbc32 + vardp[k]
      bc32b = bc32b+dp[k]
      smvbc32b = smvbc32b + vardp[k]
    }
    bb1 = 0
    bb2 = bc32
  }
  if(kc3 == (nthsum+1)){
    for(k in 1:(nthsum+1)){
      bc32 = bc32 - dp[k]
      smvbc32 = smvbc32 + vardp[k]
      bc32u = bc32u-dp[k]
      smvbc32u = smvbc32u + vardp[k]
    }
    bb1=bc32
    bb2=0
  }
  if((kc3 > 1) & (kc3 < (nthsum+1))){
    for(k in 1:(kc3-1)){
      bc32 = bc32 - dp[k]
      bb1=bb1-dp[k]
      smvbc32 = smvbc32 + vardp[k]
      bc32b = bc32b-dp[k]
      smvbc32b = smvbc32b + vardp[k]
    }
    for(k in (kc3+1):(nthsum+1)){
      bc32 = bc32 + dp[k]
      bb2=bb2+dp[k]
      smvbc32 = smvbc32 + vardp[k]
    }
  }
  if((kc3 > 1) & (kc3 < (nthsum+1))){
  for(k in (kc3):(nthsum+1)){
    bc32u = bc32u+dp[k]
    smvbc32u = smvbc32u + vardp[k]
    }
  }
  dfs <- 0
  if(smvbc32u  >0){
    dfs <- dfs +1
  }
  if(smvbc32b > 0){
    dfs <- dfs + 1
  }
  x1 <-  ifelse(smvbc32b > 0, (bc32b^2/smvbc32b),0  )
  x2 <-  ifelse(smvbc32u > 0, (bc32u^2/smvbc32u),0  )
  BCROSS2 <- x1+x2
  p_bcross2 <- pchisq(BCROSS2, dfs, lower.tail = FALSE)
  bcro32 = bc32 / ( smvbc32^0.5 )
  se32 = (smvbc32^0.5)*(netsum)
  betacroC = (bc32b+bc32u)*(netsum)
  betacro = bc32*(netsum)
  bb1=bb1*(netsum)
  bb2=bb2*(netsum)
  prop1=abs(bb1)/(abs(bb1)+abs(bb2))
  prop2=abs(bb2)/(abs(bb1)+abs(bb2))
  nprop1=round(prop1*100.0,0)
  nprop2=round(prop2*100.0,0)
  idg1=' '
  idg2=' '
  if(bb1 < 0){ idg1="R"}
  if(bb1 > 0){ idg1="F"}
  if(bb2 < 0){ idg2="F"}
  if(bb2 > 0){ idg2="R"}
  nwant = 1000
  nli=nthsum+1
  newu1 <- matrix(runif((nli*nwant),0,1),nrow = nli,byrow = T)
  crostat <- rep(9999,nwant)
  for(j in 1:nwant){
    newu <- ifelse(newu1[,j] >= .5,1,0)
    for(k in 1:(nthsum+1)){
      if(newu[k] == 0 ){
        y[k] = -1*y[k]
        dp[k] = -1*dp[k]
      }
      if(newu[k] == 1){
        y[k] = y[k]
        dp[k] = dp[k]
      }
    }
    suma <- sumb <- sumc <- sumd <- sume <- 0
    for(i in 1:(nthsum+1)){
      if ((jknmaj1[i] > killr) & (jknmin1[i] > killf)) {
        suma = suma + weight[i]
        sumb = sumb +  x[i] * weight[i]
        sumc = sumc + x[i]*x[i]*weight[i]
        sumd = sumd + y[i] * weight[i]
        sume = sume + x[i] * y[i] * weight[i]
      }
    }
    beta0new=(sumc*sumd-sumb*sume)/(suma*sumc-sumb^2.0)
    beta1new=(suma*sume-sumb*sumd)/(suma*sumc-sumb^2.0)
    if((-beta0new/beta1new) > 0){
      if(((-beta0new/beta1new) - floor((-beta0new/beta1new))) >= .5){
        knew = ceiling((-beta0new/beta1new))}
      if(((-beta0new/beta1new) - floor((-beta0new/beta1new))) < .5){
        knew = floor((-beta0new/beta1new))}
    }
    if((-beta0new/beta1new) == 0){ knew = 0}
    if((-beta0new/beta1new) < 0){
      if(((-beta0new/beta1new) - floor((-beta0new/beta1new))) >= .5){
        knew = floor((-beta0new/beta1new))}
      if(((-beta0new/beta1new) - floor((-beta0new/beta1new))) < .5){
        knew = ceiling((-beta0new/beta1new))}
    }
    if (knew <= 0){ knew=1}
    if (knew >= (nthsum+1)){ knew=(nthsum+1)}
    if ((knew > 0) & (knew < (nthsum+1))){ knew=knew}
    kc3 = knew
    bc32 = 0
    smvbc32 = 0
    if (kc3 == 1){
      for(k in 1:(nthsum+1)){
        bc32 = bc32+dp[k]
        smvbc32 = smvbc32 + vardp[k]
      }
    }
    if(kc3 == (nthsum+1)){
      for(k in 1:(nthsum+1)){
        bc32 = bc32 - dp[k]
        smvbc32 = smvbc32 + vardp[k]
      }

    }
    if((kc3 > 1) & (kc3 < (nthsum+1))){
      for(k in 1:(kc3-1)){
        bc32 = bc32 - dp[k]
        smvbc32 = smvbc32 + vardp[k]
        
      }
      for(k in (kc3+1):(nthsum+1)){
        bc32 = bc32 + dp[k]
        smvbc32 = smvbc32 + vardp[k]
      }
    }
    bcronew = bc32*(netsum)
    crostat[j] = bc32 /((smvbc32^0.5))
  }
  crostat <- sort(crostat)
  ir <-sum(abs(crostat) >= abs(bcro32),na.rm = TRUE)
  pcro = (ir)/(nwant)
  out <- data.frame(
    matrix(c(suspect_items = paste(suspect_items,collapse = " "),
                    beta = round(output[[1]][1],3),sigma_uni = round(output[[1]][2],3),  z = round(output[[1]][3],3),df = NA,
                    p_value = round(output[[1]][4],3),
             suspect_items = paste(suspect_items,collapse = " "),
                     beta_cross = round(betacro,3), cross_se = round(se32,3),
                    bcross = NA, df = NA, p_cross = round(pcro,3),
  suspect_items = paste(suspect_items,collapse = " "),
                     beta_cross = round(betacroC,3),  se = NA,
                     bcross2 = round(BCROSS2,3), df = dfs, bcross2_p = round(p_bcross2,3)), nrow = 3, byrow = T))

  colnames(out) <- c("Suspect Item","Beta","SE","Z/X2","df","p")
  rownames(out) <- c("Uniform Original","Crossing (Li & Stout)","Crossing (Chalmers)")
  return(out)
}
