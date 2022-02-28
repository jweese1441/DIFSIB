#' POLYSIBTEST
#'
#' This function assesses uniform DIF for polytomous DATA using the POLYSIBTEST procedure (Chang et al., 1996).
#' The function outputs the uniform DIF magnitudes as well as the significance testing.
#'
#' @param data_ref  The dataset that contains the reference group participants
#' @param data_foc  The dataset that contains the focal group participants
#' @param minc  The minimum number of individuals in a give total test score category. Initialized to 2 per Shealy and Stout (1993) recommendation
#' @param cusr  The guess correction factor. This is initialized to 0, but can be changed based on the items used. The cusr value is the average of guess correction across the items in the dataset.
#' @param idw  The weighting factor. This is initialized to use the reference and focal group weighting recommended by Shealy and Stout (1993)
#' @param suspect_items  the item(s) to be assess for DIF. Listing one item number will assess for DIF on that item, multiple items (e.g. c(2,3,8)) will assess DIF for the bundle of items.
#' @param matching_items  the list of items to match on.
#' @param listwise  initialized to 1, this inidicates that the procedure will not delete any rows with missing data, changing this to 2 will listwise delete.
#' @param nch  a vector indicating the number of categories for each item. Example if there are 5 items and each has 3 catgories nch <- c(3,3,3,3,3). If there are 5 items and items 1 - 3 have 4 categories and items 4 and 5 have 2 then nch <- (4,4,4,2,2)
#' @aliases psib
#' @export psib
#' @references Chang, H. H., Mazzeo, J. & Roussos, L. (1996). DIF for Polytomously Scored Items: An Adaptation of the SIBTEST Procedure. Journal of Educational Measurement, 33, 333-353.
#' @references DIF-Pack (2021) Measured Progress. https://psychometrics.onlinehelp.measuredprogress.org/tools/dif/
#' @references Shealy, R. & Stout, W. (1993). A model-based standardization approach that separates true bias/DIF from group ability differences and detect test bias/DTF as well as item bias/DIF. Psychometrika, 58, 159-194.
#' @references Weese, J. D., Turner, R. C., Liang, X., Ames, A., & Crawford, B.  (accepted). Implementing a Standardized Effect Size in the POLYSIBTEST Procedure
#' @examples
#'
#' \dontrun{
#'
#' #perform POLYSIBTEST with one suspect item. Sample data has 20 items with 5 categories each.
#' nch <- rep(5,20)
#' psib(data_ref = data_refP, data_foc = data_focP, 
#' minc = 2,suspect_items = c(20), matching_items = c(1:19),nch = nch)
#'
#' #perform POLYSIBTEST with a bundle of suspect items 
#'
#' psib(data_ref = data_refP, data_foc = data_focP, 
#' minc = 2, suspect_items = c(16,17,18,19,20), matching_items = c(1:15),nch = nch)
#'
#' }
psib <- function(data_ref, data_foc,minc=2,cusr = 0,idw = 0,suspect_items, matching_items,listwise = 1,nch){
  ############################################################################################################1
  ############################################################################################################1
  # NOTE: THIS CODE IS ADAPTED FROM THE ORIGINAL FORTRAN CODE THAT POWERED THE SIBTEST V1.7 (2005 SOFTWARE
  # RELEASE DATE). THAT CODE HAS BEEN RELEASED BY MEASURED PROGRESS AND IS NOW CONSIDERED OPEN SOURCED.
  # YOU CAN FIND THAT CODE AT THE  FOLLOWING WEBSITE
  # https://psychometrics.onlinehelp.measuredprogress.org/tools/dif/.
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
    return("POLYSIBTEST requires that both reference and focal group datasets do not have missing values. Please select TRUE to remove missing values to listwise remove missing values")}
  if(ncol(data_ref) != ncol(data_foc)){
    return("Verify that your data is correct. The number of columns in the reference group data does not equal the number of columns in the focal group data.")}
  if(any(suspect_items %in% matching_items)){
    return("POLYSIBTEST cannot run when when suspect items and matching subtest items contain the same the same items. Please check your suspect and matching subtest item selections")}
  psibuni <- function(data_ref,data_foc,xmaj,ymaj,j0maj,xmin,ymin,j0min,itdifj,
                    itdifn,n0,m0,ni,minc,cusr,iflg,idw,
                    fnamer,fnamef,hvscore,hscore,hsscore,ndiv,idiv,ipv,jknmaj,jknmin){
  estaujn <- estaunn <- estauj <-estaun <-
    adjmj <- adjmn <- obprmj<-obprmn<-scores<-iclind <- NULL
  itdfwj = itdifj
  itdfwn = itdifn
  ybarmj <- ybarmn <-ysum2j <- ysum2n <-
  spbmaj <- spbmin <- spbmj2 <- spbmn2 <- rep(0,hvscore+1)
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
  mhsscor = hsscore/n0
  ybarmj[c(2:(hvscore+1))] <- spbmaj[c(1:(hvscore))]/n0
  ybarmn[c(2:(hvscore+1))] <- spbmin[c(1:(hvscore))]/n0
  ysum2j[c(2:(hvscore+1))] <- spbmj2[c(1:(hvscore))]/(n0*n0)
  ysum2n[c(2:(hvscore+1))] <- spbmn2[c(1:(hvscore))]/(n0*n0)
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
  iclind <- rep(1,(hvscore+1))
  iclind[1] <- 0
  iclind[hvscore+1] <- 0
  iclind <- ifelse(jknmaj < (minc),0, iclind)
  iclind <- ifelse(jknmin < (minc),0, iclind)
  lowcts = round(cusr*(ndiv),0)
  iclind[c(1:(lowcts+1))] <- 0
  ipelimf = 0
  ipelimr = 0
  ntsc = 0
  ntsc <- sum(iclind,na.rm = T)
  ipelimr <- sum(jknmaj[!iclind],na.rm = T)
  ipelimf <- sum(jknmin[!iclind],na.rm = T)
  pcel = (ntsc)/(hvscore-1)
  pelimr = (ipelimr)/(j0maj)
  pelimf = (ipelimf)/(j0min)
  if(ntsc == 0 | pelimr == j0maj | pelimf == j0min){
    return("No usable data in the datasets")
  }
  scores <- c(0:hvscore)
  rm <- scores*obprmj
  rm[!iclind] <- 0
  rf <- obprmj
  rf[!iclind] <- 0
  vsavmj <- sum(rm)/sum(rf)
  vsvrmj <- (((scores-vsavmj)^2)*obprmj)
  vsvrmj[!iclind] <- 0
  vsvrmj <- sum(vsvrmj)/sum(rf)
  rm <- scores*obprmn
  rm[!iclind] <- 0
  rf <- obprmn
  rf[!iclind] <- 0
  vsavmn <- sum(rm)/sum(rf)
  vsvrmn <- (((scores-vsavmn)^2)*obprmn)
  vsvrmn[!iclind] <- 0
  vsvrmn <- sum(vsvrmn)/sum(rf)
  kk <- c(1:(hvscore+1))
  estaun <- estauj <- (kk-1)/hvscore
  mix = 1
  npv = m0 - ndiv
  if(ndiv > 0){
    for(i in 1:ndiv){
      k0=idiv[i]
      itdfwj[k0] = (itdfwj[k0]-cusr)/(1-cusr)
      if (itdfwj[k0] < 0){
        itdfwj[k0] = 0
      }
      itdfwn[k0] = (itdfwn[k0]-cusr)/(1-cusr)
      if (itdfwn[k0] < 0){
        itdfwn[k0] = 0
      }
    }
    divarj=sum(itdfwj*(1-itdfwj),na.rm=T)
    divarn=sum(itdfwn*(1-itdfwn),na.rm=T)
  }else{
    divarj=0
    divarn=0
  }
  if(ndiv != m0){
    xbar <- xvar <- rep(0,npv)
    for(i in 1:npv){
      iv = ipv[i]
      xbar[i] = sum(data_ref[,iv])
      xvar[i] = sum(data_ref[,iv]^2)
    }
    njj = nrow(data_ref)
    f1 = njj
    xbar <- xbar/f1
    xvar <- xvar/f1
    xvar = xvar-xbar*xbar
    vmj <- sum(xvar,na.rm = T)+divarj
    xbar <- xvar <- rep(0,npv)
    for(i in 1:npv){
      iv = ipv[i]
      xbar[i] = sum(data_foc[,iv])
      xvar[i] = sum(data_foc[,iv]^2)
    }
    njj = nrow(data_foc)
    f1 = njj
    xbar <- xbar/f1
    xvar <- xvar/f1
    xvar = xvar-xbar*xbar
    vmn <- sum(xvar,na.rm = T)+divarn
  }
  if(ndiv == m0){
    vmj = divarj
    vmn = divarn
  }
  kr20pj <- (m0/(m0-1))*(1-(vmj/vsvrmj))
  kr20pn <- (m0/(m0-1))*(1-(vmn/vsvrmn))
  kk <- c(1:(hvscore+1))
  estauj <- (vsavmj+kr20pj*((kk-1)-vsavmj))/hvscore
  estaun <- (vsavmn+kr20pn*((kk-1)-vsavmn))/hvscore
  adjmj  <- shnew(m=1,adjmj,ybarmj,iclind,estauj,estaun,iflg,hvscore,mhsscor)
  if(length(adjmj) == 1){
    stop("There was an error computing the adjusted mean for the reference group")
  }
  adjmn <- shnew(m=2,adjmn,ybarmn,iclind,estauj,estaun,iflg,hvscore,mhsscor)
  if(length(adjmj) == 1){
    stop("There was an error computing the adjusted mean for the focal group")
  }
  xnk <- adjmj - adjmn
  xnk[!iclind] <- 0
  dhatnk <- ((1/jknmaj)*ysdmaj^2)+((1/jknmin)*(ysdmin)^2)
  dhatnk[!iclind] <- 0
  yyy <- beta_calc(xnk,jknmaj,jknmin,dhatnk,iclind,hvscore)
  se <- yyy[2]
  buni <- yyy[1]
  r1 <- yyy[3]
  sd1 <- yyy[4]
  es1 <- yyy[5]
  se = se*n0
  buni = buni*n0
  if(r1 < 0){
    p <- (pnorm(r1,0,1))*2
  }else{
    p <- (1-pnorm(r1,0,1))*2
  }
  return(c(buni,se,r1,p,sd1,es1,pelimr,pelimf))
}
beta_calc <- function(xnk,jknmaj,jknmin,dhatnk,iclind,n){
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
        sd_denom[k] = (wgt[k]^2)*(1/jknmaj[k]+1/jknmin[k])
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
    warning("There was an error computing the standard error, it is 0 and therfore, no p-value can be calculated")
    iflag=6 
    se=0.0  
    r1 = -99 
    return(c(buni,se,r1,0,0))
  }else{ 
    se = rd^0.5
    if(sum(sd_denom) > 0){
    sd = se/sqrt(sum(sd_denom))
    es1 <- buni/sd  
    }
    if(sum(sd_denom) == 0){
      warning("There was an error computing the standard deviation, it is 0 and therfore, no effect size can be calculated")
      sd = 0
      es1 = 0
    }
    r1 = rn /se 
  }
  return(c(buni,se,r1,sd,es1))
}
shnew <- function(m,ybardt,sercor,iclind,estauj,
                  estaun,iflg,n0,mhsscor){
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
  for(i in rev(1:n0+1)){
    if(iclind[i] == 1){
      r = i
      break  }

  }
  if(r <= l){
    iflg <- -99
    return(iflg)

  }
  taumid[1] = max(estauj[l],estaun[l])
  taumid[kount] = min(estauj[r],estaun[r])

  for(kc in (2):(kount-1)){
    x <- ivec[kc]
    taumid[kc] = (estauj[x]+estaun[x])/2
  }
  if(m == 1){etprm1 <- rep(0,n0+1)}
  if(m == 2){etprm2 <- rep(0,n0+1)}
  if (m == 1){
    for(kc in (2):(kount)){
      etprm1[kc] = (sercor[ivec[kc+1]] - sercor[ivec[kc-1]]) /(estauj[ivec[kc+1]] - estauj[ivec[kc-1]] )
    }

  }
  if(m == 2){
    for(kc in (2):(kount)){
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
      if(ybardt[ivec[kc]] > mhsscor){ ybardt[ivec[kc]]=mhsscor }
    }
  }else{
    if(m == 2){
      for(kc in  (2):(kount-1)){
        estfix = (estaun[ivec[kc]] - taumid[kc])*etprm2[kc]
        ybardt[ivec[kc]] = ybardt[ivec[kc]] - estfix

        if(ybardt[ivec[kc]] < 0){ ybardt[ivec[kc]]=0 }
        if(ybardt[ivec[kc]] > mhsscor ){ ybardt[ivec[kc]]=mhsscor }
      }

    }
  }
  if(m == 1){
    ybardt[l] = ethat1(taumid[1],estauj,sercor,n0)
    ybardt[r] = ethat1(taumid[kount],estauj,sercor,n0)
  }else{

    if(m == 2){
      ybardt[l] = ethat1(taumid[1],estaun,sercor,n0)
      ybardt[r] = ethat1(taumid[kount],estaun,sercor,n0)
    }
  }
  return(ybardt)
}
ethat1 <- function( tau,estauj,ybarmj,n0){
  if (tau <= estauj[1]){
    ethat1 = ybarmj[1]
    return(ethat1)
  }else{
    if (estauj[n0+1] <= tau){
      ethat1 = ybarmj[n0+1]
      return(ethat1)
    }
  }
  kfound = -1
  for(k in (2):(n0)){
    if (estauj[k] <= tau & tau < estauj[k+1]){
      kfound = k
      alpha = (tau - estauj[kfound]) /( estauj[kfound+1] - estauj[kfound] )
      ethat1 = (alpha)*ybarmj[kfound+1] + (1-alpha)*ybarmj[kfound]
      return(ethat1)
    }
  }
}
  itdifj <- itdifn <- NULL
  nitem = ncol(data_ref)
  xnitem <- nitem
  ni = nitem
  netsum <- length(suspect_items)
  n0 = netsum
  isusp <- suspect_items
  nthsum <- length(matching_items)
  m0 = nthsum
  ivalid <- matching_items
  hscore = sum(nch)
  iden <-  idifj <- idifn <-ixtsum <- sad<-  NULL
  jknmin <- jknmaj <- rep(0,hscore)
  idifj <- idifn <- ixtsum <- rep(0,nitem)
  iden <- rep(2,nitem)
  iflag = 0
  iden[isusp] <- 0
  iden[ivalid] <-1
  hsscore = hscore
  idiv <- ipv <- rep(0,nthsum)
  hvscore = ndiv = npv = 0
  for(i in 1:nthsum){
    if(nch[ivalid[i]] == 1){
      ndiv = ndiv + 1
      idiv[ndiv] = ivalid[i]
    }else{
      npv = npv +1
      ipv[npv] = ivalid[i]
    }
    hvscore = hvscore + nch[ivalid[i]]
  }

  xmaj <- ymaj <- rep(0,nrow(data_ref))
  xmin <- ymin <- rep(0,nrow(data_foc))
  u1 <- NULL
  isum2 <- idifj <- idifn <- ixtsum <- rep(0,nitem)
  ixsum <- ix2sum <- iyjsum <- ixjsum <- ix2jsum <- 0
  for(j in 1:nrow(data_ref)){
    xmaj[j] <- 0
    ymaj[j] <- 0
    xr = 0
    u <- data_ref[j,]
    for(ii in 1:nitem){
      u1[ii] <- u[[ii]]
    }
    for(k in 1:nitem){
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
      isum2[k] = isum2[k] + u1[k]*u1[k]
    }
    xx <- xmaj[j]
    jknmaj[xx]=jknmaj[xx]+1
    iyjsum = iyjsum + ymaj[j]
    ixjsum = ixjsum + xmaj[j]
    ix2jsum = ix2jsum + xmaj[j]*xmaj[j]
    ix = xmaj[j] + ymaj[j] + xr
    ixsum = ixsum + ix
    ix2sum = ix2sum + ix*ix
    for(k in 1:nitem){
      ixtsum[k] = ixtsum[k] + ix*u1[k]
    }
  }
  jr = nrow(data_ref)
  xjr = jr
  xbr = ixsum/(xjr)
  sxr = sqrt(ix2sum/xjr - xbr*xbr)
  ix2sumt = ix2sum
  ixsum = 0
  ix2sum = 0
  iynsum = 0
  ixnsum = 0
  ix2nsum = 0
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
      isum2[k] = isum2[k] + u1[k]*u1[k]
    }
    xx <- xmin[j]
    jknmin[xx]=jknmin[xx]+1
    iynsum = iynsum + ymin[j]
    ixnsum = ixnsum + xmin[j]
    ix2nsum = ix2nsum + xmin[j]*xmin[j]
    ix = xmin[j] + ymin[j] + xr
    ixsum = ixsum + ix
    ix2sum = ix2sum + ix*ix
    for(k in 1:nitem){
      ixtsum[k] = ixtsum[k] + ix*u1[k]
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
  if(sdavg > 0){sad = ad/sdavg}
  sad1 <- sad
  term2 <- 0
  xbar = ((xjr)*xbr + (xjf)*xbf)/(fnex)
  term1 = ((xjr)*xbr + (xjf)*xbf)*((xjr)*xbr + (xjf)*xbf)
  sdx = sqrt((ix2sumt - (term1/(fnex)))/(fnex))
  for(k in 1:nitem){
    iterm4 <- idifn[k] +idifj[k]
    iterm3 <- hscore*(jr + jf) - idifj[k] - idifn[k]
    if(iterm3 > 0 & iterm4 > 0){
      if(sdx < 0){
        pbsflag = 1;
        pbs[k] = -99
      }
      if(sdx >= 0){
        sdk = sqrt((isum2[k]-(iterm4*iterm4)/fnex)/fnex)
        xnum = ixtsum[k]/fnex - xbar*iterm4/fnex
        pbs[k] = xnum/(sdx*sdk)
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
  itdifj <- idifj/xjr
  itdifn <- idifn/xjf
  jknmaj1 <- jknmin1 <-rep(0,(hvscore+1))
  jknmaj11 <- table(xmaj)
  jknmin11 <-table(xmin)

  jknmaj1[as.integer(names(jknmaj11))+1] <- jknmaj11
  jknmin1[as.integer(names(jknmin11))+1] <- jknmin11
  output <- psibuni(data_ref=data_ref,data_foc=data_foc,xmaj,ymaj,jr,xmin,ymin,
                    jf,itdifj,itdifn,
                    netsum,nthsum,nitem,minc,cusr,iflag,idw,
                    fnamer,fnamef,hvscore,hscore,hsscore,ndiv,idiv,ipv,jknmaj1,jknmin1)
  out <- data.frame(suspect_items = paste(suspect_items,collapse = " "), beta = output[1],sigma_uni = output[2], std = output[5], effect_size = output[6],z = output[3],
                    p_value = output[4])
  out <- cbind(out[1],round(out[2:7],3))
  colnames(out) <- c("Suspect Item(s)","Beta","SE","SD","Delta","z","p")
  return(out)
}
