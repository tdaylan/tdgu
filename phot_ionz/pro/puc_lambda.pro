pro puc_lambda

  common puc_common

  nptc = 100

  minrptc = 0.
  maxrptc = 1.

  minsptc = -0.5
  maxsptc = 0.5

  minzptc = -1.
  maxzptc = 1.
 
  rptc = dindgen(nptc) / (nptc - 1) * (maxrptc - minrptc) + minrptc
  sptc = dindgen(nptc) / (nptc - 1) * (maxsptc - minsptc) + minsptc
  zptc = dindgen(nptc) / (nptc - 1) * (maxzptc - minzptc) + minzptc
 
  lptc = dblarr(nptc,nptc,nptc) 
  ;for i=0, nptc-1 do $
  ;  lptc[i,*,*] = (2. / 3. / sqrt(2. * !pi) * rptc[i]^4 * exp(-rptc[i]^2 / 2.)) * (1.5 - 6 * sptc^2) # (1. / sqrt(2. * !pi) * exp(-zptc^2 / 2.))
  for i=0, nptc-1 do $
    lptc[i,*,*] = (rptc[i]^4 * exp(-rptc[i]^2 / 2.)) * (1.5 - 6 * sptc^2) # exp(-zptc^2 / 2.)
  lptc /= max(lptc)

  lambda = dblarr(npatch, 3)
  for i=0, npatch-1 do begin
    mindif = min(abs(lptc - randomu(seed)), iptc)
    iptc = array_indices(lptc, iptc) 
    rptc_ = rptc[iptc[0]]
    sptc_ = sptc[iptc[1]]
    zptc_ = zptc[iptc[2]]

    lambda[i,0] = (zptc_ + rptc_ * (sptc_ + sqrt(3. * (1. - sptc_^2))) / sqrt(5)) / 3.
    lambda[i,1] = (zptc_ - 2. * rptc_ * sptc_ / sqrt(5)) / 3.
    lambda[i,2] = (zptc_ + rptc_ * (sptc_ - sqrt(3. * (1. - sptc_^2))) / sqrt(5)) / 3.

  endfor

  ;print, rptc, sptc, zptc
  ;print, lptc
  ;print, transpose(lambda)

end
