pro puc_numdens
 
  common puc_common

; -------- get a model for ISRF

  isrfpath = '$PUC_PATH/fits/MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits'
  isrf = mrdfits(isrfpath, 0, h) * 1d-6 ; [MeV/cm^3]
  crval1 = sxpar(h, 'CRVAL1')
  cdelt1 = sxpar(h, 'CDELT1')
  risrf = dindgen(n_elements(isrf[*,0,0,0])) * cdelt1 + crval1
  crval2 = sxpar(h, 'CRVAL2')
  cdelt2 = sxpar(h, 'CDELT2')
  zisrf = dindgen(n_elements(isrf[0,*,0,0])) * cdelt2 + crval2
  crval3 = sxpar(h, 'CRVAL3')
  cdelt3 = sxpar(h, 'CDELT3')
  llisrf = dindgen(n_elements(isrf[0,0,*,0])) * cdelt3 + crval3
  lisrf = 10.d^llisrf * 1d-4 ; [cm]
  enisrf = planck * lightspeed / lisrf ; [MeV]
  enisrf = reverse(enisrf)
  nenisrf = n_elements(enisrf)
  nrisrf = n_elements(risrf)
  nzisrf = n_elements(zisrf)
  inr = max(r) * 0.1
  inz = max(z) * 0.1
  irin = where(risrf lt inr, ninr)
  izin = where(abs(zisrf) lt inz, ninz)
  ncellin = ninr + ninz
  ncell =  n_elements(isrf[*,*,0,0])


; -------- get average and spatially-dependent photon number densities

  nphopuv = dblarr(nenisrf,nrisrf,nzisrf)
  nphmfir = dblarr(nenisrf,nrisrf,nzisrf)
  nphopuvav = dblarr(nenisrf)
  nphmfirav = dblarr(nenisrf)
  nphopuvavin = dblarr(nenisrf)
  nphmfiravin = dblarr(nenisrf)
  for a=0, nenisrf-1 do begin
    nphopuv[a,*,*] = isrf[*,*,a,0] / enisrf[a]^2 ; [1/cm^3/MeV]
    nphmfir[a,*,*] = isrf[*,*,a,1] / enisrf[a]^2 ; [1/cm^3/MeV]
    nphopuvav[a] = total(isrf[*,*,a,0]) / ncell / enisrf[a]^2 ; [1/cm^3/MeV]
    nphmfirav[a] = total(isrf[*,*,a,1]) / ncell / enisrf[a]^2 ; [1/cm^3/MeV]
    nphopuvavin[a] = total(isrf[irin,izin,a,0]) / ncellin / enisrf[a]^2 ; [1/cm^3/MeV]
    nphmfiravin[a] = total(isrf[irin,izin,a,1]) / ncellin / enisrf[a]^2 ; [1/cm^3/MeV]
  endfor
    uisrfav = enisrf^2 * (nphopuvav + nphmfirav)
    uisrfavin = enisrf^2 * (nphopuvavin + nphmfiravin)
end
