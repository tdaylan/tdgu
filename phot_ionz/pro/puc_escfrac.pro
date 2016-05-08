pro puc_escfrac

  common puc_common



; -------- compute ISRF energy density
  
    puc_mw_isrf

  
; -------- get electron emissivity of MW
  
    rhomw = puc_rho(mmw, cmw)
    emismwsun =  sigmav_ * rhosun^2 ; [MeV/s/cm^3]
    emismw = dblarr(nenel,nr,nz)
    for a=0, nenel-1 do emismw[a,*,*] = sigmav_[a] * rhomw^2 ; [MeV/s/cm^3]
    lumimw = puc_lumihalo(mmw, cmw) ; [MeV/s]
  
    if not keyword_set(noplot) then puc_emismw_plot
  
  
; -------- compute the electron halo function as a function of position, injection and final electron energies
  
    puc_haloelec
    time1 = systime(1)
  
    haloelectag = '$PUC_PATH/fits/haloelec' + string(nenel,format='(I2.2)') + string(nr,format='(I2.2)') + string(nz,format='(I2.2)') + '.fits'
    if file_test(haloelectag) eq 1 then begin
      haloelec = mrdfits(haloelectag)
    endif else begin
      lambda = dblarr(nenel,nenel)
      haloelec = dblarr(nenel,nenel,nr,nz)
      for a=0, nenel-1 do begin
        for b=0, nenel-1 do begin
          if enel[a] ge enel[b] then continue
          lambda2 = qromo('puc_lambda2_integrand', enel[a], enel[b], eps=eps)
          lambda[a,b] = sqrt(lambda2)
          for i=0, nr-1 do begin
            for j=0, nz-1 do begin
              if abs(z[j]) gt lz then continue
              haloelec[a,b,i,j] = int_3d('puc_haloelec_integrand',[minr,maxr],'puc_gal_zall','puc_gal_phiall', 96)
            endfor
          endfor
        endfor
        print, 'puc_haloelec_integrand: ' + string(float(a+1)/float(nenel)*100d, format='(I3.3)') + '%'
      endfor
      mwrfits, haloelec, haloelectag, /create
    endelse
    if not keyword_set(noplot) then puc_haloelec_plot
  
    time2 = systime(1)
    print, 'puc_haloelec: ', time2 - time1, ' seconds'
  

; -------- integrate the halo function over the electron injection energies

  time1 = systime(1)

  haloinj = dblarr(nenel,nr,nz)
  for a=0, nenel-1 do begin
    for i=0, nr-1 do begin
      for j=0, nz-1 do begin
        haloinj[a,i,j] = qromo('puc_haloinj_integrand', enel[0], enel[nenel-1], eps=eps) ; [1/s/cm^3]
      endfor
    endfor
    print, 'puc_haloinj_integrand: ' + string(float(a+1)/float(nenel)*100d, format='(I3.3)') + '%'
  endfor
  
  if not keyword_set(noplot) then puc_haloinj_plot

  time2 = systime(1)
  print, 'puc_haloinj: ', time2 - time1, ' seconds'
  

; -------- get the propagated electron number density in a MW-sized halo

  elpsi = dblarr(nenel,nr,nz)
  for i=0, nr-1 do begin
    for j=0, nz-1 do begin
      elpsi[*,i,j] = haloinj[*,i,j] / puc_edot(enel) ; [1/cm^3/MeV]
    endfor
  endfor

  if not keyword_set(noplot) then puc_elpsi_plot

 
; -------- compute the power escape fraction of a MW-sized halo

  time1 = systime(1)

  numbmwsst = dblarr(nenel) ; steady-state number of electrons in a MW-sized halo
  lumimwesc = dblarr(nenel) ; electron luminosity of a MW-sized halo
  fluxmwesc = dblarr(nenel, nr) ; electron flux out of a MW-sized halo
  escfracmw = dblarr(nenel) ; fraction of annihilation power that escapes a MW-sized halo
  escfrac = dblarr(nenel, nmass) ; fraction of annihilation power that escapes a halo of mass M
  for a=0, nenel-1 do begin
    numbmwsst[a] = 50. * 2. * !pi * kpc2cm^3 * int_2d('puc_numbmwsst_integrand', [minr,maxr], 'puc_gal_zall',48) ; [1/MeV]
    izb = max(where(elpsi[a,0,*] gt 0.))
    fluxmwesc[a,*] = reform(k0 * (enel[a] / 1d3)^delta * abs((elpsi[a,*,izb] - elpsi[a,*,izb-1]) / (z[izb] - z[izb-1]))) * kpc2cm ; [1/cm^2/s/MeV]
    lumimwesc[a] = 4. * !pi * kpc2cm^2 * qromo('puc_lumimwesc_integrand', minr, maxr, eps=eps) * enel[a]^2 ; [MeV/s]
    print, 'puc_escfrac_integrand: ' + string(float(a+1)/float(nenel)*100d, format='(I3.3)') + '%'
  endfor
  lumimwesc *= (total(lumimw) - total(numbmwsst * puc_edot(enel) * enel)) / total(lumimwesc)
  numbmwesc = lumimwesc * hubbletime / enel^2 ; [1/MeV]
  escfracmw[where(lumimw gt 0)] = lumimwesc[where(lumimw gt 0)] / lumimw[where(lumimw gt 0)] ; [1]
 
  if not keyword_set(noplot) then puc_escfrac_plot

  time2 = systime(1)
  print, 'puc_escfrac: ', time2 - time1, ' seconds'



end
