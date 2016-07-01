pro puc_eloss_plot

  common puc_common


  puc_edot(enedot, )

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
  nenisrf = n_elements(enisrf)
  nphstar = dblarr(nenisrf)
  nphcmbr = dblarr(nenisrf)
  nphblac = dblarr(nenisrf)
  confinerinds = where(risrf lt escaper, nconr)
  confinezinds = where(abs(zisrf) lt escapez, nconz)
  ncellcon = nconr + nconz
  ncelltot =  n_elements(isrf[*,*,0,0])
  for a=0, nenisrf-1 do begin
    nphstar[a] = total(isrf[confinerinds,confinezinds,a,0]) / ncellcon / enisrf[a]^2 ; [1/cm^3/MeV]
    nphcmbr[a] = total(isrf[*,*,a,2]) / ncelltot / enisrf[a]^2 ; [1/cm^3/MeV]
    if enisrf[a] / kt lt 50d then nphblac[a] = 8 * !dpi * enisrf[a]^2 / (lightspeed * planck)^3 / (exp(enisrf[a] / kt) - 1) ; [1/cm^3/MeV]
  endfor
  enisrf = reverse(enisrf)
  nphcmbr = reverse(nphcmbr)
  nphstar = reverse(nphstar)
  nphblac = reverse(nphblac)

  if not keyword_set(noplot) then puc_isrf_plot

  puc_2dplot, edotgrid, r, z, xtitle=textoidl(' / 1 GeV)'), $
    ytitle=textoidl('log(E_s / 1 GeV)'), title=textoidl('\lambda(E_{el},E_s)')

; (r,z) maps for different final electron energies at a certain injection energy

  inds = [0,floor(nenel/4),floor(nenel/2),floor(3*nenel/4), nenel-3, nenel-2]
  if not keyword_set(loggrid) then $
    puc_2dplot, alog10(edotgrid[a,nenel-1,*,*]), r, z, xtitle='r [kpc]', ytitle='z [kpc]', $
      title=textoidl('log(I(E,E_s,r,z)) at E_{el} = ') + string(enel[a]/1d3,format='(G6.3)') + $
      textoidl(' GeV, E_{s} = ') + string(enel[nenel-1]/1d3,format='(G6.3)') + ' GeV' $
  else $
    foreach a, inds do puc_2dplot, alog10(reform(haloelec[a,nenel-1,*,*])), alog10(r), alog10(z), xtitle='log(r / 1 kpc)', ytitle='log(z / 1 kpc)', $
      title=textoidl('log(I(E,E_s,r,z)) at E_{el} = ') + string(enel[a]/1d3,format='(G6.3)') + $
      textoidl(' GeV, E_{s} = ') + string(enel[nenel-1]/1d3,format='(G6.3)') + ' GeV'
  restore_sysvars, state
  dfpsclose

end
