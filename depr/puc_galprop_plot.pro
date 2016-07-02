pro puc_galprop_plot

  common puc_common

  emis = mrdfits('$PUC_PATH/puc/fits/DM_emiss_50p_puc.fits',0,heademis)
  flux = mrdfits('$PUC_PATH/puc/fits/nuclei_full_50p_puc.fits',0,headflux)

  crval1 = sxpar(heademis, 'CRVAL1')
  cdelt1 = sxpar(heademis, 'CDELT1')
  radgal = dindgen(n_elements(emis[*,0,0,0])) * cdelt1 + crval1
  jradgal = where(radgal lt 1, njradgal)

  crval2 = sxpar(heademis, 'CRVAL2')
  cdelt2 = sxpar(heademis, 'CDELT2')
  zcogal = dindgen(n_elements(emis[0,*,0,0])) * cdelt2 + crval2
  jzcogal = where(abs(zcogal) lt 1, njzcogal)

  crval3 = sxpar(heademis, 'CRVAL3')
  cdelt3 = sxpar(heademis, 'CDELT3')
  neneemis = n_elements(emis[0,0,*,0])
  eneemis = 10^(dindgen(neneemis) * cdelt3 + crval3) ; [MeV]

  crval3 = sxpar(headflux, 'CRVAL3')
  cdelt3 = sxpar(headflux, 'CDELT3')
  neneflux = n_elements(flux[0,0,*,0])
  eneflux = 10^(dindgen(neneflux) * cdelt3 + crval3) ; [MeV]


;  emissivity

  emisav = dblarr(neneemis)
  for a=0, neneemis-1 do begin
    pname = 'puc_galprop_emis_' + string(a, format='(I3.3)') + '.ps'
    dfpsplot, pname, /color, xs=8, ys=6
    state = sysvars(/print)
    !p.font = -1
  
    if not (min(emis[*,*,a]) eq max(emis[*,*,a])) then puc_contour, reform(emis[*,*,a]), tit='emis [MeV/cm^3/s]'
  
    restore_sysvars, state
    dfpsclose

    emisav[a] = total(emis[min(jradgal):max(jradgal),min(jzcogal):max(jzcogal),a]) / njradgal / njzcogal 
  endfor


; flux

  tits = ['flux 0 [MeV/cm^2/s/sr]','flux 1 [MeV/cm^2/s/sr]','flux 2 [MeV/cm^2/s/sr]','flux 3 [MeV/cm^2/s/sr]','flux 4 [MeV/cm^2/s/sr]']

  nspecies = n_elements(flux[0,0,0,*])
  fluxav = dblarr(nspecies, neneflux)
  for j=0, nspecies-1 do begin
    for a=0, neneflux-1 do begin
      pname = 'puc_galprop_flux_' + string(j, format='(I1.1)') + '_' + string(a, format='(I3.3)') + '.ps'
      dfpsplot, pname, /color, xs=8, ys=6
      state = sysvars(/print)
      !p.font = -1
  
      if not (min(flux[*,*,a,j]) eq max(flux[*,*,a,j])) then puc_contour, reform(flux[*,*,a,j]), tit=tits[j]
  
      restore_sysvars, state
      dfpsclose

      fluxav[j,a] = total(flux[min(jradgal):max(jradgal),min(jzcogal):max(jzcogal),a,j]) / njradgal / njzcogal 
    endfor
  endfor

  pname = 'puc_galprop_spec.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1
 
  cgplot, eneemis, emisav, /xlog, /ylog
  cgplot, eneflux, fluxav[0,*], /xlog, /ylog
  cgplot, eneflux, fluxav[1,*], /over, line=1
  cgplot, eneflux, fluxav[2,*], /over, line=2
  cgplot, eneflux, fluxav[3,*], /over, line=3
  cgplot, eneflux, fluxav[4,*], /over, line=4

  restore_sysvars, state
  dfpsclose


  file_mkdir, '$PUC_OUTPUT_PATH/ps/galprop' 
  spawn, 'mv ~/puc_galprop_*.ps $PUC_OUTPUT_PATH/ps/galprop/'
end
