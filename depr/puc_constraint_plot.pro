pro puc_constraint_plot

  colors = ['red5', 'grn5', 'blu5']
  colors_ = ['black', colors]

  if not keyword_set(nenph) then nenph = 100
  if not keyword_set(mineph) then mineph = 1d-6 ; [MeV]
  if not keyword_set(maxeph) then maxeph = 1d-1 ; [MeV]
  enph = exp(dindgen(nenph) /(nenph-1) * (alog(maxeph) - alog(mineph)) + alog(mineph)) ; [MeV]

  if not keyword_set(nmasschi) then nmasschi = 5
  if not keyword_set(minmasschi) then minmasschi = 1d4 ; [MeV]
  if not keyword_set(maxmasschi) then maxmasschi = 1d5 ; [MeV]
  masschi = exp(dindgen(nmasschi) /(nmasschi-1) * (alog(maxmasschi) - alog(minmasschi)) + alog(minmasschi))

  if not keyword_set(nenel) then nenel = 100
  if not keyword_set(mineel) then mineel = 5d1 ; [MeV]
  if not keyword_set(maxeel) then maxeel = 1d5 ; [MeV]
  enel = exp(dindgen(nenel) /(nenel-1) * (alog(maxeel) - alog(mineel)) + alog(mineel)) ; [MeV]


  sigmav0 = 3d-26 ; [cm^3/s] velocity averaged thermal annihilation cross section
  if not keyword_set(ncsv) then ncsv = 99
  if not keyword_set(mincsv) then mincsv = 3d-29 ; [cm^3/s]
  if not keyword_set(maxcsv) then maxcsv = 3d-23 ; [cm^3/s]
  csv = exp(dindgen(ncsv) /(ncsv-1) * (alog(maxcsv) - alog(mincsv)) + alog(mincsv))
  csvboost = csv / sigmav0


; get measured X-ray background
 
  readcol, '$PUC_PATH/data/xray_background.dat', enms, phmsflux
  enms *= 1d-6 ; [MeV]
  minenms = min(enms)
  maxenms = max(enms)
  jenph = where(enph gt minenms and enph lt maxenms)
  phmsflux *= 1d-3 ; [MeV/s/sr/cm^2]
  phmsflux = interpol(phmsflux, enms, enph)


; compare photon fluxes from DM annihilations to the X-ray background

  anncharr = ['b', 'e', 'm', 't']
  nannch = n_elements(anncharr)
  phdmfluxgrid = dblarr(nenph, nmasschi, nannch)
  multgrid = dblarr(nmasschi, nenel, nannch) 
  excl = dblarr(nmasschi, nannch) 
  runtag = 1
  for i=0, nannch-1 do begin
    for d=0, nmasschi-1 do begin 
      ;puc_main, mdm=masschi[d], nred=9, /noplot, annch=anncharr[i]
      puc_main, mdm=masschi[d], nred=9, runtag=runtag, /noplot, annch=anncharr[i]
      readcol, '$PUC_OUTPUT_PATH/data/puc_phdmflux_' + runtag + '.dat', red, phdmflux
      for e=0, ncsv-1 do begin
        nbad = 0
        ibad = where(phmsflux[jenph] lt phdmflux[jenph] * csvboost[e], nbad)
        if nbad gt 0 then begin
          excl[d,i] = csv[e]
          break
        endif
      endfor
      readcol, '$PUC_OUTPUT_PATH/data/puc_mult_' + runtag + '.dat', enelmult, mult
      jenel = where(enel gt min(enelmult) and enel lt max(enelmult))
      multgrid[d,jenel,i] = interpol(mult, enelmult, enel[jenel])
      phdmfluxgrid[*,d,i] = phdmflux
    endfor
  endfor

  pname = 'puc_constraint.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1


; exclusion plot 1

  readcol, '$PUC_PATH/data/cmb_constraints.dat', masscmb0, $
    exclcmb0, masscmb1, exclcmb1, masscmb2, exclcmb2, masscmb3, exclcmb3, masscmb4, exclcmb4
  exclcmb0 = interpol(exclcmb0, masscmb0, masschi)
  exclcmb1 = interpol(exclcmb1, masscmb1, masschi)
  exclcmb2 = interpol(exclcmb2, masscmb2, masschi)
  exclcmb3 = interpol(exclcmb3, masscmb3, masschi)
  exclcmb4 = interpol(exclcmb4, masscmb4, masschi)

  cgplot, masschi * 1d-3, excl[*,0], /xlog, /ylog, ytit=textoidl('<\sigma v> [cm^3/s]'), $
    xtit=textoidl('M_{\chi} [Gev]'), /nodata, yr=[1d-29,1d-23], pos=[0.2,0.,0.9,1.]
  hline, 3d-26, color='grey', line=1, thick=5
  for i=0, nannch-1 do cgplot, masschi * 1d-3, excl[*,i], /over, thick=5, line=i
  cgplot, masschi * 1d-3, exclcmb4, /over, color=colors[0], thick=5
  cgplot, masschi * 1d-3, exclcmb3, /over, color=colors[1], thick=5
  cgplot, masschi * 1d-3, exclcmb0, /over, color=colors[2], thick=5

  cgLegend, Title=['This work (Preliminary)', 'WMAP9', 'Planck', 'CV Limited'], $
    Location=[0.4,0.2], VSpace=2, thick=5, color=colors_


; exclusion plot 2

  cgplot, masschi * 1d-3, excl[*,0], /xlog, /ylog, ytit=textoidl('<\sigma v> [cm^3/s]'), $
    xtit=textoidl('M_{\chi} [Gev]'), /nodata, yr=[1d-29,1d-23], pos=[0.2,0.,0.9,1.]
  hline, 3d-26, color='grey', line=1, thick=5
  for i=0, nannch-1 do cgplot, masschi * 1d-3, excl[*,i], /over, thick=5, color=colors_[i]
 
  cgLegend, Title=[anncharr[0], anncharr[1], anncharr[2], anncharr[3]], $
    Location=[0.4,0.2], VSpace=2, thick=5, color=colors_


; photon flux comparison plot

  cgplot, enph*1d6, phdmfluxgrid[*,0,0], /xlog, /ylog, /nodata, yr=[1d-8,1d2], $
    ytitle=textoidl('I(E) [MeV/s/cm^2/sr]'), xtitle=textoidl('Photon Energy [eV]')
  for i=0, nmasschi-1 do $
    for j=0, nannch-1 do $
      cgplot, enph*1d6, phdmfluxgrid[*,i,j], /over, color=colors_[j], thick=5, line=i
  cgplot, enph[jenph]*1d6, phmsflux[jenph], /over, thick=5, color='grey'



; annihilation energy

  cgplot, enel * 1d-3, enel^2 * multgrid[0,*,0] * 1d-3, color=colors[0], thick=5, /xlog, /ylog, $
    ytitle=textoidl('E^2dN/dE [GeV]'), xtitle=textoidl('Electron Energy [GeV]'), yr=[0.1,50], /nodata
  for i=0, nmasschi-1 do $
    for j=0, nannch-1 do $
      cgplot, enel * 1d-3, enel^2 * multgrid[i,*,j] * 1d-3, color=colors_[j], thick=5, /over, line=i

  restore_sysvars, state
  dfpsclose

  file_mkdir, '$PUC_OUTPUT_PATH/ps' 
  spawn, 'mv ~/puc_constraint.ps $PUC_OUTPUT_PATH/ps/'
end
