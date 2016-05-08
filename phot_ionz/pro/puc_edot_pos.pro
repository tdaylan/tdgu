pro puc_edot_pos
  if keyword_set(pos) then begin
    umag = dblarr(nr,nz)
    ngas = dblarr(nr,nz)
    ngasmax = 10. ; [1/cm^3]
    bmax = 10. ; [muG]
    rg = 4. ; [kpc]
    zg = 0.5 ; [kpc]
    rb = 5. ; [kpc]
    zb = 1. ; [kpc]
    for i=0, nr-1 do begin
      umag[i,*] = bedensfac * (bmax * exp(-r[i] / rb) * exp(-abs(z) / zb))^2 / 2. ; [MeV/cm^3]
      ngas[i,*] = ngasmax * exp(-r[i] / rg) * exp(-abs(z) / zg) ; [1/cm^3]
    endfor
    uisrftag = 'uisrf' + string(nr,format='(I2.2)') + string(nz,format='(I2.2)') + '.fits'
    if file_test(uisrftag) eq 1 then begin
      uisrf = mrdfits(uisrftag)
    endif else begin
      uisrf = dblarr(nr,nz)
      uopuv = dblarr(nr,nz)
      umfir = dblarr(nr,nz)
      isrfmethod = 2 ; CMB
      for i=0, nr-1 do begin
        for j=0, nz-1 do begin
          isrfmethod = 1 ; starlight
          uopuv[i,j] = qromo('puc_isrf_integrand', minephopuv, maxephopuv, eps=eps); [MeV/cm^3]
          isrfmethod = 3 ; IR
          umfir[i,j] = qromo('puc_isrf_integrand', minephmfir, maxephmfir, eps=eps); [MeV/cm^3]
          uisrf[i,j] = uopuv[i,j] + ucmbr + umfir[i,j]
        endfor
      endfor
      mwrfits, uisrf, uisrftag, /create
    endelse
    edotsync = 4d * thomcs * lightspeed / 3d / elmass^2 * umag * enel_^2 ; [MeV/s]
    edotbrem = alphaem * lightspeed * thomcs / 8d / !dpi * (4d * 4.579d1 - 4.446d1) * ngas * enel_ ; [MeV/s]
    edotincs = 4d * thomcs * lightspeed / 3d / elmass^2 * uisrf * enel_^2 ; [MeV/s]
  endif else begin

  for k=0, n_elements(edotens)-1 do begin 
    edot = puc_edot(edotens[k], edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs, /pos)
    tstot = edotens[k] / edot * tfac
    escapemap = tstot
    escapemap[*,*] = 0.
    escapemap[where(tstot gt interpol(tsdiff*tfac,enel,edotens[k]))] = 1.
    if min(escapemap) eq 1. then continue
    puc_2dplot, escapemap, r, z, xtitle=textoidl('r [kpc]'), $
    ytitle=textoidl('z [kpc]'), $
    title=textoidl('Escape Map at ') + string(edotens[k]/1d3, format='(G6.3)') + textoidl(' GeV')
  endfor

  for k=0, n_elements(edotens)-1 do begin 
    edot = puc_edot(edotens[k], edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs, /pos)
    tstot = edotens[k] / edot * tfac
    puc_2dplot, alog10(tstot), r, z, xtitle=textoidl('r [kpc]'), $
    ytitle=textoidl('z [kpc]'), $
    title=textoidl('\tau [Gyrs] at ') + string(edotens[k]/1d3, format='(G6.3)') + textoidl(' GeV')
  endfor

  for k=0, n_elements(edotens)-1 do begin 
    edot = puc_edot(edotens[k], edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs, /pos)
    tsbrem = edotens[k] / edotbrem * tfac
    puc_2dplot, alog10(tsbrem), r, z, xtitle=textoidl('r [kpc]'), $
    ytitle=textoidl('z [kpc]'), $
    title=textoidl('\tau_{brem} [Gyrs] at ') + string(edotens[k]/1d3, format='(G6.3)') + textoidl(' GeV')
  endfor

  for k=0, n_elements(edotens)-1 do begin 
    edot = puc_edot(edotens[k], edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs, /pos)
    tsincs = edotens[k] / edotincs * tfac
    puc_2dplot, alog10(tsincs), r, z, xtitle=textoidl('r [kpc]'), $
    ytitle=textoidl('z [kpc]'), $
    title=textoidl('\tau_{ICS} [Gyrs] at ') + string(edotens[k]/1d3, format='(G6.3)') + textoidl(' GeV')
  endfor

  for k=0, n_elements(edotens)-1 do begin 
    edot = puc_edot(edotens[k], edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs, /pos)
    tssync = edotens[k] / edotsync * tfac
    puc_2dplot, alog10(tssync), r, z, xtitle=textoidl('r [kpc]'), $
    ytitle=textoidl('z [kpc]'), $
    title=textoidl('\tau_{sync} [Gyrs] at ') + string(edotens[k]/1d3, format='(G6.3)') + textoidl(' GeV')
  endfor

  revenisrf = reverse(enisrf)

  for k=0, n_elements(enisrf)-1 do begin
    ibad = where(isrf[*,*,k,0] le 1d-9, nbad)
    if nbad gt 0 then continue
    
end
