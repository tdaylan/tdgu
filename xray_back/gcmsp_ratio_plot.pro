pro gcmsp_ratio_plot, ratio, fpnt, nploc, psf, evt, cir, mevt, mcir, distee, distec, nevtneig, ncirneig, $
  fdmlglcdf, fdmbglcdf, lgl, bgl, fdm, $
  uniform=uniform, flat=flat

  if not (mevt[0] eq -1) then nmevt = n_elements(mevt[*,0]) else nmevt = 0
  if not (mcir[0] eq -1) then nmcir = n_elements(mcir[*,0]) else nmcir = 0
  nevt = n_elements(evt[*,0])
  ncir = n_elements(cir[*,0])

  runtag = $
    'nevt' + string(nevt, format='(I4.4)') + $
    '_ncir' + string(ncir, format='(I4.4)') + $
    '_fpnt' + string(fpnt, format='(F3.1)') + $
    '_nploc' + string(nploc, format='(I5.5)') + $
    '_psf' + string(psf, format='(F3.1)')
  if keyword_set(uniform) then runtag += '_unif'
  if keyword_set(flat) then runtag += '_flat'

  colors = ['red5', 'grn5', 'blu5']

  dfpsplot, 'gcmsp_ratio_' + runtag + '.ps', bits=8, ys=8, xs=8
  state = sysvars(/print)
  !p.font = -1

  tit = $
    'R = ' + string(ratio, format='(F5.2)') + $
    textoidl(', f_{pnt} = ') + string(fpnt, format='(F3.1)') + $
    textoidl(', N_{pnt} = ') + string(nploc, format='(I5.5)')
;    '!C' + $
;    textoidl(', \sigma = ') + string(psf, format='(F3.1)') + ' deg'
;    textoidl(', N_{evt} = ') + string(nevt, format='(I4.4)') + $
;    textoidl(', N_{cir} = ') + string(ncir, format='(I4.4)') + $
  cgplot, 1, 1, /nodata, xr=[-10.,10], yr=[-10.,10], tit=tit, $
    xtit='l', ytit='b', pos=[0.2,0.1,0.9,0.9]

  cgplots, evt[*,0], evt[*,1], /data, psym=16

  for k=0, nevt-1 do begin
    data = circle(evt[k,0], evt[k,1], psf)
    cgplots, data[0,*], data[1,*], /data, psym=-3, color=colors[2]
  endfor

  for k=0, nmevt-1 do begin
    data = circle(mevt[k,0], mevt[k,1], psf)
    cgplots, data[0,*], data[1,*], /data, psym=-3, color=colors[2], thick=5
  endfor

  for k=0, ncir-1 do begin
    data = circle(cir[k,0], cir[k,1], psf)
    cgplots, data[0,*], data[1,*], /data, psym=-3, color=colors[1]
  endfor

  for k=0, nmcir-1 do begin
    data = circle(mcir[k,0], mcir[k,1], psf)
    cgplots, data[0,*], data[1,*], /data, psym=-3, color=colors[1], thick=5
  endfor

  !p.multi = [0,1,4]

;  if nmevt gt 1 then begin
;    if n_elements(distee[uniq(distee, sort(distee))]) then cghistoplot, distee, ytit=textoidl('#'), xtit=textoidl('\theta_{ee} [degree]')
;    if n_elements(nevtneig[uniq(nevtneig, sort(nevtneig))]) then cghistoplot, nevtneig, ytit=textoidl('#'), xtit=textoidl('N_{event neighbor}')
;  endif
;
;  if nmcir gt 1 then begin
;    if n_elements(distec[uniq(distec, sort(distec))]) then cghistoplot, distec, ytit=textoidl('#'), xtit=textoidl('\theta_{ee} [degree]'), nbins=40
;    if n_elements(ncirneig[uniq(ncirneig, sort(ncirneig))]) then cghistoplot, ncirneig, ytit=textoidl('#'), xtit=textoidl('N_{cirle neigbor}'), nbins=40
;  endif
;  
;  if (not keyword_set(uniform)) and fpnt lt 1 then begin
;    !p.multi = [0,1,2]
;    cgplot, lgl, fdmlglcdf, xtit='l [degree]', ytit='P(>l)'
;    cgplot, bgl, fdmbglcdf, xtit='b [degree]', ytit='P(>l)'
;    !p.multi = 0
;    puc_contour, fdm, lgl, bgl, tit='P(l,b)', xtit='l [degree]', ytit='b [degree]'
;  endif

  !p.multi = 0

  restore_sysvars, state
  dfpsclose

  spawn, 'mv gcmsp_ratio_' + runtag + '.ps $GCMSP_OUTPUT_PATH/ps'

end
