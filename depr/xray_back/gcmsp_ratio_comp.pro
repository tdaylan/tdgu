pro gcmsp_ratio_comp
 
  psf = 0.5
  nevt = 100
  ncir = 100

; array of fraction of total flux due to point sources

  fpnt = findgen(10) / 9.


; array of number of point sources

  nploc = floor(10^(4. * findgen(10) / 9.))


; calculate the test statistics for each variation

  ratio = fltarr(n_elements(fpnt), n_elements(nploc))
  for k=0, n_elements(fpnt)-1 do for l=0, n_elements(nploc)-1 do ratio[k,l] = gcmsp_ratio_main(nploc[l], psf, nevt, ncir, fpnt[k], /noplot)

  ratiounif = fltarr(n_elements(fpnt), n_elements(nploc))
  for k=0, n_elements(fpnt)-1 do for l=0, n_elements(nploc)-1 do ratiounif[k,l] = gcmsp_ratio_main(nploc[l], psf, nevt, ncir, fpnt[k], /uniform, /noplot)

  runtag = $
    'nevt' + string(nevt, format='(I4.4)') + $
    '_ncir' + string(ncir, format='(I4.4)') + $
    '_psf' + string(psf, format='(F3.1)')

  dfpsplot, 'gcmsp_ratio_comp_' + runtag + '.ps', bits=8, ys=8, xs=8
  state = sysvars(/print)
  !p.font = -1

  puc_contour, ratio, fpnt, nploc, /ylog, xtit=textoidl('f_{pnt}'), ytit=textoidl('N_{pnt}'), tit='R'
  puc_contour, ratio, fpnt, nploc, /ylog, xtit=textoidl('f_{pnt}'), ytit=textoidl('N_{pnt}'), tit='R', range=[0.75*max(ratio),max(ratio)]
 
  puc_contour, ratiounif, fpnt, nploc, /ylog, xtit=textoidl('f_{pnt}'), ytit=textoidl('N_{pnt}'), tit='R'
  puc_contour, ratiounif, fpnt, nploc, /ylog, xtit=textoidl('f_{pnt}'), ytit=textoidl('N_{pnt}'), tit='R', range=[0.75*max(ratio),max(ratio)]
 
  restore_sysvars, state
  dfpsclose

  file_mkdir, '$GCMSP_OUTPUT_PATH/ps' 
  spawn, 'mv gcmsp_ratio_comp_' + runtag + '.ps $GCMSP_OUTPUT_PATH/ps'

end
