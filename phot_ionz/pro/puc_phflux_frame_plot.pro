pro puc_phflux_frame_plot

  common puc_common

  pnameall = 'puc_phflux_frame_*_' + runtag + '.ps'

  
  for c=0, max(where(red lt 14)) do begin

    pname = 'puc_phflux_frame_' + string(nred - 1 - c, format='(I3.3)') + '_' + runtag + '.ps'

    dfpsplot, pname, /color, xs=8, ys=6
    state = sysvars(/print)
    !p.font = -1
  
    colors = ['red5','grn5','blu5']
   
    cgplot, enph * 1d6, phdmflux[*,c], /xlog, /ylog, pos=[0.2,0.05,0.95,0.95], $
        yrange=[1d-6,1d2], xrange=[1.,1d4], /nodata, title='z = ' + string(red[c], format='(G5.2)'), $
        thick=3, ytitle=textoidl('I(E,z) [MeV/s/cm^2/sr]'), xtitle=textoidl('Photon Energy [eV]')
    cgplot, enph * 1d6, phdmflux[*,c], /overplot, thick=3, color=colors[0]
    cgplot, enph * 1d6, phqgflux[*,c], /overplot, thick=3, color=colors[1]
    cgplot, enph * 1d6, phqgfluxhm[*,c], /overplot, thick=3, color=colors[1], line=1
    cgplot, enph * 1d6, phqsflux[*,c], /overplot, thick=3, color=colors[2]
  
    vline, ryd*1d6, line=1, color='grey'
    vline, 4.*ryd*1d6, line=1, color='grey'
   
    cgLegend, Title=[ $
      'DM Annihilation', $
      'Galaxy/Quasar', $
      'HM12 Galaxy/Quasar', $
      'Quasar' $
      ], $
      Location=[0.55,0.9], VSpace=2, thick=3, line=[0,0,1,0], $
      color=[colors[0],colors[1],colors[1],colors[2]]
  
    restore_sysvars, state
    dfpsclose
  endfor

  file_mkdir, '$PUC_OUTPUT_PATH/ps/phflux'
  spawn, 'mv puc_phflux_frame_*.ps $PUC_OUTPUT_PATH/ps/phflux'

end
