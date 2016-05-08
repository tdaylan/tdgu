pro puc_phemis_frame_plot

  common puc_common

  pnameall = 'puc_phemis_frame_*_' + runtag + '.ps'

  for c=0, max(where(red lt 14)) do begin

    pname = 'puc_phemis_frame_' + string(nred - 1 - c, format='(I3.3)') + '_' + runtag + '.ps'

    dfpsplot, pname, /color, xs=8, ys=6
    state = sysvars(/print)
    !p.font = -1
  
    colors = ['red5','grn5','blu5']
   
    cgplot, enph * 1d6, phdmemis[*,c], /xlog, /ylog, pos=[0.2,0.05,0.95,0.95], $
        yrange=[1d-33,1d-25], xrange=[1.,1d4], color=colors[0], title='z = ' + string(red[c], format='(G5.2)'), $
        thick=3, ytitle=textoidl('j(E,z) [MeV/s/cm^3]'), xtitle=textoidl('Photon Energy [eV]')
    cgplot, enph * 1d6, phqgemis[*,c], /overplot, thick=3, color=colors[1]
    cgplot, enph * 1d6, phqsemis[*,c], /overplot, thick=3, color=colors[2]
  
    vline, ryd*1d6, line=1, color='grey'
    vline, 4.*ryd*1d6, line=1, color='grey'
   
    cgLegend, Title=[ $
      'DM Annihilation', $
      'Galaxy/Quasar', $
      'Quasar' $
      ], $
      Location=[0.55,0.9], VSpace=2, thick=3, $
      color=colors
  
    restore_sysvars, state
    dfpsclose
  endfor

  file_mkdir, '$PUC_OUTPUT_PATH/ps/phemis'
  spawn, 'mv puc_phemis_frame_*.ps $PUC_OUTPUT_PATH/ps/phemis'
end
