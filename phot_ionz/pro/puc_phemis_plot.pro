pro puc_phemis_plot

  common puc_common

  pname = 'puc_phemis_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  colors = ['red5','grn5','blu5']
 
  cgplot, enph * 1d6, phdmemis[*,iredp[0]], /xlog, /ylog, pos=[0.2,0.05,0.95,0.95], $
      yrange=[1d-33,1d-25], xrange=[1.,1d4], $
      thick=3, ytitle=textoidl('j(E,z) [MeV/s/cm^3]'), xtitle=textoidl('Photon Energy [eV]'), /nodata
  cgplot, enph * 1d6, phdmemis[*,iredp[0]], color=colors[0], thick=3, /over
  cgplot, enph * 1d6, phdmemis[*,iredp[1]], color=colors[0], thick=3, /over, line=1
  cgplot, enph * 1d6, phdmemis[*,iredp[2]], color=colors[0], thick=3, /over, line=2
  cgplot, enph * 1d6, phqgemis[*,iredp[0]], /overplot, thick=3, color=colors[1]
  cgplot, enph * 1d6, phqgemis[*,iredp[1]], /overplot, thick=3, color=colors[1], line=1
  cgplot, enph * 1d6, phqgemis[*,iredp[2]], /overplot, thick=3, color=colors[1], line=2
  cgplot, enph * 1d6, phqsemis[*,iredp[0]], /overplot, thick=3, color=colors[2]
  cgplot, enph * 1d6, phqsemis[*,iredp[1]], /overplot, thick=3, color=colors[2], line=1
  cgplot, enph * 1d6, phqsemis[*,iredp[2]], /overplot, thick=3, color=colors[2], line=2

  vline, ryd*1d6, line=1, color='grey'
  vline, 4.*ryd*1d6, line=1, color='grey'

  cgLegend, Title=[ $
    'DM Annihilation, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '                 z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '                 z = ' + string(red[iredp[2]],format='(F3.1)'), $ 
    'HM12 Galaxy/Quasar, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '                       z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '                       z = ' + string(red[iredp[2]],format='(F3.1)'), $
    'HM12 Quasar, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '             z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '             z = ' + string(red[iredp[2]],format='(F3.1)') $
    ], $
    Location=[0.7,0.8], VSpace=1.2, Alignment=8, thick=3, chars=1.2, $
    color=[colors[0],colors[0],colors[0],colors[1],colors[1],colors[1],colors[2],colors[2],colors[2]], line=[0,1,2,0,1,2,0,1,2]


  cgplot, red, phdmemis1, /ylog, /xlog, thick=5, pos=[0.2,0.05,0.95,0.95], $
    ytitle=textoidl('j_{ryd}(z) [MeV/cm^3/s]'), color=colors[0], yr=[1d-32,1d-25]
  cgplot, red, phqgemis1, /over, thick=5, color=colors[1]
  cgplot, red, phqsemis1, /over, thick=5, color=colors[2]

  cgLegend, Title=[ $
    'DM Annihilation, 1 Ryd', $
    'Quasar + Galaxy, 1 Ryd', $
    'Quasar, 1 Ryd' $
    ], $
    Location=[0.6,0.7], VSpace=2.0, Alignment=8, thick=5, chars=1.2, color=colors


  restore_sysvars, state
  dfpsclose

end
