pro puc_phflux_plot

  common puc_common

  colors = ['red5','grn5','blu5']

  readcol, '$PUC_PATH/data/xray_background.dat', enrosat, rosat
  rosat /= 1d3 ; [MeV/cm^2/s/sr]

  pname = 'puc_phflux_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  multiplot, [2,2], /square, /dox, /doy
  min = min(opeff[0,*,*])
  max = max(opeff[0,*,*])
  cgplot, red, opeff[0,*,0], yr=[min,max], color=colors[0], thick=5, /ylog, /xlog ; , xtickv=[0,500,1000], xticks=3, ytickv=[-0.15,-.1,-.05,0,.05,.1], yticks=5 
  cgplot, red, opeff[0,*,nred/4], color=colors[1], /over, thick=5
  cgplot, red, opeff[0,*,nred/2], color=colors[2], /over, thick=5
;  xyouts, -500, -0.15, /data, al=0.5, '%'
  multiplot
  min = min(opeff[enel/3,*,*])
  max = max(opeff[enel/3,*,*])
  cgplot, red, opeff[nenel/3.,*,0], yr=[min,max], color=colors[0], thick=5, /ylog, /xlog; , xtickv=[0,500,1000], xticks=3, ytickv=[-0.15,-.1,-.05,0,.05,.1], yticks=5 
  cgplot, red, opeff[nenel/3.,*,nred/4], color=colors[1], /over, thick=5
  cgplot, red, opeff[nenel/3.,*,nred/2], color=colors[2], /over, thick=5
  multiplot
  min = min(opeff[2.*nenel/3,*,*])
  max = max(opeff[2.*nenel/3,*,*])
  cgplot, red, opeff[2.*nenel/3.,*,0], yr=[min,max], color=colors[0], thick=5, /ylog, /xlog; , xtickv=[0,500,1000], xticks=3, ytickv=[-0.15,-.1,-.05,0,.05,.1], yticks=5 
  cgplot, red, opeff[2.*nenel/3.,*,nred/4], color=colors[1], /over, thick=5
  cgplot, red, opeff[2.*nenel/3.,*,nred/2], color=colors[2], /over, thick=5
  multiplot
  min = min(opeff[nenel-2,*,*])
  max = max(opeff[nenel-2,*,*])
  cgplot, red, opeff[nenel-2,*,0], yr=[min,max], color=colors[0], thick=5, /ylog, /xlog ; , xtickv=[0,500,1000], xticks=3, ytickv=[-0.15,-.1,-.05,0,.05,.1], yticks=5 
  cgplot, red, opeff[nenel-2,*,nred/4], color=colors[1], /over, thick=5
  cgplot, red, opeff[nenel-2,*,nred/2], color=colors[2], /over, thick=5
;  xyouts, 0, -0.22, /data,  al=0.5, textoidl('\tau [ns]')
  multiplot, /reset


; --------- photon flux

  cgplot, enph * 1d6, phdmflux[*,iredp[0]], /xlog, /ylog, pos=[0.2,0.05,0.95,0.95], $
      yrange=[1d-6,1d2], xrange=[1.,1d4], $
      thick=3, ytitle=textoidl('I(E,z) [MeV/s/cm^2/sr]'), xtitle=textoidl('Photon Energy [eV]'), /nodata

; data

  r2hig = [9.14, 12.5, 18.8] /1d3
  r2low = [15.5, 11.7, 7.38] / 1d3
  r2 = (r2hig + r2low) / 2.
  r2x = [0.5, 1.14, 2.35] * 1d3
  cgColorFill, [r2x, Reverse(r2x), r2x[0]], [r2hig, r2low, r2hig[0]], color=colors[2]

  r1hig = 8.74d-3
  r1low = 4.96d-3
  r1 = (r1hig + r1low) / 2.
  r1x = 250
  cgplot, r1x, r1, err_yhigh=r1hig, err_ylow=r1low, psym=14, /over, color=colors[2]


; dm annihilations

  cgplot, enph * 1d6, phdmflux[*,iredp[0]], color=colors[0], thick=3, /over
  cgplot, enph * 1d6, phdmflux[*,iredp[1]], color=colors[0], thick=3, /over, line=1
  cgplot, enph * 1d6, phdmflux[*,iredp[2]], color=colors[0], thick=3, /over, line=2


; galaxy + quasar

  cgplot, enph * 1d6, phqgflux[*,iredp[0]], /overplot, thick=3, color=colors[1]
  cgplot, enph * 1d6, phqgflux[*,iredp[1]], /overplot, thick=3, color=colors[1], line=1
  cgplot, enph * 1d6, phqgflux[*,iredp[2]], /overplot, thick=3, color=colors[1], line=2


; HM12 galaxy + quasar

  cgplot, enph * 1d6, phqgfluxhm[*,iredp[0]], /overplot, thick=3, color=colors[1], psym=1
  cgplot, enph * 1d6, phqgfluxhm[*,iredp[1]], /overplot, thick=3, color=colors[1], psym=1, line=1
  cgplot, enph * 1d6, phqgfluxhm[*,iredp[2]], /overplot, thick=3, color=colors[1], psym=1, line=2


; quasar

  cgplot, enph * 1d6, phqsflux[*,iredp[0]], /overplot, thick=3, color=colors[2]
  cgplot, enph * 1d6, phqsflux[*,iredp[1]], /overplot, thick=3, color=colors[2], line=1
  cgplot, enph * 1d6, phqsflux[*,iredp[2]], /overplot, thick=3, color=colors[2], line=2


  vline, ryd*1d6, line=1, color='grey'
  vline, 4.*ryd*1d6, line=1, color='grey'


  cgLegend, Title=[ $
    'DM Annihilation, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '                 z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '                 z = ' + string(red[iredp[2]],format='(F3.1)'), $ 
    'Galaxy/Quasar, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '               z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '               z = ' + string(red[iredp[2]],format='(F3.1)'), $
    'HM12 Galaxy/Quasar, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '                    z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '                    z = ' + string(red[iredp[2]],format='(F3.1)'), $
    'Quasar, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    '        z = ' + string(red[iredp[1]],format='(F3.1)'), $
    '        z = ' + string(red[iredp[2]],format='(F3.1)') $
    ], $
    Location=[0.7,0.8], VSpace=1.2, Alignment=8, thick=3, chars=1.2, $
    color=[colors[0],colors[0],colors[0],colors[1],colors[1],colors[1],colors[1],colors[1],colors[1],colors[2],colors[2],colors[2]], $
    line=[0,1,2,0,1,2,0,1,2,0,1,2], psym=[14,14,14,14,14,14,1,1,1,14,14,14]

  cgLegend, Title=['   ROSAT Extragalactic Background'], Location=[0.50,0.72], thick=3, chars=1.2, $
    psym=14, color=colors[2], len=0
  cgLegend, Title=[''], Location=[0.46,0.72], thick=3, chars=1.2, $
    color=colors[2]

  cgplot, red, phdmflux1, /xlog, /ylog, pos=[0.2,0.05,0.95,0.95], $
    color=colors[0], yr=[1d-8,1d1], thick=5, xr=[min(red),red[iredp[2]]], $
    ytitle=textoidl('I_{ryd} [MeV/s/cm^2/sr]'), xtitle=textoidl('Redshift')
  cgplot, red, phqgflux1, /over, color=colors[1], thick=5
  cgplot, red, phqgfluxhm1, /over, color=colors[1], thick=5, line=1
  cgplot, red, phqsflux1, /over, color=colors[2], thick=5

  cgLegend, Title=[ $
    'DM Annihilation', $
    'Galaxy/Quasar', $
    'HM12 Galaxy/Quasar', $
    'Quasar' $
    ], $
    Location=[0.55,0.25], VSpace=2.0, thick=5, color=[colors[0],colors[1],colors[1],colors[2]], line=[0,0,1,0]

  restore_sysvars, state
  dfpsclose

end
