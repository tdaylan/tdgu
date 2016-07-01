pro puc_phionr_plot

  common puc_common

  colors = ['red5', 'grn5', 'blu5']

  pname = 'puc_phionr_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

;  cgplot, enph*1d6, dmionrinteg[*,iredp[0]] * enph, /ylog, yr=[1-22,1d-10], /xlog, $
;    ytitle=textoidl('d\Gamma/dE [1/s/MeV]'), xtitle=textoidl('Photon Energy'), $
;    pos=[0.2,0.05,0.95,0.95], /nodata
;  cgplot, enph*1d6, dmionrinteg[*,iredp[0]], /over, thick=5, color=colors[0], line=0, xr=[13.6,max(enph*1d6)]
;  cgplot, enph*1d6, dmionrinteg[*,iredp[1]], /over, thick=5, color=colors[0], line=1, xr=[13.6,max(enph*1d6)]
;  cgplot, enph*1d6, sigmaion, /over, thick=5, color=colors[1], xr=[13.6,max(enph*1d6)]
;  cgplot, enph*1d6, 4. * !pi * phflux[*,iredp[0]] / enph, /over, thick=5, color=colors[2]
;  cgplot, enph*1d6, 4. * !pi * phflux[*,iredp[1]] / enph, /over, thick=5, color=colors[2], line=1
;  vline, 13.6, line=1, thick=5, color='grey'
; 
;  cgLegend, Title=[ $
;    textoidl('d\Gamma/dE, z = ') + string(red[iredp[0]], format='(F3.1)'), $
;    textoidl('d\Gamma/dE, z = ') + string(red[iredp[0]], format='(F3.1)'), $
;    textoidl('\sigma(E)'), $
;    textoidl('I(E,z), z = ') + string(red[iredp[0]], format='(F3.1)'), $
;    textoidl('I(E,z), z = ') + string(red[iredp[1]], format='(F3.1)')], $
;    Location=[0.6,0.3], VSpace=2.0, Alignment=8, thick=5, chars=1.2, $
;    color=[colors[0],colors[0],colors[1],colors[2],colors[2]], line=[0,1,0,0,1]


  cgplot, red, phionrdm, /ylog, yrange=[1d-17,1d-11], /xlog, xr=[0.09,14.], $
    thick=5, ytitle=textoidl('\Gamma [1/s]'), xtitle=textoidl('Redshift'), $
    pos=[0.2,0.05,0.95,0.95], /nodata
  cgplot, red, phionrdm, color=colors[0], thick=5, /over
  cgplot, red, phionrqg, color=colors[1], thick=5, /over
  cgplot, red, phionrqghm, color=colors[1], thick=5, /over, line=1
  cgplot, red, phionrqs, color=colors[2], thick=5, /over


; Kollmeier et al.

  cgplot, 0.1, 1.8d-13, psym=16, color='orange', /over


; Fauchere-Giguerre

  readcol, '$PUC_PATH/data/fg_hig.dat', fgx, fghig
  readcol, '$PUC_PATH/data/fg_low.dat', fgx, fglow
  fg = (fghig + fglow) / 2. * 1d-12
  fghig *= 1d-12
  fglow *= 1d-12
 

; Becker et al.

  readcol, '$PUC_PATH/data/be_hig.dat', bex, behig
  readcol, '$PUC_PATH/data/be_low.dat', bex, below
  be = (behig + below) / 2. * 1d-12
  behig *= 1d-12
  below *= 1d-12

  cgplot, bex, be, psym=14, err_yhigh=behig, err_ylow=below, color='olive', /over, thick=5

  cgLegend, Title=[ $
    'DM Annihilation', $
    'Galaxy/Quasar', $
    'HM12 Galaxy/Quasar', $
    'Quasar' $
    ], $
    Location=[0.3,0.3], VSpace=2.0, thick=5, color=[colors[0],colors[1],colors[1],colors[2]], $
    line=[0,0,1,0]

  cgLegend, Title=[ $
    'Becker et al', $
    'Kollmeier et al'], $
    Location=[0.3,0.8], VSpace=2.0, thick=5, color=['olive','orange'], psym=[16,14], len=0

  restore_sysvars, state
  dfpsclose

end
