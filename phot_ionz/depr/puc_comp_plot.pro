pro puc_comp_plot

  runtag = 1
  puc_main, runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr0
  puc_main, mdm=1d4, runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr1
  puc_main, mdm=1d5, runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr2
  puc_main, conmod='hig', runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr3
  puc_main, conmod='low', runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr4
  puc_main, conmod='cut', runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr5
  puc_main, minhalomass=1d-6, runtag=runtag
  readcol, 'dmionr_' + runtag + '.dat', red, dmionr6

  colors = ['black', 'red5', 'grn5', 'blu5', 'yellow', 'orange', 'olive']

  pname = 'puc_comp.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  cgplot, red, dmionr0, /ylog, yrange=[1d-18,1d-8], /xlog, xr=[0.09,40.], $
    thick=5, ytitle=textoidl('\Gamma [1/s]'), xtitle=textoidl('Redshift'), $
    pos=[0.2,0.05,0.95,0.95], color=colors[0]
  cgplot, red, dmionr1, color=colors[1], thick=5, /over
  cgplot, red, dmionr2, color=colors[2], thick=5, /over
  cgplot, red, dmionr3, color=colors[3], thick=5, /over
  cgplot, red, dmionr4, color=colors[4], thick=5, /over
  cgplot, red, dmionr5, color=colors[5], thick=5, /over
  cgplot, red, dmionr6, color=colors[6], thick=5, /over

  cgLegend, Title=[ $
    textoidl('Baseline'), $
    textoidl('M_\chi = 10 GeV'), $
    textoidl('M_\chi = 100 GeV'), $
    textoidl('c = 6.5(M/M_0)^{-0.15}H^{-2/3}(z)'), $
    textoidl('c = 6.5(M/M_0)^{-0.05}H^{-2/3}(z)'), $
    textoidl('Saturated c'), $
    textoidl('M_{min} = 10^6 M') + sunsymbol()], $
    Location=[0.5,0.92], VSpace=1.6, thick=5, color=colors, len=0.02

  restore_sysvars, state
  dfpsclose

  file_mkdir, '$PUC_OUTPUT_PATH/ps' 
  spawn, 'mv ~/puc_comp.ps $PUC_OUTPUT_PATH/ps/'

end
