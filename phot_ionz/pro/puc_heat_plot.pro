pro puc_heat_plot

  common puc_common

  colors = ['red5', 'grn5', 'blu5']

  pname = 'puc_heat_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  cgplot, red, phheatdm, /ylog, /xlog, yr=[1d-24,1d-17], xr=[0.09,14.], $
    thick=5, ytitle=textoidl('Q(z) [MeV/s]'), xtitle=textoidl('Redshift'), $
    pos=[0.2,0.05,0.95,0.95], /nodata
  cgplot, red, phheatdm, color=colors[0], thick=5, /over
  cgplot, red, phheatqg, color=colors[1], thick=5, /over
  cgplot, red, phheatqghm, color=colors[1], thick=5, /over, line=1
  cgplot, red, phheatqs, color=colors[2], thick=5, /over
  cgplot, red, plheatbl, color='orange', thick=5, /over


  cgLegend, Title=[ $
    'DM Annihilation', $
    'Quasar/Galaxy', $
    'HM12 Quasar/Galaxy', $
    'Quasar', $
    'Blazar' $
    ], $
    Location=[0.3,0.3], VSpace=2.0, thick=5, color=[colors[0],colors[1],colors[1],colors[2],'orange'], line=[0,0,1,0,0]

  restore_sysvars, state
  dfpsclose

end
