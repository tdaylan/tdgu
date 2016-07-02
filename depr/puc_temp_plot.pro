pro puc_temp_plot

  common puc_common

  pname = 'puc_temp_' + runtag + '.ps'
  dfpsplot, pname, /color, ysize=6, xsize=8
  state = sysvars(/print)
  !p.font = -1

  colors = ['red5','grn5','blu5']

  cgplot, red, tempdm, /xlog, /ylog, thick=5, /nodata, yr=[1d3,1d5], $
    ytitle=textoidl('T [K]'), xtitle=textoidl('Redshift'), pos=[0.15,0.3,0.95,0.95]
  cgplot, red, tempdm, color=colors[0], thick=5, /over
  cgplot, red, tempqs, color=colors[1], thick=5, /over
  cgplot, red, tempbl, color=colors[2], thick=5, /over
;  cgplot, red, cmbt, color='grey', line=1, thick=5, /over

  cgLegend, Title=[ $
    'DM', $
    'Quasar', $
    'Blazar' $
;    'CMB Temperature' $
    ], $
    Location=[0.22,0.45], VSpace=2.0, thick=5, color=colors

  restore_sysvars, state
  dfpsclose

end
