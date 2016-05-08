pro puc_isrf_plot

  common puc_common

  pname = 'puc_isrf_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  colors = ['blue', 'red', 'red', 'red', 'green', 'green', 'orange', 'orange']

  cgplot, enisrf * 1d6, nphcmbrgp * enisrf^2 * 1d6, /xlog, /ylog, color=colors[0], yrange=[1d-5,1d5], $
    title='ISRF Energy Density', thick=5, ytitle=textoidl('E^2 dN/dE [eV/cm^3]'), xtitle=textoidl('Energy [eV]')
  cgplot, enisrf * 1d6, nphcmbr[*,0] * enisrf^2 * 1d6, color=colors[1], thick=5, line=2, /overplot
  cgplot, enisrf * 1d6, nphcmbr[*,nred/4] * enisrf^2 * 1d6, color=colors[2], thick=5, line=3, /overplot
  cgplot, enisrf * 1d6, nphcmbr[*,nred-1] * enisrf^2 * 1d6, color=colors[3], thick=5, line=4, /overplot
  cgplot, enisrf * 1d6, nphopuvav * enisrf^2 * 1d6, color=colors[4], thick=5, /overplot
  cgplot, enisrf * 1d6, nphopuvinav * enisrf^2 * 1d6, color=colors[5], line=1, thick=5, /overplot
  cgplot, enisrf * 1d6, nphmfirav * enisrf^2 * 1d6, color=colors[6], thick=5, /overplot
  cgplot, enisrf * 1d6, nphmfirinav * enisrf^2 * 1d6, color=colors[7], line=1, thick=5, /overplot
  
  cgLegend, Title=['GALPROP CMB at z = 0', $
    'CMB at z = ' + string(red[0],format='(G5.3)'), $
    'CMB at z = ' + string(red[nred/4],format='(G5.3)'), $
    'CMB at z = ' + string(red[nred-1],format='(G5.3)'), $
    'Optical/UV, whole galaxy', 'Optical/UV inner galaxy', $
    'IR, whole galaxy', 'IR, inner galaxy'], $
    Location=[0.73,0.75], VSpace=2.0, Alignment=8, thick=5, chars=1.2, color=colors, line=[0,2,3,4,0,1,0,1]

  restore_sysvars, state
  dfpsclose

end
