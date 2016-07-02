pro puc_haloelec_plot

  common puc_common

  pname = 'puc_haloelec_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  
  ;puc_2dplot2, lambda[0:nenel-2,*], enel[0:nenel-2]*1d-3, enel*1d-3, xtitle=textoidl('log(E_{el} / 1 GeV)'), $
  ;  ytitle=textoidl('log(E_s / 1 GeV)'), title=textoidl('\lambda(E_{el},E_s)'), /ylog, /xlog

  inds = [0,floor(nenel/4),floor(nenel/2),floor(3*nenel/4), nenel-3, nenel-2]
  foreach a, inds do puc_2dplot, reform(haloelec[a,nenel-1,*,*]), r, z, xtitle='r [kpc]', ytitle='z [kpc]'
    title=textoidl('I(E,E_s,r,z) at E_{el} = ') + string(enel[a]/1d3,format='(G6.3)') + $
    textoidl(' GeV, E_{s} = ') + string(enel[nenel-1]/1d3,format='(G6.3)') + ' GeV'
 
  restore_sysvars, state
  dfpsclose

end
