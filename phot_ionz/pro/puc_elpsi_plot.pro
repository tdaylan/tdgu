pro puc_elpsi_plot

  common puc_common

  pname = 'puc_elpsi_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  for k=0, nenel-2 do puc_2dplot, reform(elpsi[k,*,*]), r, z, xtitle='r [kpc]', ytitle='z [kpc]', $
    title=textoidl('\Psi(E,r,z) [1/cm^3/MeV] at E_{el} = ') + string(enel[k]/1d3,format='(G6.3)') + textoidl(' GeV')
    
  restore_sysvars, state
  dfpsclose

end
