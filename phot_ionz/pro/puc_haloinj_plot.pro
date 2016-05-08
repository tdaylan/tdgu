pro puc_haloinj_plot

  common puc_common

  pname = 'puc_haloinj_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  for k=0, nenel-2 do puc_2dplot, reform(haloinj[k,*,*]), r, z, xtitle='r [kpc]', ytitle='z [kpc]', $
    title=textoidl('I_{inj}(E,r,z) [1/cm^3/s] at E_{el} = ') + string(enel[k]/1d3,format='(G6.3)') + textoidl(' GeV')

  restore_sysvars, state
  dfpsclose

end
