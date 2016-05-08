pro puc_emismw_plot

  common puc_common

  pname = 'puc_emismw_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  puc_2dplot, rhomw, r, z, xtitle=textoidl('r [kpc]'), $
    ytitle=textoidl('z [kpc]'), title=textoidl('Hooperon energy density in the MW, \rho(r,z) [MeV/cm^3]')
  for a=0, nenel-1 do begin
    puc_2dplot, reform(emismw[a,*,*]), r, z, xtitle=textoidl('r [kpc]'), $
      ytitle=textoidl('z [kpc]'), title=textoidl('Hooperon annihilation e^-/e^+ emissivity [MeV/cm^3/s] in the MW, at ') $
      + string(enel[a]/1d3, format='(G6.3)') + ' GeV'
  endfor
  
  cgplot, enel*1d-3, emismwsun, /xlog, /ylog, title=textoidl('e^-+e^+ emissivity in the Solar neighborhood'), $
    thick=5, ytitle=textoidl('E^2dN/dE/dV/dt [MeV/cm^3/s]'), xtitle=textoidl('Energy [GeV]')

  restore_sysvars, state
  dfpsclose

end
