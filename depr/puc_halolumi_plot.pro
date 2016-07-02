pro puc_halolumi_plot

  common puc_common

  pname = 'puc_halolumi_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  cgplot, enel[0:nenel-2] * 1d-3, inmult[0:nenel-2] * enel[0:nenel-2]^2 * 1d-3, /xlog, /ylog, title='Electron Injection Multiplicity', $
      thick=5, ytitle=textoidl('E^2 dN/dE [GeV]'), xtitle=textoidl('Energy [GeV]')

  cgplot, enel[0:nenel-2] * 1d-3, inemissun * enel[0:nenel-2]^2, /xlog, /ylog, title='Electron Injection Emissivity near the Sun', $
      thick=5, ytitle=textoidl('E^2 dN/dE [MeV/cm^3/s]'), xtitle=textoidl('Energy [GeV]')

  puc_2dplot, alog10(rho2), r, z, xtitle='r [kpc]', ytitle='z [kpc]', title=textoidl('log(\rho^2 * cm^6/MeV^2)')
  puc_2dplot, alog10(inemis[*,*,0] * enel[0]^2), r, z, xtitle='r [kpc]', ytitle='z [kpc]', $
    title=textoidl('log(j * cm^3*s/MeV) at ') + string(enel[0]/1d3, format='(G6.3)') + ' GeV'
  puc_2dplot, alog10(inemis[*,*,nenel/2] * enel[nenel/2]^2), r, z, xtitle='r [kpc]', ytitle='z [kpc]', $
    title=textoidl('log(j * cm^3*s/MeV) at ') + string(enel[nenel/2]/1d3, format='(G6.3)') + ' GeV'
  puc_2dplot, alog10(inemis[*,*,nenel-1] * enel[nenel-1]^2), r, z, xtitle='r [kpc]', ytitle='z [kpc]', $
    title=textoidl('log(j * cm^3*s/MeV) at ') + string(enel[nenel-1]/1d3, format='(G6.3)') + ' GeV'

  restore_sysvars, state
  dfpsclose

end
