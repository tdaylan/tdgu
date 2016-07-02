pro puc_edot_plot

  common puc_common

  colors = ['red5','grn5','blu5']
  magfdplt = 5. ; [muG]
  gasdnplt = 1. ; [1/cm^3]
  uisrfplt = 3d-7 ; [MeV/cm^3]

  edot = puc_edot(ienel, magfd=magfdplt, gasdn=gasdnplt, uisrf=uisrfplt, edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs)

  magfdpltin = 50. ; [muG]
  gasdnpltin = 10. ; [1/cm^3]
  uisrfpltin = 6d-7 ; [MeV/cm^3]

  temp = puc_edot(ienel, magfd=magfdpltin, gasdn=gasdnpltin, uisrf=uisrfpltin, edotbrem=edotbremin, edotsync=edotsyncin, edotincs=edotincsin)

  tscmbr = dblarr(nenel,nred)
  for c=0, nred-1 do tscmbr[*,c] = enel / puc_edot(ienel, c, /cmbr)

  tstot = enel / edot
  tsbrem = enel / edotbrem
  tsincs = enel / edotincs
  tssync = enel / edotsync
  tsbremin = enel / edotbremin
  tsincsin = enel / edotincsin
  tssyncin = enel / edotsyncin

  tshubb = tstot
  tshubb[*] = htime ; [s]

  tfac = 3.171d-17 ; [s/Gyrs]

  minedot = min(edotbrem) < min(edotincs) < min(edotsync)
  maxedot = max(edot)

  mints = min(tstot)
  maxts = max(tsbrem) > max(tsincs) > max(tssync)

  tsdiff = lz^2 / (k0 * (enel/1d3)^delta)

  pname = 'puc_edot_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1
  
  cgplot, enel, tstot*tfac, /xlog, /ylog, $ ;ytickformat='(G7.1)', aspect=1.0
    position=[0.2,0.1,0.9,0.9], thick=5, ytitle=textoidl('\tau [Gyrs]'), xtitle=textoidl('Energy [MeV]'), $
    yrange=[1d-3,1d3], /nodata
  cgplot, enel, tsbrem*tfac, thick=5, /overplot, color='orange'
  cgplot, enel, tsbremin*tfac, thick=5, /overplot, color='orange', line=1
  cgplot, enel, tsincs*tfac, thick=5, /overplot, color=colors[0]
  cgplot, enel, tsincsin*tfac, thick=5, /overplot, color=colors[0], line=1
  cgplot, enel, tssync*tfac, thick=5, /overplot, color=colors[1]
  cgplot, enel, tssyncin*tfac, thick=5, /overplot, color=colors[1], line=1
  cgplot, enel, tshubb*tfac, thick=5, /over, line=1, color='grey'
  cgplot, enel, tsdiff*tfac, thick=5, /over
  cgplot, enel, tscmbr[*,iredp[0]]*tfac, thick=5, /over, color=colors[2]
  cgplot, enel, tscmbr[*,iredp[1]]*tfac, thick=5, /over, color=colors[2], line=1
  cgplot, enel, tscmbr[*,iredp[2]]*tfac, thick=5, /over, color=colors[2], line=2

  cgLegend, Title=[ $
    textoidl('Brem (inner galaxy), \rho_{gas} = ') + string(gasdnplt, format='(G6.3)') + textoidl(' cm^{-3}'), $
    textoidl('Brem (inner galaxy), \rho_{gas} = ') + string(gasdnpltin, format='(G6.3)') + textoidl(' cm^{-3}'), $
    textoidl('ICS  (whole galaxy), u_{isrf} = ') + string(uisrfplt*1d6, format='(G6.3)') + textoidl(' eV/cm^{3}'), $
    textoidl('ICS  (inner galaxy), u_{isrf} = ') + string(uisrfpltin*1d6, format='(G6.3)') + textoidl(' eV/cm^{3}'), $
    textoidl('Sync (whole galaxy), B = ') + string(magfdplt, format='(G6.3)') + textoidl(' \mu G'), $
    textoidl('Sync (inner galaxy), B = ') + string(magfdpltin, format='(G6.3)') + textoidl(' \mu G'), $
    'Hubble time', $
    'Diffusion Time Scale', $
    'ICS on CMB, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    'ICS on CMB, z = ' + string(red[iredp[1]],format='(F3.1)'), $
    'ICS on CMB, z = ' + string(red[iredp[2]],format='(F3.1)') $
    ], $
    Location=[0.63,0.78], VSpace=2.0, Alignment=8, thick=5, chars=1.2, $
    color=['orange', 'orange', colors[0], colors[0], colors[1], colors[1], 'grey', 'black', colors[2], colors[2], colors[2]], $
    line=[0,1,0,1,0,1,1,0,0,1,2]
 
  restore_sysvars, state
  dfpsclose

end
