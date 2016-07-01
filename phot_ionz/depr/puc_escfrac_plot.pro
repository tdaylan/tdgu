pro puc_escfrac_plot

  common puc_common

  pname = 'puc_escfrac_plot.ps'
  dfpsplot, pname, /color
  state = sysvars(/print)
  !p.font = -1

  edot = puc_edot(enel, edotbrem=edotbrem, edotsync=edotsync, edotincs=edotincs)

  cgplot, enel[0:nenel-2]*1d-3, lumimw[0:nenel-2], /ylog, /xlog, xtitle=textoidl('E_{el} [GeV]'), $
  ;  title='Power budget of a MW-sized halo', $
  ytitle=textoidl('E^2dN/dEdt [MeV/s]'), thick=5, yrange=[1d40,1d45], xrange=[0.1, 10.]
  cgplot, enel[0:nenel-2]*1d-3, lumimwesc[0:nenel-2], thick=5, /overplot, color='red'
  cgplot, enel[0:nenel-2]*1d-3, numbmwsst[0:nenel-2] * edotincs[0:nenel-2] * enel[0:nenel-2], thick=5, /overplot, color='blue', line=1
  cgplot, enel[0:nenel-2]*1d-3, numbmwsst[0:nenel-2] * edotsync[0:nenel-2] * enel[0:nenel-2], thick=5, /overplot, color='blue', line=2
  cgplot, enel[0:nenel-2]*1d-3, numbmwsst[0:nenel-2] * edotbrem[0:nenel-2] * enel[0:nenel-2], thick=5, /overplot, color='blue', line=3
  cgplot, enel[0:nenel-2]*1d-3, numbmwsst[0:nenel-2] * edot[0:nenel-2] * enel[0:nenel-2], thick=5, /overplot, color='blue'
  cgLegend, Title=['Total injection power', 'Diffusive escaping power', $
    'ICS energy loss rate', 'Synchrontron energy loss rate', 'Bremsstrahlung energy loss rate', $
    'Total energy loss rate'], $
    Location=[0.41,0.75], VSpace=2.0, Alignment=8, thick=5, chars=1.2, $
    color=['black','red','blue','blue','blue','blue'], line=[0,0,1,2,3,0]


  colors = ['red', 'blue']

  cgplot, enel[0:nenel-2]*1d-3, numbmwsst[0:nenel-2], color='red', $
    title=textoidl('Number of e^-+e^+ in the MW'), yrange=[1d46,1d57], $
    /ylog, /xlog, xtitle=textoidl('E_{el} [GeV]'), ytitle=textoidl('dN/dE [1/MeV]'), thick=5
  cgplot, enel[0:nenel-2]*1d-3, numbmwesc[0:nenel-2], color='blue', thick=5, /overplot
  cgLegend, Title=['steady-state e^-+e^+ in the MW', textoidl('e^-+e^+ that escape the MW over a Hubble time')], $
    Location=[0.48,0.8], VSpace=2.0, Alignment=8, thick=5, chars=1.2, $
    color=colors

  readcol, '$PUC_PATH/data/ams02p.dat', enams02p, ams02p
  readcol, '$PUC_PATH/data/ams02e.dat', enams02e, ams02e

  ams02p = interpol(ams02p, enams02p, enel/1d3) * 0.1 ; [MeV/cm^2/s/sr]
  ams02e = interpol(ams02e, enams02e, enel/1d3) * 0.1 ; [MeV/cm^2/s/sr]
 
  elfluxsun = dblarr(nenel)
  for a=0, nenel-1 do elfluxsun[a] = interpolate(elpsi[a,*,*],interpol(indgen(nr),r,8.5),interpol(indgen(nz),z,0.), /grid) * lightspeed / 4. / !pi
  cgplot, enel[0:nenel-2]*1d-3, elfluxsun[0:nenel-2] * enel[0:nenel-2]^2, color='red', $
    title=textoidl('e^-/e^+ fluxes in the Solar neighborhood'), yrange=[1d-5,1d2], $
    /ylog, /xlog, xtitle=textoidl('E_{el} [GeV]'), ytitle=textoidl('E^2dN/dEdAdtd\Omega [MeV/cm^2/s/sr]'), thick=5
  cgplot, enel[0:nenel-2]*1d-3, ams02e[0:nenel-2], color='blue', line=1, thick=5, /overplot
  cgplot, enel[0:nenel-2]*1d-3, ams02p[0:nenel-2], color='blue', line=2, thick=5, /overplot
  cgplot, enel[0:nenel-2]*1d-3, ams02e[0:nenel-2] + ams02p[0:nenel-2], color='blue', thick=5, /overplot
  cgLegend, Title=['Hooperon e-/e+ flux', 'AMS-02 e-', 'AMS-02 e+', 'AMS-02 e-/e+'], $
    Location=[0.36,0.8], VSpace=2.0, Alignment=8, thick=5, chars=1.2, line=[0,1,2,0], $
    color=['red','blue','blue','blue']

  cgplot, enel[0:nenel-2]*1d-3, escfracmw[0:nenel-2], title=textoidl('Diffusive escape fraction of e^-/e^+ from the MW'), $
    /xlog, xtitle=textoidl('E_{el} [GeV]'), ytitle=textoidl('f_{esc}'), thick=5, /ylog

  puc_2dplot2, alog10(reform(escfrac[0:nenel-2,*])), enel[0:nenel-2]/1d3, mass, xtitle='E [GeV]', ytitle=textoidl('M [M_{sun}]'), /xlog, /ylog, $
    title=textoidl('log(f_{esc}), diffusive escape fraction of e^-+e^+'), nlevels=5

  restore_sysvars, state
  dfpsclose

end
