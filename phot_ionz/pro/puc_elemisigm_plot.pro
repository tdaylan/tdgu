pro puc_elemisigm_plot

  common puc_common

  pname = 'puc_elemisigm_' + runtag + '.ps'
  dfpsplot, pname, /color, ysize=6, xsize=8
  state = sysvars(/print)
  !p.font = -1

  colors_ = ['red5','grn5','blu5']
  colors = ['red5','grn5','blu5','red5','grn5','blu5']
  
; -------- electron emissivity

  min = min([elemisigmint,elemisigmint6,elemisavgint]) < min(elemisigm)
  max = max([elemisigmint,elemisigmint6,elemisavgint]) > max(elemisigm)
  multiplot, [2,1], /square

  cgplot, enel[0:nenel-2]/1d3, elemisigm[0:nenel-2,iredp[0]], /xlog, /ylog, color=colors[0], thick=5, xr=[enel[0]/1d3,enel[nenel-2]/1d3], $
    ytitle=textoidl('j [MeV/cm^3/s]'), xtitle=textoidl('Electron Energy [GeV]'), yr=[min,max] ;, pos=[0.15,0.3,0.95,0.95]
  cgplot, enel[0:nenel-2]/1d3, elemisigm[0:nenel-2,iredp[1]], /over, color=colors[1], thick=5
  cgplot, enel[0:nenel-2]/1d3, elemisigm[0:nenel-2,iredp[2]], /over, color=colors[2], thick=5
  cgLegend, Title=[ $
    'z = ' + string(red[iredp[0]],format='(F3.1)'), $
    'z = ' + string(red[iredp[1]],format='(F3.1)'), $
    'z = ' + string(red[iredp[2]],format='(F3.1)')], $
    Location=[0.22,0.45], VSpace=2.0, thick=5, chars=1.2, color=colors_, length=0.03

  multiplot

  cgplot, red, elemisigmint, /ylog, yr=[min, max], /xlog, $
    thick=5, xtitle=textoidl('Redshift')
  cgplot, red, elemisigmint6, /over, thick=5, line=1
  cgplot, red, elemisavgint, /over, thick=5, line=2

  cgLegend, Title=[ $
    textoidl('M_{200}^{min} = 10^{-9} M') + sunsymbol(), $
    textoidl('M_{200}^{min} = 10^{-6} M') + sunsymbol(), $
    'No DM Clustering'], $
    Location=[0.35,0.2], VSpace=2.0, Alignment=8, thick=5, chars=1.2, line=[0,1,2], $
    length=0.03

  multiplot, /reset


; -------- halo mass function & halo e-/e+ luminosity

  cgplot, mass, mass * dndm[*,0], /xlog, xr=[1d-10,1d15], $
    pos=[0.15,0.05,0.8,0.95], ystyle=4, xtitle=textoidl('M_{200} [M') + sunsymbol() + ']', /nodata

  cgaxis, /ylog, yaxis=0, yr=[1d-20,1d10], $;yrange=[min(mass * dndm), max(mass * dndm)], $
    /save,  tit=textoidl('M_{200}dN_{halo}/dM_{200}(M_{200},z) [1/kpc^3]'), color=colors[0]
  cgplot, mass, mass * dndm[*,iredp[0]], /over, color=colors[0], thick=5
  cgplot, mass, mass * dndm[*,iredp[1]], /over, color=colors[0], thick=5, line=1
  cgplot, mass, mass * dndm[*,iredp[2]], /over, color=colors[0], thick=5, line=2

  cgaxis, /ylog, /yaxis, yrange=[min(lumihaloint), max(lumihaloint)], $
    /save, tit=textoidl('L_{halo}(M_{200},z) [MeV/s]'), color=colors[1]
  cgplot, mass, lumihaloint[*,iredp[0]], color=colors[1], /over, thick=5
  cgplot, mass, lumihaloint[*,iredp[1]], color=colors[1], /over, line=1, thick=5
  cgplot, mass, lumihaloint[*,iredp[2]], color=colors[1], /over, line=2, thick=5

  cgaxis, 0.93, 0.15, /normal, /ylog, /yaxis, yrange=[min(mass * dndm * lumihaloint / kpc2cm^3), max(mass * dndm * lumihaloint / kpc2cm^3)], $
    /save, tit=textoidl('M_{200}dj/dM_{200} [MeV/cm^3/s]'), color=colors[2]
  cgplot, mass, mass * dndm[*,iredp[0]] * lumihaloint[*,iredp[0]] / kpc2cm^3, color=colors[2], /over, thick=5
  cgplot, mass, mass * dndm[*,iredp[1]] * lumihaloint[*,iredp[1]] / kpc2cm^3, color=colors[2], /over, line=1, thick=5
  cgplot, mass, mass * dndm[*,iredp[2]] * lumihaloint[*,iredp[2]] / kpc2cm^3, color=colors[2], /over, line=2, thick=5

  vline, 1d-9, color='grey', line=1
  vline, 1d-6, color='grey', line=2
  
  cgLegend, Title=[ $
    'z = ' + string(red[iredp[0]],format='(F3.1)'), $
    'z = ' + string(red[iredp[1]],format='(F3.1)'), $
    'z = ' + string(red[iredp[2]],format='(F3.1)')], $
    Location=[0.50,0.23], VSpace=2.0, thick=5, $
    chars=1.2, line=[0,1,2], color=[colors[0],colors[0],colors[0]], len=0.03
  cgLegend, Title=[ '','',''], Location=[0.45,0.23], VSpace=2.0, $
    thick=5, chars=1.2, line=[0,1,2], color=[colors[1],colors[1],colors[1]], len=0.03
  cgLegend, Title=[ '','',''], Location=[0.40,0.23], VSpace=2.0, $
    thick=5, chars=1.2, line=[0,1,2], color=[colors[2],colors[2],colors[2]], len=0.03


; -------- halo mass function comparison

  cgplot, mass, mass * dndm[*,0], /xlog, /ylog, xr=[1d-10,1d15], pos=[0.15,0.05,0.8,0.95], $
    xtitle=textoidl('M_{200} [M') + sunsymbol() + ']', yr=[1d-20,1d10], $
    ytit=textoidl('M_{200}dN_{halo}/dM_{200}(M_{200},z) [1/kpc^3]'), color=colors[0]
  cgplot, mass, mass * dndm[*,iredp[1]], /over, color=colors[1]
  cgplot, mass, mass * dndm[*,iredp[2]], /over, color=colors[2]
  cgplot, masshmf, masshmf * dndmhmf0, /over, color=colors[0], line=1
  cgplot, masshmf, masshmf * dndmhmf3, /over, color=colors[1], line=1
  cgplot, masshmf, masshmf * dndmhmf6, /over, color=colors[2], line=1

  dndmhmf0 = interpol(dndmhmf0, masshmf, mass)
  cgplot, mass, mass * dndmhmf0, /over, color=colors[0], line=2

  restore_sysvars, state
  dfpsclose

end
