pro puc_icspower_plot

  common puc_common

  pname = 'puc_icspower_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  min = min(icspowercmbr)
  max = max(icspowercmbr)

  colors = ['red5','red5','grn5','grn5','blu5','blu5']
  colors_ = ['red5','grn5','blu5']

  cgplot, enph*1d6, icspowercmbr[*,0,iredp[0]], $
    xtit=textoidl('Photon Energy [eV]'), ytit=textoidl('P_{IC} [MeV/s]'), yr=[1d-22,1d-8], $
    thick=5, /xlog, /ylog, color=colors[0], pos=[0.25,0.,1.,1.]
  cgplot, enph*1d6, icspowercmbr[*,0,iredp[1]], thick=5, /over, color=colors[1], line=1
  cgplot, enph*1d6, icspowercmbr[*,nenel/4,iredp[0]], thick=5, /over, color=colors[2]
  cgplot, enph*1d6, icspowercmbr[*,nenel/4,iredp[1]], thick=5, /over, color=colors[3], line=1
  cgplot, enph*1d6, icspowercmbr[*,nenel/2,iredp[0]], thick=5, /over, color=colors[4]
  cgplot, enph*1d6, icspowercmbr[*,nenel/2,iredp[1]], thick=5, /over, color=colors[5], line=1

  cgLegend, Title=[ $
    string(enel[0]/1d3,format='(G5.2)') + ' GeV, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    string(enel[0]/1d3,format='(G5.2)') + ' GeV, z = ' + string(red[iredp[1]],format='(F3.1)'), $
    string(enel[nenel/4]/1d3,format='(G5.2)') + ' GeV, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    string(enel[nenel/4]/1d3,format='(G5.2)') + ' GeV, z = ' + string(red[iredp[1]],format='(F3.1)'), $
    string(enel[nenel/2]/1d3,format='(G5.2)') + ' GeV, z = ' + string(red[iredp[0]],format='(F3.1)'), $
    string(enel[nenel/2]/1d3,format='(G5.2)') + ' GeV, z = ' + string(red[iredp[1]],format='(F3.1)')], $
    Location=[0.2,0.9], VSpace=1.5, thick=5, color=colors, line=[0,1,0,1,0,1]


  cgplot, red, ucmbr, xtit='Redshift', ytit=textoidl('\epsilon_{CMB}(z) [MeV/cm^3]'), /xlog, /ylog, thick=5, pos=[0.2,0.05,0.9,0.95]

  icspowercmbr[where(icspowercmbr lt 1d-30)] = 0.
  puc_contour, icspowercmbr[*,*,0], enph*1d6, enel*1d-3, /xlog, /ylog, /log, $
    xtit=textoidl('Photon Energy [eV]'), tit=textoidl('P_{IC} [MeV/s]'), ytit=textoidl('Electron Energy [GeV]')

  cgplot, enel * 1d-3, icsedotcmbr, xtit='Electron Energy [GeV]', ytit='dE/dt [MeV]', thick=2, /xlog, /ylog, $
    pos=[0.2,0.05,0.9,0.95], yr=[1d-20, 1d1]
  cgplot, enel * 1d-3, puc_edot(ienel, 0, /cmbr), color=colors[1], /over


  restore_sysvars, state
  dfpsclose

end
