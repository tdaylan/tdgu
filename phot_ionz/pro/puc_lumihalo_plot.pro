pro puc_lumihalo_plot

  common puc_common

  pname = 'puc_lumihalo_' + runtag + '.ps'
  dfpsplot, pname, /color, ysize=6, xsize=8
  state = sysvars(/print)
  !p.font = -1

  colors = ['black','black','black','blu5','blu5']

  cgplot, mass, lumihaloint[*,0], xtitle=textoidl('M_{200} [M') + sunsymbol() + ']', $
    thick=5, pos=[0.15,0.05,0.85,0.95], /xlog, ystyle=4, /nodata

  cgaxis, yaxis=0, tit=textoidl('L_{halo}(M_{200}) [MeV/s]'), yr=[1d27,1d50], /save, /ylog
  cgplot, mass, lumihaloint[*,0], thick=5, /over
  cgplot, mass, lumihalointne[*,0], thick=5, line=1, /over
  cgplot, mass, lumihalointnsne[*,0], thick=5, /over, line=2

  jenel = nenel/2
  cgaxis, /yaxis, tit=textoidl('f_{esc}'), yr=[min(escfrac),1./min(escfrac)], /save, color=colors[3], $
    yticks=2, ytickv=[min(escfrac),sqrt(min(escfrac)),1.], /ylog
  cgplot, mass, escfrac[jenel,*,iredp[0]], thick=5, color=colors[3], /over
  cgplot, mass, escfrac[jenel,*,iredp[2]], thick=5, color=colors[4], /over, line=1

  cgLegend, Title=[ $
    textoidl('L_H'), $
    textoidl('L_H, f_{esc} = 1'), $
    textoidl('L_H, no substructure'), $
    textoidl('f_{esc} (') + string(enel[jenel]/1d3,format='(I1.1)') + ' GeV,z=' + string(red[iredp[0]], format='(I1.1)') + ')', $
    textoidl('f_{esc} (') + string(enel[jenel]/1d3,format='(I1.1)') + ' GeV,z=' + string(red[iredp[2]], format='(I1.1)') + ')'], $
    Location=[0.25,0.8], VSpace=2.0, thick=5, color=colors, line=[0,1,2,0,1]

  restore_sysvars, state
  dfpsclose

end
