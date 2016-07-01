pro puc_opeff_plot

  common puc_common

  pname = 'puc_opeff_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  colors = ['red5','grn5','blu5']

  puc_contour, (cdd # (1. + fltarr(nred))) * d2fdndz, cdd, red, /xlog, /ylog, $
    xtit=textoidl('N_{HI} [1/cm^2]'), tit=textoidl('log(d^2f/dN_{HI}dz)'), ytit=textoidl('Redshift')
  vline, 10^15, line=1, color='grey', thick=5, /log

  puc_contour, (cdd # (1. + fltarr(nred)) * optotint[*,0,*]), cdd, red, /xlog, /ylog, /log, $
    xtit=textoidl('N_{HI} [1/cm^2]'), tit=textoidl('log(N_{HI}d^2f/dlogN_{HI}dz(1-exp(N_{HI}\sigma_{ion})))'), ytit=textoidl('Redshift')

  puc_contour, optot, enph * 1d6, red, /xlog, /ylog, /log, $
    xtit=textoidl('Photon Energy [eV]'), tit=textoidl('d\tau/dz(E,z)'), ytit=textoidl('Redshift')

  puc_contour, opeff[*,iredp[0]+1:nred-1,iredp[0]], enph * 1d6, red[iredp[0]+1:nred-1], /log, /xlog, /ylog, $
    xtit=textoidl('Photon Energy [eV]'), tit=textoidl('tau_{eff}(E,z) at z_{obs} = ') + string(iredp[0], format='(G4.2)'), ytit=textoidl('Redshift')
 
  puc_contour, opeff[*,iredp[1]+1:nred-1,iredp[1]], enph * 1d6, red[iredp[1]+1:nred-1], /log, /xlog, /ylog, $
    xtit=textoidl('Photon Energy [eV]'), tit=textoidl('tau_{eff}(E,z) at z_{obs} = ') + string(iredp[1], format='(G4.2)'), ytit=textoidl('Redshift')
 
  restore_sysvars, state
  dfpsclose

end
