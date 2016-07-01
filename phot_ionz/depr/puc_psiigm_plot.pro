pro puc_psiigm_plot

  common puc_common

  pname = 'puc_psiigm_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  colors = ['red5','grn5','blu5']

  cgplot, enel * 1d-3, enel^2 * psiigm[*,iredp[0]], /xlog, /ylog, $
    color=colors[0], pos=[0.2,0.05,0.95,0.95], xr=[enel[0]*1d-3,enel[nenel-2]*1d-3], $
    thick=5, ytitle=textoidl('\Psi(E,z) [MeV/cm^3]'), xtitle=textoidl('Electron Energy [GeV]')
  cgplot, enel * 1d-3, enel^2 * psiigm[*,iredp[1]], /overplot, color=colors[1], thick=5
  cgplot, enel * 1d-3, enel^2 * psiigm[*,iredp[2]], /overplot, color=colors[2], thick=5

  cgLegend, Title=[ $
    'z = ' + string(red[iredp[0]],format='(F3.1)'), $
    'z = ' + string(red[iredp[1]],format='(F3.1)'), $
    'z = ' + string(red[iredp[2]],format='(F3.1)')], $
    Location=[0.4,0.25], VSpace=2.0, Alignment=8, thick=5, chars=1.2, color=colors

  puc_contour, psiigm * (enel^2 # (1. + dblarr(nred))), enel*1d-3, red, /xlog, /ylog, /log, $
    xtit=textoidl('Electron Energy [GeV]'), ytit=textoidl('Redshift')

  puc_contour, psiigmsst * (enel^2 # (1. + dblarr(nred))), enel*1d-3, red, /xlog, /ylog, /log, $
    xtit=textoidl('Electron Energy [GeV]'), ytit=textoidl('Redshift')

  restore_sysvars, state
  dfpsclose

end
