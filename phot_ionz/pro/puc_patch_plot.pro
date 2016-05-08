pro puc_patch_plot

  common puc_common

  colors = ['red5', 'grn5', 'blu5']

  pname = 'puc_patch_' + runtag + '.ps'
  dfpsplot, pname, /color, xs=8, ys=6
  state = sysvars(/print)
  !p.font = -1

  cgScatter2D, 1. + patchdensarr[*,iredp[0]], temppatchdm[*,iredp[0]], xtit=textoidl('1 + \delta'), ytit='T [K]', $
    color=colors[0], psym=14, fcolor=colors[0], /nodisplay, /xlog, /ylog, yr=[1d3,1d5], xr=[0.1,10.]

  cgScatter2D, 1. + patchdensarr[*,iredp[1]], temppatchdm[*,iredp[1]], $
    color=colors[1], psym=14, /over, fcolor=colors[1], /nodisplay

  cgScatter2D, 1. + patchdensarr[*,iredp[2]], temppatchdm[*,iredp[2]], $
    color=colors[2], psym=14, /over, fcolor=colors[2], /nodisplay


  cgplot, red, 1. + patchdensarr[0,*], xtit='Redshift', ytit=textoidl('1 + \delta'), /xlog
  for i=1, npatch-1 do cgplot, red, 1. + patchdensarr[i,*], /over
  cgplot, red, gwf, line=1, /over


; temperature ODE
 
  dtempdred1 = 2. / (1. + red)
  dtempdred2 = 2. / 3. / (1. + patchdensarr) * dpatchdensdzarr
  dtempdred3 = 2. * heat / 3. / boltzmann / (1. + red) / hubbf * htime

  cgplot, red, dtempdred2[0,*], xtit='Redshift', ytit=textoidl('dT/dz'), /xlog, /nodata, ystyle=4

  cgaxis, yaxis=0, tit=textoidl('dT/dz/T'), /save, yr=[min(dtempdred2) < min(dtempdred1),max(dtempdred2) > max(dtempdred1)]
  for i=0, npatch-1 do cgplot, red, dtempdred2[i,*], /over, color=colors[1]
  cgplot, red, dtempdred1, /over, color=colors[0]

  cgaxis, yaxis=1, tit=textoidl('dT/dz'), /save, yr=[min(dtempdred3),max(dtempdred3)]
  cgplot, red, dtempdred3, /over, colors=colors[2]


  restore_sysvars, state
  dfpsclose

end
