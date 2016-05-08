pro puc_sigma2_plot

  common puc_common

  colors = ['red5', 'grn5', 'blu5']

  nk = 1d6
  ds = [0, nmass/2, nmass-1]
  nd = n_elements(ds)
  dcount = 0
  sigma2int = dblarr(nk, nd)
  k = dblarr(nk, nd)
  minlimk = dblarr(nd)
  maxlimk = dblarr(nd)
  mink = dblarr(nd)
  maxk = dblarr(nd)
  foreach d, ds do begin
    minlimk[dcount] = 1d-15 / rmass[d]
    maxlimk[dcount] = 1d2 / rmass[d]
    mink[dcount] = minlimk[dcount] * 1d-5
    maxk[dcount] = maxlimk[dcount] * 1d5
    k[*,dcount] = exp(dindgen(nk) / (nk - 1d) * (alog(maxk[dcount]) - alog(mink[dcount])) + alog(mink[dcount]))
    sigma2int[*,dcount] = puc_sigma2_integrand(k[*,dcount])
    dcount++
  endforeach

  minkcom = min(k)
  maxkcom = max(k)
  kcom = exp(dindgen(nk) / (nk - 1d) * (alog(maxkcom) - alog(minkcom)) + alog(minkcom))
  k2pk = kcom^2 * puc_powerspec(kcom)
  pname = 'puc_sigma2_' + runtag + '.ps'

  dfpsplot, pname, /color, ysize=6, xsize=8
  state = sysvars(/print)
  !p.font = -1

; -------- sigma & power spectrum

  multiplot, [2,1], gap=0.05

  cgplot, massp, sigmam, /xlog, /ylog, thick=5, xr=[1d-10,1d15], $
    ytitle=textoidl('\sigma(M_{200})'), xtitle=textoidl('M_{200} [M_{sun}]'), pos=[0.15,0.05,0.45,0.95]
  cgplot, massp, fltarr(nmassp)+deltaz[iredp[0]], /overplot, thick=5, line=2, color=colors[0]
  cgplot, massp, fltarr(nmassp)+deltaz[iredp[1]], /overplot, thick=5, line=2, color=colors[1]
  cgplot, massp, fltarr(nmassp)+deltaz[iredp[2]], /overplot, thick=5, line=2, color=colors[2]

; -------- sigma2int

  multiplot

  cgplot, k[*,0], sigma2int[*,0], xtit='k [1/Mpc]', /xlog, yr=[1d-16,1d2], xr=[1d-6,1d10], $
    /ylog, ytit=textoidl('d\sigma^2/dk=k^2P(k)W^2(kR)'), /nodata, pos=[0.55,0.05,1,0.95]
  ik = where(k[*,0] gt 1d1)
  cgplot, k[ik,0], sigma2int[ik,0], color=colors[0], thick=3, /over
  ik = where(k[*,1] gt 1d-3)
  cgplot, k[ik,1], sigma2int[ik,1], color=colors[1], thick=3, /over
  cgplot, k[*,2], sigma2int[*,2], color=colors[2], thick=3, /over
  cgplot, kcom, k2pk, line=1,  thick=3, /over

  multiplot, /reset

  cgLegend, Title=[textoidl('\sigma(M)'), $
    textoidl('\delta(z = ') + string(red[iredp[0]],format='(F3.1)') + ')', $
    textoidl('\delta(z = ') + string(red[iredp[1]],format='(F3.1)') + ')', $
    textoidl('\delta(z = ') + string(red[iredp[2]],format='(F3.1)') + ')'], $
    line=[0,2,2,2], $
    Location=[0.3,0.22], VSpace=2.0, Alignment=8, thick=5, color=['black',colors]

  cgLegend, Title=[ $
    textoidl('k^2P(k)W^2(kR),R=') + string(rmass[ds[0]]*1d6,format='(F5.2)') + ' pc', $
    textoidl('k^2P(k)W^2(kR),R=') + string(rmass[ds[1]]*1d3,format='(F5.2)') + ' Kpc', $
    textoidl('k^2P(k)W^2(kR),R=') + string(rmass[ds[2]],format='(F5.2)') + ' Mpc', $
    textoidl('k^2P(k)')], $
    Location=[0.78,0.7], VSpace=2.0, Alignment=8, color=[colors, 'black'], $
    line=[0,0,0,1], thick=3, len=0.02



  restore_sysvars, state
  dfpsclose

end
