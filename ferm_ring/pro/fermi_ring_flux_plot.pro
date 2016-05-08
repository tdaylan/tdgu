;+
; NAME:
;   fermi_ring_flux_plot
;
; PURPOSE:
;   goodness of fit plot for ring and averaged NFW fluxes
;
; CALLING SEQUENCE:
;  see fermi_ring_fit
;
; INPUTS:
;
; OUTPUTS:
;
; EXAMPLES:
;
; COMMENTS:
;   
; REVISION HISTORY:
;   2014-Jan-28 - written by Tansu Daylan
;
;----------------------------------------------------------------------
pro fermi_ring_flux_plot, psiavg, erpsiavg, rinflux, errinflux, $
      nfwflux, nring, chisq, runtag=runtag, difflux=difflux, erdifflux=erdifflux

  colors=['red5', 'grn5', 'blu5']


  pname = 'fermi_ring_flux_1_' + runtag + '.ps'
  dfpsplot, pname, /color, ys=6, xs=8
  state = sysvars(/print)
  !p.font = -1

  cgplot, psiavg[1:nring-2], rinflux[1:nring-2], $
    yrange=[1d-7,8d-5], xr=[2,11], /ylog, $
    ytitle=textoidl('E^2dN/dE [GeV/cm^2/s/sr]'), $
    xtitle=textoidl('\Psi [degrees]'), $
    err_yhigh=errinflux[1:nring-2], err_ylow=errinflux[1:nring-2], $
    err_xhigh=erpsiavg[1:nring-2], err_xlow=erpsiavg[1:nring-2], $
    psym=14, $
    pos=[0.2,0.1,0.9,0.9]

  if keyword_set(difflux) then begin
    cgplot, psiavg[1:nring-2], difflux[1:nring-2], /over, $
      err_yhigh=erdifflux[1:nring-2], err_ylow=erdifflux[1:nring-2], $
      err_xhigh=erpsiavg[1:nring-2], err_xlow=erpsiavg[1:nring-2], $
      psym=15
  endif

  cgplot, psiavg, nfwflux[0,*], color=colors[0], /over, thick=3
  cgplot, psiavg, nfwflux[1,*], color=colors[1], /over, thick=3
  cgplot, psiavg, nfwflux[2,*], color=colors[2], /over, thick=3

  cgLegend, Title=[ $
    textoidl('NFW fluxes, \gamma = 1.30, \chi^2 = ') + string(chisq[0], format='(F5.2)'), $
    textoidl('NFW fluxes, \gamma = 1.35, \chi^2 = ') + string(chisq[1], format='(F5.2)'), $
    textoidl('NFW fluxes, \gamma = 1.40, \chi^2 = ') + string(chisq[2], format='(F5.2)') $
    ], $
    color=colors, Location=[0.3,0.84], VSpace=2, thick=3

  restore_sysvars, state
  dfpsclose


  pname = 'fermi_ring_flux_2_' + runtag + '.ps'
  dfpsplot, pname, /color, ys=6, xs=8
  state = sysvars(/print)
  !p.font = -1

  cgplot, psiavg, rinflux, $
    yrange=[1d-7,8d-5], xr=[2,11], /ylog, $
    ytitle=textoidl('E^2dN/dE [GeV/cm^2/s/sr]'), $
    xtitle=textoidl('\Psi [degrees]'), $
    err_yhigh=errinflux, err_ylow=errinflux, $
    err_xhigh=erpsiavg, err_xlow=erpsiavg, $
    psym=14, $
    pos=[0.2,0.1,0.9,0.9]
  cgplot, psiavg, nfwflux[0,*], /over, thick=3

  restore_sysvars, state
  dfpsclose

  pname = 'fermi_ring_flux_3_' + runtag + '.ps'
  dfpsplot, pname, /color, ys=6, xs=8
  state = sysvars(/print)
  !p.font = -1

  cgplot, psiavg[1:nring-2], rinflux[1:nring-2], $
    yrange=[1d-7,8d-5], xr=[2,11], /ylog, $
    ytitle=textoidl('E^2dN/dE [GeV/cm^2/s/sr]'), $
    xtitle=textoidl('\Psi [degrees]'), $
    err_yhigh=errinflux[1:nring-2], err_ylow=errinflux[1:nring-2], $
    err_xhigh=erpsiavg[1:nring-2], err_xlow=erpsiavg[1:nring-2], $
    psym=14, $
    pos=[0.2,0.1,0.9,0.9]
  cgplot, psiavg, nfwflux[0,*], /over, thick=3

  restore_sysvars, state
  dfpsclose

end
