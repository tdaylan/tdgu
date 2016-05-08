;+
; NAME:
;   fermi_ring_map_plot
;
; PURPOSE:
;   produce counts plots for all fit components
;
; CALLING SEQUENCE:
;   fermi_counts_plot, xx, cts, Amatrix, ex, energy, ipix, npixfit, npixcum, runtag, nring
;
; INPUTS:
;  see fermi_ring_fit.pro
;
; KEYWORDS:
;  
; OUTPUTS:
;  count (data and model) sky maps 
;
; EXAMPLES:
;
; COMMENTS:
; 
; REVISION HISTORY:
;   2014-Jan-28 - Written by Tansu Daylan
;----------------------------------------------------------------------

pro fermi_ring_map_plot, xx, cts, Amatrix, ex, energy, ipix, npixfit, npixcum, runtag, nring, $
    doubring=doubring, dmbrems=dmbrems, nslice=nslice, fdmver=fdmver, ienergy=ienergy

  nside=256L
  npix = nside * nside * 12L

  lr = 15.
  br = 15.

  nenergy = n_elements(energy)
  dfpsplot, 'fermi_ring_map_' + runtag + '.ps', bits=8
  state = sysvars(/print)
  !p.font = -1

  i = 0
  while i lt nenergy do begin
    !p.multi = [12,3,4]
    fermi_ring_map, epsmap=epsmap, datmap=datmap, allmap=allmap, rinmap=rinmap, ri2map=ri2map, isomap=isomap, difmap=difmap, $
      bubmap=bubmap, dmbmap=dmbmap, modmap=modmap, resmap=resmap, xx=xx, cts=cts, Amatrix=Amatrix, $
      ex=ex, energy=energy, ipix=ipix, npixfit=npixfit, npixcum=npixcum, nring=nring, ienergy=i, $
      doubring=doubring, dmbrems=dmbrems, nslice=nslice, fac=1d4, fdmver=fdmver

    healcart_display, rinmap, [min(rinmap),max(rinmap)], $
      title='DM Ring', $
      lrange=[-lr, lr], brange=[-br, br]
    healcart_display, ri2map, [min(ri2map),max(ri2map)], $
      title='FDM Ring ' + string(energy[i],format='(F5.1)') + ' GeV, ' + textoidl('[10^{-4} GeV/cm^2/s/sr]'), $
      lrange=[-lr, lr], brange=[-br, br]
    healcart_display, isomap, [min(isomap),max(isomap)], $
      title='Isotropic', $
      lrange=[-lr, lr], brange=[-br, br]

    healcart_display, difmap, [min(difmap),max(difmap)], $
      title='FDM', $
      lrange=[-lr, lr], brange=[-br, br]
    healcart_display, bubmap, [min(bubmap),max(bubmap)], $
      title='Bubble', $
      lrange=[-lr, lr], brange=[-br, br]
    healcart_display, dmbmap, [min(dmbmap),max(dmbmap)], $
      title='DM Brems', $
      lrange=[-lr, lr], brange=[-br, br]

    healcart_display, datmap, [min(datmap),max(datmap)], $
      title='Data', $
      lrange=[-lr, lr], brange=[-br, br]
    healcart_display, epsmap, [min(epsmap),max(epsmap)], $
      title='Exposure', $
      lrange=[-lr, lr], brange=[-br, br]
    healcart_display, resmap, [min(resmap),max(resmap)], $
      title='Residual', $
      lrange=[-lr, lr], brange=[-br, br]

    !p.multi = 0
    i++
  endwhile

  fermi_ring_map, datmap=datmap, allmap=allmap, rinmap=rinmap, ri2map=ri2map, isomap=isomap, difmap=difmap, $
    bubmap=bubmap, dmbmap=dmbmap, modmap=modmap, resmap=resmap, xx=xx, cts=cts, Amatrix=Amatrix, $
    ex=ex, energy=energy, ipix=ipix, npixfit=npixfit, npixcum=npixcum, nring=nring, ienergy=ienergy, $
    doubring=doubring, dmbrems=dmbrems, nslice=nslice, fac=1d4, fdmver=fdmver

  for i=0, nring-1 do $
    healcart_display, allmap[*,i], [min(allmap[*,i]),max(allmap[*,i])], $
      title='', lrange=[-lr, lr], brange=[-br, br]

  rings = fermi_ring_template(9,2.5,1.)
  healcart_display, rings[*,4], [min(rings[*,4]),max(rings[*,4])], title='', $
    lrange=[-10., 10.], brange=[-10., 10.]
  rings = heal_smooth(rings, 60.) 
  healcart_display, rings[*,4], [min(rings[*,4]),max(rings[*,4])], title='', $
    lrange=[-10., 10.], brange=[-10., 10.]

  nfwtemp = fermi_nfw_map(1.3)
  healcart_display, nfwtemp, [min(nfwtemp),max(nfwtemp)], title='', $
    lrange=[-10., 10.], brange=[-10., 10.]
  nfwtemp = heal_smooth(nfwtemp, 120.)
  healcart_display, nfwtemp, [min(nfwtemp),max(nfwtemp)], title='', $
    lrange=[-10., 10.], brange=[-10., 10.]

  healcart_display, nfwtemp*rings[*,4], [min(nfwtemp*rings[*,4]),max(nfwtemp*rings[*,4])], title='', $
    lrange=[-10., 10.], brange=[-10., 10.]

  restore_sysvars, state
  dfpsclose
end
