;+
; NAME:
;   fermi_ring_map
;
; PURPOSE:
;   extract maps from the fit results
;
; CALLING SEQUENCE:
;  see fermi_ring_map_plot.pro
;
; INPUTS:
;
; KEYWORDS:
;  
; OUTPUTS:
;
; EXAMPLES:
;
; COMMENTS:
; 
; REVISION HISTORY:
;   2014-Jan-28 - Written by Tansu Daylan
;----------------------------------------------------------------------
pro fermi_ring_map, epsmap=epsmap, datmap=datmap, allmap=allmap, rinmap=rinmap, ri2map=ri2map, isomap=isomap, difmap=difmap, $
      bubmap=bubmap, dmbmap=dmbmap, modmap=modmap, resmap=resmap, xx=xx, cts=cts, Amatrix=Amatrix, $
      ex=ex, energy=energy, ipix=ipix, npixfit=npixfit, npixcum=npixcum, nring=nring, ienergy=ienergy, $
      doubring=doubring, dmbrems=dmbrems, nslice=nslice, fac=fac, fdmver=fdmver

    nside = 256L
    npix = nside^2 * 12

    nenergy = n_elements(energy)

    allmap = fltarr(npix, nring)
    epsmap = fltarr(npix)
    datmap = fltarr(npix)
    rinmap = fltarr(npix)
    ri2map = fltarr(npix)
    isomap = fltarr(npix)
    difmap = fltarr(npix)
    bubmap = fltarr(npix)
    dmbmap = fltarr(npix)

    ie = npixcum[ienergy]
    ie_ = npixcum[ienergy] + npixfit[ienergy] - 1
    
    datmap[ipix[ienergy]] = cts[ie:ie_] / ex[ie:ie_] * fac
    epsmap[ipix[ienergy]] = ex[ie:ie_]
    
    strocc = 0
    while strocc lt nring do begin
      rinmap[ipix[ienergy]] += Amatrix[strocc, ie:ie_] * xx[strocc] * fac
      for i=0, strocc do allmap[ipix[ienergy],strocc] += Amatrix[i, ie:ie_] * xx[i] * fac
      strocc++
    endwhile
 
    if keyword_set(doubring) then begin
      while strocc lt 2*nring do begin
        ri2map[ipix[ienergy]] += Amatrix[strocc, ie:ie_] * xx[strocc] * fac
        strocc++
      endwhile
    endif
 
    if not keyword_set(fixdifiso) then begin 
      isomap[ipix[ienergy]] = Amatrix[strocc+ienergy, ie:ie_] * xx[strocc+ienergy] * fac
      strocc += nenergy
      difmap[ipix[ienergy]] = Amatrix[strocc+ienergy, ie:ie_] * xx[strocc+ienergy] * fac
      strocc += nenergy
    endif else begin
      isomap[ipix[ienergy]] = Amatrix[strocc, ie:ie_] * xx[strocc] * fac
      strocc++
      difmap[ipix[ienergy]] = Amatrix[strocc, ie:ie_] * xx[strocc] * fac
      strocc++
    endelse

    if not (fdmver eq 4) then begin
      for j=0, nslice-1 do $ 
        bubmap[ipix[ienergy]] += Amatrix[strocc+j, ie:ie_] * xx[strocc+j] * fac
      strocc += nslice
    endif

    if keyword_set(dmbrems) then $
      dmbmap[ipix[ienergy]] = Amatrix[strocc, ie:ie_] * xx[strocc] * fac

    datmap *= energy[ienergy]^2
    rinmap *= energy[ienergy]^2
    ri2map *= energy[ienergy]^2
    isomap *= energy[ienergy]^2
    difmap *= energy[ienergy]^2
    bubmap *= energy[ienergy]^2
    dmbmap *= energy[ienergy]^2

    modmap = rinmap + ri2map + isomap + difmap + bubmap + dmbmap
    resmap = datmap - modmap

end
