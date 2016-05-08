;+
; NAME:
;   fermi_ring_template
;
; PURPOSE:
;   Make ring templates for spatial profile fits
;
; CALLING SEQUENCE:
;   rings=fermi_ring_template(nring, rmin, dr, nside=, lcen=, bcen=, $
;                             fwhm=, bcut=)
;
; INPUTS:
;   nring      - number of rings
;   rmin       - inner ring radius [deg]
;   dr         - ring width [deg]
;
; OPTIONAL INPUTS:
;   nside      - healpix nside
;   lcen,bcen  - (l,b) of center
;   fwhm       - fwhm for smoothing [arcmin]
;   bcut       - latitude cut [deg]
;
; KEYWORDS:
;   write      - set to write the template array
;
; OUTPUTS:
;   rings      - (npix, nring) array of ring templates
;
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;   2013-Oct-24 - Written by Douglas Finkbeiner, CfA
;
;----------------------------------------------------------------------
function fermi_ring_template, nring, rmin, dr, nside=nside, lcen=lcen, $
         bcen=bcen, fwhm=fwhm, bcut=bcut, write=write

  if ~keyword_set(nside) then nside = 256L
  if n_elements(lcen) EQ 0 then lcen = 0.0
  if n_elements(bcen) EQ 0 then bcen = 0.0

  if ~keyword_set(nring) then message, 'must set nring'
  if ~keyword_set(rmin)  then message, 'must set rmin'
  if ~keyword_set(dr)    then message, 'must set dr'

; -------- angle to reference point
  healgen_lb, nside, l, b
  ang = djs_diff_angle(l, b, lcen, bcen)

; -------- allocate output array
  npix = 12L*nside*nside
  rings = dblarr(npix, nring)
  
; -------- ring index
  ringind = floor((ang-rmin)/dr)

; -------- longitude cut
  if keyword_set(bcut) then begin 
     bgood = abs(b) GT bcut
  endif else bgood = 1

; -------- define rings
  w = where(ringind GE 0 and ringind LT nring and bgood, nw)
  for i=0L, nw-1 do rings[w[i], ringind[w[i]]]++

; -------- smooth if necessary
  if keyword_set(fwhm) then begin 
     lmax = (180.*60*4/fwhm) < (3*nside)
     sm = heal_smooth(rings, fwhm, lmax=lmax)
     wsmall = where(sm LT 1E-4, nsmall)
     if nsmall GT 0 then sm[wsmall] = 0
     rings = sm
     fwhmstr = string(fwhm, format='(I3.3)')
  endif else begin
    fwhmstr = '000'
  endelse
  id = string(nring, format='(I02)') + '_' + string(rmin, format='(F3.1)') $
  + '_' + string(dr, format='(F3.1)')
  if keyword_set(write) then writefits, 'ring_template_' + id + '_fwhm' + fwhmstr + '.fits', rings


  return, rings
end


pro tryit
  
  rings = fermi_ring_template(10, 5, 2, fwhm=120)

  return
end
