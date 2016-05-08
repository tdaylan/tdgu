;+
; NAME:
;   fermi_ring_nfw
;
; PURPOSE:
;   produce squared projected NFW profile map
;
; CALLING SEQUENCE:
;   map = fermi_ring_nfw(nfwgamma, sunnorm=sunnorm, nside=nside, makeplot=makeplot)
;
; INPUTS:
;  nfwgamma  - inner slope, gamma
;
;  OPTIONAL INPUTS:
;  sunnorm   - density normalization at the position of the Solar system
;  nside     - resolution parameter for HealPix
;
; KEYWORDS:
;  makeplot  - make a plot of the squared projected NFW map and radial density profile
;  
; OUTPUTS:
;  map       - NFW map
;
; OPTIONAL OUTPUTS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;  2014-Jan-28 - Written by Tansu Daylan
;----------------------------------------------------------------------
function fermi_ring_nfw, nfwgamma, sunnorm=sunnorm, nside=nside, makeplot=makeplot

  if not keyword_set(nside) then nside = 256L
  if not keyword_set(sunnorm) then sunnorm = 0.3; dark matter energy density in the Solar neighborhood [GeV/cm^3]
  if not keyword_set(radsun) then radsun = 8.5 ; [kpc]
  if not keyword_set(radsca) then radsca = 23.1 ; [kpc]


; -------- galactocentric radius

  if not keyword_set(nrad) then nrad = 100
  if not keyword_set(minrad) then minrad = 1d-2
  if not keyword_set(maxrad) then maxrad = 1d2
  rad = exp(findgen(nrad) * (alog(maxrad) - alog(minrad)) / (nrad - 1) + alog(minrad))
  

; -------- galactocentric radius

  if not keyword_set(minsad) then minsad = 0 ; [kpc]
  if not keyword_set(maxsad) then maxsad = 2. * radsun
  if not keyword_set(nsad) then nsad = 100
  sad = findgen(nsad) * (maxsad - minsad) / (nsad - 1) + minsad


; -------- galactocenteric angle

  npix = 12L * nside * nside
  healgen_lb, nside, l, b
  cpsi = cos(!DtoR * l) * cos(!DtoR * b)
  psi = !radeg * acos(cpsi)


; -------- calculate the mass density on the radial coordinate

  rho = 1. / (rad / radsca)^nfwgamma / (1. + rad / radsca)^(3. - nfwgamma)
  rho *= sunnorm / interpol(rho, rad, radsun)

  rhosli = dblarr(nsad, npix) 
  for i=0, nsad-1 do begin
    radsli = sqrt(radsun^2 + sad[i]^2 - 2 * radsun * sad[i] * cpsi)
    rhosli[i,*] = interpol(rho, rad, radsli)
  endfor
  map = total(rhosli^2, 1)


; -------- plot

  if keyword_set(makeplot) then begin
    pname = 'fermi_ring_nfw_' + string(nfwgamma, format='(F3.1)') + '.ps'

    dfpsplot, pname, /color
    state = sysvars(/print)

    cgplot, rad, rho, thick=5, title=title, ytitle=textoidl('\rho(r) [GeV/cm^3]'), $
      xtitle='r [kpc]', /ylog, /xlog
    cgtext, 0.3, 0.3, textoidl('\gamma: ') + string(nfwgamma,format='(F3.1)')

    healcart_display, map, [0,max(map)], title='Squared Projected Density', $
      lrange=[-20., 20.], brange=[-20., 20.]

    restore_sysvars, state
    dfpsclose
  endif

  return, map

end
