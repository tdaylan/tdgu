;+
; NAME:
;   fermi_ring_design
;
; PURPOSE:
;   construct the design matrix for the template fit
;
; CALLING SEQUENCE:
;   amatrix = fermi_ring_design(str, nenergy=nenergy, ipix=ipix, npixfit=npixfit, npixcum=npixcum)
;
; INPUTS:
;   str      - template structure
;
; OPTIONAL INPUTS:
;   nenergy  - number of energy bins in the fit
;   ipix     - list of pixels to be fitted for each energy
;   npixfit  - number of pixels in each energy bin
;   npixcum  - cumulative number of pixels up to each energy bin
;
; OUTPUTS:
;   amatrix  - design matrix
;
; EXAMPLES:
;   see fermi_ring_fit.pro
;
; COMMENTS:
;
;   
; REVISION HISTORY:
;   2013-Oct-24 - Written by DPF, SKNP, TD
;   2014-Jan-28 - Modified by Tansu Daylan
;   2014-Aug-14 - Modified by Tansu Daylan
;
;---------------------------------------------------------------------
function fermi_ring_design, str, nenergy=nenergy, ipix=ipix, npixfit=npixfit, npixcum=npixcum

  if keyword_set(singener) then nenergy = 1
  if keyword_set(spectemp) then begin
    npixfit = intarr(nenergy)
    npixcum = intarr(nenergy)
    ipix = list()
    for i=0, nenergy-1 do begin
      ipix.add, 0
      npixfit[i] = 1
      npixcum[i] = i
    endfor
  endif

  nrow = total(npixfit)
  ncol = total(str.fixed) + total(str.varfixed) + total(str.float) * nenergy + total(str.varfloat) * nenergy
  Amatrix = fltarr(ncol, nrow)
  useit = str.fixed or str.float or str.varfixed or str.varfloat
  i = 0L ; column of A 
  k = 0L ; template number
  while (i LT Ncol) do begin 
    while useit[k] EQ 0 do k++
      if str[k].fixed then begin 
        spec = str[k].spectrum
        for j=0L, nenergy-1 do begin
          thistemplate = str[k].template[ipix[j]]
          Amatrix[i, npixcum[j] : npixcum[j]+npixfit[j]-1] = thistemplate * spec[j]
        endfor
        i++
      endif
      if str[k].float then begin
        for j=0L, Nenergy-1 do begin 
          thistemplate = str[k].template[ipix[j]]
          Amatrix[i, npixcum[j] : npixcum[j]+npixfit[j]-1] = thistemplate
          i++
        endfor
      endif
      if str[k].varfixed then begin
        spec = str[k].spectrum
        for j=0L, Nenergy-1 do begin
          thistemplate = str[k].templatearr[ipix[j],j]
          Amatrix[i, npixcum[j] : npixcum[j]+npixfit[j]-1] = thistemplate * spec[j]
        endfor
 	i++
      endif
      if str[k].varfloat then begin
        for j=0L, Nenergy-1 do begin
          thistemplate = str[k].templatearr[ipix[j],j]
          Amatrix[i, npixcum[j] : npixcum[j]+npixfit[j]-1] = thistemplate
 	  i++
        endfor
      endif
    k++
  endwhile
  return, Amatrix
  end
