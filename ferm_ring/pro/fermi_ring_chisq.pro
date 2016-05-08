;+
; NAME:
;   fermi_ring_chisq
;
; PURPOSE:
;   get chisq for the NFW fit to the ring coefficients
;
; CALLING SEQUENCE:
;   fermi_ring_nfw_chisq, nfwpar
;
; INPUTS:
;   nfwpar    - free parameters for the NFW model
;   nfwpar[0] - normalization
;   nfwpar[1] - constant
;  
; OUTPUTS:
;   chisq - chisq for the NFW fit to the ring coefficients
;
; EXAMPLES:
;
; COMMENTS:
;   mpfit function minimizes the chisq
;   gradient is computed numerically
;   
; REVISION HISTORY:
;   2014-Jan-28 - written by Tansu Daylan
;
;----------------------------------------------------------------------
function fermi_ring_chisq, nfwpar

  common fermi_ring_chisq_common

; -------- average the NFW map with the current parameters

  nfwflux_ = nfwflux[a,*] * nfwpar[0]
  if keyword_set(nfwconstt) then nfwflux_ += nfwpar[1]

 
; -------- compute chisq assuming correlated ring flux errors

  resflux = rinflux[jringfit] - nfwflux_[jringfit]

  if keyword_set(correrr) then begin
    resfcin = resflux # rinfcin
    chisq_ = resfcin # resflux
  endif else begin
    chisq_ = total(resflux^2/errinflux[jringfit]^2)
  endelse

  return, chisq_
end
