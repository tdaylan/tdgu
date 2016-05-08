;+
; NAME:
;   fermi_ring_spec_plot
;
; PURPOSE:
;   produce correlation spectra plots
;
; CALLING SEQUENCE:
;   fermi_ring_spec_plot, xx, xxerr, ringflux, errringflux, nring, energy, runtag, nfwgamma, bubspec
;
; INPUTS:
;  see fermi_ring_profile.pro
;
; KEYWORDS:
;  
; OUTPUTS:
;  correlation spectra plot
;
; EXAMPLES:
;
; COMMENTS:
; 
; REVISION HISTORY:
;   2014-Jan-28 - Written by Tansu Daylan
;----------------------------------------------------------------------

pro fermi_ring_spec_plot, energy, nfwspec, ernfwspec, runtag, $
    bubspec, difspec, isospec, pulspec, powspec, dmaspec

  nenergy= n_elements(energy)
  specs = fltarr(7, nenergy)
  specs[0,*] = nfwspec
  specs[1,*] = pulspec
  specs[2,*] = bubspec
  specs[3,*] = powspec
  specs[4,*] = difspec
  specs[5,*] = isospec
  specs[6,*] = dmaspec

  dfpsplot, 'fermi_ring_spectra_' + runtag + '.ps', /color, bits=8
  state = sysvars(/print)
  colors = ['red', 'orange', 'yellow', 'green', 'blue']

  cgplot, energy, specs[0,*] * energy^2, $
    color=colors[0], /ylog, /xlog, psym=-3, thick=5, $
    xrange=[0.3,300], yrange=[1d-8,1d-4], $
    xtitle='Energy [GeV]', ytitle=textoidl('E^2 dN/dE [GeV/cm^2/s/sr]'), $
    ERR_YLOW=ernfwspec * energy^2, ERR_YHIGH=ernfwspec * energy^2, $
    err_thick=5

  for i=1, 4 do begin
    cgplot, energy, specs[i,*] * energy^2, $
      color=colors[i], /ylog, /xlog, psym=-3, thick=5, /overplot
  endfor
 
  titles = ['NFW', 'Pulsar', 'Bubble', 'Power Law', 'FDM']
  cgLegend,  title=titles, color=colors, Location=[0.7,0.74], VSpace=2.0, Alignment=8, charsize=0.8, $
             thick=5

  restore_sysvars, state
  dfpsclose
end
