;+
; NAME:
;   fermi_ring_comp
;
; PURPOSE:
;   Comparison plots for different ring fit configurations
;
; CALLING SEQUENCE:
;   see fermi_ring_fit.pro
;
; OPTIONAL INPUTS:
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
;   2014-May-31 - Written by Tansu Daylan
;----------------------------------------------------------------------

pro comp_plot, titles=titles, runtags=runtags, cmints=cmints, ringpsi=ringpsi, colors=colors, dif=dif
  nvars = n_elements(runtags)
  dts = fltarr(nvars)
  for i=0, nvars-1 do begin
    d = mrdfits(runtags[i] + '.fits', 1)
    mints2 = d.mintsarr
    mints2 = mints2[0]
    dts[i] = mints2 - cmints
    color = colors[i]
    if keyword_set(dif) then begin
      single_plot, ringpsi, d.rinflux, d.errinflux, $
      difflux=d.difflux, erdifflux=d.erdifflux, energy, ienergy, color=color
      print, 'difflux: ', d.erdifflux
      print, 'rinflux: ', d.errinflux
    endif else begin
      single_plot, ringpsi, d.rinflux, d.errinflux, energy, ienergy, color=color
    endelse
    titles[i+1] += textoidl(', \Delta \chi^2 = ') + string(dts[i], format='(I10.0)')
  endfor
  colors = ['black', colors]
  line = 0
  if keyword_set(dif) then begin
    titles = [titles, 'FDM rings']
    colors = [colors, 'blue']
    line = [0,0,1]
  endif
  cgLegend, Title=titles, line=line, Location=[0.6,0.7], VSpace=2.0, Alignment=8, color=colors, thick=5, $
    chars=0.9
end

pro single_plot, ringpsi, rinflux, errinflux, energy, ienergy, $
      difflux=difflux, erdifflux=erdifflux, cplot=cplot, color=color, title=title

  ringx = reform(ringpsi[1,*])
  erringx = reform(ringpsi[1,*]-ringpsi[0,*])
  if keyword_set(cplot) then begin
    cgplot, ringx, rinflux, /ylog, psym=-14, $
      err_yhigh=errinflux, err_ylow=errinflux, $
      err_xhigh=erringx, err_xlow=erringx, $
      title=title, $
      yrange=[1d-7,8d-5], color=color, thick=5, $
      ytitle=textoidl('E^2 dN/dE [GeV/ph/cm^2/s/sr]'), xtitle=textoidl('\Psi [degrees]')
  endif else begin
  cgplot, ringx, rinflux, /ylog, psym=-14, $
    err_yhigh=errinflux, err_ylow=errinflux, thick=5, $
    err_xhigh=erringx, err_xlow=erringx, /overplot, color=color
  endelse
  if keyword_set(difflux) then begin
    color = 'blue'
    cgplot, ringx, difflux, $
      /ylog, /overplot, $
      err_yhigh=erdifflux, err_ylow=erdifflux, $
      err_xhigh=erringx, err_xlow=erringx, $
      psym=-15, color=color, line=1, thick=5
  endif
end

pro fermi_ring_comp

  ienergy = 9

  dfpsplot, 'fermi_ring_variations.ps', /color, /landscape, ysize=4
  state = sysvars(/print)
  !p.font = -1
  colors = ['red', 'orange', 'yellow', 'green', 'blue', 'violet']

  ;get control run
  fermi_ring_fit, /comfwhm, /getruntag, runtag=cruntag
  d = mrdfits(cruntag + '.fits', 1)
  cerrinflux = d.errinflux
  crinflux = d.rinflux
  cmints = d.mintsarr
  cmints = cmints[0]
  ringpsi = d.ringpsi
  energy = d.energy

  ;FDM version
  runtags = ['','']
  fermi_ring_fit, /comfwhm, /getruntag, fdmver=2, runtag=runtag
  runtags[0] = runtag  
  fermi_ring_fit, /comfwhm, /getruntag, fdmver=4, runtag=runtag
  runtags[1] = runtag  
  title = 'Variation of ring fluxes wrt the FDM version'
  titles = ['P6V11', 'P6V3', 'P7V6']
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

  ;DM brems
  fermi_ring_fit, /comfwhm, /getruntag, /dmbrems, runtag=runtags
  title = 'Variation of ring fluxes with the addition of a DM brems template'
  titles = [ 'No DM brems template', 'With DM brems template']
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

  ;high energy
  runtags = ['']
  fermi_ring_fit, /comfwhm, /getruntag, maxe=17, runtag=runtag
  runtags[0] = runtag  
  title = 'Variation of ring fluxes wrt the maximum energy of the fit'
  titles = [ $
    '0.3 GeV - ' + string(energy[19], format='(F4.1)') + ' GeV', $
    '0.3 GeV - ' + string(energy[17], format='(F4.1)') + ' GeV' $
    ]
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  print, runtags
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

  ;low energy
  runtags = ['']
  fermi_ring_fit, /comfwhm, /getruntag, mine=2, runtag=runtag
  runtags[0] = runtag  
  title = 'Variation of ring fluxes wrt the minimum energy of the fit'
  titles = [ $
    string(energy[0], format='(F4.1)') + ' GeV - 30 GeV', $
    string(energy[2], format='(F4.1)') + ' GeV - 30 GeV' $
    ]
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

  ;FDM rings
  fermi_ring_fit, /comfwhm, /getruntag, /doubring, runtag=runtags
  title = 'Variation of ring fluxes with the addition of FDM spectrum rings'
  titles = [ 'Rings with no FDM rings', 'Rings with FDM rings']
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1], /dif

  ;fixed spectrum variation
  runtags = ['','','']
  fermi_ring_fit, /comfwhm, /getruntag, fixedspec='dif', runtag=runtag
  runtags[0] = runtag  
  fermi_ring_fit, /comfwhm, /getruntag, fixedspec='pul', runtag=runtag
  runtags[1] = runtag  
  fermi_ring_fit, /comfwhm, /getruntag, fixedspec='pow', runtag=runtag
  runtags[2] = runtag  
  title = 'Variation of ring fluxes wrt the fixed spectrum on the rings'
  titles = [ $
    textoidl('NFW spectrum fixed rings'), $
    textoidl('FDM spectrum fixed rings'), $
    textoidl('Pulsar spectrum fixed rings'), $
    textoidl('Power law spectrum fixed rings') $
    ]
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

  ;spectral offset
  runtags = ['','','','']
  fermi_ring_fit, /comfwhm, /getruntag, specoffset=-5., runtag=runtag
  runtags[0] = runtag  
  fermi_ring_fit, /comfwhm, /getruntag, specoffset=-1., runtag=runtag
  runtags[1] = runtag  
  fermi_ring_fit, /comfwhm, /getruntag, specoffset=1., runtag=runtag
  runtags[2] = runtag  
  fermi_ring_fit, /comfwhm, /getruntag, specoffset=5., runtag=runtag
  runtags[3] = runtag  
  title = 'Variation of ring fluxes wrt the additive offset to the spectrum'
  titles = [ $
    textoidl('No offset'), $
    textoidl('E^2dN/dE -= 5 \times 10^{-6} [GeV/cm^2/s/str]'), $
    textoidl('E^2dN/dE -= 1 \times 10^{-6} [GeV/cm^2/s/str]'), $
    textoidl('E^2dN/dE += 1 \times 10^{-6} [GeV/cm^2/s/str]'), $
    textoidl('E^2dN/dE += 5 \times 10^{-6} [GeV/cm^2/s/str]') $
    ]
  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

  ;point source subtraction radius
;  runtags = ['','']
;  fermi_ring_fit, /comfwhm, /getruntag, pscfac=0.75, runtag=runtag
;  runtags[0] = runtag  
;  fermi_ring_fit, /comfwhm, /getruntag, pscfac=1.25, runtag=runtag
;  runtags[1] = runtag  
;  title = 'Variation of ring fluxes wrt changes in point subtraction radius'
;  titles = ['120 arcmin', '090 arcmin', '150 arcmin']
;  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
;  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]
;
;  ;number bubble slices
;  runtags = ['','','']
;  fermi_ring_fit, /comfwhm, /getruntag, nslice=1, runtag=runtag
;  runtags[0] = runtag  
;  fermi_ring_fit, /comfwhm, /getruntag, nslice=2, runtag=runtag
;  runtags[1] = runtag  
;  fermi_ring_fit, /comfwhm, /getruntag, nslice=7, runtag=runtag
;  runtags[2] = runtag  
;  title = 'Variation of ring fluxes wrt the number of bubble latitude slices'
;  titles = ['5 slices', 'no slice', '2 slices', '7 slices']
;  single_plot, ringpsi, crinflux, cerrinflux, energy, ienergy, /cplot, title=title
;  comp_plot, runtags=runtags, titles=titles, cmints=cmints, ringpsi=ringpsi, colors=colors[0:n_elements(runtags)-1]

;  fitreg       - 'no': northern sky only
;  version      - data reconstruction version
;  nohighbubspec- set to use the bubble spectrum at high latitudes (abs(b)=40-50 deg)

  restore_sysvars, state
  dfpsclose
  
end
