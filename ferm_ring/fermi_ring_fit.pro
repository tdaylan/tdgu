;+
; NAME:
;   fermi_ring_fit
;
; PURPOSE:
;   Ring-Fit Analysis Main Routine
;   generates concentric ring templates, performs a maximum Poisson likelihood fit and
;   then fits the resulting ring coefficients with NFW predictions
;
; CALLING SEQUENCE:
;    fermi_ring_fit, nring=nring, rin=rin, minpsi, $
;         $
;    mine=mine, maxe=maxe, $
;    psimax=psimax, fitreg=fitreg, bcut=bcut, nsigma=nsigma, $
;    fastrings=fastrings, ienergy=ienergy, $
;    fixedspec=fixedspec, mints=mints, nfwconst=nfwconst, $
;    version=version
;
; OPTIONAL INPUTS:
;  nring        - number of ring templates
;  rin          - inner radius of the innermost ring [degrees]
;  dr           - width of a ring [degrees]
;  mine         - minimum energy bin index
;  maxe         - maximum energy bin index
;  psimax       - maximum angle from the GC to be included in the fit
;  fitreg       - 'cc': 40X40 degree^2 ROI
;                 'al': whole sky
;  bcut         - galactic plane mask [degrees]
;  fixedspec    - the fixed spectrum for the rings
;               - 'nfw': data-driven NFW correlated emission spectrum
;               - 'pul': fit to the pulsar spectrum
;               - 'dif': diffuse model spectrum
;               - 'pow': power law spectrum
;               - 'dma': PYTHIA generated dark matter annihilation spectrum
;  version      - data reconstruction version
;               - 'p11': p11 (raw data)
;               - 'p11q2': p11 Q2 cut (best 50%)
;  fdmver       - Fermi diffuse model version
;               - 2: v2 (gll_iem_v02.fit)
;               - 3: v3 (gll_iem_v02_P6_V11_DIFFUSE.fit)
;               - 4: v4 (gal_2yearp7v6_v0.fits)
;               - 5: v5 (gll_iem_v05.fits) 
;  specoffset   - additive offset to the data-driven NFW spectrum [cts/cm^2/s/sr/KeV]
;  nslice       - number of bubble slices
;
;  ********
;  Inputs not reflected in the run tag
;  Be extra careful when varying these knobs! They are only included for diagnostic purposes.
;
;  nsigma       - maximum sigma for linear regression sigma clipping
;  ienergy      - index of the energy bin to plot the ring fluxes
;  ********
;
; KEYWORDS:
;  nfwconst     - keyword to allow an additive constant when fitting the
;                 NFW template ring fluxes to the data-driven ring fluxes 
;  nohighbubspec- set to use the bubble spectrum at high latitudes (abs(b)=40-50 deg)
;                 instead of the whole bubble spectrum
;  doubring     - set not to include a second set of rings with the diffuse model spectrum
;                 otherwise templates are smoothed *to* the PSF of the data at each energy
;  dmbrems      - set not to include a template correlated with brems emission due to electrons
;                 from DM annihilations (SFD X NFW)
;  fdmdisc      - set to interpolate Fermi diffude model at discrete energies instead of
;                 sampling over logarithmic energies
;
; OUTPUTS:
;  mints        - minimum TS
;  getruntag    - set to only retrieve the runtag
;
; EXAMPLES:
;  For default analysis configuration just run
;  fermi_ring_fit
;
; COMMENTS:
;  Data and all templates other than Fermi diffuse model are smoothed by a Gaussian kernel
;  of 0.5 degree FWHM to account for the native resolution of the Fermi diffuse model
; 
;  All templates are also smoothed by the Fermi PSF to account for the Fermi data resolution
;  as a function of energy
;
;  The sky map at each energy is fitted as a linear combination of the Fermi diffuse model, ring, bubble and isotropic
;  templates. The spectrum of rings is fixed
;
; REVISION HISTORY:
;   2013-Oct-24 - copied from fermi_template_fit.pro by DPF, SP, TD
;   2014-Jan-28 - Modified by Tansu Daylan, Harvard Physics
;   2014-Jun-04 - Modified by Tansu Daylan, Harvard Physics
;----------------------------------------------------------------------

pro fermi_ring_fit, nring=nring, rin=rin, dr=dr, mine=mine, maxe=maxe, $
    fitreg=fitreg, bcut=bcut, nsigma=nsigma, nohighbubspec=nohighbubspec, $
    fdmver=fdmver, nslice=nslice, fdmdisc=fdmdisc, jringfit_=jringfit_, $
    ienergy=ienergy, nfwconst_=nfwconst_, psimax=psimax, version=version, $
    doubring=doubring, dmbrems=dmbrems, specoffset=specoffset, $
    correrr_=correrr_, fixedspec=fixedspec, mints=mints, runtag=runtag, getruntag=getruntag

  time1 = systime(1)
  
  if not keyword_set(mdm) then mdm = 35. ; [GeV]

  if not keyword_set(mine) then mine = 0
  if not keyword_set(maxe) then maxe = 19
  if not keyword_set(nslice) then nslice = 5
  if not keyword_set(psimax) then psimax = 180. ; [deg]
  if not keyword_set(bcut) then bcut = 2. ; [deg]
  if not keyword_set(nsigma) then nsigma = 5
  if not keyword_set(fixedspec) then fixedspec = 'nfw'
  if not keyword_set(version) then version = 'p11q2'
  if not keyword_set(fitreg) then fitreg = 'al'
  if not keyword_set(fdmver) then fdmver = 3
  if not keyword_set(ienergy) then ienergy = 9
  if not keyword_set(nring) then nring = 10
  if not keyword_set(jringfit_) then jringfit_ = [1:nring-2]
  if not keyword_set(dr) then dr = 1.0 ; [deg]
  if not keyword_set(rin) then rin = 1.5 ; [deg]
   

; -------- run tag to be used for output files

  runtag = 'nring' + string(nring, format='(I2.2)') + $
    '_rin' + string(rin, format='(F3.1)') + $
    '_en' + string(mine, format='(I2.2)') + string(maxe, format='(I2.2)') + $
    '_psi' + string(psimax, format='(I3.3)') + $
    '_reg' + fitreg + $
    '_bcut' + string(bcut, format='(I2.2)') + $
    '_spec' + fixedspec + $
    '_ver' + version + $
    '_fdm' + string(fdmver, format='(I1.1)') + $
    '_nsl' + string(nslice, format='(I1.1)')
  runtag += '_jrn'
  for i=0, n_elements(jringfit_)-1 do runtag += string(jringfit_[i], format='(I1.1)')
  if keyword_set(doubring) then runtag += '_doubr'
  if keyword_set(dmbrems) then runtag += '_dmbre'
  if keyword_set(correrr_) then runtag += '_ercor'
  if keyword_set(nohighbubspec) then runtag += '_highb'
  if keyword_set(specoffset) then begin
    runtag += '_off'
    if specoffset lt 0. then runtag += '-'
    runtag += string(abs(specoffset), format='(I1)')
  endif
  if keyword_set(getruntag) then return


; -------- analysis constants

  nside = 256L
  npix = 12L * nside * nside
  healgen_lb, nside, l, b
  l = ((l + 180.) mod 360.) - 180.

  nenergy = maxe - mine + 1
  ienergy = ienergy - mine
  nstr = 100
  iring = indgen(nring)


; -------- define fit region

  cpsi = cos(!DtoR * l) * cos(!DtoR * b)
  psi = !radeg * acos(cpsi) ; [degree]
  ipix = list()
  npixfit = lonarr(nenergy)
  npixcum = lonarr(nenergy)
  if fitreg eq 'al' then mpix = where(abs(b) gt bcut and psi lt psimax)
  if fitreg eq 'cc' then mpix = where(abs(b) gt bcut and abs(b) lt 20. and abs(l) lt 20.)
  for i=0, nenergy-1 do begin
    ipix.Add, mpix
    npixfit[i] = n_elements(ipix[i])
    if i eq 0 then npixcum[i] = 0 else $
      npixcum[i] = npixfit[i-1] + npixcum[i-1]
  endfor


; -------- read Fermi maps and define energy bins


  if version eq 'p11q2' then path = '$FERMI_DATA/allsky/p11_ultraclean_Q2/specbin'
  if version eq 'p11' then path = '$FERMI_DATA/allsky/p11_ultraclean/specbin'

  flist = fermi_bestmaps(dir=path, /front)

  fermimap = fltarr(npix, nenergy)
  fermiexp = fltarr(npix, nenergy)
  energy    = fltarr(nenergy)
  dearr     = fltarr(nenergy)
  emax      = fltarr(nenergy)
  emin      = fltarr(nenergy)
  for i=0, nenergy-1 do begin 
    fermimap[*, i] = readfits(flist[i+mine], h) ; [cts/cm^2/s/sr]
    fermiexp[*, i] = readfits(flist[i+mine], ext=1, /silent) ; [cm^2*s]
    energy[i] = sxpar(h, 'EMEAN') ; [GeV]
    dearr[i]  = sxpar(h, 'EMAX') - sxpar(h, 'EMIN') ; [GeV]
    emax[i] = sxpar(h, 'EMAX') ; [GeV]
    emin[i] = sxpar(h, 'EMIN') ; [GeV]
  endfor


; -------- smooth data to the inherent resolution of the Fermi diffuse model

  fdmfwhm = 30. ; [arcmin]
  fermimap = heal_smooth(fermimap, fdmfwhm)
  fermiexp = heal_smooth(fermiexp, fdmfwhm)


; -------- get spectral templates

  tracyspec  = mrdfits('$FERMI_RING_PATH/fits/tracy_spectrum.fits', 1, hdr)
  sfac = 1d-6


; isotropic spectrum

  isospec    = tracyspec.vals[*, 0] * sfac ; isotropic spectrum
  isospec    = isospec[0:nenergy-1]


; Fermi diffuse model spectrum

  difspec    = tracyspec.vals[*, 3] * sfac
  difspec    = difspec[0:nenergy-1]


; bubble spectrum
 
  if not keyword_set(nohighbubspec) then begin
    readcol, '$FERMI_RING_PATH/dat/bubble_spectrum.dat', enbubspec, bubspec
    bubspec = bubspec[mine:maxe] ; E^2dN/dE [GeV/cm^2/s/sr]
    bubspec /= energy^2 ; dN/dE [1/cm^2/s/sr/GeV]
  endif else begin
    bubspec = tracyspec.vals[*, 1] * sfac
    bubspec = bubspec[0:nenergy-1] ; dN/dE [1/cm^2/s/sr/GeV]
  endelse


; pulsar spectrum

  pulspec = mrdfits('$FERMI_RING_PATH/fits/pulsar_spectrum.fits', 0); pulsar spectrum
  pulspec = pulspec[0:nenergy-1]


; excess spectrum

  if fitreg eq 'al' then readcol, '$FERMI_RING_PATH/dat/ig_spec_full.dat', c1, c2, c3, nfwspec, ernfwspec
  if fitreg eq 'cc' then readcol, '$FERMI_RING_PATH/dat/ig_spec_in.dat', c1, c2, c3, nfwspec, ernfwspec
  nfwspec *= 1d-6 ; [1/cm^2/s/sr/GeV]
  ernfwspec *= 1d-6 ; [1/cm^2/s/sr/GeV]
  if keyword_set(specoffset) then nfwspec += specoffset / energy^2 * 1d-6
  nfwspec = nfwspec[0:nenergy-1] 
  ernfwspec = ernfwspec[0:nenergy-1] 


; power law spectrum

  powspec = 1./energy^2.5
  powspec *= nfwspec[0] / powspec[0] ; dN/dE [abu]


; bremsstrahlung spectrum

  readcol, '$FERMI_RING_PATH/dat/brems_gammaspec.dat', en, brespec
  brespec = interpol(brespec, en, energy); E^2dN/dE [GeV/cm^2/s/sr]
  brespec /= energy^2; dN/dE [1/cm^2/s/sr/GeV]


; bbbar annihilation spectrum

  pppc4dm = read_ascii('$FERMI_RING_PATH/dat/pppc4dm/AtProduction_gammas.dat') 
  mindif = min(abs(pppc4dm.field01[0,*] - mdm), ipppc4dm)
  ipppc4dm = where(pppc4dm.field01[0,*] eq pppc4dm.field01[0,ipppc4dm])
  enmult = 10^(pppc4dm.field01[1,ipppc4dm])
  enmult *= mdm ; [GeV]
  mult = pppc4dm.field01[13,ipppc4dm]
  dmaspec = interpol(mult, enmult, energy) / energy ; dN/dE [1/cm^2/s/sr/GeV]


; -------- impose spectra on the rings

  erspec = fltarr(nenergy)
  if fixedspec eq 'nfw' then begin ; data-driven NFW correlated emission spectrum (default)
  spec = nfwspec
  erspec = ernfwspec
  endif
  if fixedspec eq 'dma' then begin ; PYTHIA generated DM annihilation spectrum
  spec = dmaspec
  endif
  if fixedspec eq 'pul' then begin ; Dan's fit to the pulsar spectrum
  spec = pulspec
  endif
  if fixedspec eq 'pow' then begin ; simple power spectrum
  spec = powspec
  endif
  if fixedspec eq 'dif' then begin ; diffuse model spectrum
  spec = difspec
  endif


; -------- data structure describing templates

  str = replicate({template: fltarr(npix), $
                   templatearr: fltarr(npix, nenergy), $
                   spectrum: fltarr(nenergy), $
                   spectrumer: fltarr(nenergy), $
                   fixed: 0, $
                   varfixed: 0, $
                   varfloat: 0, $
                   float: 0}, nstr)
 
 
; -------- produce and smooth ring templates

  rintemp = fermi_ring_template(nring, rin, dr, nside=nside)
  rintemparr = fltarr(npix,nring,nenergy)
  for i=0, nenergy-1 do rintemparr[*,*,i] = fermi_heal_smooth(rintemp, energy[i], gausfwhm=fdmfwhm)


; -------- put templates into the structure

  strocc = 0 ; structure occupancy


; ring templates fixed to the data-driven NFW correlated spectrum

  for i=0, nring-1 do begin
    str[i].varfixed = 1
    str[i].templatearr = rintemparr[*,i,*]
    str[i].spectrum = spec
    str[i].spectrumer = erspec
  endfor
  strocc += nring


; ring templates fixed to the data-driven Fermi diffuse model correlated spectrum

  if keyword_set(doubring) then begin
    for i=0, nring-1 do begin
      str[strocc+i].varfixed = 1
      str[strocc+i].templatearr = rintemparr[*,i,*]
      str[strocc+i].spectrum = difspec
      str[strocc+i].spectrumer = erdifspec
    endfor
    strocc += nring
  endif


; isotropic template

  isotemp = fltarr(npix) + 1.
  str[strocc].template = isotemp
  str[strocc].float = 1
  iscale = fltarr(2*nenergy)
  for i=0, nenergy-1 do iscale[i] = strocc + i
  strocc += nenergy


; Fermi diffuse model template

  diftemparr = fltarr(npix, nenergy)
  if not keyword_set(fdmdisc) then for i=0, nenergy-1 do diftemparr[*,i] = fermi_diffuse_model([Emin[i],Emax[i]], nside=nside, ver=fdmver) $
  else for i=0, nenergy-1 do diftemparr[*,i] = fermi_diffuse_model(Energy[i], nside=nside, ver=fdmver)
  for i=0, nenergy-1 do diftemparr[*,i] = fermi_heal_smooth(diftemparr[*,i], energy[i])
  bubtemp = bubble_fill(1) + bubble_fill(3)
  idifnorm = where(bubtemp gt 0.5 and abs(b) gt bcut, ndifnorm)
  diffac = total(diftemparr[idifnorm,*])/ndifnorm
  diftemparr = diftemparr / diffac
  str[strocc].templatearr = diftemparr
  str[strocc].varfloat = 1
  for i=0, nenergy-1 do iscale[nenergy+i] = strocc + i
  strocc += nenergy


  ;bubble templates

  sbutemp = fermi_sliced_bubbles(nslice=nslice)
  sbutemparr = fltarr(npix,nslice,nenergy)
  for i=0, nenergy-1 do sbutemparr[*,*,i] = fermi_heal_smooth(sbutemp, energy[i], gausfwhm=fdmfwhm)
  for i=0, nslice-1 do begin
    str[strocc+i].templatearr = sbutemparr[*,i,*]
    str[strocc+i].spectrum = bubspec
  endfor
  if fdmver eq 4 then begin
    str[strocc:strocc+nslice-1].varfixed = 0
  endif else begin
    str[strocc:strocc+nslice-1].varfixed = 1
    strocc += nslice
  endelse


  ;DM brems template

  if keyword_set(dmbrems) then begin
    dmbtemparr = fltarr(npix,nenergy)
    sfdtemp = dust_getval(l, b, /noloop, /interp)
    nfwtemp = fermi_nfw_map(1.3)  
    dmbtemp = sfdtemp * nfwtemp
    for i=0, nenergy-1 do dmbtemparr[*,i] = fermi_heal_smooth(dmbtemp, energy[i], gausfwhm=fdmfwhm)
    str[strocc].templatearr = dmbtemparr
    str[strocc].spectrum = brespec
    str[strocc].varfixed = 1
    strocc++
  endif


; -------- create the design and error weight matrices
  Amatrix = fermi_ring_design(str, nenergy=nenergy, ipix=ipix, npixfit=npixfit, npixcum=npixcum)


; -------- do maximum likelihood fit

  fermi_region_fit, fermimap, fermiexp, Amatrix, dEarr, energy, ipix, $
    nsigma=nsigma, cov=cov, erxx=erxx, $
    mints=mints, winpix=winpix, xx=xx, exposure=exposure, $
    counts=counts, iscale=iscale, scalefac=1d6


; -------- fit data ring fluxes to the NFW template ring fluxes

  common fermi_ring_chisq_common, $
    rinfcin, rinflux, errinflux, nfwflux, $
    nfwconst, correrr, a, jringfit

  if keyword_set(nfwconst_) then nfwconst = 1
  if keyword_set(correrr_) then correrr = 1
  jringfit = jringfit_

  rinflux = xx[0:nring-1] * spec[ienergy] * energy[ienergy]^2
  errinflux = erxx[0:nring-1] * spec[ienergy] * energy[ienergy]^2
  rinfcov = fltarr(nring, nring)
  for i=0, nring-1 do rinfcov[i,*] = cov[i,0:nring-1] * (spec[ienergy] * energy[ienergy]^2)^2
  rinfcov[iring,iring] += rinflux[iring]^2 * (erspec[ienergy] / spec[ienergy])^2
  rinfcin = invert(rinfcov)

  if keyword_set(doubring) then begin
    difflux = xx[0:nring-1] * difspec[ienergy] * energy[ienergy]^2
    erdifflux = xx[0:nring-1] * erdifspec[ienergy] * energy[ienergy]^2
  endif 


; -------- fit the ring coefficients with ring weighted NFW coefficients

  ng = [1.3, 1.35, 1.4]
  nng = n_elements(ng)
  nfwtemp = fltarr(npix, nng)
  for a=0, nng-1 do nfwtemp[*,a] = fermi_ring_nfw(ng[a])
  nfwtemp = fermi_heal_smooth(nfwtemp, energy[ienergy])

  nfwflux = fltarr(nng, nring)
  for i=0, nring-1 do begin 
    ind = where(rintemp[*,i] gt 1d-4, nind)
    for a=0, nng-1 do begin
      testmap = rintemp[*,i] * nfwtemp[*,a]
      nfwflux[a,i] = total(testmap[ind]) / nind
    endfor
  endfor

  psiavg = findgen(nring) * dr + rin + dr / 2.
  erpsiavg = fltarr(nring) + dr / 2.

  if keyword_set(nfwconst) then nfwpar = [1.,0.] else nfwpar = 1.
  ndof = n_elements(jringfit) - n_elements(nfwpar)

  for a=0, nng-1 do nfwflux[a,*] *= rinflux[0] / nfwflux[a,0]

  chisq = fltarr(nng)
  for a=0, nng-1 do begin
    nfwpar = tnmin('fermi_ring_chisq', nfwpar, bestmin=bestmin, /auto)
    chisq[a] = bestmin / ndof
    nfwflux[a,*] *= nfwpar
  endfor
  pval = 1. - chisqr_pdf(chisq, ndof)


; -------- plots

  fermi_ring_flux_plot, psiavg, erpsiavg, rinflux, errinflux, $
    nfwflux, nring, chisq, runtag=runtag, difflux=difflux, erdifflux=erdifflux

  fermi_ring_map_plot, xx, counts, Amatrix, exposure, energy, ipix, npixfit, npixcum, runtag, nring, $
    doubring=doubring, dmbrems=dmbrems, nslice=nslice, fdmver=fdmver, ienergy=ienergy

;  fermi_ring_spec_plot, energy, nfwspec, ernfwspec, runtag, $
;    bubspec, difspec, isospec, pulspec, powspec, dmaspec


; ------- write results

  writecol, 'fermi_ring_' + runtag + '.dat', psiavg, rinflux, errinflux, nfwflux

  file_mkdir, '$FERMI_RING_OUTPUT_PATH/dat' 
  spawn, 'mv fermi_ring_*.dat $FERMI_RING_OUTPUT_PATH/dat'

  file_mkdir, '$FERMI_RING_OUTPUT_PATH/ps' 
  spawn, 'mv fermi_ring_*.ps $FERMI_RING_OUTPUT_PATH/ps'


; -------- measure procedure execution time  

  time2 = systime(1)
  splog, time2 - time1, ' seconds elapsed.'

end
