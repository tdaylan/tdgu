;+
; NAME:
;   puc_main
;
; PURPOSE:
;   Main routine of the Photon Underproduction Crisis (PUC) project
;
; CALLING SEQUENCE:
;   puc_main
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mdm        - mass of the DM particle [MeV]
;   csvboost   - boost factor of the annihilation cross section
;   substfrac  - tuner for substructure
;   escmod     - string tag for the particle escape model
;   conmod     - string tag for the mass-concentration prescription
;   annch      - string tag for the annihilation channel
;   minhalomass- minimum halo mass to be considered when integrating over the halo mass function [Msun]
;
;   nr_        - number of radial bins
;   nz_        - number of azimuthal bins
;   nen_       - number of energy bins
;   nmass_     - number of mass bins
;   nred_      - number of redshift bins
;   ncdd_      - number of column density bins
;   nimage_    - number of images used in the semi-analytic calculation of the electron halo function
;   eps        - floating point precision
;
; KEYWORDS:
;   noplot     - set to suppress plots
;  
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   runtag_    - run tag
;   runtime    - run time [s]
;
; EXAMPLES:
;
; COMMENTS:
;
;   The code calculates
;   - electron and positron luminosity of a DM halo of mass M
;   - the escape fraction of electrons and positrons from the DM halo
;   - the halo mass function
;   - electron and positron emissivity in the IGM
;   - electron and positron number density in the IGM
;   - radiated power due to ICS on the CMB
;   - resulting UV/X-ray emissivity
;   - cosmic background UV/X-ray luminosity due to the los projection and radiative transfer of photons in the clumpy IGM
;   - metagalactic ionization rate
;   - the temperature evolution of the IGM
;
;   The project uses a semi-analitic solution to the diffusion-loss equation
;   assuming MIN/MED/MAX propagation model, where diffusion takes place
;   inside a cylindrical slab centered around a Milky Way sized galaxy.
;   The energy loss rate, b = dE/dt, is assumed spatially constant for the solution of the
;   propagated energy number density. In particular the code computes:
;
;     - diffusion length, \lambda(E,E_s), as a function of injection and final electron energies,
;     - diffusion halo function, I(E,E_s,r,z) as a function of diffusion length and position in the galaxy,
;     - injection-integrated halo function, I'(E,r,z) by convolving the halo function with the injection spectrum
;     - steady-state e-/e+ number density in the MW, \Psi(E,r,z) = I'(E,r,z) / b(E)
;
;   where the diffusion halo function is the convolution of the scaled heat equation propagator and the injection morphology
;   To implement the boundary conditions, method of images is used since the diffusion scale is smaller
;   than the halo size. A series of images of each source is added to make the boundaries at z = +- Lz, zero.
;
;   The code is vectorized with respect to the integration varible whenever possible. All integrals are performed by the trapezoidal rule.
;
;   The project uses a common block, puc_common, for (almost) all variables. Since elements on the common block are memory-friendly, small
;   arrays, this makes it convenient to pass variables between functions and procedures.
;
;
; REVISION HISTORY:
;  2014-Aug-23 - Written by Tansu Daylan, Harvard Physics
;  2015-Jan-26 - second major revision Tansu Daylan, Harvard Physics
;----------------------------------------------------------------------
pro puc_main, $
    promod=promod, minhalomass=minhalomass, conmod=conmod, escmod=escmod, substfrac=substfrac, mdm=mdm, csvboost=csvboost, $
    eps=eps, nimage_=nimage_, annch=annch, $
    nr_=nr_, nz_=nz_, nen_=nen_, nmass_=nmass_, nred_=nred_, ncdd_=ncdd_, $
    noplot=noplot, $
    runtag_=runtag_

  timestart = systime(1)

; the main PUC procedure is set up so that all integrand functions access variables via puc_common

  common puc_common, $
    r, z, enel, denel, ienel, $
    enph, enph_, red, dred, mass, massp, $
    masshmf, dndmhmf0, dndmhmf3, dndmhmf6, $
    ired, iredp, jred, jredtemp, gwf, $
    i, j, a, b, c, d, e, $ ; respective global indices
    nr, nz, nenph, nenel, nred, nmass, nmassp, $ ; respective sizes
    minr, haloelec, npix, mineel, lambda, lambda2, rho2ij, maxr, rho2, emismw, emismwsun, $
    phdmemis, phdmemis1, $
    phqgemis, phqgemis1, $
    phqsemis, phqsemis1, $
    phdmflux, phdmflux1, $
    phqgflux, phqgflux1, $
    phqgfluxhm, phqgfluxhm1, $
    phqsflux, phqsflux1, $
    phionrdm, phionrqg, phionrqghm, phionrqs, $
    phheatdm, phheatqg, phheatqghm, phheatqs, plheatbl, heat, $
    icsedotcmbr, icsedotthom, $
    icspowercmbrint, icspowerthomint, $
    k0, e0, delta, lz, lightspeed, htime, elmass, planck, rhosun, lcompton, thomcs, alphaem, ryd, sigma2, $ ; analysis constants
    enisrf, nsubs, boltzmann, $
    nphcmbr, nphcmbrgp, ucmbr, $
    isrfmethod, redshift, rmass, rhom, smtomev, $
    icspowercmbr, icspowerthom, yic, qic, yfac, qfac, $
    inmult, rsun, kpc2cm, cmw, mmw, ioncs, sigmam, deltaz, $
    escfracmw, escfrac, isrf, risrf, zisrf, nrisrf, nzisrf, nenisrf, kt, rho, $
    sigmaion, hrad, hubbf, $
    haloinj, fluxmwesc, minephcmbr, maxephcmbr, $
    magfd, gasdn, uisrf, $
    nimage, $
    dndm, psiigm, psiigmsst, $
    lumihalo, lumihalons, lumihalone, lumihalonsne, lumihalointns, lumihalointne, lumihalointnsne, lumihaloint, lumihaloflag, $
    elemisigm, elemisigmint6, elemisigmint, elemisigm6, elemisavg, elemisavgint, $
    numbmwsst, numbmwesc, lumimw, lumimwesc, rhomw, $
    runtag, $
    opeff, optot, optotint, opind, d2fdndz, cdd, $
    omm, oml, omk, $ 
    enpi, $
    tempdm, tempqs, tempbl, temp, $
    cmbt, patchdensarr, dpatchdensdzarr, dpatchdensdz, patchdens, temppatchdm, npatch

  if not keyword_set(eps) then eps = 1d-4
  if keyword_set(nimage_) then nimage = nimage_ else nimage = 3 ; gas density [cm^-3]

; -------- constants

  elmass = 0.510998903d ; electron mass [MeV]
  lightspeed = 2.99792458d10 ; speed of light [cm/s]
  thomcs = 6.6524d-25 ; Thompson cross section [cm^2]
  planck = 4.135667516d-21 ; Planck constant [MeV s]
  alphaem = 7.2973525698d-3 ; fine structure constant

  oml = 0.7 ; Omega Lambda
  omk = 0. ; Omega curvature
  omm = 0.3 ; Omega matter
  deltac = 1.686 ; critical linear overdensity
  sigma8 = 0.83 ; matter power spectrum normalization 
  nsubs = 0.96 ; spectral index of the primordial power spectrum
  stla = 0.707 ; Sheth-Tormen parameter a
  stba = 0.3222 ; Sheth-Tormen parameter A
  stp = 0.3 ; Sheth-Tormen parameter
  rsun = 8.5 ; radial distance from the GC to the Sun [kpc]
  zsun = 0. ; azimuthial distance from the GP to the Sun [kpc]
  rhosun = 300. ; DM energy density at the Sun [MeV/cm^3]
  rhoavg = 1.31d-3 ; average DM energy density at z = 0 [MeV/cm^3]
  ovdens = 18. * !pi^2 ; linear overdensity at virialization
  ryd = 13.5984d-6 ; Rydberg energy [MeV]
  htime = 4.35d17 ; Hubble time [s]
  hrad = htime * lightspeed ; Hubble radius [cm]
  hubb = 0.704 ; reduced Hubble constant
  boltzmann = 8.6173d-11 ; Boltzmann constant [MeV/K]
  kpc2cm = 3.086d21 ; conversion factor from kpc to cm
  erg2mev = 6.241509d5 ; conversion factor from erg to MeV
  kt0 = 2.33d-10 ; CMB thermal scale at z = 0 [MeV]
  cmbt0 = 2.725 ; CMB temperature at z = 0 [K]
  sigmav0 = 3d-26 ; [cm^3/s] velocity averaged thermal annihilation cross section

  cmw = 15 ; concentration parameter of the MW
  mmw = 1d12 ; mass of the Milky Way [Msun]
  smtomev = 1.115d60 ; mass of the Milky Way [Msun]

  sigmaion0 = 6.3d-18 ; [cm^2] neutral Hydrogen photon-ionization cross section at 13.6 eV 


; -------- photon energy axis

  if keyword_set(nen_) then nenph = nen_ else nenph = 100
  if not keyword_set(mineph) then mineph = 1d-6 ; [MeV]
  if not keyword_set(maxeph) then maxeph = 1d-1 ; [MeV]
  enph = exp(dindgen(nenph) /(nenph-1) * (alog(maxeph) - alog(mineph)) + alog(mineph)) ; [MeV]
  frph = enph / planck ; [Hz]
  sigmaion = sigmaion0 * (ryd / enph)^3


; -------- gamma-ray excess parameters

  if not keyword_set(csvboost) then csvboost = 1.
  sigmav = csvboost * sigmav0

  if not keyword_set(mdm) then mdm = 3d4 ; WIMP mass [MeV]
  if not keyword_set(annch) then annch = 'b'

  pppc4dm = read_ascii('$PUC_PATH/data/pppc4dm/AtProduction_positrons.dat') 
  mindif = min(abs(pppc4dm.field01[0,*] - mdm/1d3), ipppc4dm)
  ipppc4dm = where(pppc4dm.field01[0,*] eq pppc4dm.field01[0,ipppc4dm])
  enmult = 10^(pppc4dm.field01[1,ipppc4dm])
  enmult *= mdm ; [MeV]

  if annch eq 'e' then mult = pppc4dm.field01[4,ipppc4dm]
  if annch eq 'm' then mult = pppc4dm.field01[7,ipppc4dm]
  if annch eq 't' then mult = pppc4dm.field01[10,ipppc4dm]
  if annch eq 'b' then mult = pppc4dm.field01[13,ipppc4dm]


; -------- electron energy axis

  if keyword_set(nen_) then nenel = nen_ else nenel = 100
  if not keyword_set(mineel) then mineel = 50. ; [MeV]
  if not keyword_set(maxeel) then maxeel = max(enmult) ; [MeV]
  enel = exp(dindgen(nenel) /(nenel-1) * (alog(maxeel) - alog(mineel)) + alog(mineel)) ; [MeV]
  denel = enel[1:nenel-1] - enel[0:nenel-2]
  ienel = indgen(nenel)


; -------- redshift axis

  if keyword_set(nred_) then nred = nred_ else nred = 100
  if not keyword_set(minred) then minred = 0.1
  if not keyword_set(maxred) then maxred = 1000.
  red = exp(dindgen(nred) /(nred-1) * (alog(maxred) - alog(minred)) + alog(minred))
  dred = red[1:nred-1] - red[0:nred-2]
  ired = indgen(nred)
  ireddif1 = min(abs(red - 3), ired1)
  ireddif2 = min(abs(red - 6), ired2)
  iredp = [0,ired1,ired2]
  cmbt = cmbt0 * (1. + red) ; CMB temperature [K]
  kt = cmbt * boltzmann ; CMB thermal scale [MeV]
  scf = reverse(1. / (red + 1.)) ; scale factor
  hubbf = sqrt(oml + omm * (1 + red^3)) ; Hubble function


; -------- logarithmically spaced mass axis

  if keyword_set(nmass_) then nmass = nmass_ else nmass = 100
  if not keyword_set(minmass) then minmass = 1d-10 ; [Solar Mass]
  if not keyword_set(maxmass) then maxmass = 1d15 ; [Solar Mass]
  mass = exp(dindgen(nmass) /(nmass-1) * (alog(maxmass) - alog(minmass)) + alog(minmass))


; -------- column density distriubtion axis

  if keyword_set(ncdd_) then ncdd = ncdd_ else ncdd = 100
  if not keyword_set(mincdd) then mincdd = 10^(11.) ; [1/cm^2]
  if not keyword_set(maxcdd) then maxcdd = 10^(21.55) ; [1/cm^2]
  cdd = exp(dindgen(ncdd) /(ncdd-1) * (alog(maxcdd) - alog(mincdd)) + alog(mincdd))


; -------- radial axis

  if keyword_set(nr_) then nr = nr_ else nr = 100
  if not keyword_set(minr) then minr = 0. ; [kpc]
  if not keyword_set(maxr) then maxr = 20. ; [kpc]
  r = dindgen(nr) / (nr-1) * (maxr - minr) + minr


; -------- initial photon energy axes

  if keyword_set(nen_) then nenpi = nen_ else nenpi = 1000
  if not keyword_set(minenpi) then minenpi = 1d-12 ; [MeV]
  if not keyword_set(maxenpi) then maxenpi = 1d-3 ; [MeV]
  enpi = exp(dindgen(nenpi) / (nenpi - 1) * (alog(maxenpi) - alog(minenpi)) + alog(minenpi))
  ienpi = indgen(nenpi)


; -------- z axis

  if keyword_set(nz_) then begin
    if not (nz_ mod 2 eq 0) then begin
      print, 'nz should be even!'
      return
    endif
    nz = nz_
  endif else begin
    nz = 50
  endelse
  if not keyword_set(minz) then minz = -10. ; [kpc]
  if not keyword_set(maxz) then maxz = 10. ; [kpc]
  z = dindgen(nz) / (nz-1) * (maxz - minz) + minz


; -------- halo mass range

  if not keyword_set(minhalomass) then minhalomass = 1d-9 ; [Msun]
  maxhalomass = maxmass ; [Msun]


; -------- substructure model

  if not keyword_set(substfrac) then substfrac = 1. ; [Msun]
  subst = (1. + substfrac * 110. * (mass / mmw)^0.39) # (1. + dblarr(nred))


; -------- mass concentration model

  if not keyword_set(conmod) then conmod = 'pow'
  if conmod eq 'pra' then begin
    cona = [37.5153, -1.5093, 1.636d-2, 3.66d-4, -2.89237d-5, 5.32d-7]
    con = mass # red * 0
    for i=0, 5 do con += cona[i] * alog(mass / hubb)^i # hubbf^(-2./3.)
  endif
  if conmod eq 'low' then begin
    cmind = -0.05
    cmpiv = 3.37d12 / 0.7 
    cmnor = 6.5
    con = cmnor * (mass / cmpiv)^cmind # hubbf^(-2./3.)
  endif
  if conmod eq 'hig' then begin
    cmind = -0.15
    cmpiv = 3.37d12 / 0.7 
    cmnor = 6.5
    con = cmnor * (mass / cmpiv)^cmind # hubbf^(-2./3.)
  endif
  if conmod eq 'cut' then begin
    cmind = -0.1
    cmpiv = 3.37d12 / 0.7 
    cmnor = 6.5
    mcut = 1d-4
    con = cmnor * (mass / cmpiv)^cmind # hubbf^(-2./3.)
    conint = interpolate(con, interpol(indgen(nmass), mass, mcut), ired, /grid)
    for c=0, nred-1 do con[*,c] = con[*,c] < conint[c]
  endif
  if conmod eq 'pow' then begin
    cmind = -0.1
    cmpiv = 3.37d12 / 0.7 
    cmnor = 6.5
    con = cmnor * (mass / cmpiv)^cmind # hubbf^(-2./3.)
  endif
  if conmod eq 'mac' then begin
    cmind = -0.1
    crind = -1.
    ;con = 3.9 * (interpol(hubbf, redfor, mass) / hubbf)^(2./3.)
  endif


; -------- parametric CR diffusion model with the MED model being the best fit to B/C data

  if not keyword_set(promod) then promod = 'med'
  e0 = 1. ; [GeV]
  if promod eq 'min' then begin
    k0 = 5.074d-17 ; [kpc^2/s]
    delta = 0.85
    lz = 1. ; [kpc]
  endif
  if promod eq 'med' then begin
    k0 = 3.551d-16 ; [kpc^2/s]
    delta = 0.7
    lz = 4. ; [kpc]
  endif
  if promod eq 'max' then begin
    k0 = 2.426d-15 ; [kpc^2/s]
    delta = 0.46
    lz = 8. ; [kpc]
  endif


; -------- define model for particle escape from the DM halo

  if not keyword_set(escmod) then escmod = 'sim'
  
  if escmod eq 'sim' then begin
    escmind = -1.5
    esczind = 1.
    escfracmw = dblarr(nenel) + 0.1 
  endif else begin
    if escmod eq 'min' or escmod eq 'max' or escmod eq 'med' then begin
      escmind = -1.5
      escmind = 0.5
    endif
    if escmod eq 'min' then begin
      magfd = 0.5 ; [muG]
      gasdn = 0.1 ; [1/cm^3]
      usirf = 5d-5 ; [MeV/cm^3]
    endif
    if escmod eq 'med' then begin
      magfd = 5. ; [muG]
      gasdn = 1. ; [1/cm^3]
      usirf = 5d-6 ; [MeV/cm^3]
    endif
    if escmod eq 'max' then begin
      magfd = 50. ; [muG]
      gasdn = 10. ; [1/cm^3]
      usirf = 5d-7 ; [MeV/cm^3]
    endif
    puc_escfrac
  endelse
  escfrac = dblarr(nenel,nmass,nred)
  for d=0, nmass-1 do for c=0, nred-1 do escfrac[*,d,c] = escfracmw * (mass[d] / mmw)^escmind * (1 + red[c])^esczind < 1.


; -------- run tag

  runtag = promod + '_mhm' + string(alog10(minhalomass), format='(I3.2)') + '_' + $
    conmod + '_' + escmod + '_sbf' + string(fix(substfrac), format='(F4.2)') + $ 
    '_mdm' + string(round(mdm/1d3), format='(I3.3)') + '_csb' + string(alog10(csvboost), format='(I2.2)') + '_' + annch
  if keyword_set(runtag_) then begin
    runtag_ = runtag
    return
  endif


; -------- get number density of photons in ISRF and CMB

  puc_numdens
  nphcmbr = 8. * !pi / (lightspeed * planck)^3 * (enpi^3 # (1. + dblarr(nred))) / (exp(enpi # (1. / kt)) - 1.) ; [1/cm^3]
  ucmbr = dblarr(nred)
  for c=0, nred-1 do ucmbr[c] = integral(enpi, nphcmbr[*,c]) ; [MeV/cm^3]

  if not keyword_set(noplot) then puc_edot_plot
;  if not keyword_set(noplot) then puc_galprop_plot


; -------- get annihilation spectrum

  mult /= enmult ; [1/MeV]
  mult = interpol(mult, enmult, enel) > 0.
  sigmav_ = enel^2 * mult * sigmav / mdm^2 ; [cm^3/s/MeV]


; -------- compute the luminosity per halo of mass M

  fc = con^3 / (alog(1. + con) - con / (1. + con)) * (1. - 1. / (1. + con)^3) / 3.
  rho2int = ovdens / 3. * rhoavg * smtomev * (mass # (dblarr(nred) + 1.)) * fc ; [MeV^2/cm^3]
  lumihalo = dblarr(nenel,nmass,nred)
  lumihalone = dblarr(nenel,nmass,nred)
  lumihalons = dblarr(nenel,nmass,nred)
  lumihalonsne = dblarr(nenel,nmass,nred)
  for a=0, nenel-1 do begin
    lumihalo[a,*,*] = sigmav_[a] * subst * escfrac[a,*,*] * rho2int ; [MeV/s]
    lumihalone[a,*,*] = sigmav_[a] * subst * rho2int ; [MeV/s]
    lumihalons[a,*,*] = sigmav_[a] * escfrac[a,*,*] * rho2int ; [MeV/s]
    lumihalonsne[a,*,*] = sigmav_[a] * rho2int ; [MeV/s]
  endfor
  lumihaloint = dblarr(nmass,nred)
  lumihalointne = dblarr(nmass,nred)
  lumihalointns = dblarr(nmass,nred)
  lumihalointnsne = dblarr(nmass,nred)
  for d=0, nmass-1 do $
    for c=0, nred-1 do begin
      lumihaloint[d,c] = integral(enel, lumihalo[*,d,c] / enel) 
      lumihalointne[d,c] = integral(enel, lumihalone[*,d,c] / enel) 
      lumihalointns[d,c] = integral(enel, lumihalons[*,d,c] / enel) 
      lumihalointnsne[d,c] = integral(enel, lumihalonsne[*,d,c] / enel) 
    endfor
  if not keyword_set(noplot) then puc_lumihalo_plot


; -------- integrate over the Sheth - Tormen mass function to get the emissivity in the IGM

  nmassp = nmass + 1
  massp = [10^(alog10(minmass) - (alog10(maxmass) - alog10(minmass)) / nmass), mass]
  rmass = (3. * massp * smtomev / 4. / !pi / rhoavg)^(1./3.) / kpc2cm / 1d3 ; [Mpc]
  kmass = 2. * !pi / rmass


; RMS mass fluctuations

  sigmam = dblarr(nmassp)
  for d=0, nmassp-1 do sigmam[d] = sqrt(qromo('puc_sigma2_integrand', 1d-15 / rmass[d], 1d2 / rmass[d], eps=eps))
  sigmam *= sigma8 / interpol(sigmam, rmass, 8. / hubb)


; critical overdensity

  gwf = dblarr(nred) ; growth factor
  for c=0, nred-1 do begin
    fac = sqrt(oml * scf[c]^3 + omk * scf[c] + omm) / scf[c]^1.5
    gwfint = fac * scf^1.5 / (oml * scf^3 + omk * scf + omm)^1.5
    gwf[c] = integral(scf, gwfint, xr=[min(scf),scf[c]])
  endfor
  gwf = reverse(gwf)
  gwf /= gwf[0]
  deltaz = deltac / gwf


; halo mass function

  dndm = dblarr(nmass, nred)
  for c=0, nred-2 do begin
    dndmnu = deltaz[c] / sigmam[1:nmass]
    dndm[*,c] = -stba * sqrt(2. * stla / !pi) * (1 + (dndmnu / stla)^stp) * rhoavg * kpc2cm^3 / smtomev * dndmnu * $
      (deriv(massp, alog(sigmam)))[1:nmass] * exp(-stla * dndmnu^2 / 2.) / mass ; [1/kpc^3/Msun]
  endfor

  readcol, '$PUC_PATH/data/hmfcalc/m0.txt', masshmf, c2, c3, c4, c5, dndmhmf0
  readcol, '$PUC_PATH/data/hmfcalc/m3.txt', masshmf, c2, c3, c4, c5, dndmhmf3
  readcol, '$PUC_PATH/data/hmfcalc/m6.txt', masshmf, c2, c3, c4, c5, dndmhmf6
  dndmhmf0 *= 1d-9 * hubb^4 ; [1/kpc^3/Msun]
  dndmhmf3 *= 1d-9 * hubb^4 ; [1/kpc^3/Msun]
  dndmhmf6 *= 1d-9 * hubb^4 ; [1/kpc^3/Msun]
  masshmf /= hubb ; [Msun]

  if not keyword_set(noplot) then puc_sigma2_plot

 
; e-/e+ emissivity

  elemisavg = sigmav_ * rhoavg^2 # (1. + red)^6 ; [MeV/s]
  elemisigm = dblarr(nenel,nred)
  elemisigm6 = dblarr(nenel,nred)
  for a=0, nenel-1 do begin
    for c=0, nred-1 do begin
      elemisigm[a,c] = integral(mass, dndm[*,c] * lumihalo[a,*,c], xr=[minhalomass, maxhalomass]) $
        / kpc2cm^3 * (1. + red[c])^3 + elemisavg[a,c] ; [MeV/cm^3/s]
      elemisigm6[a,c] = integral(mass, dndm[*,c] * lumihalo[a,*,c] , xr=[1d-6, maxhalomass]) $
        / kpc2cm^3 * (1. + red[c])^3 + elemisavg[a,c] ; [MeV/cm^3/s]
    endfor
  endfor

  elemisavgint = dblarr(nred)
  elemisigmint = dblarr(nred)
  elemisigmint6 = dblarr(nred)
  for c=0, nred-1 do begin
    elemisigmint[c] = integral(enel, elemisigm[*,c] / enel) 
    elemisigmint6[c] = integral(enel, elemisigm6[*,c] / enel) 
    elemisavgint[c] = integral(enel, elemisavg[*,c] / enel) 
  endfor

  if not keyword_set(noplot) then puc_elemisigm_plot


; -------- get the electron number density in the IGM

; analytic steady state solution, assuming energy loss is much faster than the Hubble expansion

  psiigmsst = dblarr(nenel,nred)
  for a=0, nenel-2 do begin
    for c=0, nred-1 do begin
      psiigmsst[a,c] = 1. / puc_edot(a, c, /cmbr) * integral(enel, elemisigm[*,c] / enel^2, xr=[enel[a],enel[nenel-1]]) ; [1/cm^3/MeV]
    endfor
  endfor


; finite difference solution to the transport equation

  psiigm = dblarr(nenel, nred)
  a = nenel - 2
  while a ge 0 do begin
    c = nred - 2
    while c ge 0 do begin
      psiigm[a,c] = psiigm[a,c+1] + htime * dred[c] / hubbf[c+1] / (1. + red[c+1]) * $
        ((psiigm[a+1,c+1] * puc_edot(a+1, c+1, /cmbr) - psiigm[a, c+1] * puc_edot(a, c+1, /cmbr)) / denel[a] + $
        elemisigm[a+1,c+1] / enel[a+1]^2) > 0. ; [1/cm^3/MeV]
      c--
    endwhile
    a--
  endwhile
  psiigm *= htime

; minimum chisq solution to the transport equation

  ;psiigmchi = psiigm
  ;psiigm_ = mpfit('pde_fun_chisq', reform(psiigm_,(nenel-1)*(nred-1)), auto=0, status=status)

; error in the analytical approximation

  psiigmerr = psiigmsst * 0.
  dpsiigmsstdz = psiigmsst * 0.
  for a=0, nenel-1  do dpsiigmsstdz[a,*] = deriv(red, psiigmsst[a,*])
  for c=0, nred-1 do begin
    dpsiigmsstde = deriv(enel, psiigmsst[*,c] * puc_edot(ienel, c, /cmbr))
    psiigmerr[*,c] = hubbf[c] * (1 + red[c]) / htime * dpsiigmsstdz[*,c] - dpsiigmsstde - elemisigm[*,c] / enel^2
  endfor

  if not keyword_set(noplot) then puc_psiigm_plot


; -------- compute the differential ICS power as a function of final photon and electron energies

  icspowercmbr = dblarr(nenph,nenel,nred)
  icsedotcmbr = dblarr(nenel)
  nqic = 100
  maxqic = 1.
  for b=0, nenel-1 do begin
    for a=0, nenph-1 do begin
      elgamma = enel[b] / elmass
      eps = enph[a] / enel[b]
      if eps ge 1. then continue
      for c=0, nred-1 do begin
        minqic = 1. / 4. / elgamma^2
        qic = exp(dindgen(nqic) / (nqic - 1) * (alog(maxqic) - alog(minqic)) + alog(minqic))
        enpitemp = elmass^2 / 4. / enel[b] * eps / (1. - eps) / qic
        qfac = 2. * qic * alog(qic) + qic + 1. - 2. * qic^2 + eps^2 * (1. - qic) / 2. / (1. - eps)    
        icspowercmbrint = 3. * thomcs * lightspeed / 4. / elgamma^2 * enph[a] * (enph[a] - enpitemp) * puc_dndecmbr(enpitemp, c) * qfac / qic ; [MeV/s]
        icspowercmbr[a,b,c] = integral(qic, icspowercmbrint) ; [MeV/s]     
      endfor
    endfor
    icsedotcmbr[b] = integral(enph, icspowercmbr[*,b,0] / enph) ; [MeV/s]
  endfor
  if not keyword_set(noplot) then puc_icspower_plot


; -------- quasar photon emissivity

  phqsemis1 = 10^24.6 * erg2mev / (kpc2cm*1d3)^3 * (1. + red)^7.68 * exp(-0.28 * red) / (exp(1.77 * red) + 26.3) ; [MeV/cm^3/s/Hz]
  phqsemis = ((enph / ryd)^(-0.44) * frph) # phqsemis1 ; [MeV/cm^3/s]
  ienph = where(enph gt planck * lightspeed / 1.3d-5)
  phqsemis[ienph,*] = ((enph[ienph] / ryd)^(-1.57) * frph) # phqsemis1 ; [MeV/cm^3/s]


; -------- quasar + galaxy photon emissivity

  phqgemisdat = read_ascii('$PUC_PATH/data/qgemis.dat') 
  redhm = phqgemisdat.field01[0:58,0]
  wlemishm = phqgemisdat.field01[0,1:378] ; [Angstrom]
  fremishm = reverse(reform(lightspeed / (wlemishm * 1d-8))) ; [Hz]
  enemishm = planck * fremishm ; [MeV]
  phqgemis = phqgemisdat.field01[1:59,1:378] ; [erg/Mpc^3/s/Hz]
  phqgemis = transpose(phqgemis)
  phqgemis = reverse(phqgemis,1)
  phqgemis *=  2.12d-68 * fremishm # (1. + dblarr(60)); [MeV/cm^3/s]
  phqgemis = interpolate(phqgemis, interpol(indgen(378), enemishm, enph), interpol(indgen(59), redhm, red), /grid) 
  phqgemis1 = interpolate(phqgemis, interpol(indgen(378), enemishm, 1. * ryd), interpol(indgen(59), redhm, red), /grid) 


; -------- HM12 quasar + galaxy background photon flux

  phqgfluxdat = read_ascii('$PUC_PATH/data/qgflux.dat') 
  wlfluxhm = phqgfluxdat.field01[0,1:575] ; [Angstrom]
  frfluxhm = reverse(reform(lightspeed / (wlfluxhm * 1d-8))) ; [Hz]
  enfluxhm = planck * frfluxhm ; [MeV]
  phqgfluxhm = phqgfluxdat.field01[1:59,1:575] ; [erg/cm^2/s/Hz/sr]
  phqgfluxhm = transpose(phqgfluxhm)
  phqgfluxhm = reverse(phqgfluxhm,1)
  phqgfluxhm *=  6.24d5 * frfluxhm # (1. + dblarr(60)); [MeV/cm^3/s/sr]
  phqgfluxhm = interpolate(phqgfluxhm, interpol(indgen(575), enfluxhm, enph), interpol(indgen(59), redhm, red), /grid) 
  phqgfluxhm1 = interpolate(phqgfluxhm, interpol(indgen(575), enfluxhm, 1. * ryd), interpol(indgen(59), redhm, red), /grid) 


; -------- convolve the ICS power with the electron number density to get the emissivity

  phdmemis = dblarr(nenph, nred)
  for c=0, nred-1 do for a=0, nenph-1 do phdmemis[a,c] = integral(enel, icspowercmbr[a,*,c] * psiigmsst[*,c]) ; [MeV/cm^3/s]
  phdmemis1 = interpolate(phdmemis, interpol(indgen(nenph), enph, ryd), ired, /grid)

  if not keyword_set(noplot) then puc_phemis_plot
  if not keyword_set(noplot) then puc_phemis_frame_plot


; -------- calculate effective optical depth

  d2fdndz = dblarr(ncdd, nred)

; Lyman - alpha forest

  jcdd = where(cdd le 10^(15.))

  jred = where(red gt 1.56 and red le 5.5)

  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(7.079) * cdd[jcdd]^(-1.5) # (1. + red[jred])^3

  jred = where(red le 1.56)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(8.238) * cdd[jcdd]^(-1.5) # (1. + red[jred])^0.16

  jred = where(red gt 5.5)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(1.470) * cdd[jcdd]^(-1.5) # (1. + red[jred])^9.9

  jcdd = where(cdd gt 10^(15.) and cdd le 10^(17.5))

  jred = where(red gt 1.56 and red le 5.5)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(14.580) * cdd[jcdd]^(-2.) # (1. + red[jred])^3

  jred = where(red le 1.56)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(15.740) * cdd[jcdd]^(-2.) # (1. + red[jred])^0.16

  jred = where(red gt 5.5)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(8.970) * cdd[jcdd]^(-2.) # (1. + red[jred])^9.9

  jmincddlls = max(jcdd)
  d2fdndzllsmin = d2fdndz[max(jcdd),*]


; Super Lyman limit systems

  jcdd = where(cdd gt 10^(19.) and cdd le 10^(20.3))

  jred = where(red gt 1.56)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(-0.347) * cdd[jcdd]^(-1.05) # (1. + red[jred])^1.27

  jred = where(red le 1.56)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(0.107) * cdd[jcdd]^(-1.05) # (1. + red[jred])^0.16

  jmaxcddlls = min(jcdd)
  d2fdndzllsmax = d2fdndz[min(jcdd),*]


; Lyman limit system

  ncddlls = jmaxcddlls - jmincddlls + 1
  for c=0, nred-1 do $
    d2fdndz[jmincddlls:jmaxcddlls,c] = 10^(findgen(ncddlls) / (ncddlls - 1) * $
    (alog10(d2fdndzllsmax[c]) - alog10(d2fdndzllsmin[c])) + alog10(d2fdndzllsmin[c]))


; Damped Lyman alpha systems

  jcdd = where(cdd gt 10^(20.3))

  jred = where(red gt 1.56)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(18.940) * cdd[jcdd]^(-2.) # (1. + red[jred])^1.27

  jred = where(red le 1.56)
  d2fdndz[min(jcdd):max(jcdd),min(jred):max(jred)] = 10^(19.393) * cdd[jcdd]^(-2.) # (1. + red[jred])^0.16


  opind = cdd # sigmaion

  optot = dblarr(nenph, nred)
  optotint = dblarr(ncdd, nenph, nred)
  for c=0, nred-1 do $
    for a=0, nenph-1 do begin
      optotint[*,a,c] = d2fdndz[*,c] * (1. - exp(-opind[*,a]))
      optot[a,c] = integral(cdd, optotint[*,a,c])
    endfor
  opeff = dblarr(nenph, nred, nred)
  for a=0, nenph-1 do for c=0, nred-1 do for g=0, c-1 do opeff[a,c,g] = integral(red[g:c], optot[a,g:c])

  ; temp !!!!
  ;opeff[*] = 0.
  ; temp !!!!

  if not keyword_set(noplot) then puc_opeff_plot
  

; -------- project the photon emissivity over previous redshifts accounting for the effective optical depth of the clumpy IGM

  phdmflux =  dblarr(nenph, nred)
  phqgflux =  dblarr(nenph, nred)
  phqsflux =  dblarr(nenph, nred)
  for a=0, nenph-1 do begin
    for c=0, nred-2 do begin
      enph_ = enph[a] * (1. + red) / (1. + red[c])
      phdmflux[a,c] = integral(red, puc_phfluxint(phdmemis), xr=[red[c], maxred]) ; [MeV/s/cm^2/sr]
      phqgflux[a,c] = integral(red, puc_phfluxint(phqgemis), xr=[red[c], maxred]) ; [MeV/s/cm^2/sr]
      phqsflux[a,c] = integral(red, puc_phfluxint(phqsemis), xr=[red[c], maxred]) ; [MeV/s/cm^2/sr]
    endfor
  endfor


; get photon fluxes at 1 Rydberg

  phdmflux1 = interpolate(phdmflux, interpol(indgen(nenph), enph, ryd), ired, /grid)
  phqgflux1 = interpolate(phqgflux, interpol(indgen(nenph), enph, ryd), ired, /grid)
  phqsflux1 = interpolate(phqsflux, interpol(indgen(nenph), enph, ryd), ired, /grid)

  if not keyword_set(noplot) then puc_phflux_plot
  if not keyword_set(noplot) then puc_phflux_frame_plot


; -------- convolve the UV photon flux with the photoionization cross section to get the metagalactic ionization rate per ion

  phionrdm = dblarr(nred)
  phionrqg = dblarr(nred)
  phionrqghm = dblarr(nred)
  phionrqs = dblarr(nred)
  for c=0, nred-1 do begin
    phionrdm[c] = integral(enph, puc_phionrint(phdmflux[*,c]), xr=[ryd, enph[nenph-1]]) ; [1/s]
    phionrqg[c] = integral(enph, puc_phionrint(phqgflux[*,c]), xr=[ryd, enph[nenph-1]]) ; [1/s]
    phionrqghm[c] = integral(enph, puc_phionrint(phqgfluxhm[*,c]), xr=[ryd, enph[nenph-1]]) ; [1/s]
    phionrqs[c] = integral(enph, puc_phionrint(phqsflux[*,c]), xr=[ryd, enph[nenph-1]]) ; [1/s]
  endfor
 
  if not keyword_set(noplot) then puc_phionr_plot


; -------- determine the thermal state of the IGM due to photoheating

; blazar heating rate

  plheatbl = 10^(0.0315 * (1. + red)^3 - 0.512 * (1. + red)^2 + 2.27 * (1. + red) - 2.38) / 3.154d22 ; [MeV/s]
  plheatbl[where(red gt 5.7)] = 0.


; photoheating rate

  phheatdm = dblarr(nred)
  phheatqg = dblarr(nred)
  phheatqghm = dblarr(nred)
  phheatqs = dblarr(nred)
  for c=0, nred-1 do begin
    phheatdm[c] = integral(enph, puc_phheatint(phdmflux[*,c]), xr=[ryd, enph[nenph-1]]) ; [MeV/s]
    phheatqg[c] = integral(enph, puc_phheatint(phqgflux[*,c]), xr=[ryd, enph[nenph-1]]) ; [MeV/s]
    phheatqghm[c] = integral(enph, puc_phheatint(phqgfluxhm[*,c]), xr=[ryd, enph[nenph-1]]) ; [MeV/s]
    phheatqs[c] = integral(enph, puc_phheatint(phqsflux[*,c]), xr=[ryd, enph[nenph-1]]) ; [MeV/s]
  endfor

  if not keyword_set(noplot) then puc_heat_plot


; density - temperature relation

  npatch = 100

  
  initred = 20.
  mindif = min(abs(red - initred), jredtemp)
  temp = cmbt

  puc_lambda
  heat = phheatdm
  temppatchdm = dblarr(npatch, nred)
  dpatchdensdzarr = dblarr(npatch, nred)
  patchdensarr = dblarr(npatch, nred)
  for i=0, npatch-1 do begin
    patchdensarr[i,*] = 1. / ((1. + lambda[i,0] * gwf) * (1. + lambda[i,1] * gwf) * (1. + lambda[i,2] * gwf)) - 1.
    dpatchdensdzarr[i,*] = deriv(red, patchdensarr[i,*])
    patchdens = patchdensarr[i,*]
    dpatchdensdz = dpatchdensdzarr[i,*]
    temppatchdm[i,*] = puc_temp_solve()
  endfor    
  if not keyword_set(noplot) then puc_patch_plot 


; thermal evolution of the IGM at the mean overdensity patch

  patchdens = patchdensarr[0,*]
  dpatchdensdz = dpatchdensdzarr[0,*]

  ionfrac = 1d-6

  heat = phheatdm * ionfrac
  tempdm = puc_temp_solve()

  heat = phheatqs * ionfrac
  tempqs = puc_temp_solve()

  heat = plheatbl * ionfrac
  tempbl = puc_temp_solve()

  if not keyword_set(noplot) then puc_temp_plot 



; -------- write results

  writecol, 'puc_dmionr_' + runtag + '.dat', red, phionrdm
  writecol, 'puc_phdmflux_' + runtag + '.dat', enph, phdmflux[*,0]
  writecol, 'puc_mult_' + runtag + '.dat', enel, mult


; -------- calculate computation time

  timeend = systime(1)
  runtime = timeend - timestart
  print, 'puc_main: ', runtime, ' seconds'

  file_mkdir, '$PUC_OUTPUT_PATH/ps' 
  spawn, 'mv ~/puc_*.ps $PUC_OUTPUT_PATH/ps/'

  file_mkdir, '$PUC_OUTPUT_PATH/data' 
  spawn, 'mv ~/puc_*.dat $PUC_OUTPUT_PATH/data/'

end
