function gcmsp_ratio_main, nploc, psf, nevt, ncir, fpnt, $
  uniform=uniform, flat=flat, noplot=noplot


  lmax = 10. ; [degrees]
  bmax = 10. ; [degrees]
  ene = 1. ; [GeV]

  nlgl = 100
  nbgl = 100
  lgl = 2. * lmax * findgen(nlgl) / (nlgl - 1) - lmax 
  bgl = 2. * bmax * findgen(nbgl) / (nbgl - 1) - bmax 
  dlgl = 2. * lmax / nlgl
  dbgl = 2. * bmax / nbgl


; -------- generate random events

  evt = fltarr(nevt, 2)
  npnt = floor(nevt * fpnt)
  ndif = nevt - npnt


 ; -------- diffuse component

  if ndif gt 0 then begin
    if keyword_set(uniform) then begin
      evt[0:ndif-1,0] = 2. * lmax * randomu(seed, ndif) - lmax
      evt[0:ndif-1,1] = 2. * bmax * randomu(seed, ndif) - bmax
    endif else begin


; get cumulative probability distribution to sample from the Fermi diffuse model

  gcmsp_time_ranfdm()
      fdmtag = 'gcmsp_fdm_' + string(nlgl, format='(I3.3)') + '_' + string(nbgl, format='(I3.3)') + '.fits' 
      if file_test('$GCMSP_PATH/fits/' + fdmtag) eq 1 then begin
        fdm = mrdfits('$GCMSP_PATH/fits/' + fdmtag, 0)
        fdmlglcdf = mrdfits('$GCMSP_PATH/fits/' + fdmtag, 1)
        fdmbglcdf = mrdfits('$GCMSP_PATH/fits/' + fdmtag, 2)
      endif else begin
        nside = 256L
        npix = 12L * nside^2
        healgen_lb, nside, lglheal, bglheal
        lglheal = ((lglheal + 180.) mod 360.) - 180.
        fdmheal = fermi_diffuse_model(ene, nside=nside, ver=3)
        fdm = griddata(lglheal, bglheal, fdmheal, /sphere, /grid, xout=lgl, yout=bgl, /degrees)
        fdm /= total(fdm)
        mwrfits, fdm, fdmtag, /create
        mwrfits, fdmlglcdf, fdmtag
        mwrfits, fdmbglcdf, fdmtag
        file_mkdir, '$GCMSP_PATH/fits' 
        spawn, 'mv ' + fdmtag + ' $GCMSP_PATH/fits'
      endelse
      fdm /= max(fdm)
      ndif_ = ndif
      while ndif_ gt 0 do begin
        randomarr = randomu(seed, 3 * ndif_)
        lglint = interpol(lgl, findgen(nlgl) / (nlgl - 1), randomarr[0:ndif_-1])
        bglint = interpol(bgl, findgen(nbgl) / (nbgl - 1), randomarr[ndif_:2*ndif_-1])
        fdmint = interpolate(fdm, randomarr[0:ndif_-1] * nlgl, randomarr[ndif_:2*ndif_-1] * nbgl)
        jdif_ = where(fdmint gt randomarr[2*ndif_:3*ndif_-1], njdif_)
        if njdif_ gt 0 then begin
          evt[ndif-ndif_:ndif-ndif_+njdif_-1,0] = lglint[jdif_]
          evt[ndif-ndif_:ndif-ndif_+njdif_-1,1] = bglint[jdif_]
        endif
        ndif_ -= njdif_
      endwhile
    endelse
  endif
  

; -------- point sources

; generate random point source locations

  nfwgamma = 1.3
  nfw = fermi_ring_nfw(nfwgamma)
  nfw /= max(nfw)

  ploc = fltarr(nploc, 2)
  nploc_ = nploc
  while nploc_ gt 0 do begin
    randomarr = randomu(seed, 3 * nploc_)
    lglint = interpol(lgl, findgen(nlgl) / (nlgl - 1), randomarr[0:nploc_-1])
    bglint = interpol(bgl, findgen(nbgl) / (nbgl - 1), randomarr[nploc_:2*nploc_-1])
    nfwint = interpolate(nfw, randomarr[0:nploc_-1] * nlgl, randomarr[nploc_:2*nploc_-1] * nbgl)
    jploc_ = where(nfwint gt randomarr[2*nploc_:3*nploc_-1], njploc_)
    if njploc_ gt 0 then begin
      ploc[nploc-nploc_:nploc-nploc_+njploc_-1,0] = lglint[jploc_]
      ploc[nploc-nploc_:nploc-nploc_+njploc_-1,1] = bglint[jploc_]
    endif
    nploc_ -= njploc_
  endwhile


; generate random number of events per point source 

  if keyword_set(flat) then begin
    iploc = floor(randomu(seed, npnt) * nploc) 
  endif else begin
    lumiind = -2
    lfpdf = lumi^lumiind
    while not (npnt_ eq npnt) do begin
    iploc = 
  endelse


; generate random event locations

  if npnt gt 0 then begin
    evt[ndif:nevt-1,0] = ploc[iploc,0]
    evt[ndif:nevt-1,1] = ploc[iploc,1]
  endif


; generate random locations

  cir = fltarr(ncir,2)
  cir[*,0] = randomu(seed, ncir) * 20. - 10.
  cir[*,1] = randomu(seed, ncir) * 20. - 10.


; match events to circles

  t1 = systime(1)
  spherematch, cir[*,0], cir[*,1], evt[*,0], evt[*,1], psf, imcir, imevt, distec, maxmatch=0
  t2 = systime(1)
  print, 'spherematch: ', t2 - t1

  if imcir[0] eq -1 then begin
    nmcir = 0
    mcir = -1
    distec = -1
  endif else begin
    imcirunique = imcir[uniq(imcir, sort(imcir))]
    nmcir = n_elements(imcirunique)
    mcir = fltarr(nmcir, 2)
    mcir[*,0] = cir[imcirunique,0]
    mcir[*,1] = cir[imcirunique,1]
  endelse
  ncirneig = fltarr(ncir)
  for k=0, nmcir-1 do ncirneig[imcir[k]]++
 
 
; match events within a PSF

  t1 = systime(1)
  spherematch, evt[*,0], evt[*,1], evt[*,0], evt[*,1], psf, imevt, imevt, distee, maxmatch=0
  t2 = systime(1)
  print, 'spherematch: ', t2 - t1
  

  nmevt = n_elements(imevt)
  if nmevt eq nevt then begin
    nmevt = 0
    mevt = -1
    distee = -1
  endif else begin
    imevt = imevt[nevt:nmevt-1]
    distee = distee[nevt:nmevt-1]
    imevtunique = imevt[uniq(imevt, sort(imevt))]
    nmevt = n_elements(imevtunique)
    mevt = fltarr(nmevt, 2)
    mevt[*,0] = evt[imevtunique,0]
    mevt[*,1] = evt[imevtunique,1]
  endelse
  nevtneig = fltarr(nevt)
  for k=0, nmevt-1 do nevtneig[imevt[k]]++


; calculate the test statistic

  ratio = float(nevt - nmevt) / float(ncir - nmcir)


; plot

  if not keyword_set(noplot) then gcmsp_ratio_plot, ratio, fpnt, nploc, psf, evt, cir, mevt, mcir, distee, distec, nevtneig, ncirneig, $
      fdmlglcdf, fdmbglcdf, lgl, bgl, fdm, $
      uniform=uniform, flat=flat


  return, ratio

end
