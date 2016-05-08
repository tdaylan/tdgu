pro puc_haloelec_solve

  haloelecredd2r = dblarr(nenel,nenel,nr,nz)
  haloelecreddrr = dblarr(nenel,nenel,nr,nz)
  haloelecredd2z = dblarr(nenel,nenel,nr,nz)
  haloelecredde = dblarr(nenel,nenel,nr,nz)
  difeqmax = -1.
  difeq = 0.
  ntry = 1
  maxiz = max(where(abs(z) lt lz))
  miniz = min(where(abs(z) lt lz))
  pdecount = 0
  while difeq gt difeqmax do begin
    ib = floor(randomu(seed,ntry) * nenel)
    ia = floor(randomu(seed,ntry) * ib)
    ii = floor(randomu(seed,ntry) * nr)
    ij = floor(randomu(seed,ntry) * (maxiz - miniz) + miniz)
    ;print, 'haloelecred: '
    ;print, reform(haloelecred[0,1,*,*])
    for itry=0, ntry-1 do begin
      haloelecred_[ia[itry],ib[itry],ii[itry],ij[itry]] = (randomn(seed)+1.) * haloelecred[ia[itry],ib[itry],ii[itry],ij[itry]]/1d3
      ;print, 'try eezr: '
      ;print, enel[ia[itry]],enel[ib[itry]],r[ii[itry]],z[ij[itry]]
      ;print, 'haloelecred_: '
      ;print, reform(haloelecred_[0,1,*,*])
    endfor
    for a=0, nenel-1 do begin
      for b=0, nenel-1 do begin
        if enel[a] ge enel[b] then continue
        for i=0, nr-1 do begin
          haloelecredd2z[a,b,i,*] = deriv(z, deriv(z, haloelecred_[a,b,i,*]))
        endfor
        for j=0, nz-1 do begin
          if abs(z[j]) gt lz then continue
          haloelecreddr = deriv(r, haloelecred_[a,b,*,j]) 
          haloelecredd2r[a,b,*,j] = deriv(r, haloelecreddr)
          haloelecreddrr[a,b,1:nr-1,j] = haloelecreddr[1:nr-1] / r[1:nr-1]
          haloelecreddrr[a,b,0,j] = haloelecreddr[0] / 0.00001
        endfor
      endfor
    endfor
    for b=0, nenel-1 do begin
      for i=0, nr-1 do begin
        for j=0, nz-1 do begin
          if abs(z[j]) gt lz then continue
          haloelecredde[*,b,i,j] = deriv(eneps, haloelecred_[*,b,i,j] * edotpos[*,i,j] / edotsun) / k0 / tausun / eneps^(delta-2.)
        endfor
      endfor
    endfor
    if pdecount eq 0 then difeqmax = max(haloelecredd2z + haloelecredd2r + haloelecreddrr + haloelecredde) / 1d10
    difeq_ = max(haloelecredd2z + haloelecredd2r + haloelecreddrr + haloelecredde)
    ;print, 'res: '


end
