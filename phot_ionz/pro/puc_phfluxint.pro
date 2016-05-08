function puc_phfluxint, phemis

  common puc_common

  phemisint = interpolate(phemis, interpol(findgen(nenph), enph, enph_), ired, /grid) > 0. ; [MeV/cm^3/s]
  phemisint = phemisint[ired,ired]
  phfluxint = hrad / 4. / !pi * phemisint * exp(-opeff[a,c,*]) * (1. + red[c])^3 / (1. + red)^4 / hubbf  ; [MeV/cm^2/s/sr]

  return, phfluxint

end
