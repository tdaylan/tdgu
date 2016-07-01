function puc_rho, mass_, con_, rcl_, zco_

  common puc_common

  if not keyword_set(con_) then con_ = cmw * (mass_ / mmw)^cmind
 
  rad = sqrt(rcl_^2 # zco_^2) ; [kpc]

  r200 = (3 * mass_ / 4. / !pi / 200. / rhom)^(1./3.) ; [kpc]
  ovdens = 200. * con_^3 / (alog(1. + con_ ) - con_ / (1. + con_))
  rnorm = rad * con_ / r200
  rho = rhoavg * ovdens / rnorm / (1. + rnorm)^2 / kpc2cm^3 ; [MeV/cm^3]

  return, rho
end
