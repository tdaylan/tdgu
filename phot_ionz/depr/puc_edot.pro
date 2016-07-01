function puc_edot, jenel_, jred_, edotincs=edotincs, edotbrem=edotbrem, edotsync=edotsync, $
  magfd_=magfd_, gasdn_=gasdn_, uisrf_=uisrf_, $
  cmbr=cmbr, energy=energy

  common puc_common

  if not keyword_set(energy) then enel_ = enel[jenel_] else enel_ = jenel_

  if keyword_set(cmbr) then return, 4d * thomcs * lightspeed / 3d / elmass^2 * ucmbr[jred_] * enel_^2 ; [MeV/s]

  elgamma = (enel_ + elmass) / elmass
  elbeta = sqrt(1 - 1 / elgamma^2)
  elvel = elbeta * lightspeed

  magfac = 4.966835d-8 ; [(MeV/cm^3)/(muG^2/mu0)] 
  edotsync = 4d * thomcs * lightspeed / 3d / elmass^2 * magfd_ * magfac * enel_^2 ; [MeV/s]
  edotbrem = alphaem * lightspeed * thomcs / 8d / !dpi * (4d * 4.579d1 - 4.446d1) * gasdn_ * enel_ ; [MeV/s]
  edotincs = 4d * thomcs * lightspeed / 3d / elmass^2 * uisrf_ * enel_^2 ; [MeV/s]

  edot = edotsync + edotincs + edotbrem ; [MeV/s]

  return, edot
end
