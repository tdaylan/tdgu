function puc_dndecmbr, enpi_, c_

  common puc_common

  return, 8. * !pi / (lightspeed * planck)^3 * enpi_^2 / (exp(enpi_ / kt[c_]) - 1.) ; [1/cm^3/MeV]

end
