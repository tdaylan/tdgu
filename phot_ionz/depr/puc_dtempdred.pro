function puc_dtempdred, red_, temp_

  common puc_common

  dtempdred1_ = 2. * temp_ / (1. + red_)

  dtempdred2_ = 2. * temp_ / 3. / (1. + interpol(patchdens, red, red_)) * interpol(dpatchdensdz, red, red_)

  dtempdred3_ = 2. * interpol(heat, red, red_) / 3. / boltzmann / (1. + red_) / interpol(hubbf, red, red_) * htime

  dtempdred_ = dtempdred1_ + dtempdred2_ + dtempdred3_

  return, dtempdred_

end
