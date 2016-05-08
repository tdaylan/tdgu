function puc_phionrint, phflux

  common puc_common

  phionrint = 4. * !pi * sigmaion * phflux / enph^2 ; [1/MeV/s]

  return, phionrint

end
