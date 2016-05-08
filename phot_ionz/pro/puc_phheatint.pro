function puc_phheatint, phflux

  common puc_common

  phheatint = 4. * !pi * (enph - ryd) * sigmaion * phflux / enph^2 ; [1/s]

  return, phheatint

end
