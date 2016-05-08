function puc_powerspec, k

  common puc_common

  q = k / 0.15d
  cq = 14.4d + 325d / (1d + 60.5d * q^1.11d)
  lq = alog(exp(1.) + 1.84d * q)
  tk = lq / (lq + cq * q^2)
  pk = k^nsubs
  return, pk * tk^2
end
