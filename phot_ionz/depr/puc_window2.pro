function puc_window2, k

  common puc_common

  x = k * rmass[d]
  j1 = 3. * (sin(x) - x * cos(x)) / x^2
  wk = j1 / x
  return, wk^2
end
