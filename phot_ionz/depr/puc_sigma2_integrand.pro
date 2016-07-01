function puc_sigma2_integrand, k
    
  common puc_common

  pk = puc_powerspec(k)
  wk2 = puc_window2(k)
  si = k^2 * pk * wk2

  return, si

end
