pro gcmsp_ratio_submit

  ratio = gcmsp_ratio_main(10, 0.5, 100, 100, 0, /flat, /uni)
  ratio = gcmsp_ratio_main(10, 0.5, 100, 100, 0.5, /flat, /uni)
  ratio = gcmsp_ratio_main(10, 0.5, 100, 100, 1, /flat, /uni)

  ratio = gcmsp_ratio_main(10, 0.5, 100, 100, 0, /flat)
  ratio = gcmsp_ratio_main(10, 0.5, 100, 100, 0.5, /flat)
  ratio = gcmsp_ratio_main(10, 0.5, 100, 100, 1, /flat)
  
  ratio = gcmsp_ratio_main(10, 1, 100, 100, 0.5, /flat, /uni)
  ratio = gcmsp_ratio_main(10, 2, 100, 100, 0.5, /flat, /uni)
 
  ratio = gcmsp_ratio_main(100, 0.5, 100, 100, 0.5, /flat, /uni)
  ratio = gcmsp_ratio_main(1000, 0.5, 100, 100, 0.5, /flat, /uni)
  ratio = gcmsp_ratio_main(10000, 0.5, 100, 100, 0.5, /flat, /uni)

  ratio = gcmsp_ratio_main(100, 0.5, 100, 100, 0.5, /flat)
  ratio = gcmsp_ratio_main(1000, 0.5, 100, 100, 0.5, /flat)
  ratio = gcmsp_ratio_main(10000, 0.5, 100, 100, 0.5, /flat)

end
