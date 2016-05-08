function puc_temp_solve

  common puc_common

  if jredtemp eq 0 then begin
    for c=0, nred-2 do temp[c+1] = rk4(temp[c], puc_dtempdred(red[c], temp[c]), red[c], dred[c], 'puc_dtempdred')
  endif
  if jredtemp eq nred-1 then begin
    c = nred - 1
    while c ge 1 do begin
      temp[c-1] = rk4(temp[c], puc_dtempdred(red[c], temp[c]), red[c], -dred[c-1], 'puc_dtempdred')
      c--
    endwhile
  endif
  if jredtemp gt 0 and jredtemp lt nred-1 then begin
    c = jredtemp
    while c ge 1 do begin
      temp[c-1] = rk4(temp[c], puc_dtempdred(red[c], temp[c]), red[c], -dred[c-1], 'puc_dtempdred')
      c--
    endwhile
    for c=jredtemp, nred-2 do temp[c+1] = rk4(temp[c], puc_dtempdred(red[c], temp[c]), red[c], dred[c], 'puc_dtempdred')
  endif

  return, temp
end
