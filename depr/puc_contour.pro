pro puc_contour, data, x, y, nlevels=nlevels, xlog=xlog, ylog=ylog, $
   xtitle=xtitle, ytitle=ytitle, title=title, pos=pos, log=log


   if keyword_set(x) then xr = [min(x), max(x)]
   if keyword_set(y) then yr = [min(y), max(y)]

   if keyword_set(log) then begin
     ibad = where(data le 0.)
     igood = where(data gt 0.)
     data[igood] = alog10(data[igood])
     data[ibad] = !Values.F_NaN
   endif

   if not keyword_set(nlevels) then nlevels = 256
   if not keyword_set(ndivisions) then ndivisions = 5

   if not keyword_set(pos) then begin
     pos =   [0.2, 0.1, 0.9, 0.81]
     cbpos = [0.2, 0.88, 0.9, 0.95]
   endif else begin
     cbpos = pos
     pos[3] -= 0.1
     cbpos[1] = pos[3]
   endelse

   cgLoadCT, 1, ncolors=nlevels
  
   min = min(data)
   max = max(data)
   if keyword_set(range) then begin
      min = range[0] > min 
      max = range[1] < max
   endif 
   levels = findgen(nlevels) / (nlevels-1) * (max - min) + min

   cgcontour, data, x, y, levels=levels, C_Colors=Bindgen(nLevels)+1B, xr=xr, yr=yr, $
      pos=pos, xtitle=xtitle, ytitle=ytitle, xlog=xlog, ylog=ylog, missingvalue=!Values.F_NaN, /fill
   
   cgColorbar, ncolors=nlevels, pos=cbpos, $
      Range=[min, max], divisions=ndivisions, $
      title=title, TLocation='Top'

end
