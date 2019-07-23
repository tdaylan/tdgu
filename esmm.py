from tessmaps import tessmaps as tm

listindxsect = range(13)
coords = []
names = []
is_transiting = False
title = ''
savname = 'save.pdf'
for indxsect in listindxsect:
    tm.make_rect_map(indxsect, coords, names=names,
                 annotate_bools=is_transiting, title=title,
                 bkgnd_cmap='Blues', savname=savname)


