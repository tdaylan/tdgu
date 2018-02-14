from __init__ import *

pathbase = os.environ["TDGU_DATA_PATH"] + '/plot_thes/'

figr, axis = plt.subplots(figsize=(6, 6))

# data
path = pathbase + 'vdisdark_mean.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]

path = pathbase + 'vdisdark_mean.csv'
ydatdata = loadtxt(path, delimiter=',')[:, 1]

yerrdata = empty((2, xdatdata.size))

path = pathbase + 'vdisdark_lowr.csv'
yerrdata[0, :] = ydatdata - loadtxt(path, delimiter=',')[:, 1]

path = pathbase + 'vdisdark_uppr.csv'
yerrdata[1, :] = loadtxt(path, delimiter=',')[:, 1] - ydatdata

temp, listcaps, temp = axis.errorbar(xdatdata, ydatdata, yerrdata, color='black', ls='', marker='o', lw=1, capsize=5, markersize=5)
for caps in listcaps:
    caps.set_markeredgewidth(1)
            

# models
path = pathbase + 'vdisdark_0000.csv'
xdat = loadtxt(path, delimiter=',')[:, 0]
ydat = loadtxt(path, delimiter=',')[:, 1]
axis.plot(xdat, ydat, ls='-')

path = pathbase + 'vdisdark_0001.csv'
xdat = loadtxt(path, delimiter=',')[:, 0]
ydat = loadtxt(path, delimiter=',')[:, 1]
axis.plot(xdat, ydat, ls='-')

path = pathbase + 'vdisdark_0002.csv'
xdat = loadtxt(path, delimiter=',')[:, 0]
ydat = loadtxt(path, delimiter=',')[:, 1]
axis.plot(xdat, ydat, ls='-')

axis.set_xlabel('$r$ [kpc]')
axis.set_ylabel('$v$ [km s$^{-1}$]')
path = pathbase + 'vdisdark.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)




# Bhupal-Dev (2013)
figr, axis = plt.subplots(figsize=(6, 6))

# data
path = pathbase + 'Bhupal-Dev2013_fo_lowr.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1)

path = pathbase + 'Bhupal-Dev2013_fo_medi.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1)

path = pathbase + 'Bhupal-Dev2013_fo_uppr.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1)

path = pathbase + 'Bhupal-Dev2013_fo_equi.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1, color='black')
            
axis.set_xlabel('$m_{\chi}/T$')
axis.set_ylabel('$\Omega_{\chi} h^2$')
path = pathbase + 'csecdark.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)


# Markovic (2014)
figr, axis = plt.subplots(figsize=(6, 6))

path = pathbase + 'Markovic2014_lcdm.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1, color='black')

path = pathbase + 'Markovic2014_lowr.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1)

path = pathbase + 'Markovic2014_medi.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1)

path = pathbase + 'Markovic2014_uppr.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.loglog(xdatdata, ydatdata, lw=1)
            
axis.set_xlabel('$k$ [$h$/Mpc]')
axis.set_ylabel('$P(k)$ [(Mpc/$h$)$^3$]')
path = pathbase + 'psecwarm.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)


# Abazajian (2013)
figr, axis = plt.subplots(figsize=(6, 6))

path = pathbase + 'Abazajian2013_lowr.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.plot(xdatdata, ydatdata, lw=1)

path = pathbase + 'Abazajian2013_medi.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.plot(xdatdata, ydatdata, lw=1)

path = pathbase + 'Abazajian2013_uppr.csv'
xdatdata = loadtxt(path, delimiter=',')[:, 0]
ydatdata = loadtxt(path, delimiter=',')[:, 1]
axis.plot(xdatdata, ydatdata, lw=1)
            
axis.set_xscale('log')
axis.set_xlabel('$k$ [$h$/Mpc]')
axis.set_ylabel('$T(k)$')
path = pathbase + 'psecneut.pdf'
plt.tight_layout()
figr.savefig(path)
plt.close(figr)


