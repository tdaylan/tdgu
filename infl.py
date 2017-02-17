import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

path = '/Users/tansu/Desktop/infl/'
os.system('mkdir -p %s' % path)

#os.environ["FERM_IGAL_DATA_PATH"]
cambdata = camb.get_results(pars)
psec = cambdata.get_cmb_power_spectra(pars)
psectotl = psec['total']
pseculen = psec['unlensed_scalar']

ls = np.arange(psectotl.shape[0])
fig, ax = plt.subplots(2,2, figsize = (12,12))
ax[0,0].plot(ls,psectotl[:,0], color='k')
ax[0,0].plot(ls,pseculen[:,0], color='r')
ax[0,0].set_title('TT')
ax[0,1].plot(ls[2:], 1-pseculen[2:,0]/psectotl[2:,0]);
ax[0,1].set_title(r'$\Delta TT$')
ax[1,0].plot(ls,psectotl[:,1], color='k')
ax[1,0].plot(ls,pseculen[:,1], color='r')
ax[1,0].set_title(r'$EE$')
ax[1,1].plot(ls,psectotl[:,3], color='k')
ax[1,1].plot(ls,pseculen[:,3], color='r')
ax[1,1].set_title(r'$TE$');
for ax in ax.reshape(-1):
    ax.set_xlim([2,2500])
plt.savefig(path + 'psec.pdf')
plt.close()


pars.WantTensors = True
cambdata = camb.get_transfer_functions(pars)
lmax=2000
rs = np.linspace(0,0.2,6)
for r in rs:
    inflation_params = initialpower.InitialPowerParams()
    inflation_params.set_params(ns=0.96, r=r)
    cambdata.power_spectra_from_transfer(inflation_params)
    cl = cambdata.get_total_cls(lmax)
    plt.loglog(np.arange(lmax+1),cl[:,2])
plt.xlim([2,lmax])
plt.legend(rs, loc='lower right');
plt.savefig(path + 'tens.pdf')
plt.close()

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.set_dark_energy() #re-set defaults
pars.InitPower.set_params(ns=0.965)
pars.set_matter_power(redshifts=[0., 0.8], kmax=2.0)
pars.NonLinear = model.NonLinear_none
cambdata = camb.get_results(pars)
kh, z, pk = cambdata.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
s8 = np.array(cambdata.get_sigma8())

pars.NonLinear = model.NonLinear_both
cambdata.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = cambdata.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
for i, (redshift, line) in enumerate(zip(z,['-','--'])):
    plt.loglog(kh, pk[i,:], color='k', ls = line)
    plt.loglog(kh_nonlin, pk_nonlin[i,:], color='r', ls = line)
plt.xlabel('k/h Mpc');
plt.legend(['linear','non-linear'], loc='lower left');
plt.title('Matter power at z=%s and z= %s'%tuple(z));
plt.savefig(path + 'psecnonl.pdf')
plt.close()

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.InitPower.set_params(As=2e-9, ns=0.965)
pars.set_for_lmax(2000, lens_potential_accuracy=1)
ws = np.linspace(-1.5, -0.6, 5)
for w in ws:
    pars.set_dark_energy(w) 
    cambdata = camb.get_results(pars)
    cl = cambdata.get_lens_potential_cls(lmax=2000)
    plt.loglog(np.arange(2001), cl[:,0])
plt.savefig(path + 'psecdene.pdf')
plt.close()

pars.set_dark_energy() 
plt.legend(ws)
plt.ylabel('$[L(L+1)]^2C_L^{\phi\phi}/2\pi$')
plt.xlabel('$L$')
plt.xlim([2,2000]);
plt.savefig(path + 'dene.pdf')
plt.close()

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
cambdata = camb.get_background(pars)
z = np.linspace(0,5,100)
DA = cambdata.angular_diameter_distance(z)

plt.plot(z, DA)
coef = np.polyfit(z, DA, 10)
print 'coef'
print coef
plt.plot(z, np.poly1d(coef)(z))
import h5py, os
h5f = h5py.File(os.environ['PCAT_DATA_PATH'] + '/data/inpt/adiscoef.h5', 'w')
h5f.create_dataset('adiscoef', data=coef)
h5f.close()
print 'np.poly1d(coef)(z)'
print np.poly1d(coef)(z)
print 'DA'
print DA
plt.xlabel('$z$')
plt.ylabel(r'$D_A /\rm{Mpc}$')
plt.title('Angular diameter distance')
plt.ylim([0,2000]);
plt.savefig(path + 'distangl.pdf')
plt.close()

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
data = camb.get_transfer_functions(pars)
transfer = data.get_cmb_transfer_data()
fig, axs = plt.subplots(2,2, figsize=(12,8), sharex = True)
for ix, ax in zip([3, 20, 40, 60],axs.reshape(-1)):
    ax.plot(transfer.q,transfer.delta_p_l_k[0,ix,:])
    ax.set_title(r'$\ell = %s$'%transfer.l[ix])
    if ix>1: ax.set_xlabel(r'$k \rm{Mpc}$')
plt.savefig(path + 'tran.pdf')
plt.close()

trans = transfer
ix = 0
_, axs = plt.subplots(1,2, figsize=(12,6))
for source_ix, (name, ax) in enumerate(zip(['T', 'E'], axs)):
    ax.semilogx(trans.q,trans.delta_p_l_k[source_ix,ix,:])
    ax.set_xlim([1e-5, 0.05])
    ax.set_xlabel(r'$k \rm{Mpc}$')
    ax.set_title(r'%s transfer function for $\ell = %s$'%(name, trans.l[ix]))
plt.savefig(path + 'trantran.pdf')
plt.close()

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.WantScalars = False
pars.WantTensors = True
pars.set_accuracy(AccuracyBoost=2)
data = camb.get_transfer_functions(pars)
transfer = data.get_cmb_transfer_data('tensor')
plt.figure(figsize=(14,3))
ixs=[13,19,21]
ls = [transfer.l[i] for i in ixs]
cols=['b','r','c']
for ix,col in zip(ixs, cols):
    k_weight = transfer.delta_p_l_k[2,ix,:]**2
    k_weight /= np.sum(k_weight)
    plt.semilogx(transfer.q,k_weight, color=col)
plt.xlim([1e-3, 0.1])
plt.legend(ls)
plt.xlabel(r'$k \rm{Mpc}$')
plt.title(r'Contribution to B from primordial tensor power spectrum for various $\ell$')
derived = data.get_derived_params()
for l,col in zip(ls,cols):
    plt.axvline(l/(1000*derived['DAstar']), color=col, ls=':', lw=2)
plt.savefig(path + 'bmod.pdf')
plt.close()

k=10**np.linspace(-5, 1, 50)
pars.InitPower.set_params(ns=0.96, r=0.2) #this functions imposes inflation consistency relation by default
scalar_pk= pars.scalar_power(k)
tensor_pk= pars.tensor_power(k)
plt.semilogx(k,scalar_pk);
plt.semilogx(k,tensor_pk);
plt.xlabel(r'$k \rm{Mpc}$')
plt.ylabel(r'${\cal P}(k)$')
plt.legend(['scalar', 'tensor']);
plt.savefig(path + 'scaltens.pdf')
plt.close()

pars = camb.set_params(H0=67.5, ombh2=0.022, omch2=0.122, As=2e-9, ns=0.95)
data= camb.get_background(pars)
eta = 10**(np.linspace(1, 4,300))
back_ev = data.get_background_time_evolution(eta, ['x_e', 'visibility'])
fig, axs= plt.subplots(1,2, figsize=(12,5))
axs[0].semilogx(eta, back_ev['x_e'])
axs[1].loglog(eta, back_ev['visibility'])
axs[0].set_xlabel(r'$\eta/\rm{Mpc}$')
axs[0].set_ylabel('$x_e$')
axs[1].set_xlabel(r'$\eta/\rm{Mpc}$')
axs[1].set_ylabel('Visibility');
fig.suptitle('Ionization history, including both hydrogen and helium recombination and reionization');
plt.savefig(path + 'ionz.pdf')
plt.close()

z = 10**np.linspace(2, 4, 300)
back_ev = data.get_background_redshift_evolution(z, ['x_e', 'visibility'], format='array')
fig, axs= plt.subplots(1,2, figsize=(12,5))
for i, (ax, label), in enumerate(zip(axs, ['$x_e$','Visibility'])):
    ax.semilogx(z, back_ev[:,i])
    ax.set_xlabel('$z$')
    ax.set_ylabel(label)
    ax.set_xlim([500,1e4])
plt.savefig(path + 'visi.pdf')
plt.close()

eta = np.linspace(1, 400, 300)
ks = [0.02,0.1]
ev = data.get_time_evolution(ks, eta, ['delta_baryon','delta_photon'])
_, axs= plt.subplots(1,2, figsize=(12,5))
for i, ax in enumerate(axs):
    ax.plot(eta,ev[i,:, 0])
    ax.plot(eta,ev[i,:, 1])
    ax.set_title('$k= %s$'%ks[i])
    ax.set_xlabel(r'$\eta/\rm{Mpc}$');
plt.legend([r'$\Delta_b$', r'$\Delta_\gamma$'], loc = 'upper left');
plt.savefig(path + 'bary.pdf')
plt.close()

z = np.linspace(500,5000,300)
ks = [0.02,0.1]
ev = data.get_redshift_evolution(ks, z, ['delta_baryon','delta_cdm', 'delta_photon'])
_, axs= plt.subplots(1,2, figsize=(12,5))
for i, ax in enumerate(axs):
    ax.plot(z,ev[i,:, 0])
    ax.plot(z,ev[i,:, 1])
    ax.plot(z,ev[i,:, 2])
    ax.set_title(r'$k= %s/\rm{Mpc}$'%ks[i])
    ax.set_xlabel('$z$');
plt.legend([r'$\Delta_b$', r'$\Delta_c$', r'$\Delta_\gamma$'], loc = 'upper right');
plt.savefig(path + 'barybary.pdf')
plt.close()

eta = 10**(np.linspace(0, 3, 500))
def plot_ev(ev, k):
    plt.figure(figsize=(8,6))
    plt.loglog(eta,ev[:,0])
    plt.loglog(eta,np.abs(ev[:,1]))
    plt.loglog(eta,-ev[:,2])
    plt.title(r'$k= %s/\rm{Mpc}$'%k)
    plt.xlabel(r'$\eta/\rm{Mpc}$');
    plt.legend([r'$\Delta_c$', r'$|\Delta_\gamma|$', r'$-(\Phi+\Psi)/2$'], loc = 'upper left');
    plt.savefig(path + 'etaa%d.pdf' % k)
    plt.close()

k=0.3
plot_ev(data.get_time_evolution(k, eta, ['delta_cdm','delta_photon', 'Weyl']),k)
plot_ev(data.get_time_evolution(k, eta, ['delta_cdm','delta_photon', 'Weyl'],lAccuracyBoost=1),k)
plot_ev(data.get_time_evolution(k, eta, ['delta_cdm','delta_photon', 'Weyl'],lAccuracyBoost=10),k)

nz = 100 #number of steps to use for the radial/redshift integration
kmax=10  #kmax to use
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.InitPower.set_params(ns=0.965)

cambdata= camb.get_background(pars)
chistar = cambdata.conformal_time(0)- model.tau_maxvis.value
chis = np.linspace(0,chistar,nz)
zs=cambdata.redshift_at_comoving_radial_distance(chis)
dchis = (chis[2:]-chis[:-2])/2
chis = chis[1:-1]
zs = zs[1:-1]

PK = camb.get_matter_power_interpolator(pars, nonlinear=True, 
    hubble_units=False, k_hunit=False, kmax=kmax,
    var1=model.Transfer_Weyl,var2=model.Transfer_Weyl, zmax=zs[-1])

plt.figure(figsize=(8,5))
k=np.exp(np.log(10)*np.linspace(-4,2,200))
zplot = [0, 0.5, 1, 4 ,20]
for z in zplot:
    plt.loglog(k, PK.P(z,k))
plt.xlim([1e-4,kmax])
plt.xlabel('k Mpc')
plt.ylabel('$P_\Psi\, Mpc^{-3}$')
plt.legend(['z=%s'%z for z in zplot]);
plt.xlabel('$L$');
plt.savefig(path + 'psii.pdf')

 
win = ((chistar-chis)/(chis**2*chistar))**2
ls = np.arange(2,2500+1, dtype=np.float64)
cl_kappa=np.zeros(ls.shape)
w = np.ones(chis.shape) #this is just used to set to zero k values out of range of interpolation
for i, l in enumerate(ls):
    k=(l+0.5)/chis
    w[:]=1
    w[k<1e-4]=0
    w[k>=kmax]=0
    cl_kappa[i] = np.dot(dchis, w*PK.P(zs, k, grid=False)*win/k**4)
cl_kappa*= (ls*(ls+1))**2
pars.set_for_lmax(2500,lens_potential_accuracy=2)
cambdata = camb.get_results(pars)
cl_camb=cambdata.get_lens_potential_cls(2500) 

cl_limber= 4*cl_kappa/2/np.pi #convert kappa power to [l(l+1)]^2C_phi/2pi (what cl_camb is)
plt.loglog(ls,cl_limber, color='b')
plt.loglog(np.arange(2,cl_camb[:,0].size),cl_camb[2:,0], color='r')
plt.xlim([1,2000])
plt.legend(['Limber','CAMB hybrid'])
plt.ylabel('$[L(L+1)]^2C_L^{\phi}/2\pi$')
plt.xlabel('$L$');
plt.savefig(path + 'kapp.pdf')
plt.close()

camb.set_halofit_version('takahashi')
kh_nonlin, _, pk_takahashi = cambdata.get_nonlinear_matter_power_spectrum(params=pars)
camb.set_halofit_version('mead')
kh_nonlin, _, pk_mead = cambdata.get_nonlinear_matter_power_spectrum(params=pars)
fig, axs=plt.subplots(2,1, sharex=True, figsize=(8,8))
axs[0].loglog(kh_nonlin, pk_takahashi[0])
axs[0].loglog(kh_nonlin, pk_mead[0])
axs[1].semilogx(kh_nonlin, pk_mead[0]/pk_takahashi[0]-1)
axs[1].set_xlabel(r'$k/h\, \rm{Mpc}$')    
axs[1].legend(['Mead/Takahashi-1'], loc='upper left')
plt.savefig(path + 'halo.pdf')
plt.close()



