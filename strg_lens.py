from __init__ import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy import signal
from astropy.convolution import convolve, AiryDisk2DKernel

class Params(object):
    """
    The main parameter class, used to define 
    the lensmodel and source parameters. It
    contains the parameter name, its value,
    a default range in which the parameter
    should lie (used as prior for running MCMC),
    and a boolean flag determining whether the 
    parameter should be varied in an MCMC.
    """
    
    def __init__(self,name,value,prior,tv):
        """
        Initialization of a parameter.

        Args:
            name: A string containing the name of the parameter.
            value: The parameter value.
            prior: a 1D array containing the lower and upper bounds of the parameters.
            tv: a boolean operator determining whether we want to vary the parameter in MCMC.

        Returns:
            An instance of the parameter class.

        Raises:
            TypeError: if name is not a string.
            TypeError: if value is not a float.
            TypeError: if prior is not a list
            TypeError: if tv is not a boolean operator
        """
        #Check for type errors
        if not isinstance(name,str):
            raise TypeError('The parameter name must be a string.')
        if not isinstance(value,float):
            raise TypeError('The parameter value must be a float.')
        if not isinstance(prior,list):
            raise TypeError('The prior range is not a list of an upper and lower bounds.')
        if not isinstance(tv,bool):
            raise TypeError('The tv attribute must be a boolean operator.')
        
        #Assign values
        self.name  = name
        self.value = value
        self.prior = np.array(prior)
        self.tv = tv
        

class LensModel(object):
    """The main lens model class"""
    
    def __init__(self, massmodel, xgal, ygal, ellipticity, ellipt_angle, shear, shear_angle, *args):
        self.massmodel = massmodel
        #Primary lens parameters
        self.lensparams = []
        self.lensparams.append(Params('xgal',xgal,[-2.0,2.0],False))
        self.lensparams.append(Params('ygal',ygal,[-2.0,2.0],False))
        self.lensparams.append(Params('ellipticity',ellipticity,[0,0.5],False))
        self.lensparams.append(Params('ellipticity_angle',ellipt_angle,[0.0,180.0],False)) #Note that this angle is measured counterclockwise from the y-axis
        self.lensparams.append(Params('shear',shear,[0.0,0.3],False))
        self.lensparams.append(Params('shear_angle',shear_angle,[0.0,180],False)) #This angle is measured c.c. from x-axis
        self.s = 0.0001
        
        if self.massmodel == 'SIE':
            self.lensparams.append(Params('b_ein',args[0],[0.2,2.0],False))
        
        elif self.massmodel == 'alpha':
            self.lensparams.append(Params('b_ein',args[0],[0.2,2.0],False))
            self.lensparams.append(Params('alpha',args[1],[0.5,1.5],False))
            print('Note: alpha lens model has not yet been coded up')
            
        else:
            raise RuntimeError('Type of lens unrecognized!')
            
        self.numlensparams = len(self.lensparams)
            
    def potential(self,xx,yy):
        """The projected gravitational potential method"""
        
        #Translate to the center of the potential
        xxt = xx - self.lensparams[0].value
        yyt = yy - self.lensparams[1].value
        
        #Define rotation angle
        rot_angle = self.lensparams[3].value*np.pi/180.0
        
        #Rotate coordinates according to given angle
        x = -np.sin(rot_angle)*xxt + np.cos(rot_angle)*yyt
        y = -np.cos(rot_angle)*xxt - np.sin(rot_angle)*yyt
        
        #Define q
        q = 1.0 - self.lensparams[2].value
        
        #Compute the potential
        if self.massmodel == 'SIE':
            if self.lensparams[2].value > 1.0e-4:    
                psi = np.sqrt(q**2*(self.s**2+x**2)+y**2)
                phix = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                        np.arctan(np.sqrt(1.0-q**2)*x/(psi+self.s)))
                phiy = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                        np.arctanh(np.sqrt(1.0-q**2)*y/(psi+q**2*self.s)))
                phi = (x*phix + y*phiy-self.lensparams[6].value*q*self.s*np.log(np.sqrt((psi+self.s)**2+(1.0-q**2)*x**2)) 
                      +self.lensparams[6].value*q*self.s*np.log((1+q)*self.s))
            else:
                phi = self.lensparams[6].value*(np.sqrt(x**2 + y**2 + self.s**2)- self.s - 
                        self.s*np.log(0.50 + np.sqrt(x**2 + y**2 + self.s**2)/(2.0*self.s)))
           
                   
        elif self.massmodel == 'alpha':
            phi = 0.0
        
        #Add the external shear contribution
        phi += (0.5*self.lensparams[4].value*np.cos(2.0*self.lensparams[5].value*np.pi/180.0)*(xx**2 - yy**2) 
                    + self.lensparams[4].value*np.sin(2.0*self.lensparams[5].value*np.pi/180.0)*xx*yy)
        
        return phi
    
    def deflection(self,xx,yy):
        """The deflection vector at position x,y"""
        
        #Translate to the center of the potential
        xxt = xx - self.lensparams[0].value
        yyt = yy - self.lensparams[1].value
        
        #Define rotation angle
        rot_angle = self.lensparams[3].value*np.pi/180.0
        
        #Rotate coordinates according to given angle
        x = -np.sin(rot_angle)*xxt + np.cos(rot_angle)*yyt
        y = -np.cos(rot_angle)*xxt - np.sin(rot_angle)*yyt
        
        #Define q
        q = 1.0 - self.lensparams[2].value
        
        if self.massmodel == 'SIE':
            if self.lensparams[2].value > 1.0e-4:
                psi = np.sqrt(q**2*(self.s**2+x**2)+y**2)
                alphaxt = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                            np.arctan(np.sqrt(1.0-q**2)*x/(psi+self.s)))
                alphayt = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                            np.arctanh(np.sqrt(1.0-q**2)*y/(psi+q**2*self.s)))
            else:
                rint = np.sqrt(x**2 + y**2 + self.s**2)
                alphaxt = self.lensparams[6].value*(x/(rint+self.s)) 
                alphayt = self.lensparams[6].value*(y/(rint+self.s))
        
        elif self.massmodel == 'alpha':
            alphaxt = 0.0
            alphayt = 0.0
            
        #Rotate back vector to original basis
        alphax = -np.cos(rot_angle)*alphayt-np.sin(rot_angle)*alphaxt
        alphay = np.cos(rot_angle)*alphaxt-np.sin(rot_angle)*alphayt
        
        #Add the external shear contribution
        alphax += (self.lensparams[4].value*np.cos(2.0*self.lensparams[5].value*np.pi/180.0)*xx 
                + self.lensparams[4].value*np.sin(2.0*self.lensparams[5].value*np.pi/180.0)*yy)
        alphay += (-self.lensparams[4].value*np.cos(2.0*self.lensparams[5].value*np.pi/180.0)*yy 
                + self.lensparams[4].value*np.sin(2.0*self.lensparams[5].value*np.pi/180.0)*xx)
        
        return np.array([alphax,alphay])
    
def gauss_mat(size,axis_ratio,angle):
    """
    Compute the covariance matrix of a Gaussian, given a typical size,
    an axis ratio, and an angle
    """
    rot_mat = np.array([[np.cos(angle*np.pi/180.0),-np.sin(angle*np.pi/180.0)],
                        [np.sin(angle*np.pi/180.0),np.cos(angle*np.pi/180.0)]])
    pre_mat = np.array([[1.0/(axis_ratio*size)**2,0.0],[0.0,1.0/size**2]])
    
    return np.dot(np.transpose(rot_mat),np.dot(pre_mat,rot_mat))
    
           
class Source(object):
    """The main source class"""
        
    def __init__(self,model,usrc,vsrc,peak_bright,*args):
        self.model = model
        self.srcparams = []
        self.srcparams.append(Params('usrc',usrc,[-1.0,1.0],False))
        self.srcparams.append(Params('vsrc',vsrc,[-1.0,1.0],False))
        self.srcparams.append(Params('peak_brightness',peak_bright,[0.01,10],False))
        
        if self.model == 'Gaussian':
            self.srcparams.append(Params('src_size',args[0],[0.001,1.0],False))
            #How much is the Gaussian stretched along the x-axis compared to the y-axis
            self.srcparams.append(Params('src_axis_ratio',args[1],[1.0,7.0],False)) 
            self.srcparams.append(Params('src_angle',args[2],[0.0,180.0],False))
            self.cov_src = gauss_mat(args[0],args[1],args[2])
        else:
            raise RuntimeError('Type of source unrecognized!')
        
        self.numsrcparams = len(self.srcparams)
                 
    def brightness(self,u,v):
        """The source brightness profile"""
        
        if self.model == 'Gaussian':
            src_pos = np.array([u-self.srcparams[0].value,v-self.srcparams[1].value])
            S = self.srcparams[2].value*np.exp(-0.5*np.sum(src_pos*np.tensordot(self.cov_src,src_pos,(1,0)),0))
        else:
            S = 0.0
            
        return S
        
    def gradient(self,u,v):
        """The gradient of the source in the source plane"""
        
        if self.model == 'Gaussian':
            src_pos = np.array([u-self.srcparams[0].value,v-self.srcparams[1].value])
            shift = np.tensordot(self.cov_src,src_pos,(1,0))
            exp_arg = np.sum(src_pos*shift,0)
            dSdu = -self.srcparams[2].value*shift[0]*np.exp(-0.5*exp_arg)
            dSdv = -self.srcparams[2].value*shift[1]*np.exp(-0.5*exp_arg)
        else:
            dSdu = 0.0
            dSdv = 0.0
            
        return np.array([dSdu,dSdv])
    
    
class PixelMap(object):
    """A general instance of a pixelization of an image"""
    
    def __init__(self,newxmax,newdx):
        self.xmax = newxmax
        self.dx = newdx
        self.numpix1D = int(np.floor(2.0*self.xmax/self.dx) + 1)
        self.xmap = np.arange(-self.xmax,self.xmax+self.dx,self.dx)
        self.ymap = np.arange(-self.xmax,self.xmax+self.dx,self.dx)
        
def macro_only_image(pix,src,smoothlens,psf_scale):
    """
    Generate an image of a source lensed by a smooth mass model:
    pix is an instance of the PixelMap class
    src is an instance of the Source class
    smoothlens is an instance of the SIE class
    psf_scale is the size of the Airy disk PSF in pixel
    """
    
    #Set the PSF kernel
    PSF_kernel = AiryDisk2DKernel(psf_scale)
    
    #Check if we have more than one sources
    if not isinstance(src,list):
        source = [src]
    else:
        source = src

    #Generate lensed image
    xx, yy = np.meshgrid(pix.xmap, pix.ymap)
    defxmap,defymap = smoothlens.deflection(xx,yy)
    piximg0 = np.sum([source[kk].brightness(xx-defxmap,yy-defymap) for kk in range(len(source))],0)

    #Convolve with the PSF
    conv_piximg0 = convolve(piximg0,PSF_kernel)

    return conv_piximg0

def macro_only_image_no_conv(pix,src,smoothlens):
    """
    Generate an image of a source lensed by a smooth mass model. No PSF convolution.
    
    Args:
        pix is an instance of the PixelMap class.
        src is an array of instances of the Source class.
        smoothlens is an instance of the LensModel class.
        
    Returns:
        A numpy array containing the image brightness.
        
    Raises:
        TypeError if pix is not a PixelMap instance.
        TypeError if src is not a Source instance.
        TypeError if smoothlens is not a LensModel instance. 
    """

    if not isinstance(src,list):
        src = [src]
    else:
        src = src
    #Generate lensed image
    xx, yy = np.meshgrid(pix.xmap, pix.ymap)
    defxmap,defymap = smoothlens.deflection(xx,yy)
    piximg0 = np.sum([src[kk].brightness(xx-defxmap,yy-defymap) for kk in range(len(src))],0)
    
    return piximg0

def smooth_potdef_map(pix,smoothlens):
    """Generate a map of the smooth lens potential"""
    
    xx, yy = np.meshgrid(pix.xmap, pix.ymap, sparse=True)
    potmap = smoothlens.potential(xx,yy)
    defxmap, defymap = smoothlens.deflection(xx,yy)
                
    return potmap, defxmap, defymap

def smooth_def_map(pix,smoothlens):
    """Generate a deflection map for the smooth lens potential"""
    
    xx, yy = np.meshgrid(pix.xmap, pix.ymap, sparse=True)
    defxmap,defymap = smoothlens.deflection(xx,yy)
            
    return defxmap, defymap

def source_bright_map(pix,src):
    """Generate a brigthness map of the source in the source plane"""
    
    #Check if we have more than one sources
    if not isinstance(src, list):
        source = [src]
    else:
        source = src
        
    u, v = np.meshgrid(pix.xmap, pix.ymap, sparse=True)
    brightness_map = np.sum([source[kk].brightness(u,v) for kk in range(len(source))],0)
    grd_x_src, grd_y_src = np.sum([source[kk].gradient(u,v) for kk in range(len(source))],0)
                
    return brightness_map, grd_x_src, grd_y_src

def source_grad_map(pix,src,smoothlens):
    """Generate a gradient map of the source, evaluated at the image position"""
    
    xx, yy = np.meshgrid(pix.xmap, pix.ymap)
    defvecx,defvecy = smoothlens.deflection(xx,yy)
    grd_x_src,grd_y_src = np.sum([src[kk].gradient(xx - defvecx,yy - defvecy) for kk in range(len(src))],0)
                
    return grd_x_src, grd_y_src



def time_func(func, *args):
    
    numbiter = 4
    timediff = empty(numbiter)
    for k in range(numbiter):
        timeinit = time.time()
        func(*args)
        timediff[k] = time.time() - timeinit
    print 'Calling %s' % func.__name__ 
    print '%3g pm %3g ms' % (mean(timediff) * 1e3, std(timediff) * 1e3)
    print
    return timediff

#Create a pixelated image of a lensed image
Noiselevel = 0.0005
PSF_scale = 3

#Set up the structure of the pixelated image
mypix = PixelMap(49.5, 1.)

#Set up the source
mysrc = Source('Gaussian',-0.05,0.22,1.0,0.07,1.0,0.0)

#Set up the lens
mylens = LensModel('SIE',0.0,0.0,0.0,0.0,0.150,-18.435,0.89)

#Get the unperturbed image
unpert_img=macro_only_image_no_conv(mypix,mysrc,mylens)

# calculate the time performance
time_func(PixelMap, 1.9999, 0.05)
time_func(Source, 'Gaussian', -0.05, 0.22, 1.0, 0.07, 1.0, 0.0)
time_func(LensModel, 'SIE', 0., 0., 0., 0., 0.15, -18.435, 0.89)
time_func(macro_only_image, mypix,mysrc,mylens,PSF_scale)
time_func(macro_only_image_no_conv, mypix,mysrc,mylens)


print unpert_img.shape

#Plot
fig = plt.figure(figsize=(10,10)) 
im = plt.imshow(unpert_img,cmap="gray",origin='lower',extent=[-2,2,-2,2])
plt.colorbar(im)
plt.show()


