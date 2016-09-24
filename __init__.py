# plotting
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
mpl.rc('image', interpolation='nearest', origin='lower')
#mpl.rcParams['figure.figsize'] = (6., 5.)
mpl.rc('figure', figsize=(6, 5))

import seaborn as sns
sns.set(context='poster', style='ticks', color_codes=True)

# numpy
import numpy as np
from numpy import *
from numpy.random import *
from numpy.random import choice

# scipy
import scipy as sp
from scipy.interpolate import *
from scipy.special import erfinv, erf
from scipy.stats import poisson as pss

# multiprocessing
import multiprocessing as mp

# healpy
import healpy as hp
from healpy.rotator import angdist
from healpy import ang2pix

# pyfits
import pyfits as pf

# utilities
import os, time, sys, datetime, warnings, getpass, glob, fnmatch, cPickle, inspect
import functools

# tdpy
import tdpy.util
from tdpy.util import summ, show
import tdpy.mcmc

from tdpy.util import show

import pcat.main

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

np.set_printoptions(linewidth=180)

