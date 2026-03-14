
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange

from numpy import pi

ph = t3.phantom3D()
e1 =t3.ellipsoid(0.,0.,0., 1. , 5.,5.,5., 0.,0.)
e2 =t3.ellipsoid(0.,1.,0., 0.5 , 1.,2.,0.5, 0.,0.)

ph.add_component(e1)
ph.add_component(e2)

#imaging volume params.
nx0 = 32
ny0 = 32
nz0 = 32
xlen0 = 20.
ylen0 = 20.
zlen0 = 20.
x00 = -10.
y00 = -10.
z00 = -10.

#cone-beam data params.

ns0 = 24
nu0 = 128
nv0 = 128
slen0 = 2*pi
ulen0 = 10.
vlen0 = 10.
s00 = 0.
u00 = -5.
v00 = -5.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

image =phantom_image.duplicate()

# generate discrete phantom array
ph.collapse_to(phantom_image)

data = t3.sinogram3D(\
                config_name='circular_conebeam',\
                parms={"source_to_detector":30.,\
                       "radius":20.},\
                shape = (ns0,nu0,nv0),\
                slen = slen0, ulen = ulen0, vlen = vlen0,\
                s0 = s00, u0 = u00, v0 = v00,s_include_endpoints=0)
raw_input('hi')
data_test = data.duplicate()

phantom_image.project_to(data,alg = 'proj_nearest')

#turn on only one ray in the data
rf = 1.0
data.indicator *= 0
data.indicator[0,64,64] += 1
print " Is the following value non_zero ? ", data.mat[0,64,64]

image.art_step(data,reg_fac = rf,alg = 'art_nearest')

image.project_to(data_test,alg = 'proj_nearest')



print rf,data_test.mat[0,64,64], data_test.mat[0,64,64]/data.mat[0,64,64]
