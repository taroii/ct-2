
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange

from numpy import pi

ph = t3.phantom3D()
e1 =t3.ellipsoid(0.,0.,0., 1. , 5.,5.,5., 0.,0.)
e2 =t3.ellipsoid(1.,1.,1., 0.5 , 3.,7.,5., 0.,0.)

ph.add_component(e1)
ph.add_component(e2)

#imaging volume params.
nx0 = 30
ny0 = 40
nz0 = 50
xlen0 = 30.
ylen0 = 40.
zlen0 = 50.
x00 = -15.
y00 = -20.
z00 = -25.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

ph.collapse_to(phantom_image)
image =phantom_image.duplicate()

print 'OK, play with phantom_image and image'
