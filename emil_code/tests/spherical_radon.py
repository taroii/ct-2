
import sys,os
sys.path.append('../.')
import tomo3D as t3

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

#radon data params.

nsa0 = 16
nsb0 = 16
nr0 = 32
salen0 = 24.
sblen0 = 2.*pi
rlen0 = 20.
sa00 = -12.
sb00 = 0.
r00 = 2.

image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

ph.collapse_to(image)

radon_data = t3.sphericalRadon3D(\
                config_name='z_sphere',\
                parms={"radius":12.},\
                shape = (nsa0,nsb0,nr0),\
                salen = salen0, sblen = sblen0, rlen = rlen0,\
                sa0 = sa00, sb0 = sb00, r0 = r00)
