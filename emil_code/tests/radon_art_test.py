
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

#radon data params.

nsa0 = 24
nsb0 = 24
nr0 = 32
salen0 = 24.
sblen0 = 2.*pi
rlen0 = 20.
sa00 = -12.
sb00 = 0.
r00 = 2.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

image =phantom_image.duplicate()

# generate discrete phantom array
ph.collapse_to(phantom_image)

radon_data = t3.sphericalRadon3D(\
                config_name='z_sphere',\
                parms={"radius":12.},\
                shape = (nsa0,nsb0,nr0),\
                salen = salen0, sblen = sblen0, rlen = rlen0,\
                sa0 = sa00, sb0 = sb00, r0 = r00)
radon_test = radon_data.duplicate()

phantom_image.radon_project_to(radon_data,fov =  6.)

#turn on only one shell in the data
radon_data.indicator *= 0
radon_data.indicator[0,0,15] += 1
print " Is the following value non_zero ? ", radon_data.mat[0,0,15]

image.radon_art_step(radon_data,reg_fac = 1.0, fov = 6.)

image.radon_project_to(radon_test,fov = 6.)




image.mat *= 0.


# random ART

ndata = 2000

figure(1)
subplot(121)
p1 = imshow(phantom_image.mat[16])
subplot(122)
p2 = imshow(phantom_image.mat[16])
raw_input('wait')

rf= 1.0
radon_data.indicator.fill(1)
for i in range(2000):

   radon_data.indicator *= 0
   for j in range(ndata):
      i1 = randrange(nsa0)
      i2 = randrange(nsb0)
      i3 = randrange(nr0)

      radon_data.indicator[i1,i2,i3] += 1

   image.radon_art_step(radon_data,reg_fac = 1.0, fov = 6.)
   image.mat[image.mat<0.] = 0.

   imshow(image.mat[16])
   print i,rf,image.dist_to(phantom_image)
   rf *= 0.995


      
   









