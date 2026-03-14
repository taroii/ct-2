#! /usr/bin/python


import sys,os
sys.path.append('/home/emil/python/PyTomo3Diterative')
sys.path.append('/home/emil/python/PyTomo3Diterative/phantoms')

import tomo3D as t3
import time
import numpy

from head_phantom import *

from pylab import figure, subplot, imshow, cm, draw,  ion

ion()

head = build_head_phantom()


mod =numpy.mod



pi = t3.pi


nx0=256
ny0=256
nz0=16
xlen0= 25.0
ylen0= 25.0
zlen0=  25.0/16.0
x00=-xlen0/2.
y00=-ylen0/2.
z00=-zlen0/2.




#cone-beam data params.

# units are mm

ns0 = 50
#ns0 = 307
nu0 = 512
nv0 = 16
slen0 = 2.*pi
ulen0 = 60.0
vlen0 = ulen0 *nv0/nu0
s00 = 0.
u00 = -ulen0/2.
v00 = - vlen0/2.



sino = t3.sinogram3D(\
                config_name='circular_conebeam',\
                parms={"source_to_detector":80.,\
                       "radius":40.0},\
                shape = (ns0,nu0,nv0),\
                slen = slen0, ulen = ulen0, vlen = vlen0,\
                s0 = s00, u0 = u00, v0 = v00,s_include_endpoints=0)

work_sino = sino.duplicate()


image = t3.image3D(\
           shape=(nx0,ny0,nz0),\
           xlen= xlen0, ylen=ylen0,zlen=zlen0,\
	   x0 = x00, y0= y00, z0 = z00)

work_image = image.duplicate()
current_image = image.duplicate()
pre_image =image.duplicate()

print 'embedding head phantom ...'
head.embed_in(image)
print 'done'

print 'generating projection data ...'
#following line is for continuous projection
#head.project_to(sino)

#following line is for discrete projection
image.project_to(sino,alg='proj_linear')
print 'done'

reg_fac=1.00

nart = 1

red_reg=0.995


figure(1)
subplot(121)
imshow(image.mat[:,:,nz0/2],cmap = cm.gray)
subplot(122)
pl1 = imshow(image.mat[:,:,nz0/2],cmap = cm.gray)
draw()

# uncomment following line for EM
# current_image.mat.fill(1.)

done = 0
iter = 0
basetime=time.time()
while not done:
   iter+=1

   arttime = time.time()
   current_image.copy_to(pre_image)
   for i in range(0,nart):
      current_image.art_step(sino,reg_fac,alg='art_linear')
#      current_image.em_step(sino,alg='em_linear')

   arttime = time.time()-arttime
   print iter,": elapsed time (sec): ",time.time()-basetime, "POCS time: ",arttime

   pl1.set_data(current_image.mat[:,:,nz0/2])
   draw()

   current_image.project_to(work_sino)
   work_sino.mat -= sino.mat
   work_sino.mat *= sino.indicator
   print "err: ",sqrt(sum((work_sino.mat.flatten())**2))

