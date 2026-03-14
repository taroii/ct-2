
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange
import time

from numpy import pi,fromfile,float32,sqrt,sum
from pylab import imshow,hold,draw,subplot,cm
hold(False)

ph = t3.phantom3D()
e1 =t3.ellipsoid(0.,0.,0., 1. , 5.,5.,5., 0.,0.)
e2 =t3.ellipsoid(0.,1.,0., 0.5 , 1.,2.,0.5, 0.,0.)

ph.add_component(e1)
ph.add_component(e2)

#imaging volume params.
nx0 = 512
ny0 = 512
nz0 = 512
nx0 = 128
ny0 = 128
nz0 = 128
xlen0 = 20.
ylen0 = 20.
zlen0 = 20.
x00 = -10.
y00 = -10.
z00 = -10.

#cone-beam data params.

ns0 = 128
nu0 = 128
nv0 = 128
slen0 = 2*pi
ulen0 = 20.
vlen0 = 20.
s00 = 0.
u00 = -10.
v00 = -10.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

image_cpu =phantom_image.duplicate()
image_gpu =phantom_image.duplicate()

# generate discrete phantom array
print "embedding ..."

ph.collapse_to(phantom_image)
print "done"
#phantom_image.mat = fromfile('phantom512.dat',float32)
#phantom_image.mat.shape = (512,512,512)
#phantom_image.mat = fromfile('phantom256.dat',float32)
#phantom_image.mat.shape = (256,256,256)

data = t3.sinogram3D(\
                config_name='circular_conebeam',\
                parms={"source_to_detector":30.,\
                       "radius":20.},\
                shape = (ns0,nu0,nv0),\
                slen = slen0, ulen = ulen0, vlen = vlen0,\
                s0 = s00, u0 = u00, v0 = v00,s_include_endpoints=0)


work_sino = data.duplicate()

raw_input('glong')



t3.set_GPU_dev(1)
print "using GPU: #",t3.get_GPU_dev()

time0 = time.time()
#phantom_image.clin_project_to(data,1,8,8,16)
ph.project_to(data)
image_gpu.get_deriv_weights(data,2,1,8,16)

data0 = data.duplicate()
data.mat[:,:nu0-1,:] = data.mat[:,1:,:] - data.mat[:,:nu0-1,:]
#phantom_image.cnn_project_to(data,4,16,4,32)
time1 = time.time()
print "GPU projection done in", time1 - time0," seconds."
time0 = time.time()
#image_cpu.art_step(data,reg_fac = 1.0,alg='art_linear')
# CPU done in  136 seconds
time1 = time.time()
print "CPU  lin art step done in", time1 - time0," seconds."
subplot(121)
imshow(phantom_image.mat[:,:,64],cmap=cm.gray)
draw()
subplot(122)
time0 = time.time()
for giga in range(100):
   print giga
   print sqrt(sum((work_sino.mat.flatten() - data.mat.flatten())**2, dtype = 'float64'))
   image_gpu.clin_art_step(data0,2,1,8,16,reg_fac=0.1)
#   image_gpu.cderiv_art_step(data,2,1,8,16,reg_fac=0.1)
   image_gpu.make_positive()
   print image_gpu.mat.min(), image_gpu.mat.max()
   imshow(image_gpu.mat[:,:,64],cmap=cm.gray)

   image_gpu.clin_project_to(work_sino,2,1,8,16)
   work_sino.mat[:,:nu0-1,:] = work_sino.mat[:,1:,:] - work_sino.mat[:,:nu0-1,:]
   print sqrt(sum((work_sino.mat.flatten() - data.mat.flatten())**2, dtype = 'float64'))
   draw()

#image_gpu.cnn_art_step(data,1,64,128,1) # This way is a factor of 10 slower
time1 = time.time()
print "GPU  nn art step done in", time1 - time0," seconds."



