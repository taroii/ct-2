
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange
import time

from numpy import pi,fromfile,float32
from pylab import imshow

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
nu0 = 16
nv0 = 16
ns0 = 12
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


#while 1==1:
#   ii =randrange(nu0-1)
#   jj =randrange(nv0)
#   kk =randrange(ns0)
#   data.indicator.fill(0)
#   data.indicator[kk,ii,jj]=1
#   image_gpu.mat.fill(0.)
#   image_gpu.get_deriv_weights(data,1,1,4,4)
#
#   fv = data.frame_vectors[kk]
#   sp = fv[:3]
#   dp = fv[6:9]
#   eu = fv[9:12]
#   ev = fv[12:15]
#   ew = fv[15:]
#   db = dp + eu*(data.u0 + (ii+0.5)*data.get_du())+ ev*(data.v0 + (jj+0.5)*data.get_dv())
#   print db-sp

#   work_sino.mat.fill(0.)
#   image_cpu.mat.fill(0.)
#   work_sino.mat[kk,ii+1,jj] = 1.
#   work_sino.mat[kk,ii,jj] = -1.
#   work_sino.backproject_to(image_cpu)
#
#   print (image_gpu.mat - image_cpu.mat).min(),(image_gpu.mat - image_cpu.mat).max()
#   work_sino.mat.fill(0.)
#   image_cpu.project_to(work_sino)
#   res = work_sino.mat[kk,ii+1,jj] -work_sino.mat[kk,ii,jj] 
#   print kk,ii,jj,data.weights[kk,ii,jj],res
#   raw_input()
#   data.weights.fill(0.)


t3.set_GPU_dev(1)
print "using GPU: #",t3.get_GPU_dev()

time0 = time.time()
phantom_image.clin_project_to(data,1,1,16,16)
image_gpu.get_deriv_weights(data,1,1,2,2)

data.mat[:,:nu0-1,:] = data.mat[:,1:,:] - data.mat[:,:nu0-1,:]
#phantom_image.cnn_project_to(data,4,16,4,32)
time1 = time.time()
print "GPU projection done in", time1 - time0," seconds."
time0 = time.time()
#image_cpu.art_step(data,reg_fac = 1.0,alg='art_linear')
# CPU done in  136 seconds
time1 = time.time()
print "CPU  lin art step done in", time1 - time0," seconds."

time0 = time.time()
image_gpu.cderiv_art_step(data,1,1,4,4)
#image_gpu.cnn_art_step(data,1,64,128,1) # This way is a factor of 10 slower
time1 = time.time()
print "GPU  nn art step done in", time1 - time0," seconds."

ntrial=10

for i in range(ntrial):
   nsr = randrange(ns0)
   nur = randrange(nu0-1)
   nvr = randrange(nv0)
   print 'testing:',nsr,' ',nur,' ',nvr
   image_gpu.mat.fill(0.)
   data.indicator.fill(0)
   data.indicator[nsr,nur,nvr] = 1
   work_sino.indicator.fill(0)
   work_sino.indicator[nsr,nur,nvr] = 1
   work_sino.indicator[nsr,nur+1,nvr] = 1
   print data.mat[nsr,nur,nvr]
   image_gpu.cderiv_art_step(data,1,1,4,4)
   image_gpu.clin_project_to(work_sino,1,1,16,16)
   print work_sino.mat[nsr,nur+1,nvr]-work_sino.mat[nsr,nur,nvr]



