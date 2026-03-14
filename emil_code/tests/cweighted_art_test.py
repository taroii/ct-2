
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange
import time
from scipy import ndimage
convolve= ndimage.convolve


from numpy import pi,fromfile,float32,sqrt,sum,ones,array
from pylab import imshow,hold,draw,subplot,cm
#hold(False)

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
rayweights=ones([1,4,4],float32)/16.
#rayweights[0,:,:] *= -1.
#rayweights=array([[[-1.],[1.]]],float32)
nss = 1
nsu = 4
nsv = 4
nsds = (ns0+nss-1)/nss
nuds = (nu0+nsu-1)/nsu
nvds = (nv0+nsv-1)/nsv


t3.set_GPU_dev(1)
print "using GPU: #",t3.get_GPU_dev()

phantom_image.clin_project_to(data,1,8,8,16)
phantom_image.cweighted_project_to(data, rayweights,nss,nsu,nsv, 1,8,8,16)

ss= convolve(data.mat,rayweights[::-1,::-1,::-1],mode= 'constant')[::nss,::nsu,::nsv]
#ph.project_to(data)
print abs(ss-data.smallsino).max()

time0 = time.time()
#data.indicator.fill(0.)
#data.indicator[0,64,64]=1
image_gpu.put_weights(data, rayweights,nss,nsu,nsv, 1,1,8,8)
raw_input('weight...')

#for i in range(100):
#   image_gpu.mat.fill(0.)
#   data.indicator.fill(0)
#   nsr = randrange(nsds)
#   nur = randrange(nuds)
#   nvr = randrange(nvds)
#   data.indicator[nsr*nss,nur*nsu,nvr*nsv]=1

#   print nsr,nur,nvr
#   image_gpu.put_weights(data, rayweights,nss,nsu,nsv, 1,1,8,8)
#   print "weight dim. ", data.weights.shape
#   print 'test1'
#   print sum(image_gpu.mat.flatten() * phantom_image.mat.flatten(),dtype = 'float64')

#   print data.smallsino[nsr,nur,nvr]

#   print 'test2'
#   print sum(image_gpu.mat.flatten() **2,dtype = 'float64')
#   print data.weights[nsr,nur,nvr]

#   print 'test3'

#   image_cpu.mat.fill(0.)
#   image_cpu.mat = image_gpu.mat*(data.smallsino[nsr,nur,nvr] )/data.weights[nsr,nur,nvr]
#   image_cpu.cweighted_project_to(work_sino, rayweights,nss,nsu,nsv, 1,8,8,16)

#   print work_sino.smallsino[nsr,nur,nvr]
  
#   image_gpu.mat.fill(0.)
#   image_gpu.c_ss_weighted_art_step(data, rayweights,nss,nsu,nsv, 1,1,8,8,reg_fac=1.0)
#   print 'test 4'
#   print (image_gpu.mat - image_cpu.mat).max()
#   print (image_gpu.mat - image_cpu.mat).min()


#   raw_input('next')

time1 = time.time()

print "GPU get weights done in", time1 - time0," seconds."
print "weight shape: ",data.weights.shape

#data.weights.fill(0.1)
subplot(121)
imshow(phantom_image.mat[:,:,64],cmap=cm.gray)
draw()
subplot(122)

for i in range(1):
   time0 = time.time()
   image_gpu.c_ss_weighted_art_step(data, rayweights,nss,nsu,nsv, 1,1,8,8,reg_fac=1.0)
   time1 = time.time()

   imshow(image_gpu.mat[:,:,64],cmap=cm.gray)
   draw()
   image_gpu.clin_project_to(work_sino,1,1,16,16)
   print i,  work_sino.dist_to(data)

print "GPU  nn art step done in", time1 - time0," seconds."
ntrial =20;
work_sino.indicator.fill(1)
for i in range(ntrial):
   nsr = randrange(nsds)
   nur = randrange(nuds)
   nvr = randrange(nvds)
   print i,':   testing:',nsr,' ',nur,' ',nvr
   image_gpu.mat.fill(0.)
   data.indicator.fill(0)
#   convolveddata= convolve(data.mat,rayweights[::-1,::-1,::-1],mode='constant')
#   print convolveddata[nss*nsr,nsu*nur,nsv*nvr]
   print data.smallsino[nsr,nur,nvr]

   data.indicator[nss*nsr,nsu*nur,nsv*nvr] = 1
#   image_gpu.cderiv_art_step(data,1,1,4,4)

   if data.smallsino[nsr,nur,nvr] >0.0000001:
      image_gpu.c_ss_weighted_art_step(data, rayweights,nss,nsu,nsv, 1,1,8,8,reg_fac=1.0)
#   image_gpu.clin_project_to(work_sino,1,1,16,16)
      image_gpu.cweighted_project_to(work_sino, rayweights,nss,nsu,nsv, 1,1,8,8)
#   convolvedsino = convolve(work_sino.mat,rayweights[::-1,::-1,::-1],mode='constant')
#   print convolvedsino[nss*nsr,nsu*nur,nsv*nvr]
      print work_sino.smallsino[nsr,nur,nvr]
   else:
      print 'skipping ...'




