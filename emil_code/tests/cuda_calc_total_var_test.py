
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange
import time

from numpy import pi,fromfile,float32,random
randn = random.randn

ph = t3.phantom3D()
e1 =t3.ellipsoid(0.,0.,0., 1. , 5.,5.,5., 0.,0.)
e2 =t3.ellipsoid(0.,1.,0., 0.5 , 1.,2.,0.5, 0.,0.)

ph.add_component(e1)
ph.add_component(e2)

#imaging volume params.
nx0 = 512
ny0 = 512
nz0 = 512
xlen0 = 20.
ylen0 = 20.
zlen0 = 20.
x00 = -10.
y00 = -10.
z00 = -10.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

gimage_cpu =phantom_image.duplicate()
gimage_gpu =phantom_image.duplicate()

# generate discrete phantom array
#ph.collapse_to(phantom_image)
phantom_image.mat = fromfile('phantom512.dat',float32)
phantom_image.mat.shape = (512,512,512)
phantom_image.mat += randn(512,512,512)
#phantom_image.mat = fromfile('phantom256.dat',float32)
#phantom_image.mat.shape = (256,256,256)



gimage_cpu.mat = 1.*phantom_image.mat
gimage_gpu.mat = 1.*phantom_image.mat

raw_input('hi')
time0 = time.time()

# CPU tpvg took 116 sec.
res1 = gimage_cpu.calc_total_var()
time1 = time.time()
print "CPU done"
print res1

t3.set_GPU_dev(0)
print "using GPU: #",t3.get_GPU_dev()

time2 = time.time()
res2 = gimage_gpu.calc_tv_cuda(16,32)
time3 = time.time()
print "GPU done"
print res2
#phantom_image.cnn_project_to(data_test,2,2,2,2)
print "CPU time: ",time1-time0," secs."
print "GPU time: ",time3-time2," secs."

