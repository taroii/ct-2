
import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange
import time

from numpy import pi,fromfile,float32,ones

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

#cone-beam data params.

ns0=10
nu0 = 512
nv0 = 512
slen0 = 2*pi
ulen0 = 20.
vlen0 = 20.
s00 = 0.
u00 = -10.
v00 = -10.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

image_cpu =phantom_image.duplicate()
image_gpu =phantom_image.duplicate()


data = t3.sinogram3D(\
                config_name='circular_conebeam',\
                parms={"source_to_detector":30.,\
                       "radius":20.},\
                shape = (ns0,nu0,nv0),\
                slen = slen0, ulen = ulen0, vlen = vlen0,\
                s0 = s00, u0 = u00, v0 = v00,s_include_endpoints=0)
work_sino = data.duplicate()
data_fl = t3.sinogram3D(\
                config_name='circular_conebeam_flipped',\
                parms={"source_to_detector":30.,\
                       "radius":20.},\
                shape = (ns0,nv0,nu0),\
                slen = slen0, ulen = vlen0, vlen = ulen0,\
                s0 = s00, u0 = v00, v0 = u00,s_include_endpoints=0)
work_sino_fl = data_fl.duplicate()
raw_input('hi')

t3.set_GPU_dev(1)
print "using GPU: #",t3.get_GPU_dev()

time0 = time.time()
ph.project_to(data)
ph.project_to(data_fl)
time1 = time.time()
print "CPU anal. projection done in", time1 - time0," seconds."
time0 = time.time()
data.backproject_to(image_cpu,alg='adj_linear')
# CPU done in  72.36 seconds
time1 = time.time()
print "CPU  lin art step done in", time1 - time0," seconds."



nub = 1
nvb = 1
ubs = 16 
vbs = 16

image_gpu.mat.fill(0.)
print "starting GPU art..."
time0 = time.time()
#diagar = image_gpu.clin_bpcheck(data,nub,nvb,ubs,vbs)
#image_gpu.clin_art_step(data,64,1,4,2)
#image_gpu.clin_art_step(data,1,64,128,1)
#image_gpu.cnn_art_step(data,1,64,128,1) # This way is a factor of 10 slower
time1 = time.time()
print "GPU  lin bpcheck done in", time1 - time0," seconds."
raw_input()

image_gpu.mat.fill(0.)
print "starting GPU art..."
time0 = time.time()
image_gpu.clin_backproject(data,nub,nvb,ubs,vbs)
#image_gpu.clin_art_step(data,64,1,4,2)
#image_gpu.clin_art_step(data,1,64,128,1)
#image_gpu.cnn_art_step(data,1,64,128,1) # This way is a factor of 10 slower
time1 = time.time()
print "GPU  lin art step done in", time1 - time0," seconds."

temp = image_gpu.mat * 1.
image_gpu.mat.fill(0.)
print "starting GPU art on flipped data..."
time0 = time.time()
image_gpu.clin_backproject(data_fl,nvb,nub,vbs,ubs)
#image_gpu.clin_art_step(data,64,1,4,2)
#image_gpu.clin_art_step(data,1,64,128,1)
#image_gpu.cnn_art_step(data,1,64,128,1) # This way is a factor of 10 slower
time1 = time.time()
print "GPU  lin art step done in", time1 - time0," seconds."
print (temp -image_gpu.mat).min()
print (temp -image_gpu.mat).max()
