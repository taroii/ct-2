#######################################################################
#S.Lee: Fan beam filtered backprojection
#######################################################################
import sys
import numpy as n
from pylab import *

#import tomo2D as tomo
import tomo3D as t3
import Utils
from FDK import *

#choose whether sinogram with ear_right_hole is used or not.
#earhole_flag = 1 
earhole_flag = 0

#set sinogram parameters
ns0 = 360
nu0 = 256
nv0 = 256
#ns0 = 720
#nu0 = 512
#nv0 = 512
#ns0 = 1440
#nu0 = 1024
#nv0 = 1024
slen0 = 2.*pi
ulen0 = 52.
vlen0 = 52.
s00 = 0.
u00 = - ulen0/2.
v00 = - vlen0/2.
sino = t3.sinogram3D(\
             config_name='circular_conebeam',\
             parms={"source_to_detector":100.,\
                    "radius":50.},\
             shape = (ns0,nu0,nv0),\
             slen = slen0, ulen = ulen0, vlen = vlen0,\
             s0 = s00, u0 = u00, v0 = v00, s_include_endpoints=0)

#load sinogram data
print "load sinogram ..."
if earhole_flag == 1: #earholes are included
	filename = './sinograms/sino3D.%dp.%dp%dv'%(nu0,nv0,ns0)
else: #earholes are not included
	filename = './sinograms.without.earhole/sino3D.%dp.%dp.%dv'%(nu0,nv0,ns0)

sinodata = Utils.MGHread(filename)
sino.mat = sinodata

#set up coordinates of s, u in the sinogram
#s,u in sinogram
vs = s00 + n.arange(0,sino.ns) * sino.ds
#offset: midpoint in a pixel
uoffset = 0.5
voffset = 0.5
#uoffset = 0.0
#voffset = 0.0
vu = u00 + (n.arange(0,sino.nu) + uoffset) * sino.du
vv = u00 + (n.arange(0,sino.nv) + voffset) * sino.dv
#nu0 X nv0
vv, vu = n.meshgrid(vv, vu)

#weighting
w_sinomat = fdk_sino_weight3D(sino.parms["radius"], \
	sino.parms["source_to_detector"], \
	sino.ns, vu, vv, sino.mat)

#ramp filtering along u axis
print "filtering sinogram ..."
sinofilt = sino.duplicate()
opt_window = '' #none
#opt_window = 'Hann' #Hanning window
if opt_window == 'Hann':
	Lval = 0.5 #1.0, 0.75, 0.5 typically
dim_filter = 1 #filter along u axis
if opt_window == 'Hann': #Hanning window
	sinofilt.mat = fdk_sino_filter3D(sino.du, dim_filter, w_sinomat, opt_window, L=Lval)
else: #no windowing
	sinofilt.mat = fdk_sino_filter3D(sino.du, dim_filter, w_sinomat, opt_window)
del w_sinomat

#set up coordinates of x,y in the original image
#Note that meshgrid(x,y) comes as Ny X Nx array
nx0 = 256
ny0 = 256
nz0 = 256
xlen0 = 26.0
ylen0 = 26.0
zlen0 = 26.0
x00 = -13.0
y00 = -13.0
z00 = -13.0
img = t3.image3D(\
      shape=(nx0,ny0,nz0),\
      xlen= xlen0, ylen=ylen0,zlen=zlen0,\
      x0 = x00, y0= y00, z0 = z00)
#offset: midpoint in a voxel
xoffset = 0.5
yoffset = 0.5
zoffset = 0.5
#xoffset = 0.0
#yoffset = 0.0
#zoffset = 0.0
vx = x00 + (n.arange(0,nx0) + xoffset) * img.dx
vy = y00 + (n.arange(0,ny0) + yoffset) * img.dy
vz = z00 + (n.arange(0,nz0) + zoffset) * img.dz
vx,vy = n.meshgrid(vx,vy)
vx = vx.flatten()
vy = vy.flatten()
'''
#memory problem
xyz = setcoordinates(nx0, ny0, nz0, \
	x00, y00, z00, \
	img.dx, img.dy, img.dz,
	xoffset, yoffset, zoffset, mask_flag, \
	vu.max(), vv.max(), \
	sino.parms["radius"], sino.parms["source_to_detector"])
vx = xyz[0,:]
vy = xyz[1,:]
vz = xyz[2,:]
'''

#reconstruction based on linear interpolation
print "backprojection based on linear interpolation ..."
recon_linear = fdk_backproj_linear3D(sino.ns, sino.nu, sino.nv, \
	sino.du, sino.dv, sino.u0, sino.v0, \
	sino.parms["radius"], sino.parms["source_to_detector"],\
        uoffset, voffset, \
	vs, vx, vy, vz, \
	sinofilt.mat)

#back to cube
reconcube_linear=n.zeros((nx0,ny0,nz0))
tmprecon = n.zeros((ny0, nx0))
for ii in range(nz0):
	#row comes first in the output of numpy.nonzero()
	tmprecon = reshape(recon_linear[:,ii], (ny0, nx0))
	reconcube_linear[:,:,ii] = tmprecon.T
	

#save reconstructed image
if earhole_flag == 1: #ear holes are included
	filename_linear = './recon.cubes/recon3D.linear%s.%dp.%dp.%dv'%(opt_window,nu0,nv0,ns0)
else:                #ear holes are not included
	filename_linear = './recon.cubes.without.earhole/recon3D.linear%s.%dp.%dp.%dv'%(opt_window,nu0,nv0,ns0)

Utils.MGHwrite(filename_linear, reconcube_linear)
