import sys,os
sys.path.append('../.')
sys.path.append('../phantoms')
import tomo3D as t3
from phantoms import generate_3D_shepp_logan_HC
from numpy import random

#set up images with default size
nx=3
ny=4
nz=5
shape0 = (nx,ny,nz)
ranm = random.randn(nx,ny,nz)
im = t3.image3D(shape = shape0)
im.mat += ranm

wim1 = im.duplicate()
wim2 = im.duplicate()
delt = 0.000000001

print "Make sure to switch to double precision for these tests!"
p0 = float(raw_input("What p do you want for tpv test? "))

b = generate_3D_shepp_logan_HC()
b.embed_in(im)
tvgim=im.duplicate()
tpvgim = im.duplicate()
tpdvgim = im.duplicate()
tvgmag = tvgim.tv_grad()
tpvgmag = tpvgim.tpv_grad(p0)
tpdvgmag = tpdvgim.tpdv_grad(p0)

print "TV tests"
tv0 = im.calc_total_var()
print "Original tv: ",tv0
print "grad TV mag is: ", tvgmag
wim1.mat = im.mat- delt*tvgim.mat/2.
wim2.mat = im.mat+ delt*tvgim.mat/2.
tv1 = wim1.calc_total_var()
tv2 = wim2.calc_total_var()
print "tv1 = %20.15f, tv2 = %20.15f"%(tv1,tv2)
print "diff1 = %20.15f, diff2 = %20.15f"%(tv0-tv1,tv2-tv0)
print "Numerical grad is ",(tv2-tv1)/delt



print "\nTpV tests"
tpv0 = im.calc_total_var_p(p0)
print "Original tpv: ",tpv0
print "grad TpV mag is: ", tpvgmag
wim1.mat = im.mat- delt*tpvgim.mat/2.
wim2.mat = im.mat+ delt*tpvgim.mat/2.
tpv1 = wim1.calc_total_var_p(p0)
tpv2 = wim2.calc_total_var_p(p0)
print "tpv1 = %20.15f, tpv2 = %20.15f"%(tpv1,tpv2)
print "diff1 = %20.15f, diff2 = %20.15f"%(tpv0-tpv1,tpv2-tpv0)
print "Numerical grad is ",(tpv2-tpv1)/delt


print "\nTpdV tests"
print "dx = %20.15f, dy = %20.15f, dz = %20.15f"%(im.dx,im.dy,im.dz)
tpdv0 = im.calc_total_dvar_p(p0)
print "Original tpdv: ",tpdv0
print "grad TpdV mag is: ", tpdvgmag
wim1.mat = im.mat- delt*tpdvgim.mat/2.
wim2.mat = im.mat+ delt*tpdvgim.mat/2.
tpdv1 = wim1.calc_total_dvar_p(p0)
tpdv2 = wim2.calc_total_dvar_p(p0)
print "tpdv1 = %20.15f, tpdv2 = %20.15f"%(tpdv1,tpdv2)
print "diff1 = %20.15f, diff2 = %20.15f"%(tpdv0-tpdv1,tpdv2-tpdv0)
print "Numerical grad is ",(tpdv2-tpdv1)/delt




