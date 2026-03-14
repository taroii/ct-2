# linogram data example w/ random views and full detector sampling


import sys,os
sys.path.append('../.')
sys.path.append('../phantoms')

import tomo3D as t3
from phantoms import generate_3D_shepp_logan_HC
from jaw_phantom import *
import time
import numpy
import random

mod =numpy.mod
sin = numpy.sin
cos = numpy.cos


# form base directory name from name of the program
ss = sys.argv
bdtemp=''
for s in ss:
   stemp = s.split('.')
   if len(stemp) == 2:
      if stemp[1] == 'py':
         bdname = stemp[0]

if bdname =='':
   print 'basedir name missed'
   sys.exit()


baseDir='./data/'+bdname+'/'


resume=0
try:
   os.mkdir(baseDir)
except:
   ans=raw_input('Are we resuming? ')
   if ans=='n':
      print 'Data directory already exists. Choose another name.'
      sys.exit()
   else:
      resume=1



pi = t3.pi

#image size and dimesions
nx0=100
ny0=100
nz0=100
xlen0=2.
ylen0=2.
zlen0= 2.
x00=-1.
y00=-1.
z00= -1.

# projection data size and dimensions
ns0=10
nu0 = 100
ulen0 = 40.0
u00 = -20.0
delu =  ulen0/nu0
nv0 = 100
vlen0 = 40.0
v00 = -20.0
delv =  vlen0/nv0
# calc. total number of rays
nr0 = ns0*nu0*nv0

# b is the high contrast Shepp-Logan phantom
#b = generate_3D_shepp_logan_HC()

b = build_jaw_phantom()

#set up image used as numerical phantom
image = t3.image3D(\
           shape=(nx0,ny0,nz0),\
           xlen= xlen0, ylen=ylen0,zlen=zlen0,\
	   x0 = x00, y0= y00, z0 = z00)

#create necessary image variable for the TV alg.
current_image = image.duplicate()
pre_image = image.duplicate()
gimage = image.duplicate()


if resume:
   image.read(baseDir+'phantom.dat')
else:
   print "embedding ..."
# The following line makes image into the digital shepp-logan phantom
#   b.embed_in(image)

# instead we will just put an array of dots in the image
   image.mat[20:80:4,20:80:4,20:80:4] = 1.0
   filename = baseDir+"phantom.dat"
   image.write_file(filename)
   print 'wrote file'


#setup linogram3D data function
# ball_flatpanel allows the source to be any point on a sphere
# and a flat-panel detector is put opposite to the source (see linoconfigs.py)
lino = t3.linogram3D(config_name='ball_flatpanel',\
                     parms={"source_radius":      500.0,\
		            "source_detector":    700.0},\
		     nrays = nr0)


# set up ray vectors
# makes list of rays for corresponding to random source location
# on the ball and full sampling on a flat-panel detector
raycount = 0
for iv in range(ns0):
   vtheta = pi*random.random()
   thetatest = random.random()
   while thetatest> sin(vtheta):
      vtheta = pi*random.random()
      thetatest = random.random()
   vphi = 2.*pi*random.random()
   print vtheta, vphi
   for ii in range(nu0):
      u = u00 + (ii+0.5)*delu
      for jj in range(nv0):
         v = v00 + (jj+0.5)*delv
         batch = lino.frame([vtheta,vphi,u,v],lino.parms)
# the actual assignment of the ray into the lino ray list is in the next line
         lino.ray_vectors[raycount,:]=batch[:]
         raycount +=1

   


work_lino = lino.duplicate()


if resume:
   lino.read(baseDir+'lino.dat')
else:
   print "projecting ..."
   basetime=time.time()
   b.lino_project_to(lino)
   print "projection took ",time.time() -basetime

   filename = baseDir+"lino.dat"
   lino.write_file(filename)
   print 'wrote lino file'

lam1 = 0.0
truetv=image.calc_total_var(lam1)
print "the true tv is ",truetv




# set alg. parameters that vary 
# or read them in if we are resuming
if resume:

   current_image.read(baseDir+'imageTemp.dat')
   ctv=current_image.calc_total_var(lam1)
   print "starting tv rat is: ",ctv

   try:
      iterfile=open(baseDir+'iternum.txt','r')
      iter=int(iterfile.read())
      iterfile.close()
      parmfile=open(baseDir+'parms.txt','r')
      gf=float(parmfile.read())
      reg_fac=float(parmfile.read())
      meps=float(parmfile.read())
      parmfile.close()
   except:
      print "There's a problem reading iternum.txt or parms.txt ."
      sys.exit()

else:   
   iter =0
   gf=0.2
   reg_fac=1.00
# meps is the data tolerance
   meps = 0.01

#set up fixed parameters
nart = 1
ngradmax=20
red_fact=0.5
red_reg=0.999
conv = 0.9

#how often do you want an image spit out
bang_num=500
miter = mod(iter,bang_num)
# this does every 500th iteration

done = 0

#main loop
basetime=time.time()
while not done:
   iter+=1
   miter+=1

# POCS start
   arttime = time.time()
   current_image.copy_to(pre_image)
   for i in range(0,nart):
      current_image.lino_art_step(lino,reg_fac)

   arttime = time.time()-arttime
   print iter,": elapsed time (sec): ",time.time()-basetime, "POCS time: ",arttime


   current_image.zero_edges()
   current_image.make_positive()
#POCS end

# adist represents change in image and
# if it's the first iteration the gradient frac. distance is converted
# to an actual image distance
   adist=current_image.dist_to(pre_image)
   if iter == 1:
      gf*=adist


# write intermediate results
   print "writing temp image."
   current_image.write_file(baseDir+'imageTemp.dat')
   print "done writing"

   iterfile=open(baseDir+'iternum.txt','w')
   iterfile.write(str(iter))
   iterfile.close()
   parmfile=open(baseDir+'parms.txt','w')
   parmfile.write(str(gf)+"\n")
   parmfile.write(str(reg_fac)+"\n")
   parmfile.write(str(meps)+"\n")
   parmfile.close()

# should we stop?
   traf_file=open("trafficA.txt")
   traf_sign=traf_file.read()
   traf_file.close()
   traf_sign=traf_sign.strip()
   if traf_sign=='stop':
      done=1
      print "stopping ..."
      ctv=current_image.calc_total_var(lam1)
      print "current TV rat is ",ctv/truetv



   print "comp. parms."
   posttv=current_image.calc_total_var(lam1)
   current_image.lino_project_to(work_lino)
   ddist=lino.dist_to(work_lino)
   print "parms. done"


# note: we can't do the following if we don't know the actual image
   ierr=current_image.dist_to(image)

   if not done:
      derrfile=open(baseDir+'derrEvol.txt','a')
      derrfile.write(str(iter)+"\t"+str(ddist)+"\t"+str(posttv)+"\t"+str(ierr)+"\n")
      derrfile.close()




   if miter>=bang_num:
      miter=0

      filename = baseDir+"image"+str(iter)+".dat"
      current_image.write_file(filename)
      print 'wrote file'

  

   current_image.copy_to(pre_image)


# START of TV gradient descent
#if you decide to uncomment the gpu lines:
#   (1) make sure the gpu prog works!
#   (2) get rid of the gimage variable to save memory
#   (3) comment out the CPU lines

#GPU   gpu_busy = 1
#GPU   ig = 0
#GPU   while gpu_busy:
#GPU      if ig == 1:
#GPU         print "waiting for gpu ... ";
#GPU      try:
#GPU         gpuf = open("gpu.txt",'r')
#GPU         gpuf.close()
#GPU      except:
#GPU         print "entering gradient descent"
#GPU         gpu_busy = 0
#GPU         gpuf = open("gpu.txt",'w')
#GPU         gpuf.write("gpu busy")
#GPU         gpuf.close()
#GPU      ig +=1

   tvgtime = time.time()

#THIS SECTION IS CPU GRAD. DESCENT
   ngrad=0
   while (ngrad<ngradmax) and (not done):
      ngrad+=1
      gimage.mat=current_image.mat
      gimage.tv_grad()
      gimage.scale(gf)
      current_image.mat-=gimage.mat
      current_image.zero_edges()
      current_image.make_positive()
##END CPU GRAD. DESCENT

#GPU   current_image.tv_grad_descent_gpu(gf0=gf,nsteps0=ngradmax,levels_per_row=10,ilev0=1)
   tvgtime = time.time()-tvgtime
#GPU   os.system("rm gpu.txt")

   gdist=current_image.dist_to(pre_image)
   print "adist = %f, gdist = %f"%(adist,gdist)
   mag=current_image.calc_total_var(lam1)

   if gdist>=conv*adist and ddist > meps:
      gf*=red_fact

   reg_fac*=red_reg

   print 'gf = ',gf, '; reg_fac = ',reg_fac, '; meps = ',meps, ' ierr = ',ierr
   print "elapsed time (sec): ",time.time()-basetime," grad desc. time: ",tvgtime


