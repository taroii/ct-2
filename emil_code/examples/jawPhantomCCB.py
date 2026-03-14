# sinogram3D example with jaw phantom

import sys,os
sys.path.append('../.')
sys.path.append('../phantoms')

import tomo3D as t3
from jaw_phantom import *
import time
import numpy

mod =numpy.mod


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

# noise level to add to data
noise_level=0.001

#image size and dimensions
nx0=366
ny0=546
nz0=216
xlen0=12.2
ylen0=18.2
zlen0= 7.2
x00=-6.1
y00=-6.1
z00= 2.4

#projection data size and dimensions
ns0=64
slen0=pi
s00= pi
nu0=610
ulen0=61.0
u00=-30.5
nv0=93
vlen0=9.3
v00=-4.65

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
   b.embed_in(image)
   filename = baseDir+"phantom.dat"
   image.write_file(filename)
   print 'wrote file'


#set up sinogram3D data function
# see configs.py for possible configuration
sino = t3.sinogram3D(config_name='z_circular_conebeam',\
                     parms={"source_to_detector": 60.,\
		            "radius"            : 30.,\
                            "zoffset"           : 6.},
	             shape=(ns0,nu0,nv0),\
		     slen= slen0, ulen=ulen0,vlen=vlen0,\
		     s0 = s00, u0 = u00, v0 = v00,s_include_endpoints=1)

work_sino = sino.duplicate()





if resume:
   sino.read(baseDir+'sino.dat')
else:

   print "projecting ..."
   basetime=time.time()
# use analytic projection to generate data
   b.project_to(sino)
   sino.copy_to(work_sino)
   sino.add_noise(noise_level)
   print "projection took ",time.time() -basetime
   ddist = sino.dist_to(work_sino)
   print "ddist = ", ddist
   derrfile=open(baseDir+'derr.txt','w')
   derrfile.write(str(ddist))
   derrfile.close()

   filename = baseDir+"sino.dat"
   sino.write_file(filename)
   print 'wrote sino file'

   filename = baseDir+"sinoIdeal.dat"
   work_sino.write_file(filename)
   print 'wrote sino file'

lam1 = 0.0
truetv=image.calc_total_var(lam1)
print "the true tv is ",truetv




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
      print "There's a problem reading iternum.txt file."
      sys.exit()

else:   
   iter =0
   gf=0.2
   reg_fac=1.00
   meps = -1.0

nart = 1
ngradmax=20
red_fact=0.98
red_reg=0.995
conv = 0.99

#set the data tolerance to whatever the data error is at
# iteration iter_eps
iter_eps = 500
bang_num=100
miter=mod(iter,bang_num)
print "miter is ",miter

done = 0


basetime=time.time()
while not done:
   iter+=1
   miter+=1

# POCS start
   arttime = time.time()
   current_image.copy_to(pre_image)
   for i in range(0,nart):
      current_image.art_step(sino,reg_fac)

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


   traf_file=open("trafficC.txt")
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
   current_image.project_to(work_sino)
   ddist=sino.dist_to(work_sino)
   print "parms. done"

   if iter == iter_eps:
      meps = ddist


# note: we can't do the following if we don't know the actual image
# and it's not clear what this means when the data are analytic projections
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

#GPU current_image.tv_grad_descent_gpu(gf0=gf,nsteps0=ngradmax,levels_per_row=8,ilev0=9)
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


