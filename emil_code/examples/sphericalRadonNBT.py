# This program uses a new back-tracking scheme that hopefully converges better.
# the pocs step is scaled to the gradient descent dist.

import sys,os
sys.path.append('../.')
import tomo3D as t3
from random import randrange
import time

from numpy import pi

#setup data directory name, so it is the same as the program name
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



baseDir='data/'+bdname+'/'
# done with data directory


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

#setup simple phantom
ph = t3.phantom3D()
e1 =t3.ellipsoid(0.,0.,0., 1. , 5.,5.,5., 0.,0.)
e2 =t3.ellipsoid(0.,1.,0., 0.5 , 1.,2.,0.5, 0.,0.)

ph.add_component(e1)
ph.add_component(e2)

#imaging volume params.
nx0 = 32
ny0 = 32
nz0 = 32
xlen0 = 20.
ylen0 = 20.
zlen0 = 20.
x00 = -10.
y00 = -10.
z00 = -10.

fov0 = 6.
#radon data params.

nsa0 = 12
nsb0 = 12
nr0 = 32
salen0 = 24.
sblen0 = 2.*pi
rlen0 = 20.
sa00 = -12.
sb00 = 0.
r00 = 2.

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

truetv=phantom_image.calc_total_var()

current_image =phantom_image.duplicate()
pre_image =phantom_image.duplicate()
gimage =phantom_image.duplicate()
test_image =phantom_image.duplicate()

# generate discrete phantom array
ph.collapse_to(phantom_image)

radon_data = t3.sphericalRadon3D(\
                config_name='z_sphere',\
                parms={"radius":12.},\
                shape = (nsa0,nsb0,nr0),\
                salen = salen0, sblen = sblen0, rlen = rlen0,\
                sa0 = sa00, sb0 = sb00, r0 = r00)
work_data = radon_data.duplicate()

phantom_image.radon_project_to(radon_data,fov =  fov0)


if resume:

   current_image.read(baseDir+'imageTemp.dat')
   ctv=current_image.calc_total_var()
   print "starting tv is: ",ctv

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
   gf=0.5
   reg_fac=1.00
   meps = -1.0

nart = 1
ngradmax=5
red_fact=0.8
red_reg=0.995
conv = 0.9
agrat = 1.1

alpha_reduct = 0.5
cc1 = 0.001




current_image.mat *= 0.



figure(1)
subplot(121)
p1 = imshow(phantom_image.mat[:,:,nz0/2],cmap = cm.gray)
subplot(122)
p2 = imshow(phantom_image.mat[:,:,nz0/2],cmap = cm.gray)
raw_input('wait: hit enter to start TV...')



bangset=[1,2,3,4,5,10,15,20,30,40,50,100]

done = 0

#upper_bound = 1.5

basetime=time.time()
while not done:
   iter+=1

   arttime = time.time()
   current_image.copy_to(pre_image)
   for i in range(0,nart):
      current_image.radon_art_step(radon_data,reg_fac = 1.0, fov = fov0)

   arttime = time.time()-arttime
   print iter,": elapsed time (sec): ",time.time()-basetime, "POCS time: ",arttime


   current_image.zero_edges()
   current_image.make_positive()
#   current_image.mat[current_image.mat>upper_bound]=upper_bound

   adist=current_image.dist_to(pre_image)

# scale back art step
   if (iter>1) and (gdist/adist<conv):
      adistold = adist
      adist = gdist/conv
      current_image.mat = pre_image.mat + (adist/adistold)*(current_image.mat-pre_image.mat)


   print "writing temp image."
   current_image.write_file(baseDir+'imageTemp.dat')
   p2.set_data(current_image.mat[:,:,nz0/2])
   draw()
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
   ctv=current_image.calc_total_var()
   print "current TV rat is ",ctv/truetv






   if iter == 1:
      gf*=adist

   print "comp. parms."
   posttv=current_image.calc_total_var()
   current_image.radon_project_to(work_data,fov =  fov0)
   ddist=radon_data.dist_to(work_data)
   print "parms. done"



   if not done:
      derrfile=open(baseDir+'derrEvol.txt','a')
      derrfile.write(str(iter)+"\t"+str(ddist)+"\t"+str(posttv)+"\t"+str(current_image.dist_to(phantom_image))+"\n")
      derrfile.close()




   if iter in bangset:

      filename = baseDir+"image"+str(iter)+".dat"
      current_image.write_file(filename)
      print 'wrote file'

  

   current_image.copy_to(pre_image)
   ngrad=0



   tvgtime = time.time()



#THIS SECTION IS CPU GRAD. DESCENT
   rem_dist = gf
   gd_done = 0
   ncount = 0
   while (not gd_done) and (not done):
      ncount += 1
      current_tv = current_image.calc_total_var()
      gimage.mat= 1.0 * current_image.mat
      gmag = gimage.tv_grad()

      test_image.mat = current_image.mat - gf*gimage.mat
      test_tv = test_image.calc_total_var()
      alpha = 1.0
      print "%d: current_tv= %f, test_tv = %f, and rem_dist = %f"%(ncount,current_tv,test_tv, rem_dist)
      while (test_tv > current_tv - alpha*gf*gmag*cc1):
         alpha*= alpha_reduct
         test_image.mat = current_image.mat - alpha*gf*gimage.mat
         test_tv = test_image.calc_total_var()
         print "test_tv = %f, and alpha = %f"%(test_tv, alpha)
      if alpha*gf >= 0.9999*rem_dist:
         gd_done = 1
         current_image.mat -= rem_dist*gimage.mat
         rem_dist = 0.
      else:
         rem_dist -= alpha*gf
         current_image.mat = 1.0*test_image.mat

      current_image.zero_edges()

      if ncount >= ngradmax:
         gd_done = 1
         gf -= rem_dist
         print "gf now: ",gf

#END CPU GRAD. DESCENT





   tvgtime = time.time()-tvgtime


   gdist=current_image.dist_to(pre_image)
   print "adist = %f, gdist = %f"%(adist,gdist)
   mag=current_image.calc_total_var()

   if gdist>=conv*adist and ddist > meps:
      gf*=red_fact

   if adist>=agrat*gdist:
      reg_fac*=red_reg

   print 'gf = ',gf, '; reg_fac = ',reg_fac, '; meps = ',meps
   print "elapsed time (sec): ",time.time()-basetime," grad desc. time: ",tvgtime








