# inverts spherical Radon transform using basic TV algorithm (see PMB:08 paper).

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
#                center        att     half axes lengths    orientation angles
e1 =t3.ellipsoid(0.,0.,0.,      1. ,    5.,5.,5.,              0.,0.)
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
    # a is z location
    # b is azimuthal angle
    # r is shell radius
nsa0 = 12    # number of z values
nsb0 = 12    # number of phi values
nr0 = 32     # number of shells
salen0 = 24.  # range of z
sblen0 = 2.*pi # range of phi
rlen0 = 20.  # range of shell radii
sa00 = -12.  # starting z value
sb00 = 0.    # starting phi
r00 = 2.     # starting shell radius

phantom_image = t3.image3D(shape = (nx0,ny0,nz0),xlen=xlen0,ylen=ylen0,zlen=zlen0,x0=x00,y0=y00,z0=z00)

current_image =phantom_image.duplicate()
pre_image =phantom_image.duplicate()
gimage =phantom_image.duplicate()

# generate discrete phantom array
ph.collapse_to(phantom_image)
truetv=phantom_image.calc_total_var()

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
      parms=(parmfile.read()).split()
      gf=float(parms[0])
      reg_fac=float(parms[1])
      meps=float(parms[2])
      parmfile.close()
      print "gf = ",gf,"; reg_fac = ",reg_fac,"; meps = ",meps
   except:
      print "There's a problem reading iternum.txt file."
      sys.exit()

else:   
   iter =0
   gf=0.5      # alpha and dtvg in the PMB:08
   reg_fac=1.00  # beta in PMB:08
   meps = -1.0

iter_eps = 500

nart = 1
ngradmax=10
red_fact=0.8 #alpha_red from PMB:08
red_reg=0.999  #beta_red from PMB:08
red_reg =1.0

conv = 0.99   # r_max in PMB:08





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
      current_image.radon_art_step(radon_data,reg_fac = reg_fac, fov = fov0)

   arttime = time.time()-arttime
   print iter,": elapsed time (sec): ",time.time()-basetime, "POCS time: ",arttime


   current_image.zero_edges()
   current_image.make_positive()
#   current_image.mat[current_image.mat>upper_bound]=upper_bound


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


   traf_file=open("trafficD.txt")
   traf_sign=traf_file.read()
   traf_file.close()
   traf_sign=traf_sign.strip()
   if traf_sign=='stop':
      done=1
      print "stopping ..."
   ctv=current_image.calc_total_var()
   print "current TV rat is ",ctv/truetv



   adist=current_image.dist_to(pre_image)
   if iter == 1:
      gf*=adist

   print "comp. parms."
   posttv=current_image.calc_total_var()
   current_image.radon_project_to(work_data,fov =  fov0)
   ddist=radon_data.dist_to(work_data)
   if iter == iter_eps:
      meps = ddist
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
   ngrad=0
   while (ngrad<ngradmax) and (not done):
      ngrad+=1
      gimage.mat=current_image.mat*1.0
      gimage.tv_grad()
      gimage.scale(gf)
      current_image.mat-=gimage.mat
      current_image.zero_edges()
      current_image.make_positive()

#END CPU GRAD. DESCENT





   tvgtime = time.time()-tvgtime


   gdist=current_image.dist_to(pre_image)
   print "adist = %f, gdist = %f"%(adist,gdist)

   if gdist>=conv*adist and ddist > meps:
      gf*=red_fact

   reg_fac*=red_reg

   print 'gf = ',gf, '; reg_fac = ',reg_fac, '; meps = ',meps
   print "elapsed time (sec): ",time.time()-basetime," grad desc. time: ",tvgtime








