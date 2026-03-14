# TO DO: 
#  (1) see if we can rename lino_... and if
#      python can figure out from the type if its
#      sino or lino.

# (2) lino_project_to for genobj seems to be broken

import numpy
import configs3D, linoConfigs3D, sphericalRadonConfigs3D
import view_proj3D, ray_proj3D, radon3D
import view_proj3D_tomo as vpt
import totalvar as tv
import tomoscale
import tomo3Dfort
from sys import argv

from scipy import ndimage
convolve = ndimage.convolve



random = numpy.random

pi = numpy.pi
array = numpy.array
default_type = numpy.float32
#default_type = numpy.float64
indicator_type = numpy.int8
sum = numpy.sum
dot = numpy.dot
fromfile = numpy.fromfile
sqrt = numpy.sqrt
sin = numpy.sin
cos = numpy.cos
log=numpy.log
exp=numpy.exp
arctan2=numpy.arctan2
fftn= numpy.fft.fftn
ifftn= numpy.fft.ifftn
fftshift= numpy.fft.fftshift
ifftshift= numpy.fft.ifftshift
ones=numpy.ones
zeros=numpy.zeros
float32 = numpy.float32
float64 = numpy.float64

default_data_shape = 32,128,128
default_image_shape = 128,128,128 
default_num_rays = 1000000

#mat_path="/home/sidky/python/PyTomo3Diterative/materials_data/"
mat_path="/home/sidky/python/PyTomo3Diterative/materials_data/"


bp_algs_dict = {\
   "adj_nearest":view_proj3D.sinobackproject_nearest,\
   "adj_linear":view_proj3D.sinobackproject_linear,\
   "pd_nearest":view_proj3D.pd_backproj_nearest,\
   "pd_linear":view_proj3D.pd_backproj_linear\
}


proj_algs_dict = {\
   "proj_linear":view_proj3D.sinoproj_linear,\
   "proj_nearest":view_proj3D.sinoproj_nearest,\
   "zproj_linear":vpt.sinoprojz\
}


art_algs_dict = {\
   "art_linear":view_proj3D.artstep_linear,\
   "art_nearest":view_proj3D.artstep_nearest,\
   "zart_linear":vpt.artstepz\
}


em_algs_dict = {\
   "em_linear":view_proj3D.emstep_linear,\
   "em_nearest":view_proj3D.emstep_nearest,\
   "zem_linear":vpt.emstep_linearz\
}


ray_bp_algs_dict = {\
   "adj_nearest":ray_proj3D.sinobackproject_nearest\
}


ray_proj_algs_dict = {\
   "proj_nearest":ray_proj3D.sinoproj_nearest,\
   "proj_linear":ray_proj3D.sinoproj_linear\
}


ray_art_algs_dict = {\
   "art_nearest":ray_proj3D.artstep_nearest,\
   "art_linear":ray_proj3D.artstep_linear\
}


ray_em_algs_dict = {\
   "em_nearest":ray_proj3D.emstep_nearest\
}


########################################################



class sinogram3D(object):
   '''3D sinogram array for tomography'''

   def __init__(self,\
      config_name='circular_conebeam',\
      parms={"source_to_detector":10.,\
             "radius":5.},\
      shape= default_data_shape,\
      slen=pi, ulen=4.0, vlen=4.0,\
      s0=0., u0=-2.0, v0=-2.0,\
      s_include_endpoints=1):
      '''Initializes the sinogram array'''

      self.ns=shape[0]
      self.nu=shape[1]
      self.nv=shape[2]

      self.mat=zeros(shape,default_type)
      self.indicator=ones(shape,indicator_type)
      self.frame_vectors=zeros((self.ns,18) , default_type)
      self.config_name=config_name
      self.parms=parms

      self.frame=configs3D.configs[config_name]

      self.slen=slen
      self.ulen=ulen
      self.vlen=vlen
      self.s_include_endpoints=s_include_endpoints
      self.s0=s0
      self.u0=u0
      self.v0=v0
      self.ds=self.get_ds()
      self.du=self.get_du()
      self.dv=self.get_dv()
      for ip in range(self.ns):
         s=self.s0+ip*self.ds
         self.frame_vectors[ip]=self.frame(s,self.parms)


   def __str__(self):
      '''Gives sinogram matrix parameters'''


      nset=sum(sum(sum(self.indicator)))
      frac=(100.0 * nset)/(self.ns * self.nu * self.nv)
      lulu =  '''
The sinogram is %f %% full.
ns = %d  
nu = %d
nv = %d
ds = %f
du = %f
dv = %f
s0 = %f
u0 = %f
v0 = %f''' % (frac,self.ns,self.nu,self.nv,self.ds,self.du,self.dv,self.s0,self.u0,self.v0)

      return lulu

   def duplicate(self):
      
      new_sino = sinogram3D(\
         config_name=self.config_name,\
         parms=self.parms,\
         shape= (self.ns,self.nu,self.nv),\
         slen=self.slen, ulen=self.ulen, vlen=self.vlen,\
         s0=self.s0, u0=self.u0, v0=self.v0,\
         s_include_endpoints=self.s_include_endpoints)
      new_sino.mat = self.mat.copy()
      new_sino.indicator = self.indicator.copy()
      new_sino.frame_vectors = self.frame_vectors 

      return new_sino


   def get_ds(self):
      if self.s_include_endpoints:
         return self.slen/(self.ns-1.)
      else:
	 return self.slen/self.ns

   def get_du(self):
      return self.ulen/self.nu

   def get_dv(self):
      return self.vlen/self.nv

   def get_parms(self):
      return self.parms


   def read(self,file_name,data_type=default_type):
      '''reads in sinogram data.
In later versions, we need to implement selective reading of
sinogram arrays. '''
     
      sino_shape=self.mat.shape 
      self.mat=fromfile(file_name,dtype=data_type)
      self.mat.shape = sino_shape



   def flip_detector(self):
      '''Exchange u and v unit vectors and rotate detector accordingly.
Useful for parallel programs'''
      for ip in range(self.ns):
         eutemp = self.frame_vectors[ip,9:12]*1.
         evtemp = self.frame_vectors[ip,12:15]*1.
         self.frame_vectors[ip,12:15] =eutemp*1.
         self.frame_vectors[ip,9:12] =evtemp*1.
      tempmat =self.mat.transpose((0,2,1))*1.
      self.mat = tempmat

      self.shape = 1*self.mat.shape
      self.nu = self.shape[1]*1
      self.nv = self.shape[2]*1

      tempind=self.indicator.transpose((0,2,1))*1
      self.indicator = tempind
      tempval = self.u0*1.
      self.u0 = self.v0*1.
      self.v0 = tempval*1.
      tempval = self.ulen*1.
      self.ulen = self.vlen*1.
      self.vlen = tempval*1.
      tempval = self.du*1.
      self.du= self.dv*1.
      self.dv= tempval*1.

   def scramble(self,spacing):
      '''This method mixes up the view order of the
frame vector array with the hope of speeding up ART.
If it works, then come back and make this function more robust.'''

      nviews = self.ns
      scramble_indeces = array(range(nviews))
      scramble_indeces.shape = (nviews/spacing,spacing)
      scramble_indeces = scramble_indeces.transpose()
      scramble_indeces = scramble_indeces.flatten()

      self.frame_vectors = self.frame_vectors[scramble_indeces]




   def write_file(self,file_name):
      
      self.mat.tofile(file_name,sep='',format='')



   def clear_indicator(self):
      self.indicator=zeros([self.ns,self.nu,self.nv],indicator_type)



   def config_indicator(self,slow,shigh,ulow,uhigh,vlow,vhigh,ray_op='add'):
      '''Turns on rays within the specified s- ,  u-, and v-range'''


      if ray_op=='replace':
	 self.clear_indicator()



# make sure points are on the sinogram grid
      sl=max(self.s0,self.s0 + self.ds*round((slow-self.s0)/self.ds))
      sh=min(self.s0+self.ds*(self.ns-1),self.s0 + self.ds*round((shigh-self.s0)/self.ds))
      ul=max(self.u0,self.u0+self.du*round((ulow-self.u0)/self.du))
      uh=min(self.u0+self.du*(self.nu-1),self.u0 + self.du*round((uhigh-self.u0)/self.du))
      vl=max(self.v0,self.v0+self.dv*round((vlow-self.v0)/self.dv))
      vh=min(self.v0+self.dv*(self.nv-1),self.v0 + self.dv*round((vhigh-self.v0)/self.dv))
      
      s=sl
      while s<=sh:
         i=int(round((s-self.s0)/self.ds))
	 u=ul
	 while u<=uh:
            j=int(round((u-self.u0)/self.du))

            v=vl
	    while v<vh:
               k=int(round((v-self.v0)/self.dv))
	       self.indicator[i,j,k]=1
	       
	       v+=self.dv
	    u+=self.du
	 s+=self.ds


   def frame_vecs(self,s):
      '''Returns the config vectors for the trajectory parm s'''

      clist = self.frame(s,self.parms)
      tloc = array(clist[0:3])
      ttan = array(clist[3:6])
      detc = array(clist[6:9])
      euhat = array(clist[9:12])
      evhat = array(clist[12:15])
      ewhat = array(clist[15:18])

      return tloc,ttan, detc,euhat, evhat,ewhat


   def clear(self):
      self.mat*=0.


   def dist_to(self,sino):
      '''Calculates vectorial distance to sino'''
      
      dist = sqrt( sum(self.indicator.flatten()*(self.mat.flatten()-sino.mat.flatten())**2,dtype = float64) )

      return dist


   def copy_to(self,dest_sino):

      dest_sino.mat=self.mat.copy()
      dest_sino.indicator=self.indicator.copy()



   def add_noise(self,std=0.1):

      for i in range(self.ns):
         for j in range(self.nu):
            for k in range(self.nv):

               factor=max(0.,random.normal(1.,
std))
               self.mat[i,j,k]*=factor

   def add_abs_noise(self,std=0.1):

      for i in range(self.ns):
         for j in range(self.nu):
            for k in range(self.nv):

               if self.indicator[i,j,k]:
                  noise_pt = random.normal(0.,std)
                  self.mat[i,j,k]+= noise_pt



   def add_rel_noise(self,delta=0.002):

      for i in range(self.ns):
         for j in range(self.nu):
            for k in range(self.nv):

               if self.indicator[i,j,k]:
                  val = self.mat[i,j,k]
                  std = delta*(1.0+val/2.)
                  noise_pt = random.normal(0.,std)
                  self.mat[i,j,k]+= noise_pt

   def add_poisson_noise(self,photon_num_per_bin=1000):
      '''Introduces Poisson noise assuming constant number of
photons incident at each bin.
Change later to include geometry of bins.
Also make sure sinogram values are att. coef. x physical length'''

      for i in range(self.ns):
         for j in range(self.nu):
            for k in range(self.nv):

               sinoval = self.mat[i,j,k]
               phnum=exp(-sinoval)*photon_num_per_bin
               sig=(0.00001+sinoval)*sqrt(phnum)/(0.001+phnum)
               newval = numpy.random.normal(sinoval,sig)
               self.mat[i,j,k]=newval


   def generate_sigs(self,photon_num_per_bin=100):
      '''Computes sigmas according to Poisson dist. with large numbers'''

      sinovals = self.mat
      phnum=exp(-sinovals)*photon_num_per_bin
      sigs=(0.00001+sinovals)*sqrt(phnum)/(0.001+phnum)

      return sigs


   def homoscale(self):
      '''Divide each ray in the sinogram by the length of the ray'''
     
      self.mat = tomo3Dfort.sinoscale(self.mat,self.frame_vectors,\
         self.ns,self.nu,self.nv,self.du,self.dv,self.u0,self.v0)


   def backgroundsum(self,nps,npu,npv,coeffs): 

      ns = self.ns
      nu = self.nu
      nv = self.nv
      ds = self.ds
      du = self.du
      dv = self.dv
      u0 = self.u0
      v0 = self.v0

      self.mat = view_proj3D.backsum(self.indicator,\
                    ns,nu,nv,ds,du,dv,u0,v0,nps,npu,npv, coeffs)

   def backgroundsum_adj(self,nps,npu,npv): 

      ns = self.ns
      nu = self.nu
      nv = self.nv
      ds = self.ds
      du = self.du
      dv = self.dv
      u0 = self.u0
      v0 = self.v0

      coeffs = view_proj3D.backsumadj(self.mat,self.indicator,\
                    ns,nu,nv,ds,du,dv,u0,v0,nps,npu,npv)

      return coeffs

   def backproject_to(self,image,alg='adj_linear'):
      '''Backproject onto image:
adj_nearest: adjoint of the nearest neighbor projector
adj_linear:  adjoint of the linear interp. projector
pd_nearest:  pixel driven w/ nearest neighbor interp.
pd_linear:   pixel driven w/ linear interp.'''
        
      if alg not in bp_algs_dict.keys():
         print("bad algorithm choice for backprojection")
         sys.exit()

      bp_alg = bp_algs_dict[alg]

      ns=self.ns
      nu=self.nu
      nv=self.nv
      ds=self.ds
      du=self.du
      dv=self.dv
      s0=self.s0
      u0=self.u0
      v0=self.v0


      dx=image.dx
      dy=image.dy
      dz=image.dz
      x0=image.x0
      y0=image.y0
      z0=image.z0
      nx=image.nx
      ny=image.ny
      nz=image.nz


      image.mat=bp_alg(\
           self.mat,\
           self.indicator,\
           self.frame_vectors,\
           ns,nu,nv,du,dv,u0,v0,\
           dx,dy,dz,x0,y0,z0,nx,ny,nz)






### END OF SINOGRAM 3D OBJECT ###



########################################################
class raybunch(sinogram3D):

   def __init__(self,sino,nds,ndu,ndv,nws,nwu,nwv):

      self.nviews = sino.ns
      self.nds = nds
      self.ndu = ndu
      self.ndv = ndv
      self.nws = nws
      self.nwu = nwu
      self.nwv = nwv
      sinogram3D.__init__(self,\
         config_name=sino.config_name,\
         parms=sino.parms,\
         shape= (sino.ns/nds,sino.nu/ndu,sino.nv/ndv),\
         slen=sino.slen, ulen=sino.ulen, vlen=sino.vlen,\
         s0=sino.s0, u0=sino.u0, v0=sino.v0,\
         s_include_endpoints=sino.s_include_endpoints)


   def get_data(self,sino):
      ns = sino.ns
      nu = sino.nu
      nv = sino.nv
      nss = self.ns
      nus = self.nu
      nvs = self.nv
      nds = self.nds
      ndu = self.ndu
      ndv = self.ndv
      nws = self.nws
      nwu = self.nwu
      nwv = self.nwv
      
      self.mat= tomo3Dfort.sinoave(sino.mat,ns,nu,nv,nss,nus,nvs,nds,ndu,ndv,nws,nwu,nwv)

   def put_data(self,sino):
      ns = sino.ns
      nu = sino.nu
      nv = sino.nv
      nss = self.ns
      nus = self.nu
      nvs = self.nv
      nds = self.nds
      ndu = self.ndu
      ndv = self.ndv
      nws = self.nws
      nwu = self.nwu
      nwv = self.nwv
      
      sino.mat= tomo3Dfort.sinodist(self.mat,ns,nu,nv,nss,nus,nvs,nds,ndu,ndv,nws,nwu,nwv)


########################################################



class linogram3D(object):
   '''3D linogram array for tomography'''

   def __init__(self,\
      config_name='ball',\
      parms={\
        "source_radius":1.,\
        "det_radius":1.},\
      nrays = default_num_rays,\
      is_duplicate = 0):
      '''Initializes the linogram'''


      self.mat=zeros([nrays],default_type)
      self.config_name=config_name
      self.parms=parms
      self.nr = nrays

      self.frame=linoConfigs3D.configs[config_name]
      if not is_duplicate:
         self.ray_vectors=zeros((self.nr,6) , default_type)
         self.indicator=ones([nrays],indicator_type)




   def duplicate(self):

      new_lino = linogram3D(\
         config_name=self.config_name,\
         parms=self.parms.copy(),\
         nrays= self.nr,\
         is_duplicate = 1)
      new_lino.mat = self.mat.copy()
# below we copy only the pointers in order to save space
# because ray_vectors is huge.
      new_lino.indicator = self.indicator
      new_lino.ray_vectors = self.ray_vectors

      return new_lino



   def write_file(self,file_name):

      self.mat.tofile(file_name,sep='',format='')

   def frame_vecs(self,s):
      '''Returns the config vectors for the trajectory parm s'''

      clist = self.frame(s,self.parms)
      sloc = array(clist[0:3])
      bloc = array(clist[3:6])

      return sloc, bloc



   def dist_to(self,lino):
      '''Calculates vectorial distance to lino'''

      ldiff=self.mat-lino.mat
      dist = sqrt(sum(ldiff*ldiff, dtype = float64))

      return dist


   def copy_to(self,dest_lino):

      dest_lino.mat=self.mat.copy()
      dest_lino.indicator=self.indicator.copy()


   def add_noise(self,std=0.1):

      for i in range(self.nr):
         factor=max(0.,random.normal(1.,std))
         self.mat[i]*=factor


   def backproject_to(self,image):
      '''Backproject onto image:
Adjoint of ray-driven projection with nearest-neighbor interpolation.'''

      nr=self.nr


      dx=image.dx
      dy=image.dy
      dz=image.dz
      x0=image.x0
      y0=image.y0
      z0=image.z0
      nx=image.nx
      ny=image.ny
      nz=image.nz



      image.mat=ray_proj3D.sinobackproject_nearest(\
           self.mat,\
           self.indicator,\
           self.ray_vectors,\
           nr,\
           dx,dy,dz,x0,y0,z0,nx,ny,nz)





### END OF LINOGRAM 3D OBJECT ###




########################################################



class sphericalRadon3D(object):
   '''3D spherical Radon data function for ,e.g., TAT'''

   def __init__(self,\
      config_name='z_sphere',\
      parms={"radius":10.},\
      shape= default_data_shape,\
      salen=10., sblen=2.*pi, rlen=10.0,\
      sa0=0., sb0=0.0, r0=0.0\
      ):
      '''Initializes the data array'''

      self.c0 = 0.154  #speed of sound in cm/microsecond

      self.nsa=shape[0]
      self.nsb=shape[1]
      self.nr=shape[2]

      self.ns =  self.nsa * self.nsb    # compute total number of source locations

      self.mat=zeros(shape,default_type)
      self.indicator=ones(shape,indicator_type)
      self.source_locations=zeros((self.ns,3) , default_type)
      self.config_name=config_name
      self.parms=parms

      self.frame=sphericalRadonConfigs3D.configs[config_name]

      self.salen=salen
      self.sblen=sblen
      self.rlen=rlen
      self.sa0=sa0
      self.sb0=sb0
      self.r0=r0
      self.dsa=self.get_dsa()
      self.dsb=self.get_dsb()
      self.dr=self.get_dr()

      self.t0 = self.r0/self.c0
      self.tlen = self.rlen/self.c0
      self.dt = self.dr/self.c0

      ip = 0
      for jp in range(self.nsa):
         for kp in range(self.nsb):
            sa=self.sa0+jp*self.dsa
            sb=self.sb0+kp*self.dsb
            self.source_locations[ip]=self.frame(sa,sb,self.parms)
            ip += 1



   def duplicate(self):
      
      new_radon= sphericalRadon3D(\
         config_name=self.config_name,\
         parms=self.parms,\
         shape= (self.nsa,self.nsb,self.nr),\
         salen=self.salen, sblen=self.sblen, rlen=self.rlen,\
         sa0=self.sa0, sb0=self.sb0, r0=self.r0)
      new_radon.mat = self.mat.copy()
      new_radon.indicator = self.indicator.copy()

      return new_radon


   def get_dsa(self):
      return self.salen/(self.nsa-1.)

   def get_dsb(self):
      return self.sblen/(self.nsb-1.)

   def get_dr(self):
      return self.rlen/(self.nr)

   def get_parms(self):
      return self.parms


   def write_file(self,file_name):

      self.mat.tofile(file_name,sep='',format='')

   def source_locs(self,sa,sb):
      '''Returns the config vectors for the trajectory parm s'''

      clist = self.frame(sa,sb,self.parms)

      return clist



   def dist_to(self,radon):
      '''Calculates vectorial distance to radon'''

      rdiff=self.mat-radon.mat
      dist = sqrt(sum(rdiff*rdiff, dtype = float64))

      return dist


   def copy_to(self,dest_radon):

      dest_radon.mat=self.mat.copy()
      dest_radon.indicator=self.indicator.copy()





### END OF SPHERICAL RADON 3D OBJECT ###











#########################################

class image3D(object):
   '''3D image array for tomography'''

   def __init__(self,\
                shape=default_image_shape,\
                xlen=2.0,ylen=2.0,zlen=2.0,\
		x0=-1.0,y0=-1.0,z0=-1.0,
		data_type='real'):
      '''Initializes the image array'''

      if data_type=='complex':
         self.mat=zeros(shape,default_complex_type)
      else:
         self.mat=zeros(shape,default_type)
      self.indicator = ones(shape,indicator_type)
      self.nx=shape[0]
      self.ny=shape[1]
      self.nz=shape[2]
      self.x0=x0
      self.y0=y0
      self.z0=z0
      self.xlen=xlen
      self.ylen=ylen
      self.zlen=zlen
      self.dx=self.get_dx()
      self.dy=self.get_dy()
      self.dz=self.get_dz()


   def __str__(self):
      '''Gives image matrix parameters'''
      lulu =  '''
nx = %d  
ny = %d
nz = %d
dx = %f
dy = %f
dz = %f
x0 = %f
y0 = %f
z0 = %f
xlen = %f
ylen = %f
zlen = %f''' % (self.nx,self.ny,self.nz,\
              self.dx,self.dy,self.dz,\
	      self.x0,self.y0,self.z0,\
              self.xlen,self.ylen,self.zlen)

      return lulu


   def duplicate(self):
      
      new_image = image3D(\
         shape= (self.nx,self.ny,self.nz),\
         xlen=self.xlen, ylen=self.ylen, zlen=self.zlen,\
         x0=self.x0, y0=self.y0, z0=self.z0)
      new_image.mat = self.mat.copy()

      return new_image


   def get_dx(self):
      return self.xlen/self.nx

   def get_dy(self):
      return self.ylen/self.ny

   def get_dz(self):
      return self.zlen/self.nz

   def transpose(self,indexlist):
      tempmat = 1.*self.mat.transpose(indexlist)
      self.mat = tempmat
      temp = (self.nx,self.ny,self.nz)
      self.nx = temp[indexlist[0]]
      self.ny = temp[indexlist[1]]
      self.nz = temp[indexlist[2]]
      temp = (self.x0,self.y0,self.z0)
      self.x0 = temp[indexlist[0]]
      self.y0 = temp[indexlist[1]]
      self.z0 = temp[indexlist[2]]
      temp = (self.xlen,self.ylen,self.zlen)
      self.xlen = temp[indexlist[0]]
      self.ylen = temp[indexlist[1]]
      self.zlen = temp[indexlist[2]]
      temp = (self.dx,self.dy,self.dz)
      self.dx = temp[indexlist[0]]
      self.dy = temp[indexlist[1]]
      self.dz = temp[indexlist[2]]


   def read(self,file_name,data_type=default_type):

      image_shape=self.mat.shape
      self.mat=fromfile(file_name,dtype=data_type)
      self.mat.shape=image_shape


   def write_file(self,file_name):
      
      self.mat.tofile(file_name,sep='',format='')



   def copy_to(self,dest_image):

      dest_image.mat=self.mat.copy()
      dest_image.indicator=self.indicator.copy()






   def clear(self):
      self.mat*=0.

   def make_positive(self):
      self.mat[self.mat<0.]=0.

   def make_positive_softly(self,beta=1.0):
      self.mat[self.mat<0.]*=(1.0-beta)


   def clamp_edges(self):

      self.mat[0,:,:]=self.mat[2,:,:]*1.
      self.mat[1,:,:]=self.mat[2,:,:]*1.
      self.mat[-1,:,:]=self.mat[-3,:,:]*1.
      self.mat[-2,:,:]=self.mat[-3,:,:]*1.

      self.mat[:,0,:]=self.mat[:,2,:]*1.
      self.mat[:,1,:]=self.mat[:,2,:]*1.
      self.mat[:,-1,:]=self.mat[:,-3,:]*1.
      self.mat[:,-2,:]=self.mat[:,-3,:]*1.

      self.mat[:,:,0]=self.mat[:,:,2]*1.
      self.mat[:,:,1]=self.mat[:,:,2]*1.
      self.mat[:,:,-1]=self.mat[:,:,-3]*1.
      self.mat[:,:,-2]=self.mat[:,:,-3]*1.



   def zero_edges(self):

      self.mat[0,:,:]=0.
      self.mat[1,:,:]=0.
      self.mat[-1,:,:]=0.
      self.mat[-2,:,:]=0.

      self.mat[:,0,:]=0.
      self.mat[:,1,:]=0.
      self.mat[:,-1,:]=0.
      self.mat[:,-2,:]=0.

      self.mat[:,:,0]=0.
      self.mat[:,:,1]=0.
      self.mat[:,:,-1]=0.
      self.mat[:,:,-2]=0.


   def dist_to(self,image):
      '''Calculates vectorial distance to image'''
      
      dist = sqrt( sum((self.mat.flatten()-image.mat.flatten())**2,dtype = float64) )

      return dist

   def calc_total_var(self,lam1=0.):
      '''Calculates total variation of image'''

      tmag=float(tv.total_var3d(lam1,self.mat,self.nx,self.ny,self.nz))

      return tmag


   def calc_total_var_p(self,p=1.):
      '''Calculates total variation of image'''

      tmag=float(tv.total_var3d_p(p,self.mat,self.nx,self.ny,self.nz))

      return tmag

   def calc_lp(self,p=1.):
      '''Calculates total variation of image'''

      tmag=float(tv.calc_lp(p,self.mat,self.nx,self.ny,self.nz))

      return tmag


   def calc_total_dvar_p(self,p=1.):
      '''Calculates total variation of image'''

      tmag=float(tv.total_dvar3d_p(p,self.mat,self.nx,self.ny,self.nz,self.dx,self.dy,self.dz))

      return tmag

   def lp_grad(self,p=1.,normalizeq=1):
      '''Calculates total variation of image'''

      self.mat,bobo=tv.lp_grad3d(p,self.mat,self.nx,self.ny,self.nz,normalizeq)

      return bobo


   def tpv_grad(self,p=1.,normalizeq=1):
      '''Calculates total variation of image'''

      self.mat,bobo=tv.total_p_var_grad3d(p,self.mat,self.nx,self.ny,self.nz,normalizeq)

      return bobo



   def tpdv_grad(self,p=1.,normalizeq=1):
      '''Calculates total variation of image'''

      self.mat,bobo=tv.total_p_dvar_grad3d(p,self.mat,self.nx,self.ny,self.nz,self.dx,self.dy,self.dz,normalizeq)

      return bobo


   def tomoscalevol_to(self,image,xs,ys,zs):
      '''Scales tomo vol.'''

      image.mat=tomoscale.tomovolmag(self.mat,self.nx,self.ny,self.nz,\
                                             self.x0,self.y0,self.z0,\
                                             self.xlen,self.ylen,self.zlen,\
                                             xs,ys,zs)


   def tomoscalevol_trans_to(self,image,xs,ys,zs):
      '''Scales tomo vol. transpose'''

      image.mat=tomoscale.tomovolmagt(self.mat,self.nx,self.ny,self.nz,\
                                             self.x0,self.y0,self.z0,\
                                             self.xlen,self.ylen,self.zlen,\
                                             xs,ys,zs)

   def gmi(self):
      '''Calculates gradient magnitude image'''

      self.mat=tv.gmi(self.mat,self.nx,self.ny,self.nz)


   def tv_grad(self,lam1=0.,normalizeq=1):
      '''Calculates total variation of image'''

      self.mat,bobo=tv.total_var_grad3d(lam1,self.mat,self.nx,self.ny,self.nz,normalizeq)

      return bobo



   def scale(self,factor):
      '''Scales the image by factor'''

      self.mat=factor*self.mat


   def shift(self):
      '''performs scipy inverse fft shift on image'''

      self.mat=ifftshift(self.mat)

   def inv_shift(self):
      '''performs scipy fft shift on image'''

      self.mat=fftshift(self.mat)

   def fft(self):
      '''What do you think?!!'''

      self.mat=fftn(self.mat)

   def inv_fft(self):
      '''What do you think?!!'''

      self.mat=ifftn(self.mat)

   def real_part(self):

      self.mat=real(self.mat)

   def imag_part(self):
      self.mat=imag(self.mat)

   def add_shape(self,shape):
      '''Puts a 3D shape in the image array
	
The attenuation value for the shape is added to each voxel
whose center is in the shape'''

      for i in range(self.nx):
         for j in range(self.ny):
            for k in range(self.nz):
               x=self.x0+(i+0.5)*self.dx
               y=self.y0+(j+0.5)*self.dy
               z=self.z0+(k+0.5)*self.dz
               self.mat[i,j,k]=self.mat[i,j,k]+shape.voxvalf(x,y,z)






   def mat_project_to(self,sino,projmats):
      '''Use hologic projection matrix formalism'''
        
      proj_alg = view_proj3D.matproj_linear
      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0




      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz


      sino.mat=proj_alg(\
           projmats,\
           ns,nu,nv,du,dv,u0,v0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)



   def project_to(self,sino,alg='proj_linear'):
      '''project onto sinogram:
proj_nearest: nearest neighbor projector
proj_linear:  linear interp. projector
zproj_nearest:  nearest neighbor projector, going thru image slices only in z-direction
'''
        
      if alg not in proj_algs_dict.keys():
         print("bad algorithm choice for projection")
         sys.exit()

      proj_alg = proj_algs_dict[alg]

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0




      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz


      sino.mat=proj_alg(\
           sino.indicator,\
           sino.frame_vectors,\
           ns,nu,nv,du,dv,u0,v0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)



   def lino_project_to(self,lino):
      '''fill in sinogram data with projections of image

WITH FORTRAN ACCELERATION!!'''

      nr=lino.nr


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz


      lino.mat=ray_proj3D.sinoproj_linear(\
           lino.indicator,\
           lino.ray_vectors,\
           nr,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)


   def radon_project_to(self,radon,fov = -1.0 ):
      '''do spherical Radon projection

WITH FORTRAN ACCELERATION!!'''

      nsa = radon.nsa
      nsb = radon.nsb
      ns = radon.ns
      nr = radon.nr
      dr = radon.dr
      r0 = radon.r0


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz

      if fov<0.:
         fov = self.xlen/2.


      radon.mat=radon3D.spherical_radon_linear(\
           fov,\
           radon.indicator,\
           radon.source_locations,\
           nsa,nsb,ns,nr,dr,r0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)




   def art_step(self,sino,reg_fac=1.0,alg='art_linear'):
      '''art step on image with sinogram:
art_nearest: nearest neighbor ART
art_linear:  linear interp. ART
zart_nearest:  nearest neighbor ART, going thru image slices only in z-direction
'''
        
      if alg not in art_algs_dict.keys():
         print("bad algorithm choice for ART")
         sys.exit()

      art_alg = art_algs_dict[alg]

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz



      self.mat=art_alg(\
           reg_fac,\
           sino.mat,\
           sino.indicator,\
           sino.frame_vectors,\
           ns,nu,nv,du,dv,u0,v0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)




   def tomo_art_step_2vol(self,iaux,sino,reg_fac=1.0):
      '''take an ART step'''

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz

      dxa=iaux.dx
      dya=iaux.dy
      dza=iaux.dz
      x0a=iaux.x0
      y0a=iaux.y0
      z0a=iaux.z0
      nxa=iaux.nx
      nya=iaux.ny
      nza=iaux.nz



      self.mat, iaux.mat=vpt.artstepz2vol(\
           reg_fac,\
           sino.mat,\
           sino.indicator,\
           sino.frame_vectors,\
           ns,nu,nv,du,dv,u0,v0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz,\
           iaux.mat,dxa,dya,dza,x0a,y0a,z0a,nxa,nya,nza)



   def lino_art_step(self,lino,reg_fac=1.0):
      '''take an ART step'''

      nr=lino.nr


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz



      self.mat=ray_proj3D.artstep_linear(\
           reg_fac,\
           lino.mat,\
           lino.indicator,\
           lino.ray_vectors,\
           nr,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)


   def radon_art_step(self,radon,reg_fac=1.0,fov = -1.0 ):
      '''take a radon art step'''

      nsa = radon.nsa
      nsb = radon.nsb
      ns = radon.ns
      nr = radon.nr
      dr = radon.dr
      r0 = radon.r0


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz

      if fov<0.:
         fov = self.xlen/2.


      self.mat=radon3D.artstep_linear(\
           reg_fac,fov,\
           radon.mat,\
           radon.indicator,\
           radon.source_locations,\
           nsa,nsb,ns,nr,dr,r0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)






   def em_step(self,sino,alg='em_nearest'):
      '''em step on image with sinogram:
em_nearest: nearest neighbor EM
em_linear:  linear interp. EM
'''
        
      if alg not in em_algs_dict.keys():
         print("bad algorithm choice for EM")
         sys.exit()

      em_alg = em_algs_dict[alg]


      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz



      self.mat=em_alg(\
           sino.mat,\
           sino.indicator,\
           sino.frame_vectors,\
           ns,nu,nv,du,dv,u0,v0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

   def emu_step(self,sino,uimage,beta):
      '''em step linear with 
objective gradient: uimage.
'''
        

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz



      self.mat=view_proj3D.emustep_linear(\
           sino.mat,\
           sino.indicator,\
           sino.frame_vectors,\
           ns,nu,nv,du,dv,u0,v0,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz,
           uimage.mat,beta)



   def lino_em_step(self,lino):
      '''take an EM step with nearest neighbor interp.'''

      nr=lino.nr


      dx=self.dx
      dy=self.dy
      dz=self.dz
      x0=self.x0
      y0=self.y0
      z0=self.z0
      nx=self.nx
      ny=self.ny
      nz=self.nz



      self.mat=ray_proj3D.emstep_nearest(\
           lino.mat,\
           lino.indicator,\
           lino.ray_vectors,\
           nr,\
           self.mat,dx,dy,dz,x0,y0,z0,nx,ny,nz)





### END OF IMAGE 3D OBJECT ###

class phantom3D(object):
   '''A collection of 3D shapes

(later we will add images)'''
   
   def __init__(self):
      self.num_components = 0
      self.components=[]

   def __str__(self):
      lulu=""
      for i in range(1,self.num_components+1):

	      lulu+="component %d: "%(i)+ self.components[i-1].__str__()+"\n"

      return lulu

   def __add__(self,other):
      '''Creates a new phantom by combining objects in self and other'''

      new_phantom=phantom3D()

      for i in self.components:
         new_phantom.add_component(i)

      for i in other.components:
         new_phantom.add_component(i)


      return new_phantom



   def add_component(self,component):
      '''adds a shape to the phantom'''
      self.num_components+=1
      self.components.append(component)

   def collapse_to(self,image):
      '''collapse the 3D phantom to an image array or a chord.

Need to add in possibility of taking an image as part of phantom.
Note: this function is additive, so be sure to clear the image array first'''

      for i in range(0,self.num_components):
	 image.add_shape(self.components[i])


   def embed_in(self,image):
      '''embed the 3D phantom in an image array'''

      for i in self.components:
	 i.embed_in(image)


   def project_to(self,sino):
      '''Projects phantom onto a 3D sinogram.'''

      for i in range(0,self.num_components):
         print(i," : ",self.num_components)
         self.components[i].project_to(sino)


   def lino_project_to(self,lino):
      '''Projects phantom onto a 3D linogram.'''

      for i in range(0,self.num_components):
         self.components[i].lino_project_to(lino)



## WHY IS THIS HERE? IT SHOULD BE A METHOD OF TISSUE_MAP!!!
   def polychromatic_project_to(self,sino,spectrum):
      '''Projects phantom onto a 3D sinogram including the
effect of x-ray polychromaticity.

THIS method should only be called if phantom is composed of
tissue maps.'''

      num_energy=spectrum.shape[0]
      nview=sino.ns
      ndetu=sino.nu
      ndetv=sino.nv


      for i in range(nview):
         s=sino.s0+i*sino.ds
         
         exp_arg=zeros([num_energy,ndetu,ndetv],default_type)

         for itm in self.components:
            itm.project_to_single_view(sino,i)
            dens=itm.get_density()

            for en in range(num_energy):
               mac=itm.mass_attenuation_coef(spectrum[en,1])
               exp_arg[en] -= dens*mac*sino.mat[i]

       
         att_factors=exp(exp_arg) 

## TRY TO FIGURE OUT WHY I COULDn't use the dot product
         prodata=zeros([ndetu,ndetv],default_type)
         for en in range(num_energy-1):
            del_e=spectrum[en+1,1]-spectrum[en-1,1]
            specval=spectrum[en,0]
            prodata+=specval*att_factors[en]*del_e
         sino.mat[i]=prodata

      

     


################################################
class tissue_map(phantom3D):
   '''A collection of analytic shapes that represent a single tissue.

The projection algorithm is modified to take an energy spectrum
for modeling source polychromaticity'''

# parameters for energy tables
# energy tabulated in log of energy in kev
   log_e1=0.
   log_e2= log(20000.)   
   num_e = 1000

   materials=('polymethyl_methacrylate',\
              'bone')
   density={'polymethyl_methacrylate':1.211,\
            'bone':1.9}
# densities are in grams per cubic cm

   def __init__(self,material='polymethyl_methacrylate',physical_material=0):

      if physical_material:
         if material not in self.materials:
            print("I don't have ",material)
         self.material = material
         macfilename=mat_path+"mac_"+material+".dat"
         self.mac=fromfile(macfilename,sep='\t',dtype=default_type)
         self.mac.shape= tissue_map.num_e,2

      phantom3D.__init__(self)
   
     
   def mass_attenuation_coef(self,energy):
      le1=self.log_e1
      le2=self.log_e2
      dle=(le2-le1)/(self.num_e-1)
      
      ipos=int((log(energy)-le1)/dle)


      w=(log(energy)-self.mac[ipos,0])/dle
      res=(1.-w)*self.mac[ipos,1] + w*self.mac[ipos+1,1]

      return exp(res)

   def get_density(self):
      return self.density[self.material]
      
      



####### Start of analytic shape objects ###############

class shape3D(object):
   '''geometric shape class for building 3D phantoms'''

   def __init__(self,x0=0.,y0=0.,z0=0.,att=1.0):
      '''(x0,y0,z0) is the center of the shape'''
      self.x0=x0
      self.y0=y0
      self.z0=z0
      self.att=att

   def __str__(self):
      lulu="""Shape parms:
attenuation = %f
center = (%f,%f,%f)\n""" % (self.att,self.x0,self.y0,self.z0)
      return lulu

   def voxvalf(self,x,y,z):
      pass

   def project_to(self,sino):
      pass

   def project_to_single_view(self,sino,nview):
      pass

   def embed_in(self,image):
      pass

########################################################
class surface(object):
   '''Surface is used by gen_object to define a volume'''

   def __init__(self,xc=0.,yc=0.,zc=0.,alpha=0.,beta=0.,gamma=0.):
      self.xc=xc
      self.yc=yc
      self.zc=zc

# We'll keep these angles around for now, but they should be removed!
      self.alpha=alpha
      self.beta=beta
      self.gamma=gamma

#     construct frame x-unit vector

      x=1.
      y=0.
      z=0.
     
      xp=x*cos(gamma) - y*sin(gamma)
      yp=x*sin(gamma) + y*cos(gamma)
      zp= z


      xpp=xp*cos(beta) + zp*sin(beta)
      ypp=yp
      zpp=-xp*sin(beta) + zp*cos(beta)

      xppp=xpp*cos(alpha) - ypp*sin(alpha)
      yppp=ypp*cos(alpha) + xpp*sin(alpha)
      zppp=zpp

      self.xu1=xppp
      self.xu2=yppp
      self.xu3=zppp

#     construct frame y-unit vector

      x=0.
      y=1.
      z=0.
      
      xp=x*cos(gamma) - y*sin(gamma)
      yp=x*sin(gamma) + y*cos(gamma)
      zp= z


      xpp=xp*cos(beta) + zp*sin(beta)
      ypp=yp
      zpp=-xp*sin(beta) + zp*cos(beta)

      xppp=xpp*cos(alpha) - ypp*sin(alpha)
      yppp=ypp*cos(alpha) + xpp*sin(alpha)
      zppp=zpp

      self.yu1=xppp
      self.yu2=yppp
      self.yu3=zppp

#     construct frame z-unit vector

      x=0.
      y=0.
      z=1.
      
      xp=x*cos(gamma) - y*sin(gamma)
      yp=x*sin(gamma) + y*cos(gamma)
      zp= z


      xpp=xp*cos(beta) + zp*sin(beta)
      ypp=yp
      zpp=-xp*sin(beta) + zp*cos(beta)

      xppp=xpp*cos(alpha) - ypp*sin(alpha)
      yppp=ypp*cos(alpha) + xpp*sin(alpha)
      zppp=zpp

      self.zu1=xppp
      self.zu2=yppp
      self.zu3=zppp



   def __str__(self):
      lulu="""surface: %s \n
  at (%.2f,%.2f,%.2f), unit normal (%.2f,%.2f,%.2f), and x-axis along (%.2f,%.2f,%.2f)"""\
    % (self.label,self.xc,self.yc,self.zc,self.zu1,self.zu2,self.zu3,\
                                          self.xu1,self.xu2,self.xu3)
      return lulu

########################
class planar(surface):

   def __init__(self,xc=0.,yc=0.,zc=0.,alpha=0.,beta=0.,gamma=0.):

      self.label='plane'
      surface.__init__(self,xc,yc,zc,alpha,beta,gamma)

########################
class ellipsoidal(surface):

   def __init__(self,xc=0.,yc=0.,zc=0.,alpha=0.,beta=0.,gamma=0.,ax=1.,ay=1.,az=1.):

      self.label='ellipsoid'
      self.ax=ax
      self.ay=ay
      self.az=az
      surface.__init__(self,xc,yc,zc,alpha,beta,gamma)

########################
class cylindrical(surface):

   def __init__(self,xc=0.,yc=0.,zc=0.,alpha=0.,beta=0.,gamma=0.,ax=1.,ay=1.):

      self.label='cylinder'
      self.ax=ax
      self.ay=ay
      surface.__init__(self,xc,yc,zc,alpha,beta,gamma)

########################
class conical(surface):

   def __init__(self,xc=0.,yc=0.,zc=0.,alpha=0.,beta=0.,gamma=0.,ax=1.,ay=1.):

      self.label='cone'
      self.ax=ax
      self.ay=ay
      surface.__init__(self,xc,yc,zc,alpha,beta,gamma)

########################################################
class gen_object(shape3D):
   '''Analytic shape created by the intersection volume of surfaces.'''

   def __init__(self,x0=0.,y0=0.,z0=0.,att=1.0,\
                nsurf=0,\
		surfaces=[]):

      self.nsurf=nsurf
      self.surfaces=surfaces
      shape3D.__init__(self,x0,y0,z0,att)

   def __str__(self):
      lulu='''Object center at (%.2f,.2%f,%.2f)
Object surfaces:
'''%(self.x0,self.y0,self.z0)

      for i in self.surfaces:

	      lulu+="\t"+i.__str__()+"\n"

      return lulu


   def add_surface(self,surface = planar()):

      self.nsurf+= 1
      self.surfaces.append(surface)



   def rotate(self,alpha=0.,beta=0.,gamma=0.):

     for i in self.surfaces:

# rotate center of surface about object center
	xsh=i.xc
	ysh=i.yc
	zsh=i.zc

     
        xp=xsh*cos(gamma) - ysh*sin(gamma)
        yp=xsh*sin(gamma) + ysh*cos(gamma)
        zp= zsh


        xpp=xp*cos(beta) + zp*sin(beta)
        ypp=yp
        zpp=-xp*sin(beta) + zp*cos(beta)

        xppp=xpp*cos(alpha) - ypp*sin(alpha)
        yppp=ypp*cos(alpha) + xpp*sin(alpha)
        zppp=zpp


	i.xc=xppp
	i.yc=yppp
	i.zc=zppp

#rotate frame
#   x-unit vect.
	xsh=i.xu1
	ysh=i.xu2
	zsh=i.xu3

        xp=xsh*cos(gamma) - ysh*sin(gamma)
        yp=xsh*sin(gamma) + ysh*cos(gamma)
        zp= zsh


        xpp=xp*cos(beta) + zp*sin(beta)
        ypp=yp
        zpp=-xp*sin(beta) + zp*cos(beta)

        xppp=xpp*cos(alpha) - ypp*sin(alpha)
        yppp=ypp*cos(alpha) + xpp*sin(alpha)
        zppp=zpp


	i.xu1=xppp
	i.xu2=yppp
	i.xu3=zppp

#   y-unit vect.
	xsh=i.yu1
	ysh=i.yu2
	zsh=i.yu3

        xp=xsh*cos(gamma) - ysh*sin(gamma)
        yp=xsh*sin(gamma) + ysh*cos(gamma)
        zp= zsh


        xpp=xp*cos(beta) + zp*sin(beta)
        ypp=yp
        zpp=-xp*sin(beta) + zp*cos(beta)

        xppp=xpp*cos(alpha) - ypp*sin(alpha)
        yppp=ypp*cos(alpha) + xpp*sin(alpha)
        zppp=zpp


	i.yu1=xppp
	i.yu2=yppp
	i.yu3=zppp

#   z-unit vect.
	xsh=i.zu1
	ysh=i.zu2
	zsh=i.zu3

        xp=xsh*cos(gamma) - ysh*sin(gamma)
        yp=xsh*sin(gamma) + ysh*cos(gamma)
        zp= zsh


        xpp=xp*cos(beta) + zp*sin(beta)
        ypp=yp
        zpp=-xp*sin(beta) + zp*cos(beta)

        xppp=xpp*cos(alpha) - ypp*sin(alpha)
        yppp=ypp*cos(alpha) + xpp*sin(alpha)
        zppp=zpp


	i.zu1=xppp
	i.zu2=yppp
	i.zu3=zppp




      
   def embed_in(self,image):

      nx=image.nx
      ny=image.ny
      nz=image.nz
      dx=image.dx
      dy=image.dy
      dz=image.dz
      x0_image=image.x0
      y0_image=image.y0
      z0_image=image.z0

      att=self.att
      nsurf=self.nsurf

      xu1=[]
      xu2=[]
      xu3=[]

      yu1=[]
      yu2=[]
      yu3=[]
      
      zu1=[]
      zu2=[]
      zu3=[]

      x0=[]
      y0=[]
      z0=[]

      nplane=1
      
      nellipsoid=1
      ax_el=[1.]
      ay_el=[1.]
      az_el=[1.]

      ncyl=1
      ax_cyl=[1.]
      ay_cyl=[1.]
      
      ncone=1
      ax_cone=[1.]
      ay_cone=[1.]

#put surfaces in the order: planes, ellipsoid, cylinder, and cone
      for i in self.surfaces:
	 if i.label=='plane':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    nplane+=1


      for i in self.surfaces:
	 if i.label=='ellipsoid':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_el.append(i.ax)
	    ay_el.append(i.ay)
	    az_el.append(i.az)

	    nellipsoid+=1

      for i in self.surfaces:
	 if i.label=='cylinder':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cyl.append(i.ax)
	    ay_cyl.append(i.ay)

	    ncyl+=1

      for i in self.surfaces:
	 if i.label=='cone':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cone.append(i.ax)
	    ay_cone.append(i.ay)

	    ncone+=1

      image.mat = image.mat+\
           tomo3Dfort.gen_obj_embed(\
           nx,ny,nz,dx,dy,dz,x0_image,y0_image,z0_image,\
           att,\
           nsurf,\
	   xu1,xu2,xu3,yu1,yu2,yu3,zu1,zu2,zu3,\
           x0,y0,z0,\
           nplane,\
           nellipsoid,\
           ax_el,ay_el,az_el,\
           ncyl,\
           ax_cyl,ay_cyl,\
           ncone,\
           ax_cone,ay_cone)


      
   def voxvalf(self,x,y,z):
      '''Returns attenuation value if x,y,z is in the general object'''

      att=self.att
      nsurf=self.nsurf

      alpha=[]
      beta=[]
      gamma=[]
      x0=[]
      y0=[]
      z0=[]

      nplane=1
      
      nellipsoid=1
      ax_el=[1.]
      ay_el=[1.]
      az_el=[1.]

      ncyl=1
      ax_cyl=[1.]
      ay_cyl=[1.]
      
      ncone=1
      ax_cone=[1.]
      ay_cone=[1.]

#put surfaces in the order: planes, ellipsoid, cylinder, and cone
      for i in self.surfaces:
	 if i.label=='plane':
            alpha.append(i.alpha)
            beta.append(i.beta)
            gamma.append(i.gamma)
	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    nplane+=1


      for i in self.surfaces:
	 if i.label=='ellipsoid':
            alpha.append(i.alpha)
            beta.append(i.beta)
            gamma.append(i.gamma)
	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_el.append(i.ax)
	    ay_el.append(i.ay)
	    az_el.append(i.az)

	    nellipsoid+=1

      for i in self.surfaces:
	 if i.label=='cylinder':
            alpha.append(i.alpha)
            beta.append(i.beta)
            gamma.append(i.gamma)
	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cyl.append(i.ax)
	    ay_cyl.append(i.ay)

	    ncyl+=1

      for i in self.surfaces:
	 if i.label=='cone':
            alpha.append(i.alpha)
            beta.append(i.beta)
            gamma.append(i.gamma)
	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cone.append(i.ax)
	    ay_cone.append(i.ay)

	    ncone+=1


      return float(tomo3Dfort.gen_obj_voxval(x,y,z,\
           att,\
           nsurf,\
           alpha,beta,gamma,x0,y0,z0,\
           nplane,\
           nellipsoid,\
           ax_el,ay_el,az_el,\
           ncyl,\
           ax_cyl,ay_cyl,\
           ncone,\
           ax_cone,ay_cone))



   def project_to_single_view(self,sino,nview):
      '''Project object onto predefined sino object'''

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0
      
      parms=sino.parms



      att=self.att
      nsurf=self.nsurf

# the following collects info on all the surfaces to give to the projection
# routine. SEE if it makes sense to put this all in some other function.
      xu1=[]
      xu2=[]
      xu3=[]

      yu1=[]
      yu2=[]
      yu3=[]
      
      zu1=[]
      zu2=[]
      zu3=[]

      x0=[]
      y0=[]
      z0=[]

      nplane=1
      
      nellipsoid=1
      ax_el=[1.]
      ay_el=[1.]
      az_el=[1.]

      ncyl=1
      ax_cyl=[1.]
      ay_cyl=[1.]
      
      ncone=1
      ax_cone=[1.]
      ay_cone=[1.]


#put surfaces in the order: planes, ellipsoids, cylinders, and cones
      for i in self.surfaces:
	 if i.label=='plane':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    nplane+=1


      for i in self.surfaces:
	 if i.label=='ellipsoid':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_el.append(i.ax)
	    ay_el.append(i.ay)
	    az_el.append(i.az)

	    nellipsoid+=1

      for i in self.surfaces:
	 if i.label=='cylinder':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cyl.append(i.ax)
	    ay_cyl.append(i.ay)

	    ncyl+=1

      for i in self.surfaces:
	 if i.label=='cone':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cone.append(i.ax)
	    ay_cone.append(i.ay)

	    ncone+=1

      s=s0+nview*ds
      frame_vectors=sino.frame(s,parms)
      sino.mat[nview] = sino.mat[nview]+\
           view_proj3D.gen_shape_proj(\
           sino.indicator[nview],frame_vectors,\
           nu,nv,du,dv,u0,v0,\
           att,\
           nsurf,\
	   xu1,xu2,xu3,yu1,yu2,yu3,zu1,zu2,zu3,\
           x0,y0,z0,\
           nplane,\
           nellipsoid,\
           ax_el,ay_el,az_el,\
           ncyl,\
           ax_cyl,ay_cyl,\
           ncone,\
           ax_cone,ay_cone)


   def project_to(self,sino):
      '''Project object onto predefined sino object'''

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0
      
      parms=sino.parms



      att=self.att
      nsurf=self.nsurf

# the following collects info on all the surfaces to give to the projection
# routine. SEE if it makes sense to put this all in some other function.
      xu1=[]
      xu2=[]
      xu3=[]

      yu1=[]
      yu2=[]
      yu3=[]
      
      zu1=[]
      zu2=[]
      zu3=[]

      x0=[]
      y0=[]
      z0=[]

      nplane=1
      
      nellipsoid=1
      ax_el=[1.]
      ay_el=[1.]
      az_el=[1.]

      ncyl=1
      ax_cyl=[1.]
      ay_cyl=[1.]
      
      ncone=1
      ax_cone=[1.]
      ay_cone=[1.]


#put surfaces in the order: planes, ellipsoids, cylinders, and cones
      for i in self.surfaces:
	 if i.label=='plane':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    nplane+=1


      for i in self.surfaces:
	 if i.label=='ellipsoid':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_el.append(i.ax)
	    ay_el.append(i.ay)
	    az_el.append(i.az)

	    nellipsoid+=1

      for i in self.surfaces:
	 if i.label=='cylinder':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cyl.append(i.ax)
	    ay_cyl.append(i.ay)

	    ncyl+=1

      for i in self.surfaces:
	 if i.label=='cone':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cone.append(i.ax)
	    ay_cone.append(i.ay)

	    ncone+=1

      for ip in range(ns):
         s=s0+ip*ds
         frame_vectors=sino.frame(s,parms)
         sino.mat[ip] = sino.mat[ip]+\
           view_proj3D.gen_shape_proj(\
           sino.indicator[ip],frame_vectors,\
           nu,nv,du,dv,u0,v0,\
           att,\
           nsurf,\
	   xu1,xu2,xu3,yu1,yu2,yu3,zu1,zu2,zu3,\
           x0,y0,z0,\
           nplane,\
           nellipsoid,\
           ax_el,ay_el,az_el,\
           ncyl,\
           ax_cyl,ay_cyl,\
           ncone,\
           ax_cone,ay_cone)


   def lino_project_to(self,lino):
      '''Project object onto lino object'''

      nr=lino.nr
      

      att=self.att
      nsurf=self.nsurf

# the following collects info on all the surfaces to give to the projection
# routine. SEE if it makes sense to put this all in some other function.
      xu1=[]
      xu2=[]
      xu3=[]

      yu1=[]
      yu2=[]
      yu3=[]
      
      zu1=[]
      zu2=[]
      zu3=[]

      x0=[]
      y0=[]
      z0=[]

      nplane=1
      
      nellipsoid=1
      ax_el=[1.]
      ay_el=[1.]
      az_el=[1.]

      ncyl=1
      ax_cyl=[1.]
      ay_cyl=[1.]
      
      ncone=1
      ax_cone=[1.]
      ay_cone=[1.]


#put surfaces in the order: planes, ellipsoids, cylinders, and cones
      for i in self.surfaces:
	 if i.label=='plane':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    nplane+=1


      for i in self.surfaces:
	 if i.label=='ellipsoid':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_el.append(i.ax)
	    ay_el.append(i.ay)
	    az_el.append(i.az)

	    nellipsoid+=1

      for i in self.surfaces:
	 if i.label=='cylinder':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cyl.append(i.ax)
	    ay_cyl.append(i.ay)

	    ncyl+=1

      for i in self.surfaces:
	 if i.label=='cone':
            xu1.append(i.xu1)
            xu2.append(i.xu2)
            xu3.append(i.xu3)

            yu1.append(i.yu1)
            yu2.append(i.yu2)
            yu3.append(i.yu3)

            zu1.append(i.zu1)
            zu2.append(i.zu2)
            zu3.append(i.zu3)

	    x0.append(i.xc+self.x0)
	    y0.append(i.yc+self.y0)
	    z0.append(i.zc+self.z0)

	    ax_cone.append(i.ax)
	    ay_cone.append(i.ay)

	    ncone+=1

      lino.mat = lino.mat + \
         ray_proj3D.gen_shape_proj(\
           lino.indicator,lino.ray_vectors,\
           nr,\
           att,\
           nsurf,\
	   xu1,xu2,xu3,yu1,yu2,yu3,zu1,zu2,zu3,\
           x0,y0,z0,\
           nplane,\
           nellipsoid,\
           ax_el,ay_el,az_el,\
           ncyl,\
           ax_cyl,ay_cyl,\
           ncone,\
           ax_cone,ay_cone)

########################################################
class ellipsoid(shape3D):

   def __init__(self,x0=0.,y0=0.,z0=0.,\
        att=1.0,ax=1.0,ay=1.0,az=1.0,alpha=0.0,beta=0.0):
      self.ax=ax
      self.ay=ay
      self.az=az
      self.alpha=alpha
      self.beta=beta
      shape3D.__init__(self,x0,y0,z0,att)

   def __str__(self):
      lulu = "An ellipsoid!\n"+shape3D.__str__(self)+ '''ellipsoid parms:
ax = %f  
ay = %f
az = %f
alpha = %f
beta = %f\n''' % (self.ax,self.ay,self.az,self.alpha,self.beta)

      return lulu

   def voxvalf(self,x,y,z):
      '''Returns attenuation value if x,y,z is in the ellipsoid'''

      return float(tomo3Dfort.ellipse_voxval(x,y,z,\
      self.att,self.alpha,self.beta,self.x0,self.y0,self.z0,self.ax,self.ay,self.az))

      
   def embed_in(self,image):

      nx=image.nx
      ny=image.ny
      nz=image.nz
      dx=image.dx
      dy=image.dy
      dz=image.dz
      x0_image=image.x0
      y0_image=image.y0
      z0_image=image.z0

      att=self.att
      alpha=self.alpha
      beta=self.beta

      x0=self.x0
      y0=self.y0
      z0=self.z0

      ax=self.ax
      ay=self.ay
      az=self.az


      image.mat = image.mat+\
           tomo3Dfort.ellipsoid_embed(\
           nx,ny,nz,dx,dy,dz,x0_image,y0_image,z0_image,\
           att,\
           alpha,beta,x0,y0,z0,\
           ax,ay,az)


   def voxval(self,x,y,z):
      '''Returns attenuation value if x,y,z is in the ellipsoid'''

      mu=self.att

      ca = cos(self.alpha)
      sa = sin(self.alpha)
      cb = cos(self.beta)
      sb = sin(self.beta)
      ax=self.ax
      ay=self.ay
      az=self.az

      rel_x = x - self.x0
      rel_y = y - self.y0
      rel_z = z - self.z0

      rot1x= cb*rel_x + sb*rel_y
      rot1y=-sb*rel_x + cb*rel_y
      rot1z= rel_z

      rot2x= ca*rot1x + sa*rot1z
      rot2y= rot1y
      rot2z=-sa*rot1x + ca*rot1z


      ellipsoid_lhs=(rot2x/ax)**2+(rot2y/ay)**2+(rot2z/az)**2

      if ellipsoid_lhs<=1.0:
         return mu
      else:
         return 0.0



   def project_to_single_view(self,sino,nview):
      '''Project ellipsoid onto predefined sino object'''

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0

      parms=sino.parms

      ax=self.ax
      ay=self.ay
      az=self.az
      x0=self.x0
      y0=self.y0
      z0=self.z0
      alpha=self.alpha
      beta=self.beta
      att=self.att


      s=s0+nview*ds

      sino.mat[nview] = sino.mat[nview]+\
           view_proj3D.ellipsoid_proj(\
           sino.indicator[nview],frame_factors,\
           nu,nv,du,dv,u0,v0,\
           ax,ay,az,x0,y0,z0,alpha,beta,att)



   def project_to(self,sino):
      '''Project ellipsoid onto predefined sino object'''

      ns=sino.ns
      nu=sino.nu
      nv=sino.nv
      ds=sino.ds
      du=sino.du
      dv=sino.dv
      s0=sino.s0
      u0=sino.u0
      v0=sino.v0

      parms=sino.parms

      ax=self.ax
      ay=self.ay
      az=self.az
      x0=self.x0
      y0=self.y0
      z0=self.z0
      alpha=self.alpha
      beta=self.beta
      att=self.att


      for ip in range(ns):
         s=s0+ip*ds
         frame_vectors=sino.frame(s,parms)
         sino.mat[ip] = sino.mat[ip]+\
           view_proj3D.ellipsoid_proj(\
           sino.indicator[ip],frame_vectors,\
           nu,nv,du,dv,u0,v0,\
           ax,ay,az,x0,y0,z0,alpha,beta,att)



   def lino_project_to(self,lino):
      '''Project ellipsoid onto lino object'''

      nr=lino.nr


      ax=self.ax
      ay=self.ay
      az=self.az
      x0=self.x0
      y0=self.y0
      z0=self.z0
      alpha=self.alpha
      beta=self.beta
      att=self.att


      lino.mat = lino.mat+\
           ray_proj3D.ellipsoid_proj(\
           lino.indicator,lino.ray_vectors,\
           nr,\
           ax,ay,az,x0,y0,z0,alpha,beta,att)






########################################################

########################################################

# to test tomo3D run: ipython -pylab tomo3D.py
if __name__ == "__main__":

   print("hi, let's run through some tests.")


   sys.path.append('phantoms/.')
   import phantoms


## get shepp-logan phantom

   sl=phantoms.generate_3D_shepp_logan_HC()


## put it in the default image array and display
   im1=image3D()
   print('Embedding phantom in an array of size: ',im1.mat.shape)
   nx2=im1.mat.shape[0]/2
   ny2=im1.mat.shape[1]/2
   nz2=im1.mat.shape[2]/2
   sl.embed_in(im1)
   imshow(im1.mat[nx2,:,:])
   axis('equal')
   raw_input('displaying x=0 plane. ')
   imshow(im1.mat[:,ny2,:])
   axis('equal')
   raw_input('displaying y=0 plane. ')
   imshow(im1.mat[:,:,nz2])
   axis('equal')
   raw_input('displaying z=0 plane. ')

## make default sinogram

   sino_an=sinogram3D()
   sino_nn=sinogram3D()
   sino_lin=sinogram3D()
   
## analytic projection

   print('analytic projection to sino_an ...')

   sl.project_to(sino_an)
   figure(1)
   imshow(sino_an.mat[0])

   print('nearest neighbor projection to sino_nn ...')
   im1.project_to(sino_nn)
   figure(2)
   imshow(sino_nn.mat[0])



   print('linear interp. projection to sino_lin ...')
   im1.project_to(sino_lin,alg='proj_linear')
   figure(3)
   imshow(sino_lin.mat[0])




   
