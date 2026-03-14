
import numpy
import random
import sys

sys.path.append('../.')
import tomo3D as t3

pi = numpy.pi


def right_elliptical_cylinder(\
       x0=0.,y0=0.,z0=0.,\
       att=1.,\
       ax=1.,ay=1.,\
       length=1.,\
       alpha=0.,beta=0.,gamma=0.):

   nsurf=3
   half_len=length/2.
   ell_cyl=t3.cylindrical(ax=ax,ay=ay)
   ell_cyl_top=t3.planar(zc=half_len, beta=pi)
   ell_cyl_bot=t3.planar(zc=-half_len)
   rec=t3.gen_object(\
      x0=x0,y0=y0,z0=z0,\
      att=att,\
      nsurf=nsurf,\
      surfaces=[ell_cyl,ell_cyl_top,ell_cyl_bot])
   rec.rotate(alpha=alpha,beta=beta,gamma=gamma)

   return rec

def build_lower_half_defrise_phantom():
   '''Puts together objects for the head phantom.'''

# parameters for setting up component objects.

# attenuation coefficients
   bkgd_att =   1.0
   low_att=-0.7
   high_att = 0.7


   ax_cyl = 0.85
   ay_cyl = 0.85
   len = 0.85
   bkgd_cyl=right_elliptical_cylinder(\
       z0=-len/2.,\
       att=bkgd_att,\
       ax=ax_cyl,ay=ay_cyl,\
       length=len)


# ellipsoid params.
   ax_ell = 0.6
   ay_ell = 0.6
   az_ell = 0.03


   nsurf=2
   ell_surf=t3.ellipsoidal(ax=ax_ell,ay=ay_ell,az=az_ell)
   midplane=t3.planar(beta=pi)
   half_ell=t3.gen_object(\
      att=low_att,\
      nsurf=nsurf,\
      surfaces=[ell_surf,midplane])


   ell_low=[]
   for i in range(3):
      z0 = -0.2-i*0.2
      ell_low.append(t3.ellipsoid(att=low_att,\
               z0=z0,\
	       ax=ax_ell,ay=ay_ell,az=az_ell))


   ell_high=[]
   for i in range(4):
      z0 = -0.1-i*0.2
      ell_high.append(t3.ellipsoid(att=high_att,\
               z0=z0,\
	       ax=ax_ell,ay=ay_ell,az=az_ell))




#put together head phantom
   defrise_phantom=t3.phantom3D()
   defrise_phantom.add_component(bkgd_cyl)
   defrise_phantom.add_component(half_ell)
   for i in ell_low:
      defrise_phantom.add_component(i)
   for i in ell_high:
      defrise_phantom.add_component(i)
   

   return defrise_phantom
   
def build_upper_half_defrise_phantom():
   '''Puts together objects for the head phantom.'''

# parameters for setting up component objects.

# attenuation coefficients
   bkgd_att =   1.0
   low_att=-0.7
   high_att = 0.7


   ax_cyl = 0.85
   ay_cyl = 0.85
   len = 0.85
   bkgd_cyl=right_elliptical_cylinder(\
       z0=len/2.,\
       att=bkgd_att,\
       ax=ax_cyl,ay=ay_cyl,\
       length=len)


# ellipsoid params.
   ax_ell = 0.6
   ay_ell = 0.6
   az_ell = 0.03


   nsurf=2
   ell_surf=t3.ellipsoidal(ax=ax_ell,ay=ay_ell,az=az_ell)
   midplane=t3.planar()
   half_ell=t3.gen_object(\
      att=low_att,\
      nsurf=nsurf,\
      surfaces=[ell_surf,midplane])


   ell_low=[]
   for i in range(3):
      z0 = 0.2+i*0.2
      ell_low.append(t3.ellipsoid(att=low_att,\
               z0=z0,\
	       ax=ax_ell,ay=ay_ell,az=az_ell))


   ell_high=[]
   for i in range(4):
      z0 = 0.1+i*0.2
      ell_high.append(t3.ellipsoid(att=high_att,\
               z0=z0,\
	       ax=ax_ell,ay=ay_ell,az=az_ell))




#put together head phantom
   defrise_phantom=t3.phantom3D()
   defrise_phantom.add_component(bkgd_cyl)
   defrise_phantom.add_component(half_ell)
   for i in ell_low:
      defrise_phantom.add_component(i)
   for i in ell_high:
      defrise_phantom.add_component(i)
   

   return defrise_phantom
   

if __name__=="__main__":
   h=build_lower_half_defrise_phantom()
   print h

   image = t3.image3D(\
           shape=(256,256,256),\
           xlen= 2.0, ylen=2.0,zlen=1.0,\
	   x0 = -1. , y0= -1., z0 = -0.9)



   print "embedding in image  ..."

   h.embed_in(image)

