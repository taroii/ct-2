import numpy
import random
import sys

sys.path.append('../.')
import tomo3D as t3

pi = numpy.pi
sqrt = numpy.sqrt
sin = numpy.sin
cos = numpy.cos


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

def make_ear_holes(bone_att=1.8):

   eh=[]
   rad=0.15
   sha=0.2*t3.sqrt(3.)

# mid level
   for i in range(8):
      xc=5.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=5.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2*sha,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(5):
      xc=6.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3*sha,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=5.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2*sha,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(5):
      xc=6.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3*sha,z0=0.,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

# first floor

   for i in range(7):
      xc=5.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2.*sha,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3.*sha,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2.*sha,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3.*sha,z0=sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))



# first floor - lower level

   for i in range(7):
      xc=5.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2.*sha,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3.*sha,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2.*sha,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))
   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3.*sha,z0=-sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

# second floor

   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.4+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2.*sha,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3.*sha,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.4+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2.*sha,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3.*sha,z0=2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

# second floor lower level

   for i in range(7):
      xc=6.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.4+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2.*sha,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3.*sha,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.2+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(6):
      xc=6.4+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2.*sha,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3.*sha,z0=-2.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

# third floor

   for i in range(5):
      xc=6.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2.*sha,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(3):
      xc=7.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3.*sha,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2.*sha,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(3):
      xc=7.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3.*sha,z0=3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

# third floor - lower level

   for i in range(5):
      xc=6.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=0.,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=sha,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=2.*sha,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(3):
      xc=7.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=3.*sha,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(5):
      xc=6.8+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-sha,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(4):
      xc=7.0+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-2.*sha,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))

   for i in range(3):
      xc=7.6+0.4*i
      eh.append(\
         t3.ellipsoid(\
	    x0=xc,y0=-3.*sha,z0=-3.*sha,\
            att= - bone_att,\
            ax=rad,ay=rad,az=rad))




   return eh



def build_head_phantom():
   '''Puts together objects for the head phantom.'''

# parameters for setting up component objects.

# tissue attenuation coefficients
   bone_att =     1.8
   eye_att=       1.06
   brain_att =    1.05
   ventricle_att= 1.045
   spot1_att    = 1.0525
   spot2_att    = 1.0475
   hematoma_att = 1.055

#   bone_att =     1.8
#   eye_att=       1.1
#   brain_att =    1.0
#   ventricle_att= 1.05
#   spot1_att    = 1.1
#   spot2_att    = 0.9
#   hematoma_att = 1.05


# shape params.

# Outer ellipsoid params.
   ax_oe = 9.6
   ay_oe = 12.0
   az_oe = 12.5

# Inner skull ellipsoid
   ax_ie = 9.0
   ay_ie = 11.4
   az_ie = 11.9


# bone objects
   outer_ell=t3.ellipsoid(\
      att=bone_att,\
      ax=ax_oe,ay=ay_oe,az=az_oe)
   inner_ell=t3.ellipsoid(\
      att= - bone_att,\
      ax=ax_ie,ay=ay_ie,az=az_ie)

   sinus_bone1=t3.ellipsoid(\
      x0 = -1.9, y0 = 5.4, z0 = 0.0,\
      att=bone_att,\
      ax = 1.165, ay = 0.4060, az = 3.0,\
      alpha = 0.,\
      beta = -pi/4.)
   sinus_bone1_hole=t3.ellipsoid(\
      x0 = -1.9, y0 = 5.4, z0 = 0.0,\
      att=-brain_att,\
      ax = 1.165, ay = 0.4060, az = 3.0,\
      alpha = 0.,\
      beta = -pi/4.)

   sinus_bone2=t3.ellipsoid(\
      x0 = 1.9, y0 = 5.4, z0 = 0.0,\
      att=bone_att,\
      ax = 1.165, ay = 0.4060, az = 3.0,\
      alpha = 0.,\
      beta = pi/4.)
   sinus_bone2_hole=t3.ellipsoid(\
      x0 = 1.9, y0 = 5.4, z0 = 0.0,\
      att=-brain_att,\
      ax = 1.165, ay = 0.4060, az = 3.0,\
      alpha = 0.,\
      beta = pi/4.)

   x0_sb3 = 0.
   y0_sb3 = 3.6
   z0_sb3 = 0.
   ax_sb3 = 4.0
   ay_sb3 = 1.2
   len = 0.483
   sinus_bone3=right_elliptical_cylinder(\
       x0=x0_sb3,y0=y0_sb3,z0=z0_sb3,\
       att=bone_att,\
       ax=ax_sb3,ay=ay_sb3,\
       length=len,\
       alpha=pi/2.,beta=pi/3.)
   sinus_bone3_hole=right_elliptical_cylinder(\
       x0=x0_sb3,y0=y0_sb3,z0=z0_sb3,\
       att=-brain_att,\
       ax=ax_sb3,ay=ay_sb3,\
       length=len,\
       alpha=pi/2.,beta=pi/3.)


   x0_sb4 = 0.
   y0_sb4 = 9.6
   z0_sb4 = 0.
   ax_sb4 = 2.0
   ay_sb4 = 0.525561
   len = 0.4
   sinus_bone4=right_elliptical_cylinder(\
       x0=x0_sb4,y0=y0_sb4,z0=z0_sb4,\
       att=bone_att,\
       ax=ax_sb4,ay=ay_sb4,\
       length=len,\
       alpha=0.,beta=pi/2.,gamma=-pi/6.)


   x0_sb5 = -4.3
   y0_sb5 = 6.8
   z0_sb5 = -1.0
   ax_sb5 = 1.8
   ay_sb5 = 0.24
   len = 4.0
   sinus_bone5=right_elliptical_cylinder(\
       x0=x0_sb5,y0=y0_sb5,z0=z0_sb5,\
       att=bone_att,\
       ax=ax_sb5,ay=ay_sb5,\
       length=len,\
       alpha=0.,beta=0.,gamma=-pi/6.)
   sinus_bone5_hole=right_elliptical_cylinder(\
       x0=x0_sb5,y0=y0_sb5,z0=z0_sb5,\
       att=-brain_att,\
       ax=ax_sb5,ay=ay_sb5,\
       length=len,\
       alpha=0.,beta=0.,gamma=-pi/6.)


   x0_sb6 = 4.3
   y0_sb6 = 6.8
   z0_sb6 = -1.0
   ax_sb6 = 1.8
   ay_sb6 = 0.24
   len = 4.0
   sinus_bone6=right_elliptical_cylinder(\
       x0=x0_sb6,y0=y0_sb6,z0=z0_sb6,\
       att=bone_att,\
       ax=ax_sb6,ay=ay_sb6,\
       length=len,\
       alpha=0.,beta=0.,gamma=pi/6.)
   sinus_bone6_hole=right_elliptical_cylinder(\
       x0=x0_sb6,y0=y0_sb6,z0=z0_sb6,\
       att=-brain_att,\
       ax=ax_sb6,ay=ay_sb6,\
       length=len,\
       alpha=0.,beta=0.,gamma=pi/6.)


   nsurf=2
   ear_ell=t3.ellipsoidal(ax=4.2,ay=1.8,az=1.8)
   skull_ell=t3.ellipsoidal(xc=-9.1,ax=ax_ie,ay=ay_ie,az=az_ie)
   ear_right=t3.gen_object(\
      x0=9.1,y0=0.,z0=0.,\
      att=bone_att,\
      nsurf=nsurf,\
      surfaces=[ear_ell,skull_ell])

   dot1 = t3.ellipsoid(\
      x0=-8.0, y0=0.0, z0=0.0,\
      att=bone_att,\
      ax=0.1,ay=0.1,az=0.1)
   dot2 = t3.ellipsoid(\
      x0=-8.5, y0=0.5, z0=0.0,\
      att=bone_att,\
      ax=0.05,ay=0.05,az=0.05)
   dot3 = t3.ellipsoid(\
      x0=-8.5, y0=-0.5, z0=0.0,\
      att=bone_att,\
      ax=0.02 ,ay=0.02 ,az=0.02 )
   dot1_hole = t3.ellipsoid(\
      x0=-8.0, y0=0.0, z0=0.0,\
      att=-brain_att,\
      ax=0.1,ay=0.1,az=0.1)
   dot2_hole = t3.ellipsoid(\
      x0=-8.5, y0=0.5, z0=0.0,\
      att=-brain_att,\
      ax=0.05,ay=0.05,az=0.05)
   dot3_hole = t3.ellipsoid(\
      x0=-8.5, y0=-0.5, z0=0.0,\
      att=-brain_att,\
      ax=0.02 ,ay=0.02 ,az=0.02 )


   ear_holes=make_ear_holes()

   ear_right_hole=t3.gen_object(\
      x0=9.1,y0=0.,z0=0.,\
      att=-brain_att,\
      nsurf=nsurf,\
      surfaces=[ear_ell,skull_ell])



   nsurf=2
   skull_ell=t3.ellipsoidal(ax=ax_ie,ay=ay_ie,az=az_ie)
#   pro_cone=t3.conical(xc=0.,yc=-11.15,zc=0.2,\
   pro_cone=t3.conical(xc=0.,yc=-10.0,zc=0.2,\
                       alpha=pi/2.,beta=-pi/2.,\
                       ax=0.5,ay=0.2)
   protuberance=\
      t3.gen_object(\
         att=bone_att,\
         nsurf=nsurf,\
         surfaces=[pro_cone,skull_ell])

   protuberance_hole=\
      t3.gen_object(\
         att=-brain_att,\
         nsurf=nsurf,\
         surfaces=[pro_cone,skull_ell])





#put together bone map
   bone_map=t3.tissue_map(physical_material=False)
   bone_map.add_component(outer_ell)
   bone_map.add_component(inner_ell)
   bone_map.add_component(sinus_bone1)
   bone_map.add_component(sinus_bone2)
   bone_map.add_component(sinus_bone3)
   bone_map.add_component(sinus_bone4)
   bone_map.add_component(sinus_bone5)
   bone_map.add_component(sinus_bone6)

   bone_map.add_component(protuberance)

# EAR STRUCTURE IS HERE
   bone_map.add_component(ear_right)
   for i in ear_holes:
      bone_map.add_component(i)
   bone_map.add_component(dot1)
   bone_map.add_component(dot2)
   bone_map.add_component(dot3)

# eyeballs
   eye_rad=2.0

   left_eye=t3.ellipsoid(\
      x0=-4.7, y0=4.3, z0=0.872,\
      att=eye_att,\
      ax=eye_rad,ay=eye_rad,az=eye_rad)
   left_eye_hole=t3.ellipsoid(\
      x0=-4.7, y0=4.3, z0=0.872,\
      att=-brain_att,\
      ax=eye_rad,ay=eye_rad,az=eye_rad)

   right_eye=t3.ellipsoid(\
      x0=4.7, y0=4.3, z0=0.872,\
      att=eye_att,\
      ax=eye_rad,ay=eye_rad,az=eye_rad)
   right_eye_hole=t3.ellipsoid(\
      x0=4.7, y0=4.3, z0=0.872,\
      att=-brain_att,\
      ax=eye_rad,ay=eye_rad,az=eye_rad)

   eye_map=t3.tissue_map(physical_material=False)
   eye_map.add_component(left_eye)
   eye_map.add_component(right_eye)

# ventricle

   ventricle=t3.ellipsoid(\
      y0=-3.6,\
      att=ventricle_att,\
      ax=1.8, ay=3.6, az=3.6)
   ventricle_hole=t3.ellipsoid(\
      y0=-3.6,\
      att=-brain_att,\
      ax=1.8, ay=3.6, az=3.6)


   ventricle_map=t3.tissue_map(physical_material=False) 
   ventricle_map.add_component(ventricle)

# spot 1  (object 3 in FORBILD)

   spot1=t3.ellipsoid(\
      x0=-1.08, y0=-9.,\
      att=spot1_att,\
      ax=0.4, ay=0.4, az=0.4)
   spot1_hole=t3.ellipsoid(\
      x0=-1.08, y0=-9.,\
      att=-brain_att,\
      ax=0.4, ay=0.4, az=0.4)


   spot1_map=t3.tissue_map(physical_material=False) 
   spot1_map.add_component(spot1)


# spot 2  (object 4 in FORBILD)

   spot2=t3.ellipsoid(\
      x0=1.08, y0=-9.,\
      att=spot2_att,\
      ax=0.4, ay=0.4, az=0.4)
   spot2_hole=t3.ellipsoid(\
      x0=1.08, y0=-9.,\
      att=-brain_att,\
      ax=0.4, ay=0.4, az=0.4)


   spot2_map=t3.tissue_map(physical_material=False) 
   spot2_map.add_component(spot2)

# hematoma
   hematoma=t3.ellipsoid(\
      x0=6.393945, y0=-6.393945,\
      att=hematoma_att,\
      ax=1.2, ay=0.42, az=1.4,\
      beta=pi*58.1/180.)
   hematoma_hole=t3.ellipsoid(\
      x0=6.393945, y0=-6.393945,\
      att=-brain_att,\
      ax=1.2, ay=0.42, az=1.4,\
      beta=pi*58.1/180.)

   hematoma_map=t3.tissue_map(physical_material=False) 
   hematoma_map.add_component(hematoma)



#brain map
   brain=t3.ellipsoid(\
      att= brain_att,\
      ax=ax_ie,ay=ay_ie,az=az_ie)

   sinus_hole=t3.ellipsoid(\
      x0=0.,y0=8.4,z0=0.,\
      att= -brain_att,\
      ax=1.8,ay=3.0,az=3.0)


   brain_map=t3.tissue_map(physical_material=False)
   brain_map.add_component(brain)
   brain_map.add_component(ear_right_hole)
   brain_map.add_component(dot1_hole)
   brain_map.add_component(dot2_hole)
   brain_map.add_component(dot3_hole)

   brain_map.add_component(protuberance_hole)
   brain_map.add_component(left_eye_hole)
   brain_map.add_component(right_eye_hole)
   brain_map.add_component(sinus_bone1_hole)
   brain_map.add_component(sinus_bone2_hole)
   brain_map.add_component(sinus_bone3_hole)
   brain_map.add_component(sinus_bone5_hole)
   brain_map.add_component(sinus_bone6_hole)
   brain_map.add_component(sinus_hole)
   brain_map.add_component(ventricle_hole)
   brain_map.add_component(spot1_hole)
   brain_map.add_component(spot2_hole)
   brain_map.add_component(hematoma_hole)

   



#put together head phantom
   head_phantom=t3.phantom3D()
   head_phantom.add_component(bone_map)
   head_phantom.add_component(eye_map)
   head_phantom.add_component(ventricle_map)
   head_phantom.add_component(spot1_map)
   head_phantom.add_component(spot2_map)
   head_phantom.add_component(hematoma_map)
   head_phantom.add_component(brain_map)

   return head_phantom
   

if __name__=="__main__":
   h=build_head_phantom()

   print h

   image = t3.image3D(\
        shape=(512,512,32),\
        xlen= 26.0, ylen=26.0,zlen=26.0/16.0,\
        x0 = -13., y0= -13., z0 = -13./16.0)



   print "embedding in image  ..."

   h.embed_in(image)

