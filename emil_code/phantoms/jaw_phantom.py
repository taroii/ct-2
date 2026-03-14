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


def box(\
       x0=0.,y0=0.,z0=0.,\
       att=1.,\
       lenx=1.,leny=1.,lenz=1.,\
       alpha=0.,beta=0.,gamma=0.):

   nsurf=6
   half_lenx=lenx/2.
   half_leny=leny/2.
   half_lenz=lenz/2.
   xside1=t3.planar(xc=half_lenx, beta=pi/2, alpha=pi)
   xside2=t3.planar(xc=-half_lenx, beta=pi/2)
   yside1=t3.planar(yc=half_leny, alpha=-pi/2, beta=pi/2.)
   yside2=t3.planar(yc=-half_leny, alpha=pi/2, beta=pi/2.)
   zside1=t3.planar(zc=half_lenz, beta=pi)
   zside2=t3.planar(zc=-half_lenz)
   box=t3.gen_object(\
      x0=x0,y0=y0,z0=z0,\
      att=att,\
      nsurf=nsurf,\
      surfaces=[xside1,\
                xside2,\
                yside1,\
                yside2,\
                zside1,\
                zside2])
   box.rotate(alpha=alpha,beta=beta,gamma=gamma)

   return box





def build_jaw_phantom():
   '''Puts together objects for the jaw phantom.'''

# parameters for setting up component objects.

# tissue attenuation coefficients
   body_att =       1.0
   bone_att =       2.0
   marrow_att =     1.1
   teeth_att =      2.2
   gold_crown_att= 19.3
   cancer_att=      1.01


#marrow map
   spine_marrow = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 5.0,\
      att = marrow_att, \
      ax = 0.5, ay = 0.5, length = 10.0)

   marrow_map=t3.tissue_map()
   marrow_map.add_component(spine_marrow)


#cancer map
   tumor= t3.ellipsoid(\
      x0 = -1.0, y0 = 7.0, z0 = 4.0,\
      att = cancer_att, \
      ax = 0.5, ay = 0.5, az = 0.5)

   cancer_map=t3.tissue_map()
   cancer_map.add_component(tumor)


#body map
   neck1=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 0.0, z0 = 1.5,\
      att = body_att, \
      ax = 6.0, ay = 6.0, length = 3.0)

   neck2=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 0.0, z0 = 6.5,\
      att = body_att, \
      ax = 6.0, ay = 6.0, length = 7.0)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = -pi/2, beta=pi/2)
   neck2.add_surface(ell_cyl_bisector)

   throat = box(\
      x0 = 0.0, y0 = 4.0, z0 = 6.0,\
      att = body_att, \
      lenx = 5.0, leny = 2.0, lenz = 5.0)

   tongue1=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 4.25,\
      att = body_att, \
      ax = 2.0, ay = 2.0, length = 1.5)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   tongue1.add_surface(ell_cyl_bisector)

   tongue2 = box(\
      x0 = 0.0, y0 = 6.25, z0 = 4.25,\
      att = body_att, \
      lenx = 4.0, leny = 2.5, lenz = 1.5)

   head1=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 6.0, z0 = 6.5,\
      att = body_att, \
      ax = 6.0, ay = 6.0, length = 7.0)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   head1.add_surface(ell_cyl_bisector)
   

   head2=box(\
      x0 = 0.0, y0 = 3.0, z0 = 6.5,\
      att = body_att, \
      lenx = 12.0, leny = 6.0, lenz = 7.0)

   spine_hole = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 5.0,\
      att = -body_att, \
      ax = 2.5, ay = 1.5, length = 10.0)

   trachea_hole = right_elliptical_cylinder(\
      x0 = 0.0, y0 = 3.75, z0 = 2.5,\
      att = -body_att, \
      ax = 0.75, ay = 0.75, length = 5.0)

# radius is sligthly larger than phantom spec
   tumor_hole= t3.ellipsoid(\
      x0 = -1.0, y0 = 7.0, z0 = 4.0,\
      att = -body_att, \
      ax = 0.5, ay = 0.5, az = 0.5)

   jaw_hole1= right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 6.0,\
      att = -body_att, \
      ax = 3.55, ay = 3.55, length = 5.0)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   jaw_hole1.add_surface(ell_cyl_bisector)

   jaw_hole2= box(\
      x0 = 0.0, y0 = 5.25, z0 = 6.0,\
      att = -body_att, \
      lenx = 7.1, leny = 4.5, lenz = 5.0)

   jaw_hole3= right_elliptical_cylinder(\
      x0 = 0.0, y0 = 5.0, z0 = 7.5,\
      att = -body_att, \
      ax = 2.0, ay = 2.0, length = 5.0)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = -pi/2, beta=pi/2)
   jaw_hole3.add_surface(ell_cyl_bisector)

   sinus_hole1=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 9.25,\
      att = -body_att, \
      ax = 2.0, ay = 2.0, length = 1.5)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   sinus_hole1.add_surface(ell_cyl_bisector)
   sinus_hole2 = box(\
      x0 = 0.0, y0 = 6.25, z0 = 9.25,\
      att = -body_att, \
      lenx = 4.0, leny = 2.5, lenz = 1.5)

   chord_hole1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 1.5,\
      att = -body_att,\
      ax = 0.5, ay = 0.5, length = 1.0)
   disc1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 1.5,\
      att = body_att,\
      ax = 2.5, ay = 1.5, length = 1.0)

   chord_hole2 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 4.5,\
      att = -body_att,\
      ax = 0.5, ay = 0.5, length = 1.0)
   disc2 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 4.5,\
      att = body_att,\
      ax = 2.5, ay = 1.5, length = 1.0)

   chord_hole3 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 7.5,\
      att = -body_att,\
      ax = 0.5, ay = 0.5, length = 1.0)
   disc3 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 7.5,\
      att = body_att,\
      ax = 2.5, ay = 1.5, length = 1.0)


   body_map=t3.tissue_map()
   body_map.add_component(neck1)
   body_map.add_component(neck2)
   body_map.add_component(head1)
   body_map.add_component(head2)
   body_map.add_component(throat)
   body_map.add_component(tongue1)
   body_map.add_component(tongue2)
   body_map.add_component(disc1)
   body_map.add_component(disc2)
   body_map.add_component(disc3)

   body_map.add_component(trachea_hole)
   body_map.add_component(spine_hole)
   body_map.add_component(jaw_hole1)
   body_map.add_component(jaw_hole2)
   body_map.add_component(jaw_hole3)
   body_map.add_component(sinus_hole1)
   body_map.add_component(sinus_hole2)
   body_map.add_component(chord_hole1)
   body_map.add_component(chord_hole2)
   body_map.add_component(chord_hole3)
   body_map.add_component(tumor_hole)

# bone map

   spine1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 5.0,\
      att = bone_att, \
      ax = 2.5, ay = 1.5, length = 10.0)


   lower_jaw1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 4.25,\
      att = bone_att, \
      ax = 3.55, ay = 3.55, length = 1.5)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   lower_jaw1.add_surface(ell_cyl_bisector)

   lower_jaw2 = box(\
      x0 = 0.0, y0 = 6.25, z0 = 4.25,\
      att = bone_att, \
      lenx = 7.1, leny = 2.5, lenz = 1.5)

   back_jaw1 = box(\
      x0 = 3.025, y0 = 4.0, z0 = 6.0,\
      att = bone_att, \
      lenx = 1.05, leny = 2.0, lenz = 5.0)
   back_jaw2 = box(\
      x0 = -3.025, y0 = 4.0, z0 = 6.0,\
      att = bone_att, \
      lenx = 1.05, leny = 2.0, lenz = 5.0)


   upper_jaw1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 7.75,\
      att = bone_att, \
      ax = 3.55, ay = 3.55, length = 1.5)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   upper_jaw1.add_surface(ell_cyl_bisector)

   upper_jaw2 = box(\
      x0 = 0.0, y0 = 6.25, z0 = 7.75,\
      att = bone_att, \
      lenx = 7.1, leny = 2.5, lenz = 1.5)

#   throat = box(\
#      x0 = 0.0, y0 = 4.0, z0 = 6.0,\
#      att = body_att, \
#      lenx = 5.0, leny = 2.0, lenz = 5.0)

   chord_hole1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 0.5,\
      att = -bone_att,\
      ax = 0.5, ay = 0.5, length = 1.0)
   disc_hole1 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 1.5,\
      att = -bone_att,\
      ax = 2.5, ay = 1.5, length = 1.0)

   chord_hole2 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 3.0,\
      att = -bone_att,\
      ax = 0.5, ay = 0.5, length = 2.0)
   disc_hole2 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 4.5,\
      att = -bone_att,\
      ax = 2.5, ay = 1.5, length = 1.0)

   chord_hole3 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 6.0,\
      att = -bone_att,\
      ax = 0.5, ay = 0.5, length = 2.0)
   disc_hole3 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.0, z0 = 7.5,\
      att = -bone_att,\
      ax = 2.5, ay = 1.5, length = 1.0)

   chord_hole4 = right_elliptical_cylinder(\
      x0 = 0.0, y0 = -1.5, z0 = 9.0,\
      att = -bone_att,\
      ax = 0.5, ay = 0.5, length = 2.0)

   tongue_hole1=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 4.25,\
      att = -bone_att, \
      ax = 2.0, ay = 2.0, length = 1.5)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   tongue_hole1.add_surface(ell_cyl_bisector)

   tongue_hole2 = box(\
      x0 = 0.0, y0 = 6.25, z0 = 4.25,\
      att = -bone_att, \
      lenx = 4.0, leny = 2.5, lenz = 1.5)

   roof_hole1=right_elliptical_cylinder(\
      x0 = 0.0, y0 = 7.5, z0 = 7.5,\
      att = -bone_att, \
      ax = 2.0, ay = 2.0, length = 1.0)
   ell_cyl_bisector = t3.planar(yc=0.0, alpha = pi/2, beta=pi/2)
   roof_hole1.add_surface(ell_cyl_bisector)
   roof_hole2 = box(\
      x0 = 0.0, y0 = 6.25, z0 = 7.5,\
      att = -bone_att, \
      lenx = 4.0, leny = 2.5, lenz = 1.0)

   bone_map=t3.tissue_map()
   bone_map.add_component(spine1)
   bone_map.add_component(lower_jaw1)
   bone_map.add_component(lower_jaw2)
   bone_map.add_component(back_jaw1)
   bone_map.add_component(back_jaw2)
   bone_map.add_component(upper_jaw1)
   bone_map.add_component(upper_jaw2)
   bone_map.add_component(tongue_hole1)
   bone_map.add_component(tongue_hole2)
   bone_map.add_component(roof_hole1)
   bone_map.add_component(roof_hole2)
   bone_map.add_component(disc_hole1)
   bone_map.add_component(disc_hole2)
   bone_map.add_component(disc_hole3)
   bone_map.add_component(chord_hole1)
   bone_map.add_component(chord_hole2)
   bone_map.add_component(chord_hole3)
   bone_map.add_component(chord_hole4)

# crown map
   crown1=right_elliptical_cylinder(\
      x0 = 3.0, y0 = 6.5, z0 = 5.475,\
      att = gold_crown_att, \
      ax = 0.5, ay = 0.5, length = 0.95)
   crown2=right_elliptical_cylinder(\
      x0 = 3.0, y0 = 7.5, z0 = 5.475,\
      att = gold_crown_att, \
      ax = 0.5, ay = 0.5, length = 0.95)
   crown3=right_elliptical_cylinder(\
      x0 = -3.0, y0 = 6.5, z0 = 5.475,\
      att = gold_crown_att, \
      ax = 0.5, ay = 0.5, length = 0.95)

   crown_map=t3.tissue_map()
   crown_map.add_component(crown1)
   crown_map.add_component(crown2)
   crown_map.add_component(crown3)




# teeth map
   teeth_loc=[\
      (3.0, 5.5, 5.475),\
#      (3.0, 6.5, 5.475),\
#      (3.0, 7.5, 5.475),\
      (2.9, 8.4, 5.475),\
      (2.5, 9.2, 5.475),\
      (1.9, 9.9, 5.475),\
      (1.2, 10.4, 5.475),\
      (0.4, 10.6, 5.475),\
      (3.0, 5.5, 6.525),\
      (3.0, 6.5, 6.525),\
      (3.0, 7.5, 6.525),\
      (2.9, 8.4, 6.525),\
      (2.5, 9.2, 6.525),\
      (1.9, 9.9, 6.525),\
      (1.2, 10.4, 6.525),\
      (0.4, 10.6, 6.525),\
      (-3.0, 5.5, 5.475),\
      (-3.0, 7.5, 5.475),\
      (-2.9, 8.4, 5.475),\
      (-2.5, 9.2, 5.475),\
      (-1.9, 9.9, 5.475),\
      (-1.2, 10.4, 5.475),\
      (-0.4, 10.6, 5.475),\
      (-3.0, 5.5, 6.525),\
      (-3.0, 6.5, 6.525),\
#      (-3.0, 7.5, 6.525),\
      (-2.9, 8.4, 6.525),\
      (-2.5, 9.2, 6.525),\
      (-1.9, 9.9, 6.525),\
      (-1.2, 10.4, 6.525),\
      (-0.4, 10.6, 6.525)  ]

   teeth_map=t3.tissue_map()
   for tloc in teeth_loc:
      tooth_rad = 0.4
# all teeth at abs(x0) =3 are molars
      if abs(tloc[0]) > 2.99:
         tooth_rad = 0.5

      tooth=right_elliptical_cylinder(\
         x0 = tloc[0], y0 = tloc[1], z0 = tloc[2],\
         att = teeth_att, \
         ax = tooth_rad, ay = tooth_rad, length = 0.95)

      teeth_map.add_component(tooth)



#put together jaw phantom
   jaw_phantom=t3.phantom3D()
   jaw_phantom.add_component(crown_map)
   jaw_phantom.add_component(teeth_map)
   jaw_phantom.add_component(body_map)
   jaw_phantom.add_component(bone_map)
   jaw_phantom.add_component(marrow_map)
   jaw_phantom.add_component(cancer_map)

   return jaw_phantom
   

if __name__=="__main__":
   jaw=build_jaw_phantom()

   print jaw

   image = t3.image3D(\
        shape=(200,200,100),\
        xlen= 20.0, ylen=20.0,zlen=12.0,\
        x0 = -10., y0= -6., z0 = 0.)



   print "embedding in image  ..."

   jaw.embed_in(image)

