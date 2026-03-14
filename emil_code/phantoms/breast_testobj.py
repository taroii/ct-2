import sys
#sys.path.append('/home/ireiser/python/')
sys.path.append('../.')
import tomo3D as t3
from numpy import sin
from numpy import cos
pi = t3.pi
#import pycode

def build_breast_phantom(breast_xc = 1.2, breast_yc = 0, breast_zc = -15):
   '''Puts together gen_objects for the breast phantom.
   Current list of objects:
   truncated ellipsoid  (compressed breast)
   rectangular slab     (muscle wall)
   truncated cylinder   (nipple)
   crescent             (ductal tissue)
   a small sphere       (only for test purposes)'''
   
   
   # parameters for setting up component objects.
   # These should be constructed from command line parms later.
   
   # tissue attenuation coefficients (30keV)
   breast_att = 0.245
   muscle_att = 0.3972
   duct_att = 0.3931
   mass_att = 0.40768
   dense_mass_att = 0.40768*1.3
   dense_att = 0.3931
   calc_att = 4.31
   metal_att = 22.6   # for testing
   
   # breast center
 ##   breast_xc=1.
##    breast_xc=1.2
##    breast_yc=0.
##    breast_zc=-15.0
   #breast_zc= 0.0    # this is for Dexela
   #   breast_zc=0.0
   
   # mass centers relative to breast center
   # ellipsoids in pect muscle with varying az
   mass11_xc = -0.65
   mass11_yc = -2.0
   mass11_zc = 0.0
   mass12_xc = -0.65
   mass12_yc = 0.0
   mass12_zc = 0.0
   mass13_xc = -0.65
   mass13_yc = 2.0
   mass13_zc = 0.0
   
   # spheres with varying diameter
   mass21_xc = 1.25
   mass21_yc = -3.75
   mass21_zc = 0.0
   mass22_xc = 1.25
   mass22_yc = -2.25
   mass22_zc = 0.0
   mass23_xc = 1.25
   mass23_yc = -0.75
   mass23_zc = 0.0
   mass24_xc = 1.25
   mass24_yc = 0.75
   mass24_zc = 0.0
   mass25_xc = 1.25
   mass25_yc = 2.25
   mass25_zc = 0.0
   mass26_xc = 1.25
   mass26_yc = 3.75
   mass26_zc = 0.0

   # stacked spheres with varying spacing
   mass31a_xc = 3.0
   mass31a_yc = -2.0
   mass31a_zc = -1.5
   mass32a_xc = 3.0
   mass32a_yc = 0.0
   mass32a_zc = -1.5
   mass33a_xc = 3.0
   mass33a_yc = 2.0
   mass33a_zc = -1.5

   mass31b_xc = 3.0
   mass31b_yc = -2.0
   mass31b_zc = -0.5
   mass32b_xc = 3.0
   mass32b_yc = 0.0
   mass32b_zc = 0.5
   mass33b_xc = 3.0
   mass33b_yc = 2.0
   mass33b_zc = 1.5

   # ellipsoidal mass on pectoralis boundary
   mass4_xc = 0.0
   mass4_yc = -4.0
   mass4_zc = 0.0

   # mass in dense tissue
   mass5_xc = 4.75
   mass5_yc = 0.0
   mass5_zc = 0.0

   # calcification cluster centers relative to breast center
   calc1_xc = 2.5
   calc1_yc = -3.5
   calc1_zc = 0
   calc2_xc = 4.15
   calc2_yc = -3.
   calc2_zc = 0
   calc3_xc = 2.5
   calc3_yc = 3.5
   calc3_zc = 0
   calc4_xc = 4.15
   calc4_yc = 3.
   calc4_zc = 0

   # metal artifact
   metal_xc = 0.25
   metal_yc = 4.5
   metal_zc = 0.0

   # breast is constructed point toward positive x
   breast_length = 6.0
   breast_width = 6.5
   breast_thickness = 5.0
   
   breast_back = 0.
   wall_length = 19.0
   wall_thickness = 5.0
   wall_width = 2.0
   
   crescent_diff = 0.5
   crx_outer_center = 0
   cr_outer_radius = 6
   crx_sphere = -8.
   cr_sphere_rad = 12.0

   nipple_length = 0.5
   nipple_width = 0.5


   # create mass shapes. Need spheres of 6 radii, and 2 ellipsoids. 
   r1 = 0.125
   r2 = 0.25
   r3 = 0.375
   r4 = 0.5
   r5 = 0.625
   r6 = 0.75
   
   n_surfaces=1
   r = r1
   mass_r1= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)
   r = r2
   mass_r2= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)
   r = r3
   mass_r3= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)
   r = r4
   mass_r4= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)
   r = r5
   mass_r5= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)
   r = r6
   mass_r6= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)
   mass_r6= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=r,\
      ay=r,\
      az=r)

   # two elliptical shaped masses
   a_x1 = 0.5
   a_y1 = 0.5
   a_z1 = 0.375
   a_z2 = 0.25
   
   mass_ell1= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=a_x1,\
      ay=a_y1,\
      az=a_z1)
   
   mass_ell2= t3.ellipsoidal(\
      xc=0.,\
      yc=0.,\
      zc=0.,\
      alpha=0.,\
      beta=0.,\
      gamma=0.,\
      ax=a_x1,\
      ay=a_y1,\
      az=a_z2)

   # generate the masses at the appropriate locations
   # 3 masses in pectoralis muscle
   mass11=t3.gen_object(\
     x0=breast_xc+mass11_xc,y0=breast_yc+mass11_yc,z0=breast_zc+mass11_zc,\
     att=dense_mass_att-muscle_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_ell1])
   
   mass12=t3.gen_object(\
     x0=breast_xc+mass12_xc,y0=breast_yc+mass12_yc,z0=breast_zc+mass12_zc,\
     att=dense_mass_att-muscle_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])

   mass13=t3.gen_object(\
     x0=breast_xc+mass13_xc,y0=breast_yc+mass13_yc,z0=breast_zc+mass13_zc,\
     att=dense_mass_att-muscle_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_ell2])

   # 6 masses with varying sizes, in fatty tissue
   mass21=t3.gen_object(\
     x0=breast_xc+mass21_xc,y0=breast_yc+mass21_yc,z0=breast_zc+mass21_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r1])

   mass22=t3.gen_object(\
     x0=breast_xc+mass22_xc,y0=breast_yc+mass22_yc,z0=breast_zc+mass22_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r6])

   mass23=t3.gen_object(\
     x0=breast_xc+mass23_xc,y0=breast_yc+mass23_yc,z0=breast_zc+mass23_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r2])

   mass24=t3.gen_object(\
     x0=breast_xc+mass24_xc,y0=breast_yc+mass24_yc,z0=breast_zc+mass24_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r5])

   mass25=t3.gen_object(\
     x0=breast_xc+mass25_xc,y0=breast_yc+mass25_yc,z0=breast_zc+mass25_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r3])
   
   mass26=t3.gen_object(\
     x0=breast_xc+mass26_xc,y0=breast_yc+mass26_yc,z0=breast_zc+mass26_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   # 3 pairs of stacked masses, in fatty tissue
   mass31a=t3.gen_object(\
     x0=breast_xc+mass31a_xc,y0=breast_yc+mass31a_yc,z0=breast_zc+mass31a_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   mass32a=t3.gen_object(\
     x0=breast_xc+mass32a_xc,y0=breast_yc+mass32a_yc,z0=breast_zc+mass32a_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   mass33a=t3.gen_object(\
     x0=breast_xc+mass33a_xc,y0=breast_yc+mass33a_yc,z0=breast_zc+mass33a_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   mass31b=t3.gen_object(\
     x0=breast_xc+mass31b_xc,y0=breast_yc+mass31b_yc,z0=breast_zc+mass31b_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   mass32b=t3.gen_object(\
     x0=breast_xc+mass32b_xc,y0=breast_yc+mass32b_yc,z0=breast_zc+mass32b_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   mass33b=t3.gen_object(\
     x0=breast_xc+mass33b_xc,y0=breast_yc+mass33b_yc,z0=breast_zc+mass33b_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])
   
   
   # set up mass at pectoralis border
   # need two surfaces
   n_surfaces = 2
   
   sidexl =  t3.planar(xc=0.,yc=0.,zc=0.,\
                       alpha=0.,beta=pi/2.,gamma=0.)
   
   mass4A=t3.gen_object(\
     x0=breast_xc+mass4_xc,y0=breast_yc+mass4_yc,z0=breast_zc+mass4_zc,\
     att=mass_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4,sidexl])
   
   sidexl =  t3.planar(xc=0.,yc=0.,zc=0.,\
                       alpha=pi,beta=pi/2.,gamma=0.)
   
   mass4B=t3.gen_object(\
     x0=breast_xc+mass4_xc,y0=breast_yc+mass4_yc,z0=breast_zc+mass4_zc,\
     att=mass_att-muscle_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4,sidexl])
   
# mass 5 in ductal tissue
   n_surfaces = 1
   mass5=t3.gen_object(\
     x0=breast_xc+mass5_xc,y0=breast_yc+mass5_yc,z0=breast_zc+mass5_zc,\
     att=mass_att-duct_att,\
     nsurf=n_surfaces,\
     surfaces=[mass_r4])

# set up calc clusters
   n_surfaces = 1
   calc_rad = 0.015
   calc_rad_small = 0.0075
   #calc_rad = 0.03 #use this to see calcs
   #calc_rad_small = 0.02 #use this to see calcs
   calc_sphere = t3.ellipsoidal(ax=calc_rad,ay=calc_rad,az=calc_rad)
   calc_sphere_small = t3.ellipsoidal(ax=calc_rad,ay=calc_rad_small,az=calc_rad)
   calc_individ=t3.gen_object(\
     x0= breast_xc+calc1_xc,y0=breast_yc+calc1_yc,z0=breast_zc+calc1_zc,\
     att=calc_att-breast_att,\
     nsurf=n_surfaces,\
     surfaces=[calc_sphere])

   cc_rad = 0.3
   cluster1 = []
   cluster1 = []
   offs = 20./180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=-cc_rad
      cluster1.append(\
               t3.gen_object(\
               x0=breast_xc+calc1_xc+x1,\
               y0=breast_yc+calc1_yc+y1,\
               z0=breast_zc+calc1_zc+z1,\
               att=calc_att-breast_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere]))
   offs = (20.+60.)/180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=cc_rad
      cluster1.append(\
               t3.gen_object(\
               x0=breast_xc+calc1_xc+x1,\
               y0=breast_yc+calc1_yc+y1,\
               z0=breast_zc+calc1_zc+z1,\
               att=calc_att-breast_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere]))
   cluster2 = []
   offs = 20./180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=-cc_rad
      cluster2.append(\
               t3.gen_object(\
               x0=breast_xc+calc2_xc+x1,\
               y0=breast_yc+calc2_yc+y1,\
               z0=breast_zc+calc2_zc+z1,\
               att=calc_att-dense_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere]))
   offs = (20.+60.)/180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=cc_rad
      cluster2.append(\
               t3.gen_object(\
               x0=breast_xc+calc2_xc+x1,\
               y0=breast_yc+calc2_yc+y1,\
               z0=breast_zc+calc2_zc+z1,\
               att=calc_att-dense_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere]))
   cluster3 = []
   offs = 20./180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=-cc_rad
      cluster3.append(\
               t3.gen_object(\
               x0=breast_xc+calc3_xc+x1,\
               y0=breast_yc+calc3_yc+y1,\
               z0=breast_zc+calc3_zc+z1,\
               att=calc_att-breast_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere_small]))
   offs = (20.+60.)/180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=cc_rad
      cluster3.append(\
               t3.gen_object(\
               x0=breast_xc+calc3_xc+x1,\
               y0=breast_yc+calc3_yc+y1,\
               z0=breast_zc+calc3_zc+z1,\
               att=calc_att-breast_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere_small]))
   cluster4 = []
   offs = 20./180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=-cc_rad
      cluster4.append(\
               t3.gen_object(\
               x0=breast_xc+calc4_xc+x1,\
               y0=breast_yc+calc4_yc+y1,\
               z0=breast_zc+calc4_zc+z1,\
               att=calc_att-dense_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere_small]))
   offs = (20.+60.)/180.*pi
   for i in range(3):
      c_angle = offs+i*120./180.*pi
      print c_angle
      x1=cc_rad*sin(c_angle)
      y1=cc_rad*cos(c_angle)
      z1=cc_rad
      cluster4.append(\
               t3.gen_object(\
               x0=breast_xc+calc4_xc+x1,\
               y0=breast_yc+calc4_yc+y1,\
               z0=breast_zc+calc4_zc+z1,\
               att=calc_att-dense_att,\
               nsurf=n_surfaces,\
               surfaces=[calc_sphere_small]))

# metal artifact
   n_surfaces=3
   metal_length = 0.3
   metal_rad = 0.015
   #metal_rad = 0.03  # use this to visualize in object
   ang1 = 0.
   ang2 = pi/2.
   ang1 = pi/4.
   ang2 = pi/2.
   
   cyl  = t3.cylindrical(xc=0.,yc=0.,zc=0.,\
                         alpha = ang1, beta = ang2, gamma = 0.,\
			 ax=metal_rad,\
			 ay=metal_rad)
   sidexu =  t3.planar(xc = metal_length, yc = metal_length, zc = metal_length,\
                       alpha = ang1 + pi, beta = ang2, gamma = 0.)
   sidexl =  t3.planar(xc=0., yc=0., zc=0.,\
                       alpha = ang1, beta = ang2, gamma =0.)
   
   metal = t3.gen_object(\
      x0 = breast_xc + metal_xc, y0 = breast_yc + metal_yc, z0 = breast_zc + metal_zc,\
      att = metal_att-breast_att,\
      nsurf = n_surfaces,\
      surfaces = [sidexl,cyl,sidexu])


# set up truncated ellipsoid (compressed breast)
   zl=-breast_thickness/2.
   zu=breast_thickness/2.
   xl=breast_back

   n_surfaces=4

   sidezl =  t3.planar(xc=0.,yc=0.,zc=zl,\
                       alpha=0.,beta=0.,gamma=0.)
   sidezu =  t3.planar(xc=0.,yc=0.,zc=zu,\
                       alpha=0.,beta=pi,gamma=0.)
   sidexl =  t3.planar(xc=xl,yc=0.,zc=0.,\
                       alpha=0.,beta=pi/2.,gamma=0.)

   ell_surf =t3.ellipsoidal(xc=0.,\
                            yc=0.,\
			    zc=0.,\
                            alpha=0.,\
			    beta=0.,\
			    gamma=0.,\
			    ax=breast_length,\
			    ay=breast_width,\
			    az=breast_width/2.)


   breast_main = t3.gen_object(\
                    x0=breast_xc,y0=breast_yc,z0=breast_zc,\
                    att=breast_att,\
		    nsurf=n_surfaces,\
		    surfaces=[sidezl,sidezu,sidexl,ell_surf])


# set up rectangular box for chest wall

   xl=breast_back-wall_width
   xu=breast_back
   yl=-wall_length/2.
   yu=wall_length/2.
   zl=-wall_thickness/2.
   zu=wall_thickness/2.

   n_surfaces=6

   sidexl =  t3.planar(xc=xl,yc=0.,zc=0.,\
                   alpha=0.,beta=pi/2.,gamma=0.)
   sidexu =  t3.planar(xc=xu,yc=0.,zc=0.,\
                   alpha=pi,beta=pi/2.,gamma=0.)
   sideyl =  t3.planar(xc=0.,yc=yl,zc=0.,\
                   alpha=pi/2.,beta=pi/2.,gamma=0.)
   sideyu =  t3.planar(xc=0.,yc=yu,zc=0.,\
                   alpha=-pi/2.,beta=pi/2.,gamma=0.)
   sidezl =  t3.planar(xc=0.,yc=0.,zc=zl,\
                   alpha=0.,beta=0.,gamma=0.)
   sidezu =  t3.planar(xc=0.,yc=0.,zc=zu,\
                   alpha=0.,beta=pi,gamma=0.)

   chest_wall= t3.gen_object(\
                   x0=breast_xc,y0=breast_yc,z0=breast_zc,\
                   att=muscle_att,\
		   nsurf=n_surfaces,\
		   surfaces=[sidexl,sidexu,sideyl,sideyu,sidezl,sidezu])

# set up truncated cylinder (nipple)
# actually need TWO objects: cylinder - cylinder,truncated by ell_surf
   n_surfaces=3

   xu= breast_length+nipple_length

   cyl  = t3.cylindrical(xc=0.,yc=0.,zc=0.,\
                         alpha=0.,beta=pi/2.,gamma=0.,\
			 ax=nipple_width,\
			 ay=nipple_width)


   sidexu =  t3.planar(xc=xu,yc=0.,zc=0.,\
                       alpha=pi,beta=pi/2.,gamma=0.)

   sidexl =  t3.planar(xc=0.,yc=0.,zc=0.,\
                       alpha=0.,beta=pi/2.,gamma=0.)
   
   nipple_pos = t3.gen_object(\
      x0=breast_xc,y0=breast_yc,z0=breast_zc,\
      att=breast_att,\
      nsurf=n_surfaces,\
      surfaces=[sidexl,cyl,sidexu])

# second nipple object to eliminate cyl. inside breast

   n_surfaces=3

   
   nipple_neg = t3.gen_object(\
      x0=breast_xc,y0=breast_yc,z0=breast_zc,\
      att=-breast_att,\
      nsurf=n_surfaces,\
      surfaces=[sidexl,cyl,ell_surf])


# set up crescent 
# again we need to objects
   n_surfaces=2



   sidexl =  t3.planar(xc=0.,yc=0.,zc=0.,\
                       alpha=0.,beta=pi/2.,gamma=0.)

   ell_surf_crescent_old =t3.ellipsoidal(\
                      xc=0.,\
                      yc=0.,\
                      zc=0.,\
                      alpha=0.,\
                      beta=0.,\
                      gamma=0.,\
                      ax=breast_length-crescent_diff,\
                      ay=breast_width-crescent_diff,\
	              az=breast_width/2.-crescent_diff)
   
   ell_surf_crescent =t3.ellipsoidal(\
                      xc=0.,\
                      yc=0.,\
                      zc=0.,\
                      alpha=0.,\
                      beta=0.,\
                      gamma=0.,\
                      ax=breast_length-crescent_diff,\
                      ay=breast_width-crescent_diff,\
	              az=breast_width/2.-crescent_diff)

   
   crescent_pos = t3.gen_object(\
      x0=breast_xc,y0=breast_yc,z0=breast_zc,\
      att=duct_att-breast_att,\
      nsurf=n_surfaces,\
      surfaces=[sidexl,ell_surf_crescent])

# set up subtraction object
   n_surfaces=3

   sphere_crescent =t3.ellipsoidal(\
                      xc=crx_sphere,\
                      yc=0.,\
                      zc=0.,\
                      alpha=0.,\
                      beta=0.,\
                      gamma=0.,\
                      ax=cr_sphere_rad,\
                      ay=cr_sphere_rad,\
	              az=cr_sphere_rad)

   
   crescent_neg = t3.gen_object(\
      x0=breast_xc,y0=breast_yc,z0=breast_zc,\
      att=-(duct_att-breast_att),\
      nsurf=n_surfaces,\
      surfaces=[sidexl,ell_surf_crescent,sphere_crescent])


           
#put together breast phantom
   breast_phantom=t3.phantom3D()
   breast_phantom.add_component(breast_main)
   breast_phantom.add_component(nipple_pos)
   breast_phantom.add_component(nipple_neg)
   breast_phantom.add_component(crescent_pos)
   breast_phantom.add_component(crescent_neg)
   breast_phantom.add_component(chest_wall)
   breast_phantom.add_component(mass11)
   breast_phantom.add_component(mass12)
   breast_phantom.add_component(mass13)
   breast_phantom.add_component(mass21)
   breast_phantom.add_component(mass22)
   breast_phantom.add_component(mass23)
   breast_phantom.add_component(mass24)
   breast_phantom.add_component(mass25)
   breast_phantom.add_component(mass26)
   breast_phantom.add_component(mass31a)
   breast_phantom.add_component(mass31b)
   breast_phantom.add_component(mass32a)
   breast_phantom.add_component(mass32b)
   breast_phantom.add_component(mass33a)
   breast_phantom.add_component(mass33b)
   breast_phantom.add_component(mass4A)
   breast_phantom.add_component(mass4B)
   breast_phantom.add_component(mass5)
#   breast_phantom.add_component(calc_individ)
   breast_phantom.add_component(metal)
   for i in cluster1:
      breast_phantom.add_component(i)
      print i
   for i in cluster2:
      breast_phantom.add_component(i)
      print i
   for i in cluster3:
       breast_phantom.add_component(i)
   for i in cluster4:
       breast_phantom.add_component(i)

   return breast_phantom
   

if __name__=="__main__":
   b=build_breast_phantom()
   print b
   
   
   #sino = t3.sinogram3D(config_name='tomosynthesis_sliding_detector',\
   sino = t3.sinogram3D(config_name='tomosynthesis',\
                        parms={"center_to_detector": 20.,\
                               "radius"            : 46.},\
                        shape=(11,850,1800),\
                        #shape=(11,170,300),\
                        slen = pi/2., s0 = pi/4.,\
                        ulen = 8.5, u0 = 0.,\
                        #vlen = 19.0, v0 = -9.5)
                        vlen = 18.0, v0 = -9.)
   
   ###
   #  to read this into matlab:
   #  fid = fopen('file.dat');
   #  proj = fread(fid,ns*nu*nv,'float');
   #  fclose(fid);
   #  proj = reshape(proj,[nv nu ns]);
   ###


   image = t3.image3D(\
      #           nx=450,ny=1000,nz=25,\
      #           xlen= 11.25, ylen=25.,zlen=6.5,\
      #	   x0 = -0.5, y0= -12.5, z0 = -18.0)
      #            nx=450,ny=1000,nz=50,\
      #            xlen= 11.25, ylen=25.,zlen=6.5,\
      #            x0 = -0.5, y0= -12.5, z0 = -18.0)
            #shape=(425,1060,30),\
            shape=(425,1000,300),\
            #            nx=425,ny=1060,nz=120,\
            #xlen= 8.5, ylen=21.2,zlen=6.,\
            xlen= 8.5, ylen=20.0,zlen=6.,\
            x0 = -0.5, y0= -10., z0 = -18.0)

   ### to read the IMAGE into matlab:
   #  fid = fopen('breast_phantomSMALL.dat');
   #  img = fread(fid,450*1000*25,'float');
   #  fclose(fid);
   #  img=permute(reshape(img,[25 1000 450]),[3 2 1]);
   ###
 
#   print "projecting ..."
#   b.project_to(sino)
#   sino.write_file("breast_testobj_sino.dat")
   #noisy = pycode.add_noise(sino)
   #noisy.write_file("breast_testobj_sino_noisy.dat")



   print "embedding in image  ..."
   b.embed_in(image)
   image.write_file("breast_testobj.dat")

