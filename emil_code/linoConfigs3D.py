import numpy

sin=numpy.sin
cos=numpy.cos
sqrt=numpy.sqrt
pi=numpy.pi
arctan2 = numpy.arctan2


def config_vectors_ccbn(s,parms):


   srad=parms["source_radius"]
   rlen=parms["ray_len"]
  
   slambda= s[0]
   sphi1= s[1]
   sphi2= s[2]
   
   x_source=srad * cos(slambda)
   y_source=srad *sin(slambda)
   z_source= 0.

   rayx = sin(sphi2)*cos(sphi1)
   rayy = sin(sphi2)*sin(sphi1)
   rayz = cos(sphi2)

   x_bin=x_source + rlen*rayx
   y_bin=y_source + rlen*rayy
   z_bin=z_source + rlen*rayz

   return\
      x_source,y_source,z_source,\
      x_bin,y_bin,z_bin


def config_vectors_spherical(s,parms):


   srad=parms["source_radius"]
   brad=parms["det_radius"]
  
   stheta = s[0]
   sphi = s[1]
   btheta = s[2]
   bphi = s[3] 
   
   x_source=srad * cos(sphi)*sin(stheta)
   y_source=srad *sin(sphi)*sin(stheta)
   z_source=srad *cos(stheta)

   x_bin=brad * cos(bphi)*sin(btheta)
   y_bin=brad *sin(bphi)*sin(btheta)
   z_bin=brad *cos(btheta)

   return\
      x_source,y_source,z_source,\
      x_bin,y_bin,z_bin


def config_vectors_helical_conebeam_64(s,parms):


   srad=parms["source_radius"]
   sddist=parms["source_to_detector"]
   pitch = parms["pitch"]
   lshift = parms["lshift"]
  
   slam = s[0] +lshift*pi/180.
   bbeta = s[1]
   bv= s[2]
   
   x_source=srad * cos(slam)
   y_source=srad *sin(slam)
   z_source=pitch *slam/(2.*pi)

   x_bin=srad * cos(slam) - sddist*cos(slam-bbeta)
   y_bin=srad *sin(slam) -  sddist*sin(slam-bbeta)
   z_bin=pitch *slam/(2.*pi) + bv

   return\
      x_source,y_source,z_source,\
      x_bin,y_bin,z_bin





def config_vectors_spherical_flatpanel(s,parms):


   srad=parms["source_radius"]
   sddist=parms["source_detector"]
  
   stheta = s[0]
   sphi = s[1]
   u = s[2]
   v = s[3] 
   
   x_source=srad * cos(sphi)*sin(stheta)
   y_source=srad *sin(sphi)*sin(stheta)
   z_source=srad *cos(stheta)

   ewx =  cos(sphi)*sin(stheta)
   ewy = sin(sphi)*sin(stheta)
   ewz = cos(stheta)

   if stheta==0.:
      evx = 0.
      evy = 1.
      evz = 0.
      eux = 0.
      euy = 1.
      euz = 0.
   else:
      euxp = ewy
      euyp = - ewx
      euzp = 0.
      mag = sqrt(euxp*euxp + euyp*euyp)
      eux = euxp/mag
      euy = euyp/mag
      euz = 0.
      evx = ewy*euz - ewz*euy
      evy = eux*ewz - euz*ewx
      evz = ewx*euy - ewy*eux



   x_bin=(srad-sddist) * cos(sphi)*sin(stheta) + u*eux + v*evx
   y_bin=(srad-sddist) *sin(sphi)*sin(stheta) + u*euy + v*evy
   z_bin=(srad-sddist) *cos(stheta) + u*euz + v*evz

   return\
      x_source,y_source,z_source,\
      x_bin,y_bin,z_bin





configs={\
"ccbn"                             : config_vectors_ccbn,\
"ball"                             : config_vectors_spherical,\
"64slice"                          : config_vectors_helical_conebeam_64,\
"ball_flatpanel"                   : config_vectors_spherical_flatpanel\
}
