import numpy

sin=numpy.sin
sqrt=numpy.sqrt
cos=numpy.cos
tan=numpy.tan
pi=numpy.pi
array = numpy.array


def config_vectors_circular_conebeam(s,parms):


   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   
   x_source=radius * cos(s)
   y_source=radius *sin(s)
   z_source=0.

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= 0.

   eux = -sin(s)
   euy =  cos(s)
   euz = 0.

   evx = 0.
   evy = 0.
   evz = 1.

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_giotto(s,parms):


   radius=parms["radius"]
   center_to_detector=parms["center_to_detector"]
   
   x_source=radius * sin(s)
   y_source=0.
   z_source=radius * cos(s)

   x_traj_tan= radius *cos(s)
   y_traj_tan= 0.
   z_traj_tan= -radius *sin(s)



#   x_det_center= - center_to_detector*sin(s)
   x_det_center= 0.
   y_det_center=0.
   z_det_center=-center_to_detector


#   eux = cos(s)
#   euy = 0.
#   euz = -sin(s)
   eux = 1.
   euy = 0.
   euz = 0.

   evx = 0.
   evy = 1.
   evz = 0.

#   ewx = sin(s)
#   ewy =0.
#   ewz = cos(s)
   ewx = 0.
   ewy = 0.
   ewz = 1.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_zxcircularCB(s,parms):


   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   
   x_source=radius * sin(s)
   y_source=0.
   z_source=radius * cos(s)

   x_traj_tan= radius *cos(s)
   y_traj_tan= 0.
   z_traj_tan= -radius *sin(s)



   x_det_center=(radius - source_to_detector )*sin(s)
   y_det_center=0.
   z_det_center=(radius - source_to_detector )*cos(s) 


#   eux = cos(s)
#   euy = 0.
#   euz = -sin(s)
   eux = 1.
   euy = 0.
   euz = 0.

   evx = 0.
   evy = 1.
   evz = 0.

#   ewx = sin(s)
#   ewy =0.
#   ewz = cos(s)
   ewx = 0.
   ewy = 0.
   ewz = 1.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz



def config_vectors_array_circular_conebeam(si,parms):


   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   source_angles = parms["angles"]
   ucenters = parms["ucenters"]
   vcenters = parms["vcenters"]
 
   isi = int(round(si)) 
   s = source_angles[isi] 
   x_source=radius * cos(s)
   y_source=radius *sin(s)
   z_source=0.

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s) - ucenters[isi]*sin(s)
   y_det_center=(radius - source_to_detector )*sin(s) + ucenters[isi]*cos(s)
   z_det_center= vcenters[isi]

   eux = -sin(s)
   euy =  cos(s)
   euz = 0.

   evx = 0.
   evy = 0.
   evz = 1.

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_circular_conebeam_flipped(s,parms):


   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   
   x_source=radius * cos(s)
   y_source=radius *sin(s)
   z_source=0.

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= 0.

   evx = -sin(s)
   evy =  cos(s)
   evz = 0.

   eux = 0.
   euy = 0.
   euz = 1.

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz

def config_vectors_angled_circular_conebeam(s,parms):


   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   soffset = parms["soffset"]
   
   x_source=radius * cos(s)
   y_source=radius *sin(s)
   z_source=0.

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= 0.

   eux = -sin(s-soffset)
   euy =  cos(s-soffset)
   euz = 0.

   evx = 0.
   evy = 0.
   evz = 1.

   ewx = cos(s-soffset)
   ewy = sin(s-soffset)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_rot_det_circular_conebeam(s,parms):
   """ the detector is rotated by an angle tilt about its normal
   """

   radius = parms["radius"]
   source_to_detector=parms["source_to_detector"]
   ta = parms['tilt']

   x_source = radius * cos(s)
   y_source = radius *sin(s)
   z_source = 0.

   x_traj_tan = -radius *sin(s)
   y_traj_tan = radius *cos(s)
   z_traj_tan = 0.


   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= 0.

   euxt = -sin(s)
   euyt =  cos(s)
   euzt = 0.

   evxt = 0.
   evyt = 0.
   evzt = 1.

   eux = cos(ta)*euxt + sin(ta)*evxt
   euy = cos(ta)*euyt + sin(ta)*evyt
   euz = cos(ta)*euzt + sin(ta)*evzt

   evx = -sin(ta)*euxt + cos(ta)*evxt
   evy = -sin(ta)*euyt + cos(ta)*evyt
   evz = -sin(ta)*euzt + cos(ta)*evzt

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz

def config_vectors_offcenter_circular_conebeam(s,parms):

   xc = parms["xc"]
   yc = parms["yc"]
   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   
   x_source=radius * cos(s) + xc
   y_source=radius *sin(s)  + yc
   z_source=0.

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s) +xc
   y_det_center=(radius - source_to_detector )*sin(s) +yc
   z_det_center= 0.

   eux = -sin(s)
   euy =  cos(s)
   euz = 0.

   evx = 0.
   evy = 0.
   evz = 1.

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_detector_tipped_circular_conebeam(s,parms):


   sdirection=parms["sdirection"]
   tangle=parms["tangle"]
   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
  
   sp =s*sdirection
   x_source=radius * cos(sp)
   y_source=radius *sin(sp)
   z_source= 0.

   x_traj_tan=- radius *sin(sp)
   y_traj_tan=radius *cos(sp)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(sp)
   y_det_center=(radius - source_to_detector )*sin(sp)
   z_det_center= 0.

   eu0 = array([-sin(sp),cos(sp),0.])
   ev0 = array([0.,0.,1.])
   ew0 = array([cos(sp),sin(sp),0.])


   eut = eu0*cos(tangle) - ev0*sin(tangle)
   evt = eu0*sin(tangle) + ev0*cos(tangle)
   ewt = ew0

   eux =eut[0]
   euy =eut[1]
   euz =eut[2]
   evx =evt[0]
   evy =evt[1]
   evz =evt[2]
   ewx =ewt[0]
   ewy =ewt[1]
   ewz =ewt[2]



   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz



def config_vectors_tipped_circular_conebeam(s,parms):


   wangle=parms["wangle"]
   rangle=parms["rangle"]
   tangle=parms["tangle"]
   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   
   x_source=radius * cos(s)*cos(rangle)
   y_source=radius *sin(s)*cos(rangle)
   z_source=-radius*sin(rangle)

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s)*cos(rangle)
   y_det_center=(radius - source_to_detector )*sin(s)*cos(rangle)
   z_det_center= (radius - source_to_detector )*sin(rangle)

   eu0 = array([-sin(s),cos(s),0.])
   ev0 = array([0.,0.,1.])
   ew0 = array([cos(s),sin(s),0.])


   eut = eu0*cos(tangle) - ev0*sin(tangle)
   evt = eu0*sin(tangle) + ev0*cos(tangle)
   ewt = ew0
   eur = eut
   evr = evt*cos(rangle) - ewt*sin(rangle)
   ewr = evt*sin(rangle) + ewt*cos(rangle)
   euw = eur*cos(wangle) - ewr*sin(wangle)
   evw = evr
   eww = eur*sin(wangle) + ewr*cos(wangle)

   eux =euw[0]
   euy =euw[1]
   euz =euw[2]
   evx =evw[0]
   evy =evw[1]
   evz =evw[2]
   ewx =eww[0]
   ewy =eww[1]
   ewz =eww[2]



   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz




def config_vectors_z_circular_conebeam(s,parms):

   zoffset=parms["zoffset"]
   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   
   x_source=radius * cos(s)
   y_source=radius *sin(s)
   z_source=zoffset

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=0.



   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= zoffset

   eux = -sin(s)
   euy =  cos(s)
   euz = 0.

   evx = 0.
   evy = 0.
   evz = 1.

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz




def config_vectors_tipped_helical_conebeam(s,parms):


   rangle=parms["rangle"]
   tangle=parms["tangle"]
   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   pitch=parms["pitch"]

   x_source=radius *cos(s)
   y_source=radius *sin(s)
   z_source=pitch * s/(2.*pi)

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=pitch /(2.*pi)



   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= pitch * s/(2.*pi)

   eux = -sin(s)*cos(tangle)
   euy =  cos(s)*cos(tangle) 
   euz = sin(tangle)

   evx = sin(s)*sin(tangle)   *cos(rangle) + cos(s)*sin(rangle)
   evy = - cos(s)*sin(tangle) *cos(rangle) + sin(s)*sin(rangle)
   evz =  cos(tangle)         *cos(rangle)

   ewx = -sin(s)*sin(tangle)  *sin(rangle) + cos(s) *cos(rangle)
   ewy =  cos(s)*sin(tangle)  *sin(rangle) + sin(s) *cos(rangle)
   ewz = -cos(tangle)         *sin(rangle)


   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_helical_conebeam(s,parms):


   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]
   pitch=parms["pitch"]
  
 
   x_source=radius *cos(s)
   y_source=radius *sin(s)
   z_source=pitch * s/(2.*pi)

   x_traj_tan=- radius *sin(s)
   y_traj_tan=radius *cos(s)
   z_traj_tan=pitch /(2.*pi)



   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)
   z_det_center= pitch * s/(2.*pi)

   eux = -sin(s)
   euy =  cos(s)
   euz = 0.

   evx = 0.
   evy = 0.
   evz = 1.

   ewx = cos(s)
   ewy = sin(s)
   ewz = 0.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_tomosynthesis(s,parms):


   radius=parms["radius"]
   center_to_detector=parms["center_to_detector"]
  
 
   x_source=0.
   y_source=radius *cos(s)
   z_source=radius *sin(s)

   x_traj_tan= 0.
   y_traj_tan=- radius *sin(s)
   z_traj_tan=radius *cos(s)



   x_det_center= 0.
   y_det_center= 0.
   z_det_center= -center_to_detector

   eux = 1.
   euy = 0.
   euz = 0.

   evx = 0.
   evy = 1.
   evz = 0.

   ewx = 0.
   ewy = 0.
   ewz = 1.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_tomosynthesis_sliding_detector(s,parms):


   radius=parms["radius"]
   center_to_detector=parms["center_to_detector"]
  
 
   x_source=0.
   y_source=radius *cos(s)
   z_source=radius *sin(s)

   x_traj_tan= 0.
   y_traj_tan=- radius *sin(s)
   z_traj_tan=radius *cos(s)



   x_det_center= 0.
   y_det_center= -center_to_detector * cos(s)/sin(s)
   z_det_center= -center_to_detector

   eux = 1.
   euy = 0.
   euz = 0.

   evx = 0.
   evy = 1.
   evz = 0.

   ewx = 0.
   ewy = 0.
   ewz = 1.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz



def config_vectors_two_panel_pet(s,parms):
   '''configuration set up so that s rasters top detector.
(0,0,0) is the middle of the top corner pixel'''


   separation=parms["separation"]
   detector_width=parms["detector_width"]
   detector_height=parms["detector_height"]
   raster_spacing=parms["raster_spacing"]

   nrast=int(s/detector_width)
   w=s-detector_width*nrast
   h= raster_spacing*nrast
  
 
   x_source=w
   y_source=h
   z_source=0.

   x_traj_tan= 1.
   y_traj_tan=0.
   z_traj_tan=0.



   x_det_center= 0.
   y_det_center= 0.
   z_det_center= separation

   eux = 1.
   euy = 0.
   euz = 0.

   evx = 0.
   evy = 1.
   evz = 0.

   ewx = 0.
   ewy = 0.
   ewz = -1.

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz

def config_vectors_yline_scan(s,parms):
   '''linear scan along ev in x-y plane'''


   center_to_detector=parms["center_to_detector"]
   center_to_line=parms["center_to_line"]

   x_source= center_to_line
   y_source= s
   z_source= 0.

   gam = arctan2(s,center_to_line)
   
   uvx =  center_to_line/sqrt(s*s + center_to_line**2)
   uvy = s/sqrt(s*s + center_to_line**2)

   x_traj_tan= 0.0
   y_traj_tan= 1.0
   z_traj_tan= 0.0



#   x_det_center= -center_to_detector*sin(gam)
   x_det_center= -center_to_detector*uvx
   y_det_center= -center_to_detector*uvy
#   z_det_center= -center_to_detector*cos(gam)
   z_det_center= 0.

#   eux = cos(gam)
   eux = -uvy
   euy = uvx
#   euz = -sin(gam)
   euz = 0.0

   evx = 0.0
   evy = 0.0
   evz = 1.0

#   ewx = sin(gam)
   ewx = uvx
   ewy = uvy
#   ewz = cos(gam) 
   ewz = 0.
   


   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz


def config_vectors_xline_const_angle_scan(s,parms):
   '''linear scan along eu in x-y plane'''


   center_to_detector=parms["center_to_detector"]
   center_to_line=parms["center_to_line"]

   sp = center_to_line*tan(s)

   x_source= sp
   y_source= - center_to_line
   z_source= 0.

   
   uvx = sp/sqrt(sp*sp + center_to_line**2)
   uvy = - center_to_line/sqrt(sp*sp + center_to_line**2)

   x_traj_tan= center_to_line/(cos(s)**2.)
   y_traj_tan= 0.0
   z_traj_tan= 0.0



   x_det_center= -center_to_detector*uvx
   y_det_center= -center_to_detector*uvy
   z_det_center= 0.

   eux = -uvy
   euy = uvx
   euz = 0.0
#   eux = 1.0
#   euy = 0.0
#   euz = 0.0

   evx = 0.0
   evy = 0.0
   evz = 1.0

   ewx = uvx
   ewy = uvy
   ewz = 0.
#   ewx = 0.0
#   ewy = -1.0
#   ewz = 0.
   

   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz





def config_vectors_xline_scan(s,parms):
   '''linear scan along eu in x-y plane'''


   center_to_detector=parms["center_to_detector"]
   center_to_line=parms["center_to_line"]

   x_source= s
   y_source= - center_to_line
   z_source= 0.

   gam = arctan2(s,center_to_line)
   
   uvx = s/sqrt(s*s + center_to_line**2)
   uvy = - center_to_line/sqrt(s*s + center_to_line**2)

   x_traj_tan= 1.0
   y_traj_tan= 0.0
   z_traj_tan= 0.0



   x_det_center= -center_to_detector*uvx
   y_det_center= -center_to_detector*uvy
#   x_det_center = -s*center_to_detector/center_to_line
#   y_det_center = center_to_detector
   z_det_center= 0.

   eux = -uvy
   euy = uvx
   euz = 0.0
#   eux = 1.0
#   euy = 0.0
#   euz = 0.0

   evx = 0.0
   evy = 0.0
   evz = 1.0

   ewx = uvx
   ewy = uvy
   ewz = 0.
#   ewx = 0.0
#   ewy = -1.0
#   ewz = 0.
   


   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz



def config_vectors_line_scan(s,parms):
   '''linear scan along eu'''


   center_to_detector=parms["center_to_detector"]
   center_to_line=parms["center_to_line"]

   x_source= s
   y_source= 0.
   z_source= center_to_line

   gam = arctan2(s,center_to_line)
   
   uvx = s/sqrt(s*s + center_to_line**2)
   uvz = center_to_line/sqrt(s*s + center_to_line**2)

   x_traj_tan= 1.0
   y_traj_tan= 0.0
   z_traj_tan= 0.0



#   x_det_center= -center_to_detector*sin(gam)
   x_det_center= -center_to_detector*uvx
   y_det_center= 0.
#   z_det_center= -center_to_detector*cos(gam)
   z_det_center= -center_to_detector*uvz

#   eux = cos(gam)
   eux = uvz
   euy = 0.0
#   euz = -sin(gam)
   euz = -uvx

   evx = 0.0
   evy = 1.0
   evz = 0.0

#   ewx = sin(gam)
   ewx = uvx
   ewy = 0.0
#   ewz = cos(gam) 
   ewz = uvz
   


   return\
      x_source,y_source,z_source,\
      x_traj_tan,y_traj_tan,z_traj_tan,\
      x_det_center,y_det_center,z_det_center,\
      eux,euy,euz,\
      evx,evy,evz,\
      ewx,ewy,ewz




configs={\
"circular_conebeam"                : config_vectors_circular_conebeam,\
"giotto"                           : config_vectors_giotto,\
"zxcircularCB"                     : config_vectors_zxcircularCB,\
"array_circular_conebeam"          : config_vectors_array_circular_conebeam,\
"circular_conebeam_flipped"        : config_vectors_circular_conebeam_flipped,\
"angled_circular_conebeam"         : config_vectors_angled_circular_conebeam,\
"rot_det_circular_conebeam"        : config_vectors_rot_det_circular_conebeam,\
"offcenter_circular_conebeam"      : config_vectors_offcenter_circular_conebeam,\
"detector_tipped_circular_conebeam": config_vectors_detector_tipped_circular_conebeam,\
"tipped_circular_conebeam"         : config_vectors_tipped_circular_conebeam,\
"z_circular_conebeam"              : config_vectors_z_circular_conebeam,\
"helical_conebeam"                 : config_vectors_helical_conebeam,\
"tipped_helical_conebeam"          : config_vectors_tipped_helical_conebeam,\
"tomosynthesis"                    : config_vectors_tomosynthesis,\
"tomosynthesis_sliding_detector"   : config_vectors_tomosynthesis_sliding_detector,\
"line_scan"                        : config_vectors_line_scan,\
"xline_scan_const_angle"           : config_vectors_xline_const_angle_scan,\
"xline_scan"                       : config_vectors_xline_scan,\
"yline_scan"                       : config_vectors_yline_scan,\
"two_panel_pet"                    : config_vectors_two_panel_pet\
}
