import numpy
import sys
sys.path.append('/home/emil/python/PyTomo3Diterative/')
import tomo3D as t3


def MGHread(filename):
  ''' read MGH formatted file into ndarray
  '''
  # NOTE:
  # no byte swapping
  fid = open(filename,'r');
  h_string = fid.read(512);
  fid.close()
  h_lines = h_string.splitlines()
  
  header = {'HEADER_BYTES=' : 0,\
            'DIM='          : 0,\
            'SIZE1='        : 0,\
            'SIZE2='        : 0,\
            'SIZE3='        : 0,\
            'BYTE_ORDER='   :'bla',\
            'TYPE='         :'bla'}


  l = header.keys()
  for i in range(len(h_lines)):
    for k in range(len(l)):
      s = h_lines[i].find(l[k])
      if s == 0:
      #if s >= 0:
        line = h_lines[i]
        num = line[len(l[k]):len(line)-1]
        num = num.lstrip()
        if num.isdigit():
          header[l[k]]=int(num)
        else:
          header[l[k]]=num
          
          
  if header["TYPE="] == 'unsigned_short':
    datatype = 'uint16'
  elif header["TYPE="] == 'short':
    datatype = 'int16'
  elif header["TYPE="] == 'float':
    datatype = 'float32'
  elif header["TYPE="] == 'double':
    datatype = 'float64'
    
    
    
  fid = open(filename,'r');
  fid.seek(header["HEADER_BYTES="])
  a = numpy.fromfile(fid,dtype = datatype);
  fid.close()

## array layout from matlab or IDL is fortran
## python default layout is C

  if header["DIM="] == 3:
    a = a.reshape(header["SIZE3="],header["SIZE2="],header["SIZE1="])
    a = a.transpose()
    a = a.copy()  # change from fortran memory layout to C memory layout
  elif header["DIM="] == 2:
    a = a.reshape(header["SIZE2="],header["SIZE1="])
    a = a.transpose()
    a = a.copy()  # change from fortran memory layout to C memory layout

  return a


def MGHwrite(filename,a,d=dict(),write_type=''):
  '''  write ndarray a into MGH formatted file.
  Datatypes: float32, float64, int16, uint16
  '''
  if write_type != '':
    a = a.astype(write_type)
  
  if a.dtype == 'int64':
    a = a.astype('int16')
  if a.dtype == 'float64':
    a = a.astype('float32')
      
      
  header = {'HEADER_BYTES=' : 512,\
            'DIM='          : a.ndim,\
            'SIZE1='        : 0,\
            'SIZE2='        : 0,\
            'BYTE_ORDER='   :'little_endian',\
            'TYPE='         :'bla1'}
  
  
  if a.dtype == 'float32':
    header["TYPE="] = 'float'
  elif a.dtype == 'float64':
    header["TYPE="] = 'double'
  elif a.dtype == 'uint16':
    header["TYPE="] = 'unsigned_short'
  elif a.dtype == 'int16':
    header["TYPE="] = 'short'
    
  size = a.shape
  #print "%i" %(size.__len__())
  
  header["SIZE1="] = size[0]
  header["SIZE2="] = size[1]
  if size.__len__() == 3:
    header["SIZE3="] = size[2]

  
  # added 2/9/07
  header.update(d)
  
  hstr = '{\n'
  for key,value in header.items():
    hstr = hstr +  '%s %s;\n' % (key, value)
  hstr = hstr + '}\n\n'
  
  fill = 512 - len(hstr)

  fid = open(filename,'w');
  fid.write(hstr)
  for i in range(fill):
    fid.write(' ')
  #print fid.tell()
  # done writing header
  fid.close()

  fid = open(filename,'a');
  fid.seek(header["HEADER_BYTES="])
  a = a.transpose()  # change to fortran memory layout 
  a.tofile(fid);
  fid.close()

  return 

 
def WriteSino(filename,sino,comment_line = ''):
  ''' write sinogram3D object into MGH-formatted file
  '''
  # make dictionary of sinogram parameters
  sinoparms = dict(img_type = 'sinogram3D',\
                   config_name = sino.config_name,\
                   shape = sino.mat.shape,\
                   slen = sino.slen,\
                   s0 = sino.s0,\
                   ulen = sino.ulen,\
                   u0 = sino.u0,\
                   vlen = sino.vlen,\
                   v0 = sino.v0,\
                   s_include_endpoints = sino.s_include_endpoints,\
                   comment = comment_line)
  sinoparms.update(sino.parms)
  MGHwrite(filename,sino.mat,sinoparms)
    
  return 
 
def WriteImage(filename,image,comment_line = ''):
  ''' write image3D object into MGH-formatted file
  '''
  
  # make dictionary of image parameters
  imgparms = dict(img_type = 'image3D',\
                   shape = image.mat.shape,\
                   xlen = image.xlen,\
                   ylen = image.ylen,\
                   zlen = image.zlen,\
                   x0 = image.x0,\
                   y0 = image.y0,\
                   z0 = image.z0,\
                   comment = comment_line)
  
  MGHwrite(filename,image.mat,imgparms)
    
  return 
 
def GEdata2sino(filename,pixsize = 0.01,Laterality='n'):
    ''' GEdata2sino
    reads a sequence of GE projection images and sets up 
    the sinogram for reconstruction (i.e., non-equidistant angles  ..)
    '''
    img = MGHread(filename);
    #img = m.MGHread(filename);
    #img = m.MGHread('/nfs/nemo/data1/ireiser/CalcRecon/5385Rproj.en');
    #img = m.MGHread('/nfs/nemo/data1/ireiser/CalcRecon/5370Rproj.en');
    #img = m.MGHread('/nfs/nemo/data1/ireiser/CalcRecon/5447Lproj.en');

########################################################
    #   line up the image, set angles
########################################################
    if Laterality == 'R':
        img = img.transpose(2,0,1)
        img1 = numpy.zeros(img.shape,dtype='float32')
        for i in range(img.shape[1]):
            img1[:,i,:] = img[:,i,:]
        img = img1
        del(img1)
        ang = numpy.array([-24.98,-19.52,-14.19,-8.99,-3.90,1.09,5.99,10.82,15.59,20.32,25.00])
    elif Laterality == 'L':
        img = img.transpose(2,0,1)
        #img1 = numpy.ndarray(img.shape)
        img1 = numpy.zeros(img.shape,dtype='float32')
        for i in range(img.shape[1]):
            img1[:,i,:] = img[:,img.shape[1]-i-1,:]
        img = img1
        del(img1)
        ang = numpy.array([-25.00,-20.32,-15.59,-10.82,-5.99,-1.09,3.90,8.99,14.19,19.52,24.98])
    ang = (ang+90)/180.*numpy.pi
            
##########################################################
    #    set the geometry: GE. Initialize sinogram & angles
##########################################################
    configName0 = 'tomosynthesis'
    cent_to_det = 23.7
    source_radius = 44.3
    
    shape0 = img.shape
    slen0 = img.shape[0] 
    ulen0 = img.shape[1]*pixsize
    vlen0 = img.shape[2]*pixsize
    s00 = numpy.pi/4.
    u00 = 0.
    v00 = -vlen0/2.

    sino = t3.sinogram3D(config_name=configName0,\
                         parms={"center_to_detector": cent_to_det,\
                                "radius"            : source_radius},\
                         shape = img.shape,\
                         slen = slen0, s0 = s00,\
                         ulen = ulen0, u0 = u00,\
                         vlen = vlen0, v0 = v00)

    img[img < 0.0] = 0.0
    sino.mat = img
    for ip in range(sino.mat.shape[0]):
        sino.frame_vectors[ip] = sino.frame(ang[ip],sino.parms)
    return sino




def str2numid(strid):
    import string
    ''' convert a string ID to numeric ID.
    example:
    str2numid('5220L') -> 522076
    '''

    num_part = string.atoi(strid[0:4])
    str_part = strid[4]
    if str_part == 'L':
        numid = num_part*100+76
    elif str_part == 'R':
        numid = num_part*100+82

    return numid


def num2strid(numid):
    ''' convert a numeric ID to a string ID.
    example:
    str2numid(522076) -> '5220L'
    '''

    #numid = 522076
    num_part = int(numid/100)
    str_part = int(numid - 100 * num_part)
    s1 = str(num_part)
    s2 = chr(str_part)
    strid = "%s%s"%(s1,s2)

    
    return strid


def proj_point(sino,point):
    '''
    project points onto detector plane
    point: [npoints x 3] array of [x,y,z] coordinates

    returns : [sino.ns x npoints x 2] array of [u,v] coordinates on detector
    '''
    
    point_proj = numpy.zeros([sino.ns,point.shape[0],point.shape[1]-1])
    print point_proj.shape


    for ip in range(sino.ns):
        RSource = sino.frame_vectors[ip,0:3]
        DetCenter = sino.frame_vectors[ip,6:9]
        eu = sino.frame_vectors[ip,9:12]
        ev = sino.frame_vectors[ip,12:15]
        ew = sino.frame_vectors[ip,15:18]
        udr = numpy.dot((DetCenter-RSource),eu)
        vdr = numpy.dot((DetCenter-RSource),ev)
        TCent = DetCenter - udr*eu - vdr * ev
        STdist = numpy.sqrt(sum((RSource-TCent)**2))
        
        for npoint in range(point.shape[0]):
            diff = point[npoint,:] - RSource
            zoom = abs(numpy.dot(diff,ew))/STdist
            DetBin = diff/zoom + RSource
            DetPt = DetBin - DetCenter
            DetU = numpy.dot(DetPt,eu)
            DetV = numpy.dot(DetPt,ev)
            point_proj[ip,npoint,0] = DetU
            point_proj[ip,npoint,1] = DetV
            
    return point_proj




if __name__ == "__main__":
#  a = numpy.array([range(24)])
#  a = a.reshape(4,3,2)
#  
#  print a[2,0,1]
#  
#  print a.dtype
#
#  MGHwrite("out1.en",a)
#  #b = MGHread("out1.en")
#  sino = t3.sinogram3D(config_name='tomosynthesis',\
#                       parms={"center_to_detector": 20.,\
#                              "radius"            : 46.},\
#                       shape=(11,90,190),\
#                       slen= pi/2., ulen=9.0,vlen=19.0,\
#                       s0 = -pi/4., u0 =0., v0 = 9.5,s_include_endpoints=1)
  
#   a = MGHread('~/recon/data4testing/test3L.proj')
    a=32
