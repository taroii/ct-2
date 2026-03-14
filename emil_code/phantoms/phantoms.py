import numpy
import random
import sys

sys.path.append('../.')
from tomo3D import *

pi = numpy.pi


def simple_sphere():

   ph=phantom3D()
   e=ellipsoid(0.,0.,0.0,\
               1.2,\
               0.8,0.8,0.8,\
               0.0,0.0)
   ph.add_component(e)

   return ph

def simple_task(contrast=0.01,present=1,cx=0.,cy=0.,cz=0.):

   ph=phantom3D()
   e1=ellipsoid(0.,0.,0.,\
                1.0,\
                0.8,0.8,0.8,\
                0.0,0.0)
   ph.add_component(e1)

   if present:
      e2=ellipsoid(cx,cy,cz,\
                   contrast,\
                   0.2,0.2,0.2,\
                   0.0,0.0)
      ph.add_component(e2)


   return ph


def generate_defrise_with_balls():

   ph=phantom3D()
   e=ellipsoid(0.,0.,-0.7,\
               1.0,\
               0.8,0.8,0.1,\
	       0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.30,-0.7,\
               0.15,\
               0.04,0.04,0.03,\
               0.0,0.0)
   ph.add_component(e)



   e=ellipsoid(0.,0.,-0.35,\
               1.0,\
               0.8,0.8,0.1,\
	       0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.1,-0.35,\
               0.1,\
               0.04,0.04,0.03,\
               0.0,0.0)
   ph.add_component(e)


   e=ellipsoid(0.,0.,0.0,\
                1.0,\
                0.8,0.8,0.1,\
                0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,-0.10,0.0,\
               0.075,\
               0.04,0.04,0.03,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.,0.35,\
               1.0,\
		  0.8,0.8,0.1,\
		  0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,-0.3,0.35,\
               0.05,\
               0.04,0.04,0.03,\
               0.0,0.0)
   ph.add_component(e)


   e=ellipsoid(0.,0.,0.7,\
               1.0,\
               0.8,0.8,0.1,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,-0.5,0.7,\
               0.025,\
               0.04,0.04,0.03,\
               0.0,0.0)
   ph.add_component(e)

   return ph



def generate_defrise():

   ph=phantom3D()
   e=ellipsoid(0.,0.,-0.8,\
               1.0,\
               0.8,0.8,0.05,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.,-0.5,\
               1.0,\
               0.8,0.8,0.05,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.,-0.2,\
               1.0,\
               0.8,0.8,0.05,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.,0.1,\
               1.0,\
               0.8,0.8,0.05,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.,0.4,\
               1.0,\
               0.8,0.8,0.05,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,0.,0.7,\
               1.0,\
               0.8,0.8,0.05,\
               0.0,0.0)
   ph.add_component(e)
  
   return ph




def generate_3D_shepp_logan_HC():


   ph=phantom3D()

   e=ellipsoid(0.,0.,0.,\
               2.0,\
               0.69,0.90,0.92,\
               0.0,0.0)
   ph.add_component(e)
      
   e=ellipsoid(0.,0.,-0.0184,\
               -0.98,\
               0.6624,0.88,0.874,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(-0.22, -0.25, 0.0,\
               -0.2,\
               0.41, 0.21, 0.16,\
               -pi*72.0/180.,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.22, -0.25, 0.0,\
               -0.2,\
               0.31, 0.22, 0.11,\
               pi*72.0/180.,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,-0.25,0.35,\
               0.1,\
               0.21,0.35,0.25,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,-0.25,0.10,\
               0.1,\
               0.046,0.046,0.046,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(-0.08,-0.25,-0.605,\
               0.1,\
               0.046,0.02,0.023,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.06,-0.25,-0.605,\
               0.1,\
               0.046,0.02,0.023,\
               -pi/2.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.06,0.625,-0.105,\
               0.2,\
               0.056,0.1,0.04,\
               -pi/2.,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,0.625,0.10,\
               -0.2,\
               0.056,0.10,0.056,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,-0.25,-0.10,\
               0.1,\
               0.046,0.046,0.046,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,-0.25,-0.605,\
               0.1,\
               0.023,0.023,0.023,\
               0.0,0.0)
   ph.add_component(e)
  
   return ph 
   


def generate_3D_shepp_logan():


   ph=phantom3D()
   e=ellipsoid(0.,0.,0.,\
               2.0,\
               0.69,0.90,0.92,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,0.,-0.0184,\
               -0.98,\
               0.6624,0.88,0.874,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(-0.22, -0.25, 0.0,\
               -0.02,\
               0.41, 0.21, 0.16,\
               -pi*72.0/180.,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.22, -0.25, 0.0,\
               -0.02,\
               0.31, 0.22, 0.11,\
               pi*72.0/180.,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,-0.25,0.35,\
               0.01,\
               0.21,0.35,0.25,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,-0.25,0.10,\
               0.01,\
               0.046,0.046,0.046,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(-0.08,-0.25,-0.605,\
               0.01,\
               0.046,0.02,0.023,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.06,-0.25,-0.605,\
               0.01,\
               0.046,0.02,0.023,\
               -pi/2.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.06,0.625,-0.105,\
               0.02,\
               0.056,0.1,0.04,\
               -pi/2.,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,0.625,0.10,\
               -0.02,\
               0.056,0.10,0.056,\
               0.0,0.0)
   ph.add_component(e)

   e=ellipsoid(0.,-0.25,-0.10,\
               0.01,\
               0.046,0.046,0.046,\
               0.0,0.0)
   ph.add_component(e)
   
   e=ellipsoid(0.,-0.25,-0.605,\
               0.01,\
               0.023,0.023,0.023,\
               0.0,0.0)
   ph.add_component(e)

   return ph

if __name__=="__main__":

   nsize = 200
   i = image3D(shape=(nsize,nsize,nsize))
   ph  = generate_3D_shepp_logan()
   print "embedding ...."
   ph.embed_in(i) 
   figure(1)
   hold(False)
   ion()
   ans = raw_input("Do you want x-slices? ")

   if ans=="y":
      cp = imshow(i.mat[nsize/2],vmin = 0.95, vmax =1.05)
      draw()
      raw_input("hit enter")
      for j in range(nsize):
         cp.set_data(i.mat[j])
         draw()
         print("slice ",j)
         raw_input() 

   ans = raw_input("Do you want y-slices? ")

   if ans=="y":
      cp = imshow(i.mat[:,nsize/2],vmin = 0.95, vmax =1.05)
      draw()
      raw_input("hit enter")
      for j in range(nsize):
         cp.set_data(i.mat[:,j])
         draw()
         print("slice ",j)
         raw_input() 

   ans = raw_input("Do you want z-slices? ")

   if ans=="y":
      cp = imshow(i.mat[:,:,nsize/2],vmin = 0.95, vmax =1.05)
      draw()
      raw_input("hit enter")
      for j in range(nsize):
         cp.set_data(i.mat[:,:,j])
         draw()
         print("slice ",j)
         raw_input() 
