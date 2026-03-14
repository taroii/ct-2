from numpy import sin, cos, pi, sqrt



def config_source_locations_z_sphere(sa,sb,parms):


   radius=parms["radius"]
  
   rho = sqrt(radius*radius - sa*sa) 
   x_source=rho* cos(sb)
   y_source=rho*sin(sb)
   z_source=sa


   return x_source,y_source,z_source




configs={\
"z_sphere"                : config_source_locations_z_sphere\
}
