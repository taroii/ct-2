
      subroutine spherical_radon_linear(fov,radonmat,indradon,
     +       source_locations,
     +       nsa,nsb,ns,nr,dr,r0,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      Dimension source_locations(ns,3)
      Dimension indradon(nsa,nsb,nr)
      Dimension radonmat(nsa,nsb,nr),smat(nx,ny,nz)

Cf2py intent(in) nsa
Cf2py intent(in) nsb
Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) source_locations
cf2py depend(ns) source_locations
Cf2py intent(out) radonmat
cf2py depend(nsa) radonmat
cf2py depend(nsb) radonmat
cf2py depend(nr) radonmat
Cf2py intent(in) indradon
cf2py depend(nsa) indradon
cf2py depend(nsb) indradon
cf2py depend(nr) indradon
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat


      ip = 0
      pi = 4.*asin(1.0)

      Do ipa=1,nsa
      Do ipb = 1,nsb
         ip = ip +1

         xSource=source_locations(ip,1)
         ySource=source_locations(ip,2)
         zSource=source_locations(ip,3)

         dd = sqrt(xSource*xSource+ySource*ySource+zSource*zSource)
         gamSource = asin(zSource/dd)
         thetaSource = atan2(ySource,xSource)
c         print *, dd,gamSource,thetaSource
         cg = cos(gamSource)
         sg = sin(gamSource)
         ct = cos(thetaSource)
         st = sin(thetaSource)
         


      Do jp=1,nr



         r = r0+(jp-0.5)*dr 
         rho = r

         rat = (-fov*fov + dd*dd + rho*rho)/(2.*dd*rho)
         thetamax = acos(rat)
c         thetamax = pi/2.
         dtheta  = dx/rho
         ntheta = int(thetamax/dtheta)

         if (indradon(ipa,ipb,jp).eq.1) then

            Do ith = 1,ntheta
               theta = (ith - 0.5)*dtheta
               r1 = rho*sin(theta)
               x1 = dd - rho*cos(theta)
               dalpha = dx/r1
               nalpha = int(2.*pi/dalpha)

            Do ia = 1,nalpha
               alpha = (ia-0.5)*dalpha
               y1 = r1*cos(alpha)
               z1 = r1*sin(alpha)
              
               
               x =x1*ct*cg - y1*st - z1*ct*sg
               y = (x1*st*cg + y1*ct - z1*st*sg)
               z =x1*sg + z1*cg


             
               px = (x - x0 - 0.5*dx)/dx
               py = (y - y0 - 0.5*dy)/dy
               pz = (z - z0 - 0.5*dz)/dz

               ipx = floor(px)
               ipy = floor(py)
               ipz = floor(pz)
c               print *,ipx,ipy,ipz, x1,y1,z1


c               if (((ipx.ge.1).and.(ipx.lt.nx)).and.
c     +             ((ipy.ge.1).and.(ipy.lt.ny)).and.
c     +             ((ipz.ge.1).and.(ipz.lt.nz))) then
                  wx = px - ipx
                  wy = py - ipy
                  wz = pz - ipz

                  t1 = smat(ipx,ipy,ipz)*(1.-wy)*(1.-wx)
                  t2 = smat(ipx+1,ipy,ipz)*(1.-wy)*wx
                  t3 = smat(ipx,ipy+1,ipz)*(wy)*(1.-wx)
                  t4 = smat(ipx+1,ipy+1,ipz)*wy*wx
                  t5 = smat(ipx,ipy,ipz+1)*(1.-wy)*(1.-wx)
                  t6 = smat(ipx+1,ipy,ipz+1)*(1.-wy)*wx
                  t7 = smat(ipx,ipy+1,ipz+1)*(wy)*(1.-wx)
                  t8 = smat(ipx+1,ipy+1,ipz+1)*wy*wx

                  res =dx*dx*((1.-wz)*(t1+t2+t3+t4)+wz*(t5+t6+t7+t8))
                  radonmat(ipa,ipb,jp) = radonmat(ipa,ipb,jp) + res
c               end if
           end do  !alpha loop
           end do  !theta loop
         end if  ! indicator conditional

      end do ! t loop
      end do ! sb loop
      end do ! sa loop

      return
      end




      subroutine artstep_linear(reg_fac,fov,radonmat,indradon,
     +       source_locations,
     +       nsa,nsb,ns,nr,dr,r0,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      Real*8 total

      Dimension source_locations(ns,3)
      Dimension indradon(nsa,nsb,nr)
      Dimension radonmat(nsa,nsb,nr),smat(nx,ny,nz), wmat(nx,ny,nz)

Cf2py intent(in) nsa
Cf2py intent(in) nsb
Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) source_locations
cf2py depend(ns) source_locations
Cf2py intent(in) radonmat
cf2py depend(nsa) radonmat
cf2py depend(nsb) radonmat
cf2py depend(nr) radonmat
Cf2py intent(in) indradon
cf2py depend(nsa) indradon
cf2py depend(nsb) indradon
cf2py depend(nr) indradon
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat


      ip = 0
      pi = 4.*asin(1.0)

      Do i = 1, nx
         Do j = 1, ny
            Do k = 1, nz
               wmat(i,j,k) = 0.
            End Do
         End Do
      End Do

      Do ipa=1,nsa
      Do ipb = 1,nsb
         ip = ip +1

         xSource=source_locations(ip,1)
         ySource=source_locations(ip,2)
         zSource=source_locations(ip,3)

         dd = sqrt(xSource*xSource+ySource*ySource+zSource*zSource)
         gamSource = asin(zSource/dd)
         thetaSource = atan2(ySource,xSource)
c         print *, dd,gamSource,thetaSource
         cg = cos(gamSource)
         sg = sin(gamSource)
         ct = cos(thetaSource)
         st = sin(thetaSource)
         


      Do jp=1,nr


         dataval = radonmat(ipa,ipb,jp)

         r = r0+(jp-0.5)*dr 
         rho = r

         rat = (-fov*fov + dd*dd + rho*rho)/(2.*dd*rho)
         thetamax = acos(rat)
c         thetamax = pi/2.
         dtheta  = dx/rho
         ntheta = int(thetamax/dtheta)

         if (indradon(ipa,ipb,jp).eq.1) then

            irp = 0
            surfsum = 0.

            Do ith = 1,ntheta
               theta = (ith - 0.5)*dtheta
               r1 = rho*sin(theta)
               x1 = dd - rho*cos(theta)
               dalpha = dx/r1
               nalpha = int(2.*pi/dalpha)

            Do ia = 1,nalpha
               alpha = (ia-0.5)*dalpha
               y1 = r1*cos(alpha)
               z1 = r1*sin(alpha)
              
               
               x =x1*ct*cg - y1*st - z1*ct*sg
               y = (x1*st*cg + y1*ct - z1*st*sg)
               z =x1*sg + z1*cg


             
               px = (x - x0 - 0.5*dx)/dx
               py = (y - y0 - 0.5*dy)/dy
               pz = (z - z0 - 0.5*dz)/dz

               ipx = floor(px)
               ipy = floor(py)
               ipz = floor(pz)
c               print *,ipx,ipy,ipz, x1,y1,z1


c               if (((ipx.ge.1).and.(ipx.lt.nx)).and.
c     +             ((ipy.ge.1).and.(ipy.lt.ny)).and.
c     +             ((ipz.ge.1).and.(ipz.lt.nz))) then
                  wx = px - ipx
                  wy = py - ipy
                  wz = pz - ipz

                  t1 = smat(ipx,ipy,ipz)*(1.-wy)*(1.-wx)
                  t2 = smat(ipx+1,ipy,ipz)*(1.-wy)*wx
                  t3 = smat(ipx,ipy+1,ipz)*(wy)*(1.-wx)
                  t4 = smat(ipx+1,ipy+1,ipz)*wy*wx
                  t5 = smat(ipx,ipy,ipz+1)*(1.-wy)*(1.-wx)
                  t6 = smat(ipx+1,ipy,ipz+1)*(1.-wy)*wx
                  t7 = smat(ipx,ipy+1,ipz+1)*(wy)*(1.-wx)
                  t8 = smat(ipx+1,ipy+1,ipz+1)*wy*wx

                  res =dx*dx*((1.-wz)*(t1+t2+t3+t4)+wz*(t5+t6+t7+t8))

                  wmat(ipx,ipy,ipz) = wmat(ipx,ipy,ipz)+
     +                dx*dx*(1.-wz)*(1.-wy)*(1.-wx)

                  wmat(ipx+1,ipy,ipz) = wmat(ipx+1,ipy,ipz)+
     +                dx*dx*(1.-wz)*(1.-wy)*wx

                  wmat(ipx,ipy+1,ipz) = wmat(ipx,ipy+1,ipz)+
     +               dx*dx*(1.-wz)*(wy)*(1.-wx)

                  wmat(ipx+1,ipy+1,ipz) = wmat(ipx+1,ipy+1,ipz)+
     +               dx*dx*(1.-wz)*(wy)*(wx)

                  wmat(ipx,ipy,ipz+1) = wmat(ipx,ipy,ipz+1)+
     +                dx*dx*(wz)*(1.-wy)*(1.-wx)

                  wmat(ipx+1,ipy,ipz+1) = wmat(ipx+1,ipy,ipz+1)+
     +                dx*dx*(wz)*(1.-wy)*wx

                  wmat(ipx,ipy+1,ipz+1) = wmat(ipx,ipy+1,ipz+1)+
     +               dx*dx*(wz)*(wy)*(1.-wx)

                  wmat(ipx+1,ipy+1,ipz+1) = wmat(ipx+1,ipy+1,ipz+1)+
     +               dx*dx*(wz)*(wy)*(wx)


                  surfsum = surfsum + res
c               end if
            end do  !alpha loop
            end do  !theta loop

            total = 0.
            Do iw = 1,nx
               Do jw = 1,ny
                  Do kw = 1,nz
                     total = total + wmat(iw,jw,kw)**2
                  End Do
               End Do
            End do

c            print *,ipa, ipb, jp, total

            if (total.lt.0.0000001) then
               total = 1.
            end if

            Do iw = 1,nx
               Do jw = 1,ny
                  Do kw = 1,nz
                    smat(iw,jw,kw)= smat(iw,jw,kw)+
     +                 reg_fac*(dataval-surfsum)*wmat(iw,jw,kw)/total
                    wmat(iw,jw,kw) = 0.
                  End Do
               End Do
            End do

         end if  ! indicator conditional

      end do ! t loop
      end do ! sb loop
      end do ! sa loop

      return
      end



