c This set of functions are written for proj. and back_proj.
c using a list of rays.

c the tomo.py library attaches a configuration specific
c function that converts a sino index to beam properties used in
c these routines.



c Change poor man's floor function and nraypix
      subroutine root_plane(troot,zs,zb)

      Implicit Real*4(a-h,o-z)
      Parameter (alarge=100000000.0,eps=0.00000001)
c alarge is supposed to imitate infinity

      Dimension troot(2)

      
      if (abs(zs-zb).lt.eps) then
         if (zs.ge.0.) then
            troot(1)=-alarge
            troot(2)=alarge
         else
            troot(1)=1.0
            troot(2)=-1.0
         end if
      else
         root1=zs/(zs-zb)
         if (zs.lt.zb) then
            troot(1)=root1
            troot(2)=alarge
         else
            troot(1)=-alarge
            troot(2)=root1
         endif
      end if

      return
      end

      subroutine root_ellipsoid(troot,xs,ys,zs,xb,yb,zb,
     +                          ax,ay,az)

      Implicit Real*4(a-h,o-z)

      Dimension troot(2)


      dist2=(4*ax**2*ay**2*az**2*
     -    (-(az**2*(xs*yb - xb*ys)**2) +  ay**2*(az**2*(xb - xs)**2 + 
     -         ax**2*(zb - zs)**2 - (xs*zb - xb*zs)**2)
     -       + ax**2*(az**2*(yb - ys)**2 -  (ys*zb - yb*zs)**2)))/
     -  (ax**2*az**2*(yb - ys)**2 + 
     -     ay**2*(az**2*(xb - xs)**2 +  ax**2*(zb - zs)**2))**2


      avet=(ax**2*az**2*ys*(-yb + ys) + 
     -    ay**2*(az**2*xs*(-xb + xs) +  ax**2*zs*(-zb + zs)))/
     -  (ax**2*az**2*(yb - ys)**2 + 
     -    ay**2*(az**2*(xb - xs)**2 + ax**2*(zb - zs)**2) )

      if (dist2.lt.0.) then
         troot(1)=1.0
         troot(2)=-1.0
      else
         troot(1)=avet - 0.5*sqrt(dist2)
         troot(2)=avet + 0.5*sqrt(dist2)
      end if

      return
      end

      
      subroutine root_cylinder(troot,xs,ys,xb,yb,ax,ay)

      Implicit Real*4(a-h,o-z)

      Dimension troot(2)


      dist2=(4*ax**2*ay**2*(ay**2*(xb-xs)**2 + ax**2*(yb-ys)**2 -
     + (xs*yb - xb*ys)**2))/(ay**2*(xb-xs)**2 + ax**2*(yb-ys)**2)**2


      avet=(ay**2*xs*(xs-xb) + ax**2*ys*(ys-yb))/
     + (ay**2*(xb-xs)**2 + ax**2*(yb-ys)**2)

      if (dist2.lt.0.) then
         troot(1)=1.0
         troot(2)=-1.0
      else
         troot(1)=avet - 0.5*sqrt(dist2)
         troot(2)=avet + 0.5*sqrt(dist2)
      end if

      return
      end

      subroutine root_cone(troot,xs,ys,zs,xb,yb,zb,ax,ay)

      Implicit Real*4(a-h,o-z)
      Parameter (alarge=100000000.0,eps=0.00000001)

      Dimension troot(2)


      dist2=(4*ax**2*ay**2*(xs**2*(-yb**2 + ay**2*zb**2) + ax**2*(ys*zb-
     - yb*zs)**2 + 
     -   2*xb*xs*(yb*ys - ay**2*zb*zs) + xb**2*(-ys**2 + ay**2*zs**2)))/
     -(ax**2*(yb-ys)**2+ay**2*(xb**2-2*xb*xs+xs**2-ax**2*(zb-zs)**2))**2


      avet=(ax**2*ys*(ys-yb)+ay**2*(-(xb*xs)+xs**2+ax**2*(zb-zs)*zs))/
     -  (ax**2*(yb-ys)**2+ay**2*(xb**2-2*xb*xs+xs**2-ax**2*(zb-zs)**2))

      sd2=sqrt(dist2)*0.5

      if (abs(zs-zb).lt.eps) then
         plane_root= -alarge
      else
         plane_root=zs/(zs-zb)
      end if

      if (dist2.lt.0.) then
         troot(1)=1.0
         troot(2)=-1.0
      else
         tr1 = avet - sd2
         tr2 = avet + sd2
c find correspond z-coords to these roots
         zr1=zs*(1.-tr1)+zb*tr1
         zr2=zs*(1.-tr2)+zb*tr2

c Eliminate negative cone
         If ((zr1.le.0.).and.(zr2.le.0.)) then
            troot(1)=1.0
            troot(2)=-1.0
         else

c case where both roots are on pos. cone
            if ((zr1.gt.0.).and.(zr2.gt.0.)) then
               troot(1)=tr1
               troot(2)=tr2
            else
c case where only root1 is on pos. cone
               if (zr1.ge.0.) then
                  if (zs.lt.zb) then
                     troot(1)=tr1
                     troot(2)=alarge
                  else
                     troot(1)=-alarge
                     troot(2)=tr1
                  endif
               endif
c case where only root2 is on pos. cone
               if (zr2.ge.0.) then
                  if (zs.lt.zb) then
                     troot(1)=tr2
                     troot(2)=alarge
                  else
                     troot(1)=-alarge
                     troot(2)=tr2
                  endif
               endif
            end if
         end if
      end if

      return
      end



      subroutine gen_shape_proj(sinomat,indsino,
     +        ray_vectors,
     +        nr,
     +        att,
     +        nsurf,
     +        xu1,xu2,xu3,yu1,yu2,yu3,zu1,zu2,zu3,
     +        x0,y0,z0,
     +        nplane,
     +        nellipsoid,
     +        ax_ell,ay_ell,az_ell,
     +        ncyl,
     +        ax_cyl,ay_cyl,
     +        ncone,
     +        ax_cone,ay_cone)

c nsurf = nplane + nellipsoid + ncyl + ncone

      Implicit Real*4(a-h,o-z)

      Parameter (alarge=100000000.0)
c alarge is supposed to imitate infinity

      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr)
      Dimension troot(2), froot(2)
      Dimension xu1(nsurf),xu2(nsurf),xu3(nsurf)
      Dimension yu1(nsurf),yu2(nsurf),yu3(nsurf)
      Dimension zu1(nsurf),zu2(nsurf),zu3(nsurf)
      Dimension x0(nsurf),y0(nsurf),z0(nsurf)
      Dimension ax_ell(nellipsoid),ay_ell(nellipsoid),az_ell(nellipsoid)
      Dimension ax_cyl(ncyl),ay_cyl(ncyl)
      Dimension ax_cone(ncone),ay_cone(ncone)

Cf2py intent(in) nr
cf2py intent(in) ray_vectors
cf2py depend(nr)  ray_vectors
Cf2py intent(out) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino

cf2py intent(in) nsurf
cf2py intent(in) xu1, xu2, xu3, yu1, yu2, yu3, zu1, zu2, zu3, x0, y0, z0
cf2py depend(nsurf) xu1, xu2, xu3, yu1, yu2, yu3, zu1, zu2, zu3, x0, y0, z0

cf2py intent(in) nplane

cf2py intent(in) nellipsoid
cf2py intent(in) ax_ell,ay_ell,az_ell
cf2py depend(nellipsoid) ax_ell,ay_ell,az_ell

cf2py intent(in) ncyl
cf2py intent(in) ax_cyl,ay_cyl
cf2py depend(ncyl) ax_cyl,ay_cyl

cf2py intent(in) ncone
cf2py intent(in) ax_cone,ay_cone
cf2py depend(ncone) ax_cone,ay_cone


c         s = frame_vectors(19)


      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)

         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)





      if (indsino(ip).eq.1) then




c calc t scale factor
         sbdist=sqrt((xSource-xbin)**2+
     +               (ySource-ybin)**2+(zSource-zbin)**2)

c         print *,"sb_dist:",sbdist

c initialize transmission path segment to infinite line
      froot(1)=-alarge
      froot(2)=alarge

c loop through all surfaces
      isurf=0

c loop through planes

c      print*,"surface: ",isurf
      Do iplane=2,nplane
c         print*,"   plane: ",iplane
         isurf=isurf+1
c put source and bin points in current plane frame

         xt=xSource-x0(isurf)
         yt=ySource-y0(isurf)
         zt=zSource-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt

         xt=xbin-x0(isurf)
         yt=ybin-y0(isurf)
         zt=zbin-z0(isurf)

         xb=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         yb=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zb=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt


         Call root_plane(troot,zs,zb)


         froot(1)=max(froot(1),troot(1))
         froot(2)=min(froot(2),troot(2))

         if (froot(1).ge.froot(2)) then
            tlen= 0.
            goto 111
         end if

      End Do

c loop through ellipsoids

      Do iell=2,nellipsoid
         isurf=isurf+1
c put source and bin points in current plane frame


         xt=xSource-x0(isurf)
         yt=ySource-y0(isurf)
         zt=zSource-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt

         xt=xbin-x0(isurf)
         yt=ybin-y0(isurf)
         zt=zbin-z0(isurf)

         xb=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         yb=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zb=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt



         Call root_ellipsoid(troot,xs,ys,zs,xb,yb,zb,
     +      ax_ell(iell),ay_ell(iell),az_ell(iell))

         froot(1)=max(froot(1),troot(1))
         froot(2)=min(froot(2),troot(2))
         if (froot(1).ge.froot(2)) then
            tlen= 0.
            goto 111
         end if
      End Do

c loop through cylinders

      Do icyl=2,ncyl
         isurf=isurf+1
c put source and bin points in current plane frame

         xt=xSource-x0(isurf)
         yt=ySource-y0(isurf)
         zt=zSource-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt

         xt=xbin-x0(isurf)
         yt=ybin-y0(isurf)
         zt=zbin-z0(isurf)

         xb=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         yb=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zb=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt



         Call root_cylinder(troot,xs,ys,xb,yb,
     +      ax_cyl(icyl),ay_cyl(icyl))

         froot(1)=max(froot(1),troot(1))
         froot(2)=min(froot(2),troot(2))
         if (froot(1).ge.froot(2)) then
            tlen= 0.
            goto 111
         end if

      End Do

c loop through cones

      Do icone=2,ncone
         isurf=isurf+1
c put source and bin points in current plane frame

         xt=xSource-x0(isurf)
         yt=ySource-y0(isurf)
         zt=zSource-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt

         xt=xbin-x0(isurf)
         yt=ybin-y0(isurf)
         zt=zbin-z0(isurf)

         xb=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         yb=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zb=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt



         Call root_cone(troot,xs,ys,zs,xb,yb,zb,
     +      ax_cone(icone),ay_cone(icone))

         froot(1)=max(froot(1),troot(1))
         froot(2)=min(froot(2),troot(2))
         if (froot(1).ge.froot(2)) then
            tlen= 0.
            goto 111
         end if

      End Do


      tlen = froot(2)-froot(1)



 111  Continue      

      sinomat(ip)=sbdist*tlen*att

      endif
      enddo


      return
      end



      subroutine ellipsoid_proj(sinomat,indsino,
     +        ray_vectors,
     +        nr,
     +        ax,ay,az,x0,y0,z0,alpha,beta,att)

      Implicit Real*4(a-h,o-z)


      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr)

Cf2py intent(in) nr
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(out) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino


      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)

         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)



      if (indsino(ip).eq.1) then



c put source and detector points in ellipsoid frame
      xspp=(-x0 + xSource)*Cos(beta) - (y0 - ySource)*Sin(beta)
      yspp=(-y0 + ySource)*Cos(beta) + (x0 - xSource)*Sin(beta)
      zspp=zSource-z0
      xsp=xspp*Cos(alpha) + zspp*Sin(alpha)
      ysp=yspp
      zsp=-xspp*Sin(alpha) + zspp*Cos(alpha)

      xbpp=(-x0 + xbin)*Cos(beta) - (y0 - ybin)*Sin(beta)
      ybpp=(-y0 + ybin)*Cos(beta) + (x0 - xbin)*Sin(beta)
      zbpp=zbin-z0
      xbp=xbpp*Cos(alpha) + zbpp*Sin(alpha)
      ybp=ybpp
      zbp=-xbpp*Sin(alpha) + zbpp*Cos(alpha)


      pathLength2= (4*ax**2*ay**2*az**2*((xbp - xsp)**2 +
     - (ybp - ysp)**2 + (zbp - zsp)**2)*
     -     (-(az**2*(xsp*ybp - xbp*ysp)**2) + 
     -   ay**2*(az**2*(xbp - xsp)**2 + ax**2*(zbp -
     -          zsp)**2 - (xsp*zbp - xbp*zsp)**2) + 
     -       ax**2*(az**2*(ybp - ysp)**2 - (ysp*zbp - ybp*zsp)**2)))/
     -   (ax**2*az**2*(ybp - ysp)**2 +
     -  ay**2*(az**2*(xbp - xsp)**2 + ax**2*(zbp - zsp)**2))**2


      if (pathLength2.lt.0.) pathLength2=0.

      sinomat(ip)=sqrt(pathLength2)*att

      endif
      enddo

      return
      end







      subroutine sinoproj_linear(sinomat,indsino,
     +       ray_vectors,
     +       nr,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)


      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr),smat(nx,ny,nz)

Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(out) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat


      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)

         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)

c compute source-detector distance
         sd_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)


      if (indsino(ip).eq.1) then


      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

c scan through x direction
      if ((xad.gt.yad).and.(xad.gt.zad)) then
            yox=ydiff/xdiff
            zox=zdiff/xdiff
            travVoxlen=dx*sqrt(1.0+yox**2+zox**2)
            total=0.
            do ix=0, nx-1
               x=xl+dx*(ix+0.5)
               y=ySource+yox*(x-xSource)
               z=zSource+zox*(x-xSource)
               iy=floor((y-y0)/dy - 0.5)
               iz=floor((z-z0)/dz - 0.5)
               if (((iy.ge.0).and.(iy.lt.ny-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then
                  wy2=(y-(dy*(iy+0.5)+y0))/dy
                  wy1=1.0-wy2
                  wz2=(z-(dz*(iz+0.5)+z0))/dz
                  wz1=1.0-wz2
                  val0=smat(ix+1,iy+1,iz+1)
                  valy1=smat(ix+1,iy+2,iz+1)
                  valz1=smat(ix+1,iy+1,iz+2)
                  val2=smat(ix+1,iy+2,iz+2)
                  vi1=wy1*val0+wy2*valy1
                  vi2=wy1*valz1+wy2*val2
                  val=wz1*vi1+wz2*vi2
                  total=total+travVoxlen*val
               endif
            enddo
            goto 10
      endif

c scan through y direction

         if (yad.gt.zad) then
            xoy=xdiff/ydiff
            zoy=zdiff/ydiff
            travVoxlen=dy*sqrt(1.0+xoy**2+zoy**2)
            total=0.
            do iy=0, ny-1
               y=yl+dy*(iy+0.5)
               x=xSource+xoy*(y-ySource)
               z=zSource+zoy*(y-ySource)
               ix=floor((x-x0)/dx - 0.5)
               iz=floor((z-z0)/dz - 0.5)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then
                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
                  wz2=(z-(dz*(iz+0.5)+z0))/dz
                  wz1=1.0-wz2
                  val0=smat(ix+1,iy+1,iz+1)
                  valx1=smat(ix+2,iy+1,iz+1)
                  valz1=smat(ix+1,iy+1,iz+2)
                  val2=smat(ix+2,iy+1,iz+2)
                  vi1=wx1*val0+wx2*valx1
                  vi2=wx1*valz1+wx2*val2
                  val=wz1*vi1+wz2*vi2
                  total=total+travVoxlen*val
               endif
            enddo
            goto 10
         endif
c scan through z direction

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
            total=0.
            do iz=0, nz-1
               z=zl+dz*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0)/dx - 0.5)
               iy=floor((y-y0)/dy - 0.5)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iy.ge.0).and.(iy.lt.ny-1))) then
                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
                  wy2=(y-(dy*(iy+0.5)+y0))/dy
                  wy1=1.0-wy2
                  val0=smat(ix+1,iy+1,iz+1)
                  valx1=smat(ix+2,iy+1,iz+1)
                  valy1=smat(ix+1,iy+2,iz+1)
                  val2=smat(ix+2,iy+2,iz+1)
                  vi1=wx1*val0+wx2*valx1
                  vi2=wx1*valy1+wx2*val2
                  val=wy1*vi1+wy2*vi2
                  total=total+travVoxlen*val
               endif
            enddo

 10   continue 
      sinomat(ip)=total

      endif
      end do

      return
      end





      subroutine sinoproj_nearest(sinomat,indsino,
     +       ray_vectors,
     +       nr,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)


      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr),smat(nx,ny,nz)

Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(out) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat


      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)

         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)

c compute source-detector distance
         sd_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)


      if (indsino(ip).eq.1) then


      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

c scan through x direction
      if ((xad.gt.yad).and.(xad.gt.zad)) then
            yox=ydiff/xdiff
            zox=zdiff/xdiff
            travVoxlen=dx*sqrt(1.0+yox**2+zox**2)
            total=0.
            do ix=0, nx-1
               x=xl+dx*(ix+0.5)
               y=ySource+yox*(x-xSource)
               z=zSource+zox*(x-xSource)
               iy=floor((y-y0)/dy)
               iz=floor((z-z0)/dz)
               if (((iy.ge.0).and.(iy.le.ny-1)).and.
     +             ((iz.ge.0).and.(iz.le.nz-1))) then
                  val0=smat(ix+1,iy+1,iz+1)
                  total=total+travVoxlen*val0
               endif
            enddo
            goto 10
      endif

c scan through y direction

         if (yad.gt.zad) then
            xoy=xdiff/ydiff
            zoy=zdiff/ydiff
            travVoxlen=dy*sqrt(1.0+xoy**2+zoy**2)
            total=0.
            do iy=0, ny-1
               y=yl+dy*(iy+0.5)
               x=xSource+xoy*(y-ySource)
               z=zSource+zoy*(y-ySource)
               ix=floor((x-x0)/dx)
               iz=floor((z-z0)/dz)
               if (((ix.ge.0).and.(ix.le.nx-1)).and.
     +             ((iz.ge.0).and.(iz.le.nz-1))) then
                  val0=smat(ix+1,iy+1,iz+1)
                  total=total+travVoxlen*val0
               endif
            enddo
            goto 10
         endif
c scan through z direction

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
            total=0.
            do iz=0, nz-1
               z=zl+dz*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0)/dx)
               iy=floor((y-y0)/dy)
               if (((ix.ge.0).and.(ix.le.nx-1)).and.
     +             ((iy.ge.0).and.(iy.le.ny-1))) then
                  val0=smat(ix+1,iy+1,iz+1)
                  total=total+travVoxlen*val0
               endif
            enddo

 10   continue 
      sinomat(ip)=total

      endif
      end do

      return
      end




      subroutine sinobackproject_nearest(sinomat,indsino,
     +       ray_vectors,
     +       nr,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000)

      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr),smat(nx,ny,nz)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)

Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(in) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino
Cf2py intent(out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat


      Do ix = 1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               smat(ix,iy,iz)=0.
            End Do
         End Do
      End Do

      
      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)


         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)

c compute source-detector distance
         sd_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)


      if (indsino(ip).eq.1) then



      dataval=sinomat(ip)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

c scan through x direction
      if ((xad.gt.yad).and.(xad.gt.zad)) then
            yox=ydiff/xdiff
            zox=zdiff/xdiff
            travVoxlen=dx*sqrt(1.0+yox**2+zox**2)
            irp=0
            do ix=0, nx-1
               x=xl+dx*(ix+0.5)
               y=ySource+yox*(x-xSource)
               z=zSource+zox*(x-xSource)
               iy=floor((y-y0)/dy)
               iz=floor((z-z0)/dz)
               if (((iy.ge.0).and.(iy.lt.ny-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then


                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
               end if


            enddo
            goto 10
      endif

c scan through y direction

         if (yad.gt.zad) then
            xoy=xdiff/ydiff
            zoy=zdiff/ydiff
            travVoxlen=dy*sqrt(1.0+xoy**2+zoy**2)
            irp=0
            do iy=0, ny-1
               y=yl+dy*(iy+0.5)
               x=xSource+xoy*(y-ySource)
               z=zSource+zoy*(y-ySource)
               ix=floor((x-x0)/dx)
               iz=floor((z-z0)/dz)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then




                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
               end if

            enddo
            goto 10
         endif
c scan through z direction

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
            irp=0
            do iz=0, nz-1
               z=zl+dz*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0)/dx)
               iy=floor((y-y0)/dy)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iy.ge.0).and.(iy.lt.ny-1))) then



                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
               end if
            enddo

 10   continue 

      nrp=irp
      do irp=1,nrp
         smat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +      smat(iraypix(irp),jraypix(irp),kraypix(irp))+
     +      dataval*raypix(irp)
      enddo


      endif
      end do


      return
      end






      subroutine artstep_linear(reg_fac,sinomat,indsino,
     +       ray_vectors,
     +       nr,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000)

      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr),smat(nx,ny,nz)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)

Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(in) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      
      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)


         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)

c compute source-detector distance
         sd_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)



      if (indsino(ip).eq.1) then


      dataval=sinomat(ip)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

c scan through x direction
      if ((xad.gt.yad).and.(xad.gt.zad)) then
            yox=ydiff/xdiff
            zox=zdiff/xdiff
            travVoxlen=dx*sqrt(1.0+yox**2+zox**2)
            irp=0
            total=0.
            raysum=0.
            do ix=0, nx-1
               x=xl+dx*(ix+0.5)
               y=ySource+yox*(x-xSource)
               z=zSource+zox*(x-xSource)
               iy=floor((y-y0)/dy - 0.5)
               iz=floor((z-z0)/dz - 0.5)
               if (((iy.ge.0).and.(iy.lt.ny-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then
                  wy2=(y-(dy*(iy+0.5)+y0))/dy
                  wy1=1.0-wy2
                  wz2=(z-(dz*(iz+0.5)+z0))/dz
                  wz1=1.0-wz2
                  val0=smat(ix+1,iy+1,iz+1)
                  valy1=smat(ix+1,iy+2,iz+1)
                  valz1=smat(ix+1,iy+1,iz+2)
                  val2=smat(ix+1,iy+2,iz+2)
                  vi1=wy1*val0+wy2*valy1
                  vi2=wy1*valz1+wy2*val2
                  val=wz1*vi1+wz2*vi2

                  raysum=raysum+travVoxlen*val
c                  total=total+travVoxlen*travVoxlen


                  irp=irp+1
                  raypix(irp)=travVoxlen*wy1*wz1
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wy2*wz1
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wy1*wz2
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+2
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wy2*wz2
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+2
                  total = total + raypix(irp)**2
               end if


            enddo
            goto 10
      endif

c scan through y direction

         if (yad.gt.zad) then
            xoy=xdiff/ydiff
            zoy=zdiff/ydiff
            travVoxlen=dy*sqrt(1.0+xoy**2+zoy**2)
            irp=0
            total=0.
            raysum=0.
            do iy=0, ny-1
               y=yl+dy*(iy+0.5)
               x=xSource+xoy*(y-ySource)
               z=zSource+zoy*(y-ySource)
               ix=floor((x-x0)/dx - 0.5)
               iz=floor((z-z0)/dz - 0.5)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then

                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
                  wz2=(z-(dz*(iz+0.5)+z0))/dz
                  wz1=1.0-wz2
                  val0=smat(ix+1,iy+1,iz+1)
                  valx1=smat(ix+2,iy+1,iz+1)
                  valz1=smat(ix+1,iy+1,iz+2)
                  val2=smat(ix+2,iy+1,iz+2)
                  vi1=wx1*val0+wx2*valx1
                  vi2=wx1*valz1+wx2*val2
                  val=wz1*vi1+wz2*vi2



                  raysum=raysum+travVoxlen*val
c                  total=total+travVoxlen*travVoxlen


                  irp=irp+1
                  raypix(irp)=travVoxlen*wx1*wz1
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wz1
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx1*wz2
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+2
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wz2
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+2
                  total = total + raypix(irp)**2

               end if

            enddo
            goto 10
         endif
c scan through z direction

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
            irp=0
            total=0.
            raysum=0.
            do iz=0, nz-1
               z=zl+dz*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0)/dx - 0.5)
               iy=floor((y-y0)/dy - 0.5)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iy.ge.0).and.(iy.lt.ny-1))) then

                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
                  wy2=(y-(dy*(iy+0.5)+y0))/dy
                  wy1=1.0-wy2
                  val0=smat(ix+1,iy+1,iz+1)
                  valx1=smat(ix+2,iy+1,iz+1)
                  valy1=smat(ix+1,iy+2,iz+1)
                  val2=smat(ix+2,iy+2,iz+1)
                  vi1=wx1*val0+wx2*valx1
                  vi2=wx1*valy1+wx2*val2
                  val=wy1*vi1+wy2*vi2


                  raysum=raysum+travVoxlen*val
c                  total=total+travVoxlen*travVoxlen


                  irp=irp+1
                  raypix(irp)=travVoxlen*wx1*wy1
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wy1
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx1*wy2
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wy2
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total + raypix(irp)**2

               end if
            enddo

 10   continue 

      nrp=irp
      do irp=1,nrp
         smat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +      smat(iraypix(irp),jraypix(irp),kraypix(irp))+
     +      reg_fac*(dataval-raysum)*raypix(irp)/total
      enddo


      endif
      end do


      return
      end




      subroutine artstep_nearest(reg_fac,sinomat,indsino,
     +       ray_vectors,
     +       nr,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000)

      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr),smat(nx,ny,nz)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)

Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(in) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      
      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)


         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)

c compute source-detector distance
         sd_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)



      if (indsino(ip).eq.1) then


      dataval=sinomat(ip)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

c scan through x direction
      if ((xad.gt.yad).and.(xad.gt.zad)) then
            yox=ydiff/xdiff
            zox=zdiff/xdiff
            travVoxlen=dx*sqrt(1.0+yox**2+zox**2)
            irp=0
            total=0.
            raysum=0.
            do ix=0, nx-1
               x=xl+dx*(ix+0.5)
               y=ySource+yox*(x-xSource)
               z=zSource+zox*(x-xSource)
               iy=floor((y-y0)/dy)
               iz=floor((z-z0)/dz)
               if (((iy.ge.0).and.(iy.lt.ny-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then

                  val0=smat(ix+1,iy+1,iz+1)
                  raysum=raysum+travVoxlen*val0
                  total=total+travVoxlen*travVoxlen

                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
               end if


            enddo
            goto 10
      endif

c scan through y direction

         if (yad.gt.zad) then
            xoy=xdiff/ydiff
            zoy=zdiff/ydiff
            travVoxlen=dy*sqrt(1.0+xoy**2+zoy**2)
            irp=0
            total=0.
            raysum=0.
            do iy=0, ny-1
               y=yl+dy*(iy+0.5)
               x=xSource+xoy*(y-ySource)
               z=zSource+zoy*(y-ySource)
               ix=floor((x-x0)/dx)
               iz=floor((z-z0)/dz)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then


                  val0=smat(ix+1,iy+1,iz+1)

                  raysum=raysum+travVoxlen*val0
                  total=total+travVoxlen*travVoxlen

                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
               end if

            enddo
            goto 10
         endif
c scan through z direction

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
            irp=0
            total=0.
            raysum=0.
            do iz=0, nz-1
               z=zl+dz*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0)/dx)
               iy=floor((y-y0)/dy)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iy.ge.0).and.(iy.lt.ny-1))) then

                  val0=smat(ix+1,iy+1,iz+1)

                  raysum=raysum+travVoxlen*val0
                  total=total+travVoxlen*travVoxlen

                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
               end if
            enddo

 10   continue 

      nrp=irp
      do irp=1,nrp
         smat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +      smat(iraypix(irp),jraypix(irp),kraypix(irp))+
     +      reg_fac*(dataval-raysum)*raypix(irp)/total
      enddo


      endif
      end do


      return
      end




      subroutine emstep_nearest(sinomat,indsino,
     +       ray_vectors,
     +       nr,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000)

      Dimension ray_vectors(nr,6)
      Dimension indsino(nr)
      Dimension sinomat(nr),smat(nx,ny,nz)
      Dimension wsmat(nx,ny,nz), onemat(nx,ny,nz)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)

Cf2py intent(in) nr
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) ray_vectors
cf2py depend(nr) ray_vectors
Cf2py intent(in) sinomat
cf2py depend(nr) sinomat
Cf2py intent(in) indsino
cf2py depend(nr) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      Do ip=1,nx
         Do jp=1,ny
            Do kp=1,nz
               wsmat(ip,jp,kp) = 0.
               onemat(ip,jp,kp) = 0.
            End Do
         End Do
      End Do


      Do ip=1,nr

         xSource=ray_vectors(ip,1)
         ySource=ray_vectors(ip,2)
         zSource=ray_vectors(ip,3)


         xbin =ray_vectors(ip,4)
         ybin =ray_vectors(ip,5)
         zbin =ray_vectors(ip,6)

c compute source-detector distance
         sd_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)



      if (indsino(ip).eq.1) then


      dataval=sinomat(ip)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

c scan through x direction
      if ((xad.gt.yad).and.(xad.gt.zad)) then
            yox=ydiff/xdiff
            zox=zdiff/xdiff
            travVoxlen=dx*sqrt(1.0+yox**2+zox**2)
            irp=0
            total=0.
            raysum=0.
            do ix=0, nx-1
               x=xl+dx*(ix+0.5)
               y=ySource+yox*(x-xSource)
               z=zSource+zox*(x-xSource)
               iy=floor((y-y0)/dy)
               iz=floor((z-z0)/dz)
               if (((iy.ge.0).and.(iy.lt.ny-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then
                  val0=smat(ix+1,iy+1,iz+1)
                  raysum=raysum+travVoxlen*val0
                  total=total+travVoxlen*travVoxlen

                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1

               endif
            enddo
            goto 10
      endif

c scan through y direction

         if (yad.gt.zad) then
            xoy=xdiff/ydiff
            zoy=zdiff/ydiff
            travVoxlen=dy*sqrt(1.0+xoy**2+zoy**2)
            irp=0
            total=0.
            raysum=0.
            do iy=0, ny-1
               y=yl+dy*(iy+0.5)
               x=xSource+xoy*(y-ySource)
               z=zSource+zoy*(y-ySource)
               ix=floor((x-x0)/dx)
               iz=floor((z-z0)/dz)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iz.ge.0).and.(iz.lt.nz-1))) then
                  val0=smat(ix+1,iy+1,iz+1)

                  raysum=raysum+travVoxlen*val0
                  total=total+travVoxlen*travVoxlen

                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1

               endif
            enddo
            goto 10
         endif
c scan through z direction

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
            irp=0
            total=0.
            raysum=0.
            do iz=0, nz-1
               z=zl+dz*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0)/dx)
               iy=floor((y-y0)/dy)
               if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +             ((iy.ge.0).and.(iy.lt.ny-1))) then
                  val0=smat(ix+1,iy+1,iz+1)

                  raysum=raysum+travVoxlen*val0
                  total=total+travVoxlen*travVoxlen

                  irp=irp+1
                  raypix(irp)=travVoxlen
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1


               endif
            enddo

 10   continue 

      if (raysum.lt.0.00000001) raysum=1.0
      nrp=irp
      do irp=1,nrp
         wsmat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +    wsmat(iraypix(irp),jraypix(irp),kraypix(irp)) +
     +      (dataval/raysum)*raypix(irp)
         onemat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +    onemat(iraypix(irp),jraypix(irp),kraypix(irp)) +
     +      raypix(irp)
      enddo

      endif
      end do

      Do ip=1,nx
         Do jp=1,ny
            Do kp=1,nz
            if (onemat(ip,jp,kp).lt.0.00000001) onemat(ip,jp,kp)=1.0
            smat(ip,jp,kp) = smat(ip,jp,kp)*wsmat(ip,jp,kp)/
     +                                      onemat(ip,jp,kp)
            End Do
         End Do
      End Do


      return
      end

