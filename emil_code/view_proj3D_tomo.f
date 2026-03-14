c projection and step routines for tomosynthesis

c this projection routine goes through image array only in z-direction
      subroutine sinoprojz(sinomat,indsino,
     +       frame_vectors,
     +       ns,nu,nv,du,dv,u0,v0,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)


      Dimension frame_vectors(ns,18)
      Dimension indsino(ns,nu,nv)
      Dimension sinomat(ns,nu,nv),smat(nx,ny,nz)

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(out) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
cf2py depend(nv) indsino
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat


      Do ip=1,ns

         xSource=frame_vectors(ip,1)
         ySource=frame_vectors(ip,2)
         zSource=frame_vectors(ip,3)

         xTrajTan=frame_vectors(ip,4)
         yTrajTan=frame_vectors(ip,5)
         zTrajTan=frame_vectors(ip,6)

         xDetCenter=frame_vectors(ip,7)
         yDetCenter=frame_vectors(ip,8)
         zDetCenter=frame_vectors(ip,9)

         eux=frame_vectors(ip,10)
         euy=frame_vectors(ip,11)
         euz=frame_vectors(ip,12)

         evx=frame_vectors(ip,13)
         evy=frame_vectors(ip,14)
         evz=frame_vectors(ip,15)

         ewx=frame_vectors(ip,16)
         ewy=frame_vectors(ip,17)
         ewz=frame_vectors(ip,18)


c Calculate detector true coordinate center

         udr =(xDetCenter-xSource)*eux+
     +        (yDetCenter-ySource)*euy+
     +        (zDetCenter-zSource)*euz
         vdr =(xDetCenter-xSource)*evx+
     +        (yDetCenter-ySource)*evy+
     +        (zDetCenter-zSource)*evz
         tcx = xDetCenter - eux*udr - evx*vdr 
         tcy = yDetCenter - euy*udr - evy*vdr 
         tcz = zDetCenter - euz*udr - evz*vdr 

c compute source-detector distance
         sd_dist=Sqrt((xSource-tcx)**2+
     +                (ySource-tcy)**2+
     +                (zSource-tcz)**2)

      Do jp=1,nu
      Do kp=1,nv

         u = u0+(jp-0.5)*du 
         v = v0+(kp-0.5)*dv 

      if (indsino(ip,jp,kp).eq.1) then

c find bin location in space-fixed coords.
c         xbin = tcx + eux*u+evx*v
c         ybin = tcy + euy*u+evy*v
c         zbin = tcz + euz*u+evz*v

         xbin = xDetCenter + eux*u+evx*v
         ybin = yDetCenter + euy*u+evy*v
         zbin = zDetCenter + euz*u+evz*v

      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

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
c                  wx2=(x-(dx*ix+x0))/dx
                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
c                  wy2=(y-(dy*iy+y0))/dy
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

      sinomat(ip,jp,kp)=total

      endif
      end do
      end do
      end do

      return
      end





      subroutine artstepz(reg_fac,sinomat,indsino,
     +       frame_vectors,
     +       ns,nu,nv,du,dv,u0,v0,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000)

      Dimension frame_vectors(ns,18)
      Dimension indsino(ns,nu,nv)
      Dimension sinomat(ns,nu,nv),smat(nx,ny,nz)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
cf2py depend(nv) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      
      Do ip=1,ns

         xSource=frame_vectors(ip,1)
         ySource=frame_vectors(ip,2)
         zSource=frame_vectors(ip,3)

         xTrajTan=frame_vectors(ip,4)
         yTrajTan=frame_vectors(ip,5)
         zTrajTan=frame_vectors(ip,6)

         xDetCenter=frame_vectors(ip,7)
         yDetCenter=frame_vectors(ip,8)
         zDetCenter=frame_vectors(ip,9)

         eux=frame_vectors(ip,10)
         euy=frame_vectors(ip,11)
         euz=frame_vectors(ip,12)

         evx=frame_vectors(ip,13)
         evy=frame_vectors(ip,14)
         evz=frame_vectors(ip,15)

         ewx=frame_vectors(ip,16)
         ewy=frame_vectors(ip,17)
         ewz=frame_vectors(ip,18)



c Calculate detector true coordinate center

         udr =(xDetCenter-xSource)*eux+
     +        (yDetCenter-ySource)*euy+
     +        (zDetCenter-zSource)*euz
         vdr =(xDetCenter-xSource)*evx+
     +        (yDetCenter-ySource)*evy+
     +        (zDetCenter-zSource)*evz
         tcx = xDetCenter - eux*udr - evx*vdr 
         tcy = yDetCenter - euy*udr - evy*vdr 
         tcz = zDetCenter - euz*udr - evz*vdr 

c compute source-detector distance
         sd_dist=Sqrt((xSource-tcx)**2+
     +                (ySource-tcy)**2+
     +                (zSource-tcz)**2)

      Do jp=1,nu
      Do kp=1,nv

         u = u0+(jp-0.5)*du 
         v = v0+(kp-0.5)*dv 

      if (indsino(ip,jp,kp).eq.1) then

c find bin location in space-fixed coords.
c         xbin = tcx + eux*u+evx*v
c         ybin = tcy + euy*u+evy*v
c         zbin = tcz + euz*u+evz*v

         xbin = xDetCenter + eux*u+evx*v
         ybin = yDetCenter + euy*u+evy*v
         zbin = zDetCenter + euz*u+evz*v


      dataval=sinomat(ip,jp,kp)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz

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
c                  wx2=(x-(dx*ix+x0))/dx
                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
c                  wy2=(y-(dy*iy+y0))/dy
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
                  total = total+raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wy1
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total+raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx1*wy2
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total+raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wy2
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total+raypix(irp)**2

               endif
            enddo


      nrp=irp
      do irp=1,nrp
         smat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +      smat(iraypix(irp),jraypix(irp),kraypix(irp))+
     +      reg_fac*(dataval-raysum)*raypix(irp)/total
      enddo


      endif
      end do
      end do
      end do


      return
      end




      subroutine artstepz2vol(reg_fac,sinomat,indsino,
     +       frame_vectors,
     +       ns,nu,nv,du,dv,u0,v0,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz,
     +       smata,dxa,dya,dza,x0a,y0a,z0a,nxa,nya,nza)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000)

      Dimension frame_vectors(ns,18)
      Dimension indsino(ns,nu,nv)
      Dimension sinomat(ns,nu,nv),smat(nx,ny,nz),smata(nxa,nya,nza)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)
      Dimension raypixa(nraypix)
      Dimension iraypixa(nraypix),jraypixa(nraypix),kraypixa(nraypix)

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
cf2py depend(nv) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat
Cf2py intent(in,out) smata
cf2py depend(nxa) smata
cf2py depend(nya) smata
cf2py depend(nza) smata

      
      Do ip=1,ns

         xSource=frame_vectors(ip,1)
         ySource=frame_vectors(ip,2)
         zSource=frame_vectors(ip,3)

         xTrajTan=frame_vectors(ip,4)
         yTrajTan=frame_vectors(ip,5)
         zTrajTan=frame_vectors(ip,6)

         xDetCenter=frame_vectors(ip,7)
         yDetCenter=frame_vectors(ip,8)
         zDetCenter=frame_vectors(ip,9)

         eux=frame_vectors(ip,10)
         euy=frame_vectors(ip,11)
         euz=frame_vectors(ip,12)

         evx=frame_vectors(ip,13)
         evy=frame_vectors(ip,14)
         evz=frame_vectors(ip,15)

         ewx=frame_vectors(ip,16)
         ewy=frame_vectors(ip,17)
         ewz=frame_vectors(ip,18)



c Calculate detector true coordinate center

         udr =(xDetCenter-xSource)*eux+
     +        (yDetCenter-ySource)*euy+
     +        (zDetCenter-zSource)*euz
         vdr =(xDetCenter-xSource)*evx+
     +        (yDetCenter-ySource)*evy+
     +        (zDetCenter-zSource)*evz
         tcx = xDetCenter - eux*udr - evx*vdr 
         tcy = yDetCenter - euy*udr - evy*vdr 
         tcz = zDetCenter - euz*udr - evz*vdr 

c compute source-detector distance
         sd_dist=Sqrt((xSource-tcx)**2+
     +                (ySource-tcy)**2+
     +                (zSource-tcz)**2)

      Do jp=1,nu
      Do kp=1,nv

         u = u0+(jp-0.5)*du 
         v = v0+(kp-0.5)*dv 

      if (indsino(ip,jp,kp).eq.1) then


         xbin = xDetCenter + eux*u+evx*v
         ybin = yDetCenter + euy*u+evy*v
         zbin = zDetCenter + euz*u+evz*v


      dataval=sinomat(ip,jp,kp)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource

c scan through z direction (main volume)

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
c                  wx2=(x-(dx*ix+x0))/dx
                  wx2=(x-(dx*(ix+0.5)+x0))/dx
                  wx1=1.0-wx2
c                  wy2=(y-(dy*iy+y0))/dy
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
                  total = total+raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wy1
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+1
                  kraypix(irp)=iz+1
                  total = total+raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx1*wy2
                  iraypix(irp)=ix+1
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total+raypix(irp)**2

                  irp=irp+1
                  raypix(irp)=travVoxlen*wx2*wy2
                  iraypix(irp)=ix+2
                  jraypix(irp)=iy+2
                  kraypix(irp)=iz+1
                  total = total+raypix(irp)**2

               endif
            enddo

      nrp=irp
c scan through z direction (auxiliary volume)

            xoz=xdiff/zdiff
            yoz=ydiff/zdiff
            travVoxlen=dza*sqrt(1.0+xoz**2+yoz**2)
            irpa=0
            do iz=0, nza-1
               z=z0a+dza*(iz+0.5)
               x=xSource+xoz*(z-zSource)
               y=ySource+yoz*(z-zSource)
               ix=floor((x-x0a)/dxa-0.5)
               iy=floor((y-y0a)/dya-0.5)
               if (((ix.ge.0).and.(ix.lt.nxa-1)).and.
     +             ((iy.ge.0).and.(iy.lt.nya-1))) then
                  wx2=(x-(dxa*(ix+0.5)+x0a))/dxa
                  wx1=1.0-wx2
                  wy2=(y-(dya*(iy+0.5)+y0a))/dya
                  wy1=1.0-wy2
                  val0=smata(ix+1,iy+1,iz+1)
                  valx1=smata(ix+2,iy+1,iz+1)
                  valy1=smata(ix+1,iy+2,iz+1)
                  val2=smata(ix+2,iy+2,iz+1)
                  vi1=wx1*val0+wx2*valx1
                  vi2=wx1*valy1+wx2*val2
                  val=wy1*vi1+wy2*vi2

                  raysum=raysum+travVoxlen*val
c                  total=total+travVoxlen*travVoxlen

                  irpa=irpa+1
                  raypixa(irpa)=travVoxlen*wx1*wy1
                  iraypixa(irpa)=ix+1
                  jraypixa(irpa)=iy+1
                  kraypixa(irpa)=iz+1
                  total = total+raypixa(irpa)**2

                  irpa=irpa+1
                  raypixa(irpa)=travVoxlen*wx2*wy1
                  iraypixa(irpa)=ix+2
                  jraypixa(irpa)=iy+1
                  kraypixa(irpa)=iz+1
                  total = total+raypixa(irpa)**2

                  irpa=irpa+1
                  raypixa(irpa)=travVoxlen*wx1*wy2
                  iraypixa(irpa)=ix+1
                  jraypixa(irpa)=iy+2
                  kraypixa(irpa)=iz+1
                  total = total+raypixa(irpa)**2

                  irpa=irpa+1
                  raypixa(irpa)=travVoxlen*wx2*wy2
                  iraypixa(irpa)=ix+2
                  jraypixa(irpa)=iy+2
                  kraypixa(irpa)=iz+1
                  total = total+raypixa(irpa)**2

               endif
            enddo

      nrpa=irpa

      do irp=1,nrp
         smat(iraypix(irp),jraypix(irp),kraypix(irp))=
     +      smat(iraypix(irp),jraypix(irp),kraypix(irp))+
     +      reg_fac*(dataval-raysum)*raypix(irp)/total
      enddo

      do irpa=1,nrpa
         smata(iraypixa(irpa),jraypixa(irpa),kraypixa(irpa))=
     +      smata(iraypixa(irpa),jraypixa(irpa),kraypixa(irpa))+
     +      reg_fac*(dataval-raysum)*raypixa(irpa)/total
      enddo


      endif
      end do
      end do
      end do


      return
      end


      subroutine emstep_linearz(sinomat,indsino,
     +       frame_vectors,
     +       ns,nu,nv,du,dv,u0,v0,
     +       smat,dx,dy,dz,x0,y0,z0,nx,ny,nz)

      Implicit Real*4(a-h,o-z)

      parameter (nraypix=5000,small = 0.00001)

      Dimension frame_vectors(ns,18)
      Dimension indsino(ns,nu,nv)
      Dimension sinomat(ns,nu,nv),smat(nx,ny,nz)
      Dimension wsmat(nx,ny,nz), onemat(nx,ny,nz)
      Dimension raypix(nraypix)
      Dimension iraypix(nraypix),jraypix(nraypix),kraypix(nraypix)

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
Cf2py intent(in) nx
Cf2py intent(in) ny
Cf2py intent(in) nz
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
cf2py depend(nv) indsino
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


      Do ip=1,ns

         xSource=frame_vectors(ip,1)
         ySource=frame_vectors(ip,2)
         zSource=frame_vectors(ip,3)

         xTrajTan=frame_vectors(ip,4)
         yTrajTan=frame_vectors(ip,5)
         zTrajTan=frame_vectors(ip,6)

         xDetCenter=frame_vectors(ip,7)
         yDetCenter=frame_vectors(ip,8)
         zDetCenter=frame_vectors(ip,9)

         eux=frame_vectors(ip,10)
         euy=frame_vectors(ip,11)
         euz=frame_vectors(ip,12)

         evx=frame_vectors(ip,13)
         evy=frame_vectors(ip,14)
         evz=frame_vectors(ip,15)

         ewx=frame_vectors(ip,16)
         ewy=frame_vectors(ip,17)
         ewz=frame_vectors(ip,18)



c Calculate detector true coordinate center

         udr =(xDetCenter-xSource)*eux+
     +        (yDetCenter-ySource)*euy+
     +        (zDetCenter-zSource)*euz
         vdr =(xDetCenter-xSource)*evx+
     +        (yDetCenter-ySource)*evy+
     +        (zDetCenter-zSource)*evz
         tcx = xDetCenter - eux*udr - evx*vdr 
         tcy = yDetCenter - euy*udr - evy*vdr 
         tcz = zDetCenter - euz*udr - evz*vdr 

c compute source-detector distance
         sd_dist=Sqrt((xSource-tcx)**2+
     +                (ySource-tcy)**2+
     +                (zSource-tcz)**2)

      Do jp=1,nu
      Do kp=1,nv

         u = u0+(jp-0.5)*du 
         v = v0+(kp-0.5)*dv 

      if (indsino(ip,jp,kp).eq.1) then

c find bin location in space-fixed coords.
c         xbin = tcx + eux*u+evx*v
c         ybin = tcy + euy*u+evy*v
c         zbin = tcz + euz*u+evz*v

         xbin = xDetCenter + eux*u+evx*v
         ybin = yDetCenter + euy*u+evy*v
         zbin = zDetCenter + euz*u+evz*v


      dataval=sinomat(ip,jp,kp)
      
      xl=x0
      yl=y0
      zl=z0

      xdiff=xbin-xSource
      ydiff=ybin-ySource
      zdiff=zbin-zSource
      xad=abs(xdiff)/dx
      yad=abs(ydiff)/dy
      zad=abs(zdiff)/dz


      xoz=xdiff/zdiff
      yoz=ydiff/zdiff
      travVoxlen=dz*sqrt(1.0+xoz**2+yoz**2)
      irp=0
      raysum=0.
      do iz=0, nz-1
         z=zl+dz*(iz+0.5)
         x=xSource+xoz*(z-zSource)
         y=ySource+yoz*(z-zSource)
         ix=floor((x-x0)/dx - 0.5)
         iy=floor((y-y0)/dy - 0.5)
         if (((ix.ge.0).and.(ix.lt.nx-1)).and.
     +       ((iy.ge.0).and.(iy.lt.ny-1))) then
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

            irp=irp+1
            raypix(irp)=travVoxlen*wx1*wy1
            iraypix(irp)=ix+1
            jraypix(irp)=iy+1
            kraypix(irp)=iz+1

            irp=irp+1
            raypix(irp)=travVoxlen*wx2*wy1
            iraypix(irp)=ix+2
            jraypix(irp)=iy+1
            kraypix(irp)=iz+1

            irp=irp+1
            raypix(irp)=travVoxlen*wx1*wy2
            iraypix(irp)=ix+1
            jraypix(irp)=iy+2
            kraypix(irp)=iz+1

            irp=irp+1
            raypix(irp)=travVoxlen*wx2*wy2
            iraypix(irp)=ix+2
            jraypix(irp)=iy+2
            kraypix(irp)=iz+1

         endif
      enddo


      if (raysum.lt.small) raysum=1.0
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
      end do
      end do

      Do ip=1,nx
         Do jp=1,ny
            Do kp=1,nz
            if (onemat(ip,jp,kp).lt.small) onemat(ip,jp,kp)=1.0
            smat(ip,jp,kp) = smat(ip,jp,kp)*wsmat(ip,jp,kp)/
     +                                      onemat(ip,jp,kp)
            End Do
         End Do
      End Do


      return
      end


