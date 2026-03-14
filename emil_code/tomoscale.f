

      subroutine tomovolmag(smatout,smat,nx,ny,nz,
     +                      x0,y0,z0,xlen,ylen,zlen,xs,ys,zs)

      Implicit Real*4(a-h,o-z)

      Dimension smat(nx,ny,nz),smatout(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat
Cf2py intent(out) smatout
cf2py depend(nx) smatout
cf2py depend(ny) smatout
cf2py depend(nz) smatout

      dx = xlen/nx
      dy = ylen/ny
      dz = zlen/nz
      Do ix = 1,nx
         Do iy = 1,ny
            Do iz = 1,nz
               smatout(ix,iy,iz) = 0.0
            End Do
         End Do
      End Do

      zdiff = zs-z0
      Do iz = 1,nz
         zh = (iz-0.5)*dz 
         zfrac = zh/zdiff
         xp0 = x0 + zfrac*(xs-x0)
         yp0 = y0 + zfrac*(ys-y0)
         xlenp = xlen*(1.-zfrac)
         ylenp = ylen*(1.-zfrac)
         dxp = xlenp/nx
         dyp = ylenp/ny
         Do ix=1,nx
            xp = xp0 + (ix-0.5)*dxp
            if ((xp.gt.(x0 +0.5*dx)).and.
     +         (xp.lt.(x0+xlen-0.5*dx))) then
               ixr = int((xp - (x0+0.5*dx))/dx)+1
               fxr = (xp - (x0 + (ixr-0.5)*dx))/dx
               Do iy = 1,ny
                  yp = yp0 + (iy-0.5)*dyp
                  if ((yp.gt.(y0 +0.5*dy)).and.
     +               (yp.lt.(y0+ylen-0.5*dy))) then
                     iyr = int((yp - (y0+0.5*dy))/dy)+1
                     fyr = (yp - (y0 + (iyr-0.5)*dy))/dy
                     w11 = (1. - fxr)*(1.-fyr)
                     w21 = (fxr)*(1.-fyr)
                     w12 = (1. - fxr)*(fyr)
                     w22 = fxr*fyr
                     smatout(ix,iy,iz)=
     +                   smat(ixr,iyr,iz)*w11 +
     +                   smat(ixr+1,iyr,iz)*w21 +
     +                   smat(ixr,iyr+1,iz)*w12 +
     +                   smat(ixr+1,iyr+1,iz)*w22 
                  Endif
               End Do
            Endif
         End Do
      End Do
      
      return
      end


      subroutine tomovolmagt(smatout,smat,nx,ny,nz,
     +                      x0,y0,z0,xlen,ylen,zlen,xs,ys,zs)

      Implicit Real*4(a-h,o-z)

      Dimension smat(nx,ny,nz),smatout(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat
Cf2py intent(out) smatout
cf2py depend(nx) smatout
cf2py depend(ny) smatout
cf2py depend(nz) smatout

      dx = xlen/nx
      dy = ylen/ny
      dz = zlen/nz
      Do ix = 1,nx
         Do iy = 1,ny
            Do iz = 1,nz
               smatout(ix,iy,iz) = 0.0
            End Do
         End Do
      End Do

      zdiff = zs-z0
      Do iz = 1,nz
         zh = (iz-0.5)*dz 
         zfrac = zh/zdiff
         xp0 = x0 + zfrac*(xs-x0)
         yp0 = y0 + zfrac*(ys-y0)
         xlenp = xlen*(1.-zfrac)
         ylenp = ylen*(1.-zfrac)
         dxp = xlenp/nx
         dyp = ylenp/ny
         Do ix=1,nx
            xp = xp0 + (ix-0.5)*dxp
            if ((xp.gt.(x0 +0.5*dx)).and.
     +         (xp.lt.(x0+xlen-0.5*dx))) then
               ixr = int((xp - (x0+0.5*dx))/dx)+1
               fxr = (xp - (x0 + (ixr-0.5)*dx))/dx
               Do iy = 1,ny
                  yp = yp0 + (iy-0.5)*dyp
                  if ((yp.gt.(y0 +0.5*dy)).and.
     +               (yp.lt.(y0+ylen-0.5*dy))) then
                     iyr = int((yp - (y0+0.5*dy))/dy)+1
                     fyr = (yp - (y0 + (iyr-0.5)*dy))/dy
                     w11 = (1. - fxr)*(1.-fyr)
                     w21 = (fxr)*(1.-fyr)
                     w12 = (1. - fxr)*(fyr)
                     w22 = fxr*fyr
                     smatout(ixr,iyr,iz)=
     +                  smatout(ixr,iyr,iz)+smat(ix,iy,iz)*w11
                     smatout(ixr+1,iyr,iz)=
     +                  smatout(ixr+1,iyr,iz)+smat(ix,iy,iz)*w21
                     smatout(ixr,iyr+1,iz)=
     +                  smatout(ixr,iyr+1,iz)+smat(ix,iy,iz)*w12
                     smatout(ixr+1,iyr+1,iz)=
     +                  smatout(ixr+1,iyr+1,iz)+smat(ix,iy,iz)*w22 
                  Endif
               End Do
            Endif
         End Do
      End Do
      
      return
      end
