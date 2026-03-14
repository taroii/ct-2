      subroutine sinoscale(sinomat,
     +       frame_vectors,
     +       ns,nu,nv,du,dv,u0,v0)

      Implicit Real*4(a-h,o-z)


      Dimension frame_vectors(ns,18)
      Dimension sinomat(ns,nu,nv)


Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in,out) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat

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

      Do jp=1,nu
      Do kp=1,nv

         u = u0+(jp-0.5)*du 
         v = v0+(kp-0.5)*dv 


         xbin = xDetCenter + eux*u+evx*v
         ybin = yDetCenter + euy*u+evy*v
         zbin = zDetCenter + euz*u+evz*v

         sb_dist=Sqrt((xSource-xbin)**2+
     +                (ySource-ybin)**2+
     +                (zSource-zbin)**2)
         sinomat(ip,jp,kp) = sinomat(ip,jp,kp)/sb_dist
      End do
      end do
      end do

      return
      end


      subroutine sinoave(sinomat,sinoout,
     +       ns,nu,nv,nss,nus,nvs,nds,ndu,ndv,nws,nwu,nwv)

      Implicit Real*4(a-h,o-z)


      Dimension sinomat(ns,nu,nv),sinoout(nss,nus,nvs)


Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
Cf2py intent(in) nss
Cf2py intent(in) nus
Cf2py intent(in) nvs
Cf2py intent(in) nds
Cf2py intent(in) ndu
Cf2py intent(in) ndv
Cf2py intent(in) nws
Cf2py intent(in) nwu
Cf2py intent(in) nwv
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat
Cf2py intent(out) sinoout
cf2py depend(nss) sinoout
cf2py depend(nus) sinoout
cf2py depend(nvs) sinoout


      do i = 1,nss
         do j = 1,nus
            do k = 1, nvs

               do ii = 1, nws
                  do jj = 1, nwu
                     do kk = 1, nwv
                        i0 = (i-1)*nds + ii
                        j0 = (j-1)*ndu + jj
                        k0 = (k-1)*ndv + kk
                        sinoout(i,j,k)=sinoout(i,j,k)+sinomat(i0,j0,k0)
                     end do
                  end do
               end do
               sinoout(i,j,k) = sinoout(i,j,k)/(nws*nwu*nwv)

            end do
         end do
      end do

      return
      end

      subroutine sinodist(sinomat,sinoin,
     +       ns,nu,nv,nss,nus,nvs,nds,ndu,ndv,nws,nwu,nwv)

      Implicit Real*4(a-h,o-z)


      Dimension sinomat(ns,nu,nv),sinoin(nss,nus,nvs)


Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nv
Cf2py intent(in) nss
Cf2py intent(in) nus
Cf2py intent(in) nvs
Cf2py intent(in) nds
Cf2py intent(in) ndu
Cf2py intent(in) ndv
Cf2py intent(in) nws
Cf2py intent(in) nwu
Cf2py intent(in) nwv
Cf2py intent(out) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
cf2py depend(nv) sinomat
Cf2py intent(in) sinoin
cf2py depend(nss) sinoin
cf2py depend(nus) sinoin
cf2py depend(nvs) sinoin


      do i = 1,nss
         do j = 1,nus
            do k = 1, nvs

               do ii = 1, nws
                  do jj = 1, nwu
                     do kk = 1, nwv
                        i0 = (i-1)*nds + ii
                        j0 = (j-1)*ndu + jj
                        k0 = (k-1)*ndv + kk
                        sinomat(i0,j0,k0) = sinoin(i,j,k)
                     end do
                  end do
               end do

            end do
         end do
      end do

      return
      end


      function ellipse_voxval(x,y,z,
     +   att,alpha,beta,x0,y0,z0,ax,ay,az) 

      Implicit Real*4(a-h,o-z)


      ca = cos(alpha)
      sa = sin(alpha)
      cb = cos(beta)
      sb = sin(beta)

      rel_x = x - x0
      rel_y = y - y0
      rel_z = z - z0

      rot1x= cb*rel_x + sb*rel_y
      rot1y=-sb*rel_x + cb*rel_y
      rot1z= rel_z

      rot2x= ca*rot1x + sa*rot1z
      rot2y= rot1y
      rot2z=-sa*rot1x + ca*rot1z


      ellipsoid_lhs=(rot2x/ax)**2+(rot2y/ay)**2+(rot2z/az)**2

      ellipse_voxval=0.
      if (ellipsoid_lhs.lt.1.) ellipse_voxval=att

      return
      end


      subroutine ellipsoid_embed(smat,
     +        nx,ny,nz,dx,dy,dz,x0i,y0i,z0i,
     +        att,
     +        alpha,beta,x0,y0,z0,
     +        ax,ay,az)


      Implicit Real*4(a-h,o-z)

      Dimension smat(nx,ny,nz)


cf2py intent(in) nx,ny,nz,dx,dy,dz,x0i,y0i,z0i
cf2py intent(out) smat
cf2py depend(nx,ny,nz) smat


      do ip = 1, nx
      do jp = 1, ny
      do kp = 1, nz

      x=x0i+(ip-0.5)*dx
      y=y0i+(jp-0.5)*dy
      z=z0i+(kp-0.5)*dz

      ca = cos(alpha)
      sa = sin(alpha)
      cb = cos(beta)
      sb = sin(beta)

      rel_x = x - x0
      rel_y = y - y0
      rel_z = z - z0

      rot1x= cb*rel_x + sb*rel_y
      rot1y=-sb*rel_x + cb*rel_y
      rot1z= rel_z

      rot2x= ca*rot1x + sa*rot1z
      rot2y= rot1y
      rot2z=-sa*rot1x + ca*rot1z


      ellipsoid_lhs=(rot2x/ax)**2+(rot2y/ay)**2+(rot2z/az)**2

      ellipse_voxval=0.
      if (ellipsoid_lhs.lt.1.) smat(ip,jp,kp)=att
      end do
      end do
      end do

      return
      end





      function step(x)

      Real*4 x, step

      step=0.0
      if (x.ge.0) step=1.0

      return
      end


      subroutine framer(xppp,yppp,zppp,x,y,z,
     +      alpha,beta,gam,x0,y0,z0)

      Implicit Real*4(a-h,o-z)



      xp=(x-x0)*Cos(alpha) - (y0 - y)*Sin(alpha)
      yp=(y-y0)*Cos(alpha) + (x0 - x)*Sin(alpha)
      zp=z-z0

c      xpp=xp
c      ypp=yp*Cos(beta) + zp*Sin(beta)
c      zpp=-yp*Sin(beta) + zp*Cos(beta)

      xpp=xp*Cos(beta) - zp*Sin(beta)
      ypp=yp
      zpp=xp*Sin(beta) + zp*Cos(beta)


      xppp=xpp*Cos(gam) + ypp*Sin(gam)
      yppp=-xpp*Sin(gam) + ypp*Cos(gam)
      zppp= zpp

      return
      end




      function gen_obj_voxval(x,y,z,
     +        att,
     +        nsurf,
     +        alpha,beta,gam,x0,y0,z0,
     +        nplane,
     +        nellipsoid,
     +        ax_ell,ay_ell,az_ell,
     +        ncyl,
     +        ax_cyl,ay_cyl,
     +        ncone,
     +        ax_cone,ay_cone)

c nsurf = nplane + nellipsoid + ncyl + ncone

      Implicit Real*4(a-h,o-z)


      Dimension alpha(nsurf),beta(nsurf),gam(nsurf)
      Dimension x0(nsurf),y0(nsurf),z0(nsurf)
      Dimension ax_ell(nellipsoid),ay_ell(nellipsoid),az_ell(nellipsoid)
      Dimension ax_cyl(ncyl),ay_cyl(ncyl)
      Dimension ax_cone(ncone),ay_cone(ncone)

cf2py intent(in) nsurf
cf2py intent(in) alpha, beta, gam, x0, y0, z0
cf2py depend(nsurf) alpha, beta, gam, x0, y0, z0

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



c loop through all surfaces
      isurf=0

      surf_interior = 1.0

c loop through planes

      Do iplane=2,nplane
         isurf=isurf+1
         Call framer(xs,ys,zs,
     +      x,y,z,
     +      alpha(isurf),beta(isurf),gam(isurf),
     +      x0(isurf),y0(isurf),z0(isurf))

         surf_interior=surf_interior*step(zs)
      End Do

c loop through ellipsoids

      Do iell=2,nellipsoid
         isurf=isurf+1
         Call framer(xs,ys,zs,
     +      x,y,z,
     +      alpha(isurf),beta(isurf),gam(isurf),
     +      x0(isurf),y0(isurf),z0(isurf))

         eq_ell= (xs/ax_ell(iell))**2+(ys/ay_ell(iell))**2+
     +          (zs/az_ell(iell))**2
         surf_interior=surf_interior*step(1.0-eq_ell)
      End Do

c loop through cylinders

      Do icyl=2,ncyl
         isurf=isurf+1
c put source and bin points in current plane frame
         Call framer(xs,ys,zs,
     +      x,y,z,
     +      alpha(isurf),beta(isurf),gam(isurf),
     +      x0(isurf),y0(isurf),z0(isurf))

         eq_cyl= (xs/ax_cyl(icyl))**2+(ys/ay_cyl(icyl))**2
         surf_interior=surf_interior*step(1.0-eq_cyl)
      End Do

c loop through cones

      Do icone=2,ncone
         isurf=isurf+1
c put source and bin points in current plane frame
         Call framer(xs,ys,zs,
     +      x,y,z,
     +      alpha(isurf),beta(isurf),gam(isurf),
     +      x0(isurf),y0(isurf),z0(isurf))

         eq_cone= (xs/ax_cone(icone))**2+(ys/ay_cone(icone))**2
         surf_interior=surf_interior*step(zs**2-eq_cone)*step(zs)
      End Do

      gen_obj_voxval=att*surf_interior

      return
      end



      subroutine gen_obj_embed(smat,
     +        nx,ny,nz,dx,dy,dz,x0i,y0i,z0i,
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

      Dimension smat(nx,ny,nz)

      Dimension  xu1(nsurf),xu2(nsurf),xu3(nsurf)
      Dimension  yu1(nsurf),yu2(nsurf),yu3(nsurf)
      Dimension  zu1(nsurf),zu2(nsurf),zu3(nsurf)
      Dimension  x0(nsurf),y0(nsurf),z0(nsurf)
      Dimension  ax_ell(nellipsoid),ay_ell(nellipsoid)
      Dimension az_ell(nellipsoid)
      Dimension  ax_cyl(ncyl),ay_cyl(ncyl)
      Dimension  ax_cone(ncone),ay_cone(ncone)

cf2py intent(in) nx,ny,nz,dx,dy,dz,x0i,y0i,z0i
cf2py intent(out) smat
cf2py depend(nx,ny,nz) smat

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

      do ip = 1, nx
      do jp = 1, ny
      do kp = 1, nz

      x=x0i+(ip-0.5)*dx
      y=y0i+(jp-0.5)*dy
      z=z0i+(kp-0.5)*dz


c loop through all surfaces
      isurf=0

      surf_interior = 1.0

c loop through planes

      Do iplane=2,nplane
         isurf=isurf+1

         xt=x-x0(isurf)
         yt=y-y0(isurf)
         zt=z-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt


c         Call framer(xs,ys,zs,
c     +      x,y,z,
c     +      alpha(isurf),beta(isurf),gam(isurf),
c     +      x0(isurf),y0(isurf),z0(isurf))

         surf_interior=surf_interior*step(zs)
      End Do

c loop through ellipsoids

      Do iell=2,nellipsoid
         isurf=isurf+1

         xt=x-x0(isurf)
         yt=y-y0(isurf)
         zt=z-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt


c         Call framer(xs,ys,zs,
c     +      x,y,z,
c     +      alpha(isurf),beta(isurf),gam(isurf),
c     +      x0(isurf),y0(isurf),z0(isurf))

         eq_ell= (xs/ax_ell(iell))**2+(ys/ay_ell(iell))**2+
     +          (zs/az_ell(iell))**2
         surf_interior=surf_interior*step(1.0-eq_ell)
      End Do

c loop through cylinders

      Do icyl=2,ncyl
         isurf=isurf+1

         xt=x-x0(isurf)
         yt=y-y0(isurf)
         zt=z-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt


c         Call framer(xs,ys,zs,
c     +      x,y,z,
c     +      alpha(isurf),beta(isurf),gam(isurf),
c     +      x0(isurf),y0(isurf),z0(isurf))

         eq_cyl= (xs/ax_cyl(icyl))**2+(ys/ay_cyl(icyl))**2
         surf_interior=surf_interior*step(1.0-eq_cyl)
      End Do

c loop through cones

      Do icone=2,ncone
         isurf=isurf+1

         xt=x-x0(isurf)
         yt=y-y0(isurf)
         zt=z-z0(isurf)

         xs=xu1(isurf)*xt+xu2(isurf)*yt+xu3(isurf)*zt
         ys=yu1(isurf)*xt+yu2(isurf)*yt+yu3(isurf)*zt
         zs=zu1(isurf)*xt+zu2(isurf)*yt+zu3(isurf)*zt

c         Call framer(xs,ys,zs,
c     +      x,y,z,
c     +      alpha(isurf),beta(isurf),gam(isurf),
c     +      x0(isurf),y0(isurf),z0(isurf))

         eq_cone= (xs/ax_cone(icone))**2+(ys/ay_cone(icone))**2
         surf_interior=surf_interior*step(zs**2-eq_cone)*step(zs)
      End Do

      smat(ip,jp,kp)=att*surf_interior
      end do
      end do
      end do

      return
      end
