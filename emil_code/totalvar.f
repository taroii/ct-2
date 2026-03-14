


      subroutine gmi(smatout,smat,nx,ny,nz)

      Implicit Real*4(a-h,o-z)
      Parameter(eps=0.000000001)

      Dimension smat(nx,ny,nz),smatout(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat
Cf2py intent(out) smatout
cf2py depend(nx) smatout
cf2py depend(ny) smatout
cf2py depend(nz) smatout


      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               smatout(ix,iy,iz)=
     +              sqrt(eps+ (smat(ix,iy,iz)-smat(ix-1,iy,iz))**2+
     +                      (smat(ix,iy,iz)-smat(ix,iy-1,iz))**2+
     +                      (smat(ix,iy,iz)-smat(ix,iy,iz-1))**2)
            End Do
         End Do
      End Do
      
      return
      end


      function total_var3d(ambdala,smat,nx,ny,nz)

      Implicit Real*4(a-h,o-z)
      Parameter(eps=0.000000001)

      Real*8 tot

      Dimension smat(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      tot=0.0
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               tot=tot+sqrt(eps+ (smat(ix,iy,iz)-smat(ix-1,iy,iz))**2+
     +                      (smat(ix,iy,iz)-smat(ix,iy-1,iz))**2+
     +                      (smat(ix,iy,iz)-smat(ix,iy,iz-1))**2)+
     +                  ambdala*abs(smat(ix,iy,iz))
            End Do
         End Do
      End Do
      shift = (nx-1.)*(ny-1.)*(nz-1.)*sqrt(eps)
      total_var3d=tot-shift
      return
      end


      function calc_lp(p,smat,nx,ny,nz)

      Implicit Real*4(a-h,o-z)
      Parameter(eps=0.000000001)

      Real*8 tot
      Dimension smat(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      tot=0.0
      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1,nz
               tot=tot+sqrt(eps + smat(ix,iy,iz)**2)**p
            End Do
         End Do
      End Do
      calc_lp=tot
      return
      end



      function total_var3d_p(p,smat,nx,ny,nz)

      Implicit Real*4(a-h,o-z)
      Parameter(eps=0.000000001)

      Real*8 tot
      Dimension smat(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      tot=0.0
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               tot=tot+sqrt(eps +
     +                      (smat(ix,iy,iz)-smat(ix-1,iy,iz))**2+
     +                      (smat(ix,iy,iz)-smat(ix,iy-1,iz))**2+
     +                      (smat(ix,iy,iz)-smat(ix,iy,iz-1))**2 )**p
            End Do
         End Do
      End Do
      shift = (nx-1.)*(ny-1.)*(nz-1.)*sqrt(eps)**p
      total_var3d_p=tot-shift
      return
      end



      function total_dvar3d_p(p,smat,nx,ny,nz,dx,dy,dz)

      Implicit Real*4(a-h,o-z)
      Parameter(eps=0.000000001)

      Real*8 tot
      Dimension smat(nx,ny,nz)
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      tot=0.0
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               tot=tot+sqrt(eps +
     +                   ((smat(ix,iy,iz)-smat(ix-1,iy,iz))/dx)**2+
     +                   ((smat(ix,iy,iz)-smat(ix,iy-1,iz))/dy)**2+
     +                   ((smat(ix,iy,iz)-smat(ix,iy,iz-1))/dz)**2 )**p
            End Do
         End Do
      End Do
      total_dvar3d_p=tot
      return
      end



      subroutine lp_grad3d(ppp,tvgmat,smat,nx,ny,nz,
     +                            normalizeq,bobo)

      Implicit Real*4(a-h,o-z)

      Parameter(eps=0.000000001)

      Real*8 beebee

      Dimension smat(nx,ny,nz),tvgmat(nx,ny,nz)
Cf2py intent(out) tvgmat
Cf2py intent(out) bobo
cf2py depend(nx) tvgmat
cf2py depend(ny) tvgmat
cf2py depend(nz) tvgmat
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               tvgmat(ix,iy,iz)=0.
            End Do
         End Do
      End Do

c compute main grad term
      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1,nz
               w14=smat(ix  ,iy  ,iz  )

               tvgmat(ix,iy,iz)= ppp*w14 *sqrt(w14*w14+eps)**(ppp-2)

            End Do
         End Do
      End Do

      beebee = 0.
      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               beebee = beebee + tvgmat(ix,iy,iz)**2
            End Do
         End Do
      End Do

      bobo=sqrt(beebee)
      if (abs(bobo).lt.0.00001) bobo=1.0

      if (normalizeq.eq.1) then
         Do ix=1,nx
            Do iy = 1,ny
               Do iz =1, nz
                  tvgmat(ix,iy,iz)=tvgmat(ix,iy,iz)/bobo
               End Do
            End Do
         End Do
      end if

      return
      end





      subroutine total_p_var_grad3d(ppp,tvgmat,smat,nx,ny,nz,
     +                            normalizeq,bobo)

      Implicit Real*4(a-h,o-z)

      Parameter(eps=0.000000001)

      Real*8 beebee

      Dimension smat(nx,ny,nz),tvgmat(nx,ny,nz)
Cf2py intent(out) tvgmat
Cf2py intent(out) bobo
cf2py depend(nx) tvgmat
cf2py depend(ny) tvgmat
cf2py depend(nz) tvgmat
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               tvgmat(ix,iy,iz)=0.
            End Do
         End Do
      End Do

c compute main grad term
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               w5=smat(ix  ,iy  ,iz-1)
               w11=smat(ix  ,iy-1,iz  )
               w13=smat(ix-1,iy  ,iz  )
               w14=smat(ix  ,iy  ,iz  )

               term1 = (-w11 - w13 + 3*w14 - w5)*
     -  (Sqrt((w11 - w14)**2 + (w13 - w14)**2 + (-w14 + w5)**2 + eps) 
     -   )**(ppp-2)

               tvgmat(ix,iy,iz)= ppp*(term1)

            End Do
         End Do
      End Do

c compute term one-ahead in x
      Do ix=1,nx-1
         Do iy = 2,ny
            Do iz = 2,nz
               w6=smat(ix+1,iy  ,iz-1)
               w12=smat(ix+1,iy-1,iz  )
               w14=smat(ix  ,iy  ,iz  )
               w15=smat(ix+1,iy  ,iz  )
               termx  =  (w14 - w15)*
     -  (Sqrt((w12 - w15)**2 + (w14 - w15)**2 + (-w15 + w6)**2 + eps)
     -   )**(ppp-2)

               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + ppp*termx

            End Do
         End Do
      End Do

c compute term one-ahead in y
      Do ix=2,nx
         Do iy = 1,ny-1
            Do iz = 2,nz
               w8=smat(ix  ,iy+1,iz-1)
               w14=smat(ix  ,iy  ,iz  )
               w16=smat(ix-1,iy+1,iz  )
               w17=smat(ix  ,iy+1,iz  )
               termy = (w14 - w17)*
     -  (Sqrt((w14 - w17)**2 + (w16 - w17)**2 + (-w17 + w8)**2 + eps)
     -   )**(ppp-2)

               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + ppp*termy

            End Do
         End Do
      End Do

c compute term one-ahead in z
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 1,nz-1
               w14=smat(ix  ,iy  ,iz  )
               w20=smat(ix  ,iy-1,iz+1)
               w22=smat(ix-1,iy  ,iz+1)
               w23=smat(ix  ,iy  ,iz+1)
               termz = (w14 - w23)*
     -  (Sqrt((w14 - w23)**2 + (w20 - w23)**2 + (w22 - w23)**2 + eps)
     -   )**(ppp-2)

               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + ppp*termz

            End Do
         End Do
      End Do

      beebee = 0.
      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               beebee = beebee + tvgmat(ix,iy,iz)**2
            End Do
         End Do
      End Do

      bobo=sqrt(beebee)
      if (abs(bobo).lt.0.00001) bobo=1.0

      if (normalizeq.eq.1) then
         Do ix=1,nx
            Do iy = 1,ny
               Do iz =1, nz
                  tvgmat(ix,iy,iz)=tvgmat(ix,iy,iz)/bobo
               End Do
            End Do
         End Do
      end if

      return
      end



      subroutine total_p_dvar_grad3d(ppp,tvgmat,smat,nx,ny,nz,
     +                            dx,dy,dz,normalizeq,bobo)

      Implicit Real*4(a-h,o-z)

      Parameter(eps=0.000000001)

      Real*8 beebee

      Dimension smat(nx,ny,nz),tvgmat(nx,ny,nz)
Cf2py intent(out) tvgmat
Cf2py intent(out) bobo
cf2py depend(nx) tvgmat
cf2py depend(ny) tvgmat
cf2py depend(nz) tvgmat
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               tvgmat(ix,iy,iz)=0.
            End Do
         End Do
      End Do


c main term
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               w5=smat(ix  ,iy  ,iz-1)
               w11=smat(ix  ,iy-1,iz  )
               w13=smat(ix-1,iy  ,iz  )
               w14=smat(ix  ,iy  ,iz  )
        term1 = ((w14-w11)/(dy*dy)+(w14-w13)/(dx*dx)+(w14-w5)/(dz*dz))*
     - (Sqrt(((w11-w14)/dy)**2+((w13-w14)/dx)**2+((-w14+w5)/dz)**2+eps) 
     -  )**(ppp-2)

               tvgmat(ix,iy,iz)= ppp*term1

            End Do
         End Do
      End Do

c one-ahead in x
      Do ix=1,nx-1
         Do iy = 2,ny
            Do iz = 2,nz
               w6=smat(ix+1,iy  ,iz-1)
               w12=smat(ix+1,iy-1,iz  )
               w14=smat(ix  ,iy  ,iz  )
               w15=smat(ix+1,iy  ,iz  )
               termx =  ((w14 - w15)/(dx*dx))*
     - (Sqrt(((w12-w15)/dy)**2+((w14-w15)/dx)**2+((-w15+w6)/dz)**2+eps)
     -  )**(ppp-2)

               tvgmat(ix,iy,iz)=tvgmat(ix,iy,iz) + ppp*termx

            End Do
         End Do
      End Do

c one-ahead in y
      Do ix=2,nx
         Do iy = 1,ny-1
            Do iz = 2,nz
               w8=smat(ix  ,iy+1,iz-1)
               w14=smat(ix  ,iy  ,iz  )
               w16=smat(ix-1,iy+1,iz  )
               w17=smat(ix  ,iy+1,iz  )
               termy = ((w14 - w17)/(dy*dy))*
     -  (Sqrt(((w14-w17)/dy)**2+((w16-w17)/dx)**2+((-w17+w8)/dz)**2+eps)
     -   )**(ppp-2)

               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + ppp*termy

            End Do
         End Do
      End Do

c one-ahead in z
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 1,nz-1
               w14=smat(ix  ,iy  ,iz  )
               w20=smat(ix  ,iy-1,iz+1)
               w22=smat(ix-1,iy  ,iz+1)
               w23=smat(ix  ,iy  ,iz+1)
               termz = ((w14 - w23)/(dz*dz))*
     - (Sqrt(((w14-w23)/dz)**2+((w20-w23)/dy)**2+((w22-w23)/dx)**2+eps)
     -  )**(ppp-2)

               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + ppp*termz

            End Do
         End Do
      End Do

      beebee = 0.
      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               beebee = beebee + tvgmat(ix,iy,iz)**2
            End Do
         End Do
      End Do

      bobo=sqrt(beebee)
      if (abs(bobo).lt.0.00001) bobo=1.0

      if (normalizeq.eq.1) then
         Do ix=1,nx
            Do iy = 1,ny
               Do iz =1, nz
                  tvgmat(ix,iy,iz)=tvgmat(ix,iy,iz)/bobo
               End Do
            End Do
         End Do
      end if

      return
      end




      subroutine total_var_grad3d(ambdala,tvgmat,smat,nx,ny,nz,
     +                            normalizeq,bobo)

      Implicit Real*4(a-h,o-z)

      Parameter(eps=0.000000001)

      Real*8 beebee

      Dimension smat(nx,ny,nz),tvgmat(nx,ny,nz)
Cf2py intent(out) tvgmat
Cf2py intent(out) bobo
cf2py depend(nx) tvgmat
cf2py depend(ny) tvgmat
cf2py depend(nz) tvgmat
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
cf2py depend(nz) smat

      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               tvgmat(ix,iy,iz)=0.
            End Do
         End Do
      End Do

c main TV term and l1 term
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 2,nz
               w5=smat(ix  ,iy  ,iz-1)
               w11=smat(ix  ,iy-1,iz  )
               w13=smat(ix-1,iy  ,iz  )
               w14=smat(ix  ,iy  ,iz  )
               term1 = (-w11 - w13 + 3*w14 - w5)/
     -  Sqrt((w11 - w14)**2 + (w13 - w14)**2 + (-w14 + w5)**2 + eps) 
               term2 = ambdala*w14/(eps+abs(w14))


               tvgmat(ix,iy,iz)= term1 + term2

            End Do
         End Do
      End Do

c one-ahead in x
      Do ix=1,nx-1
         Do iy = 2,ny
            Do iz = 2,nz
               w6=smat(ix+1,iy  ,iz-1)
               w12=smat(ix+1,iy-1,iz  )
               w14=smat(ix  ,iy  ,iz  )
               w15=smat(ix+1,iy  ,iz  )
               termx =  (w14 - w15)/
     -  Sqrt((w12 - w15)**2 + (w14 - w15)**2 + (-w15 + w6)**2 + eps)


               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + termx

            End Do
         End Do
      End Do

c one-ahead in y
      Do ix=2,nx
         Do iy = 1,ny-1
            Do iz = 2,nz
               w8=smat(ix  ,iy+1,iz-1)
               w14=smat(ix  ,iy  ,iz  )
               w16=smat(ix-1,iy+1,iz  )
               w17=smat(ix  ,iy+1,iz  )
               termy = (w14 - w17)/
     -  Sqrt((w14 - w17)**2 + (w16 - w17)**2 + (-w17 + w8)**2 + eps)

               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + termy

            End Do
         End Do
      End Do

c one-ahead in z
      Do ix=2,nx
         Do iy = 2,ny
            Do iz = 1,nz-1
               w14=smat(ix  ,iy  ,iz  )
               w20=smat(ix  ,iy-1,iz+1)
               w22=smat(ix-1,iy  ,iz+1)
               w23=smat(ix  ,iy  ,iz+1)
               termz = (w14 - w23)/
     -  Sqrt((w14 - w23)**2 + (w20 - w23)**2 + (w22 - w23)**2 + eps)


               tvgmat(ix,iy,iz)= tvgmat(ix,iy,iz) + termz

            End Do
         End Do
      End Do

      beebee = 0.
      Do ix=1,nx
         Do iy = 1,ny
            Do iz = 1, nz
               beebee = beebee + tvgmat(ix,iy,iz)**2
            End Do
         End Do
      End Do

      bobo=sqrt(beebee)
      if (abs(bobo).lt.0.00001) bobo=1.0


      if (normalizeq.eq.1) then
         Do ix=1,nx
            Do iy = 1,ny
               Do iz =1, nz
                  tvgmat(ix,iy,iz)=tvgmat(ix,iy,iz)/bobo
               End Do
            End Do
         End Do
      end if



      return
      end

