      subroutine  write_results(lx,ly,obst,rho,upx,upy,p, n)
      implicit none
      include "head.inc"
      integer slicex2,lx,ly
      integer  x,y,i,j,n,obsval,jcen,k1,k,obst(lx,ly),slicex,slicey
      real*8  rho1,h(lx),rho(lx,ly),upx(lx,ly),upy(lx,ly)
     & , p(lx,ly)

      character*40 V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12
      REAL*4 ZONEMARKER, EOHMARKER
      character filename*16, B8*8
      integer*1 obob
      integer*2 ii, jj
      integer, parameter:: kmax=1
      character*40 Title,var, Zonename1

      write(B8,'(i8.8)') n

      open(41,file='./out/2D'//B8//'.plt',form="BINARY")


c--------------------------------------
      ZONEMARKER=           299.0
      EOHMARKER =           357.0

c       I. The header section.
      write(41) "#!TDV101"

c-----Integer value of 1.------------------
      write(41) 1


      Title="Rayleigh-Taylor"
      call dumpstring(Title)

c-----Number of variables (NumVar) in the datafile.
      write(41) 7
c------Variable names.
      V1='X'
        call dumpstring(V1)
      V2='Y'
        call dumpstring(V2)
      V3='u'
        call dumpstring(V3)
      V4='v'
        call dumpstring(V4)
      V5='rho1'
        call dumpstring(V5)
      V7='p'
        call dumpstring(V7)
      V8='obst'
        call dumpstring(V8)

      write(41) ZONEMARKER
      Zonename1="ZONE 001"
      call dumpstring(Zonename1)
c---------Zone Color
      write(41) -1
c---------ZoneType  !
      write(41) 0
c---------DataPacking 0=Block, 1=Point
      write(41) 1
c-Specify Var Location. 0 = Don't specify, all data is located at the nodes
      write(41) 0
c---------Number of user defined face neighbor connections (value >= 0)
      write(41) 0
 
      write(41) lx
      write(41) ly
      write(41) kmax

      write(41) 0
      write(41) EOHMARKER


c-----2.1 zone ----------------------------------
      write(41) Zonemarker

C--------variable data format, 1=Float, 2=Double,
c---3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
       write(41) 4
       write(41) 4
       write(41) 1
       write(41) 1
       write(41) 1
       write(41) 1
       write(41) 5
       
       write(41) 0
       write(41) -1

      do k=1,kmax
       do j=1,ly
        do i=1,lx
         obob = obst(i,j)
         ii = i
         jj = j
         write(41) ii
         write(41) jj
         write(41) real(upx(i,j))
         write(41) real(upy(i,j))
         write(41) real(rho(i,j))
         write(41) real(p(i,j))
         write(41) obob
         end do
       end do
      end do

       close(41)
       end

c----------------------------------------------------------

      subroutine dumpstring(instring)
       character(40) instring
       integer len

       len=LEN_TRIM(instring)

       do ii=1,len
         I=ICHAR(instring(ii:ii))
         write(41) I
       end do
       write(41) 0
      return
      end
c--------------------------------
      subroutine track_interface(lx,ly,rho,n)
      implicit none
      include "head.inc"
      integer lx,ly, n, x,y
      real*8  rho(lx,ly), rho_d, pos_d,rho_u, pos_u, pos, pos2

      open(45,file='./out/posi.dat', access='append')
      x=lx/2
      do y=1, ly
       if(rho(x,y) .lt. (rho_h+rho_l)/2.d0 ) then
        rho_d=rho(x,y)
        pos_d=y
       endif
      enddo

      do y=ly, 1, -1
       if(rho(x,y) .gt. (rho_h+rho_l)/2.d0 ) then
        rho_u=rho(x,y)
        pos_u=y
       endif
      enddo

      pos=pos_d+ (pos_u-pos_d)/(rho_u-rho_d)*((rho_h+rho_l)/2.-rho_d)

      x=1
      do y=1, ly
       if(rho(x,y) .lt. (rho_h+rho_l)/2.d0 ) then
        rho_d=rho(x,y)
        pos_d=y
       endif
      enddo

      do y=ly, 1, -1
      if(rho(x,y) .gt. (rho_h+rho_l)/2.d0 ) then
       rho_u=rho(x,y)
       pos_u=y
      endif
      enddo
      pos2=pos_d+ (pos_u-pos_d)/(rho_u-rho_d)*((rho_h+rho_l)/2.-rho_d)
      write(45,'(i8, 2f15.4)') n, pos, pos2

      end