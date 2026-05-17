
      subroutine  write_results2(obst,rho,p,upx,upy,upz, n)
      implicit none
      include "head.inc"
      integer  x,y,z,i,n
      real*8  rho(lx,ly,lz),upx(lx,ly,lz),upy(lx,ly,lz)
      real*8  upz(lx,ly,lz), p(lx,ly,lz)
      integer  obst(lx,ly,lz)

      character filename*16, B*6

      write(B,'(i6.6)') n
      filename='out/3D'//B//'.plt'

      open(41,file=filename)

      write(41,*) 'variables = x, y, z, rho, upx, upy, upz, p, obst'
      write(41,*) 'zone i=', lx, ', j=', ly, ', k=', lz, ', f=point'
      do 10 z = 1, lz
       do 10 y = 1, ly
        do 10 x = 1, lx
          write(41,9) x, y, z, rho(x,y,z),
     & upx(x,y,z), upy(x,y,z), upz(x,y,z),
     &  p(x,y,z), obst(x,y,z)
   10 continue

   9  format(3i4, 5f15.8, i4)

      close(41)
      end
