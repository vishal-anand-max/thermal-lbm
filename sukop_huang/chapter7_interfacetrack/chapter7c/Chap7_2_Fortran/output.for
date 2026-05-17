      subroutine  write_results2(obst,rho,upx,upy,upz,p, n)
      implicit none
      include "head.inc"
      integer  x,y,z,i,n,obsval,k1,k,obst(lx,ly,lz)
      real*8  rho1,h(lx),rho(lx,ly,lz),upx(lx,ly,lz),upy(lx,ly,lz)
      real*8  upz(lx,ly,lz),p(lx,ly,lz)
      character filename*16, B*7, C*1

      write(B,'(i7.7)') n
      filename='3D'//B//'.plt'

      open(41,file=filename)

      write(41,*) 'variables = x, y, z, rho, p, upx, upy, upz,  obst'
      write(41,*) 'zone i=', lx, ', j=', ly, ', k=', lz, ', f=point'
      do 10 z = 1, lz
         do 10 y = 1, ly
          do 10 x = 1, lx
c.........write results to file
          write(41,9) x, y, z,
     &     rho(x,y,z), p(x,y,z),
     &     upx(x,y,z), upy(x,y,z), upz(x,y,z), obst(x,y,z)
   10 continue

   9  format(3i4, 5f15.8, i4)


      close(41)
c----------------------------------
      end
