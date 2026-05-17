        real*8 h_ref
        integer lx,ly,lz
        PARAMETER(lx=70,ly=70,lz=70)

        common/AA/  BUOYANCY,BD
        integer  BUOYANCY, BD

        common/AA/ ex(0:18),ey(0:18),ez(0:18)
        real*8 ex,ey,ez

        common/b/ error,vel,xc(0:18),yc(0:18),zc(0:18),t_k(0:18)
        real*8 error,vel,xc,yc,zc,t_k

        common/vel/ c_squ,cc,Nwri, opp(18)
        real*8 c_squ,cc
        integer Nwri,opp

        common/app/ psi_max,psi_min,rho_h, rho_l, t_0,t_1,t_2,xsten
        real*8 psi_max,psi_min,rho_h, rho_l, t_0,t_1,t_2,xsten

        common /prho/ dprho_x(lx,ly,lz), dprho_y(lx,ly,lz)
        real*8 dprho_x, dprho_y
        common /fai/ dfai_x(lx,ly,lz),dfai_y(lx,ly,lz),dfai_z(lx,ly,lz)
        real*8 dfai_x, dfai_y,dfai_z
        common /fai2/ dprho_z(lx,ly,lz)
        real*8 dprho_z
        common /velo/ Ra,velS, Kappa, gforce, tau_f, tau_g, UU,DD
        real*8  Ra,velS, Kappa, gforce, tau_f, tau_g, UU,DD
