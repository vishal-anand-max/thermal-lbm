! Domain size is lx \times ly \times \lz
        integer lx,ly,lz
        PARAMETER(lx=51,ly=51,lz=41)
        common/AA/ G,tau,Gads(2),ex(0:18),ey(0:18),ez(0:18),opp(18)
        real*8 G,tau,Gads
        integer ex,ey,ez, opp

        common/b/ error,vel,xc(0:18),yc(0:18),zc(0:18),t_k(0:18)
        real*8 error,vel,xc,yc,zc,t_k

        common/vel/ c_squ,cc,TT0,rho_w,Ca,RR, Nwri
        real*8 c_squ,cc,TT0, rho_w, Ca,RR
        integer Nwri

        common/app/ t_0,t_1,t_2, rho_h, rho_l
        real*8 t_0,t_1,t_2, rho_h, rho_l
