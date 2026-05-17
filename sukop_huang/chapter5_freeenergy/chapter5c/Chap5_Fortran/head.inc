       integer lx,ly
c Define the dimension of the computational domain
        PARAMETER(lx=100,ly=100)

c Array xc, yc are $e_{ix}$ and $e_{iy}$, respectively.
c Array t_k is the weighting coefficient in f_i^{eq}
        common/b/ xc(0:8),yc(0:8),t_k(0:8)
        real*8 xc,yc,t_k

c 'cc' is the lattice velocity. 'RR' is the droplet's initial radius.
c a, b, TT are parameters in van der Waals EOS.
        common/vel/ c_squ,cc, RR,df,  a,b, TT
        real*8 c_squ,cc, RR, df,  a,b,  TT

c \tau is the relaxation parameter and 'nu' is the kinematic viscosity.
c \kappa is a key parameter to adjust the surface tension strength
        common/gg/   kappa, tau, nu
        real*8  kappa, tau, nu

c 'Nwri' is the frequency of output (dump data every 'Nwri' time step).
c t_max is the desired maximum time step.
        common/oo/  Nwri, t_max
        integer Nwri, t_max