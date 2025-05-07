program test_gammaf

      implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      INTEGER           JZ, N, JX, ID, i, NPTS
      real*8            r, p, r1, r2, p1, p2, dr, dp
      real*8            dens1, dens2
      character*32      archivo


      COMMON /energy/   EN,PION
      COMMON /capas/    NLVAL
      COMMON /koga/     AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/      ABCIS,WEIGHT


      DIMENSION ABCIS(48),WEIGHT(48)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION EN(103,-1:20),PION(103,-1:20),NLVAL(103,-1:20)

      call transfer
      call gauss(48, abc, weight)


      NPTS = 100
      dr = 0.1d0
      dp = 0.1d0


      write(*,*) 'Introduce Z (carga nuclear) y N (electrones):'
      read(*,*) JZ, N
      JX = JZ - N


      if (JZ < 1 .or. JZ > 103) stop 'Z fuera de rango'
      if (JX < -1 .or. JX > 20) stop 'JX fuera de rango'


      archivo = 'densidades.dat'
      open(10, file=archivo, status='unknown')
      write(10,'(A)') '# x    rho(x)    rho_2(x,x+dx)    gamma(p)   gamma_2(p,p+dp)'


      do i = 0, NPTS
         r  = i * dr
         r1 = r
         r2 = r 
         p  = i * dp
         p1 = p
         p2 = p 

         dens1 = RO(JZ, JX, r)
         dens2 = RO_DOB(JZ, JX, r1, r2)
         gamma1 = GAMMAF(JZ, JX, p)
         gamma2 = GAMMA_DOB(JZ, JX, p1, p2)

         write(10,'(f10.4,2x,e20.10,2x,e20.10,2x,e20.10,2x,e20.10)') &
              r, dens1, dens2, gamma1, gamma2
      end do

      close(10)
      write(*,*) 'Calculo completado. Resultados en densidades.dat'

      end
