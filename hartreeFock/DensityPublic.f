C#############################################################
C
C
C
C                   #######################################
C                   #######################################
C                   ##                                   ##
C                   ##  INSTRUCCIONES SOBRE LA LIBRERIA  ##
C                   ##                                   ##
C                   #######################################
C                   #######################################
C
C
C
C#############################################################
C
C  EN ESTA LIBRERIA SE UTILIZAN COMO DENSIDADES LAS SIGUIENTES:
C
C   0 = rho(r)   (Norma = N)
C   1 = gammaf(p) (Norma = N)
C   
C   rho(r)   = Densidad de carga
C   gammaf(p) = Densidad de momento
C   
C  LOS DOS PRIMEROS PARAMETROS (JZ,JX) CARACTERIZAN 
C  EL SISTEMA ATOMICO (ATOMOS NEUTROS, IONES, SERIES ISOELECTRONICAS 
C  E HIDROGENO EXCITADO; CALCULO H-F O BCF)
C
C  
C
C#########################################################################
C#########################################################################
C
C
C
C                   ###################################
C                   ###################################
C                   ##                               ##
C                   ##  LAS FUNCIONES DENSIDAD TOTAL ##
C                   ##                               ##
C                   ###################################
C                   ###################################
C
C
C
C#########################################################################

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	DEN(JZ,JX,R,ID)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD 'ID'

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)



      
      IF(ID.EQ.0)DEN=RO(JZ,JX,R)
      IF(ID.EQ.1)DEN=GAMMAF(JZ,JX,R)
      
     
      RETURN

      END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	DENDOB(JZ,JX,R,R1,ID)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD 'ID'

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)



      
      IF(ID.EQ.0)DENDOB=RO_DOB(JZ,JX,R,R1)
      IF(ID.EQ.1)DENDOB=GAMMA_DOB(JZ,JX,R,R1)
      IF(ID.EQ.2)DENDOB=RO(JZ,JX,R)*RO(JZ,JX,R1)
      IF(ID.EQ.3)DENDOB=GAMMAF(JZ,JX,R)*GAMMAF(JZ,JX,R1)
      IF(ID.EQ.4)DENDOB=0.5d0*(RO_DOB(JZ,JX,R,R1)+
     &  RO(JZ,JX,R)*RO(JZ,JX,R1))
      IF(ID.EQ.5)DENDOB=0.5d0*(GAMMA_DOB(JZ,JX,R,R1)+
     &  GAMMAF(JZ,JX,R)*GAMMAF(JZ,JX,R1))    
      RETURN

      END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	RO(JZ,JX,R)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD DE CARGA ATOMICA (NORMALIZADA A 1)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)

        pi=4.d0*datan(1.d0)

      z=abs(jz)

      if(jx.gt.1)then
         ne=jx
         q=z-jx
         else
            ne=z-jx
            q=jx
            endif


      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111


        RO=0.0D0
      
C***************SUMATORIA EN LOS ORBITALES **************************
      
      DO 50 l=0,3
         do 51 nn=l+1,l+norb(z,q,l)
         SUMind=0.0
         
C***************SUMATORIA EN LAS FUNCIONES DE SLATER ****************
         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0
            IF(R.EQ.0.0D0)THEN
               IF (N.EQ.1)THEN
                  ERRE=1.0D0
               ELSE
                  ERRE=0.0D0
               ENDIF
            ELSE
               ERRE=pot(R,N-1)
            ENDIF
            
            SUMind=SUMind+C*((2.d0*A)**(N+.5d0))*ERRE*dEXP(-A*R)/
     #           ((FACT(2*N))**(.5d0))
 60      END DO
         
C***********FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ****************
         
         RO=RO+SUMind*SUMind
 51   END DO
 50   END DO
      
C***********FIN DE LA SUMA EN LOS ORBITALES **************************

      RO=RO/(4.d0*pi*ne)    

      
      RETURN


 111  continue


	      bcf=0.d0

	      do 10 l=0,3
                
                 do 11 nn=l+1,l+norb(z,q,l)
                    nenl=mzi(z,q,nn,l)

		    dens=rho_nl(z,nn,l,r)


		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    ro=bcf

	      return


 222          n=(-jx)/10
              l=-jx-10*n
              ro=rho_nl(Z,n,l,r)

              return

              end

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	RO_DOB(JZ,JX,R,R1)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD DE CARGA ATOMICA DE DOS PARTICULAS (NORMALIZADA A 1)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,fil,CONT,I,J,KONT,NEN
      real*8            SUMIND,SUMIND1,MATRIZ(4,100)
      real*8            mult,noci,nocj,li,lj,h
     
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)



        z=abs(jz)

      if(jx.gt.1)then
         ne=jx
         q=z-jx
         else
            ne=z-jx
            q=jx
            endif


      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111


        RO_DOB=0.0D0
	CONT=0

      
C***************CALCULO DE LOS ORBITALES **************************
   
      DO 50 l=0,3
	 
         do 51 nn=l+1,l+norb(z,q,l)
         SUMind=0.0d0
         SUMind1=0.0d0

         CONT=CONT+1
C***************SUMATORIA EN LAS FUNCIONES DE SLATER ****************         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0

c            IF(R.EQ.0.0D0)THEN
c               IF (N.EQ.1)THEN
c                 ERRE=1.0D0	!SUMIND-VALOR DE LA SUMA DE FUNCIONES DE SLATER EN LA PRIMERA VARIABLE
c               ELSE
c                 ERRE=0.0D0
c               ENDIF
c           ELSE
c               ERRE=R**(N-1)
c            ENDIF

            ERRE=POT(R,N-1)
            
            SUMind=SUMind+C*((2.d0*A)**(N+.5d0))*ERRE*DEXP(-A*R)/
     &           ((FACT(2*N))**(.5d0))






c           IF(R1.EQ.0.0D0)THEN
c               IF (N.EQ.1)THEN
c                  ERRE1=1.0D0
c               ELSE                       !SUMIND1-VALOR DE LA SUMA DE FUNCIONES DE SLATER EN LA PRIMERA VARIABLE
c                  ERRE1=0.0D0
c              ENDIF
c            ELSE
c               ERRE1=R1**(N-1)
c            ENDIF
 
            ERRE1=POT(R1,N-1)


            SUMind1=SUMind1+C*((2*A)**(N+.5d0))*ERRE1*DEXP(-A*R1)/
     &           ((FACT(2*N))**(.5d0))




 60      END DO         
C**************FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ****************

	 MATRIZ(3,CONT)=(NN*10.D0+L*1.D0)*1.0D0   !!ATOMOS NEUTROS
	 MATRIZ(4,CONT)=nesub(z,q,int(matriz(3,cont)))*1.d0
	 MATRIZ(1,CONT)=SUMind/matriz(4,cont)**0.5d0
	 MATRIZ(2,CONT)=SUMind1/matriz(4,cont)**0.5d0

 51   END DO
 50   END DO


      DO I=1,cont
	DO J=I,cont

	noci=matriz(4,i)
	nocj=matriz(4,j)
	li=mod(matriz(3,i),10.d0*int(matriz(3,i)/10))
	lj=mod(matriz(3,j),10.d0*int(matriz(3,j)/10))

	mult=noci*nocj
	if(i==j) mult=noci*(noci-1.d0)/2.d0

	ro_dob=ro_dob+mult*((matriz(1,i)*matriz(2,j))**2.d0+
     &                      (matriz(1,j)*matriz(2,i))**2.d0)

	if((li.eq.lj).and.(i.ne.j))then
	   h=noci
	   if(h.gt.nocj) h=nocj
	  
	   ro_dob=ro_dob-2.d0*h*matriz(1,j)*matriz(2,i)
     &                      *matriz(1,i)*matriz(2,j)
	endif
       enddo
	enddo



      RO_DOB=RO_DOB/(12.5663706143591729539**2.d0*ne*(ne-1))      !!! RO_DOB/(N(N-1)(4PI)^2)

      RETURN


 111  continue                      !ESTO ES PARA IONES Y ESTA SIN TOCAR


	      bcf=0.d0

	      do 10 l=0,3
                
                 do 11 nn=l+1,l+norb(z,q,l)
                    nenl=mzi(z,q,nn,l)

		    dens=rho_nl(z,nn,l,r)


		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    ro=bcf

	      return


 222          n=(-jx)/10
              l=-jx-10*n
              ro=rho_nl(Z,n,l,r)

              return

              end
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	GAMMA_DOB(JZ,JX,p,p1)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD DE CARGA ATOMICA DE DOS PARTICULAS (NORMALIZADA A 1)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,fil,CONT,I,J,KONT,NEN,li,lj,s
      real*8            SUMIND,SUMIND1,MATRIZ(4,100),pepe,pepe1,sums
      real*8            mult,noci,nocj,h,ax,ay

     
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3),d(7)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)


        pi=4.d0*datan(1.d0)
        z=abs(jz)

      if(jx.gt.1)then
         ne=jx
         q=z-jx
         else
            ne=z-jx
            q=jx
            endif


      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111


        GAMMA_DOB=0.0D0
	CONT=0

      
C***************CALCULO DE LOS ORBITALES **************************
   
      DO 50 l=0,3
	 
         do 51 nn=l+1,l+norb(z,q,l)
         SUMind=0.0d0
         SUMind1=0.0d0

         CONT=CONT+1
C***************SUMATORIA EN LAS FUNCIONES DE SLATER ****************         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0

            x=p/a
            y=1.d0+x*x

            tk=c*pot(x,l)*pot(4.d0,n)/(FACT(2*N)*a*a*a*0.25d0)**(.5d0)

            sfk=fact(n)/pot(y,n+1)

        do k=0,int((n-l)/2)
        sumind=sumind+tk*sfk
       sfk=-sfk*(n-l-2*k)*(n-l-2*k-1)*y/((n-k)*(k+1.d0)*4.d0)
        enddo

            x=p1/a
            y=1.d0+x*x

            tk=c*pot(x,l)*pot(4.d0,n)/(FACT(2*N)*a*a*a*0.25d0)**(.5d0)

            sfk=fact(n)/pot(y,n+1)

             do k=0,int((n-l)/2)
       sumind1=sumind1+tk*sfk
       sfk=-sfk*(n-l-2*k)*(n-l-2*k-1)*y/((n-k)*(k+1.d0)*4.d0)
        enddo



      !SUMIND-VALOR DE LA SUMA DE FUNCIONES DE SLATER EN LA PRIMERA VARIABLE
c       do k=0,int((n-l)/2)
c       sumind=sumind+c*(2*A)**(N+.5)/FACT(2*N)**(.5d0)            
c     &  *fact(n-l)*(2*a)**n
c     &  *(p/a)**l*fact(n-k)/(fact(k)*fact(n-l-2*k)*(p**2+a**2)**(n+1-k))
c     &  *(-1/(4*a**2))**k*(2/3.1415D0)**0.5  
c     end do


      !SUMIND1-VALOR DE LA SUMA DE FUNCIONES DE SLATER EN LA PRIMERA VARIABLE                
          
c       do k=0,int((n-l)/2)
c      sumind1=sumind1+c*(2*A)**(N+.5)/FACT(2*N)**(.5d0)
c     &  *fact(n-l)*(2*a)**n                              
c     &  *(p1/a)**l*fact(n-k)/(fact(k)*fact(n-l-2*k)*
c     &  (p1**2+a**2)**(n+1-k))*(-1/(4*a**2))**k*(2/3.1415D0)**0.5
c       end do

 






 60      END DO         
C**************FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ****************
	 MATRIZ(3,CONT)=(NN*10.D0+L*1.D0)*1.0D0   !!ATOMOS NEUTROS
	 MATRIZ(4,CONT)=nesub(z,q,int(matriz(3,cont)))*1.d0
	 MATRIZ(1,CONT)=SUMind/matriz(4,cont)**0.5
	 MATRIZ(2,CONT)=SUMind1/matriz(4,cont)**0.5

 51   END DO
 50   END DO


      DO I=1,cont
	DO J=I,cont

	noci=matriz(4,i)
	nocj=matriz(4,j)
	li=int(mod(matriz(3,i),10.d0*int(matriz(3,i)/10)))
	lj=int(mod(matriz(3,j),10.d0*int(matriz(3,j)/10)))

	mult=noci*nocj
	if(i==j) mult=noci*(noci-1)/2

	gamma_dob=gamma_dob+mult*((matriz(1,i)*matriz(2,j))**2+
     &                      (matriz(1,j)*matriz(2,i))**2)


	if((li.eq.lj).and.(i.ne.j))then
	   h=noci
	   if(h.gt.nocj) h=nocj
	   gamma_dob=gamma_dob-2*h*matriz(1,j)*matriz(2,i)
     &  *matriz(1,i)*matriz(2,j)

	   endif
       enddo
	   enddo



      GAMMA_DOB=GAMMA_DOB/(16.d0*pi*pi*pi*pi*ne*(ne-1))      !!! GAMMA_DOB/(N(N-1)(4PI)^2)

      RETURN


 111  continue                      !ESTO ES PARA IONES Y ESTA SIN TOCAR


	      bcf=0.d0

	      do 10 l=0,3
                
                 do 11 nn=l+1,l+norb(z,q,l)
                    nenl=mzi(z,q,nn,l)

		    dens=rho_nl(z,nn,l,r)


		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    ro=bcf

	      return


 222          n=(-jx)/10
              l=-jx-10*n
              ro=rho_nl(Z,n,l,r)

              return

              end
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function RHO_NL(Z,n,l,r)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c densidad rho(r) para hidrogenoide (n,l) carga Z
c

      implicit real*8 (a-h,o-z)

      integer z

      z=abs(z)

      x=2.d0*z*r/n

      pi=4.d0*datan(1.d0)

      rhonl=0.d0

      do 10 i=0,n-l-1

         do 20 m=0,n-l-1

            term=pot(-1.d0,m+i)*pot(x,m+i+2*l)/(fact(m)*fact(i))

            term=term/(fact(n-l-1-m)*fact(2*l+m+1))

            term=term/(fact(n-l-1-i)*fact(2*l+i+1))

            rhonl=rhonl+term

 20         continue

 10         continue

         rho_nl=rhonl*fact(n-l-1)*fact(n+l)*dexp(-x)*z*z*z/(pi*n*n*n*n)

            return

            end



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	GAMMAF(JZ,JX,PE)
      
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD DE CARGA ATOMICA (NORMALIZADA A N)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,s
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3),d(7)


        pi=4.d0*datan(1.d0)

        z=abs(jz)



      if(jx.gt.1)then
         ne=jx
         q=z-jx
         else
            ne=z-jx
            q=jx
            endif


      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111


      GammaF=0.0D0
      
C**********SUMA EN LOS ORBITALES *********************************
      
             DO 50 l=0,3
             do 52 nn=l+1,l+norb(z,q,l) 


         sumind=0.d0

C**********SUMA EN LAS FUNCIONES DE SLATER ***********************
         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0

            x=pe/a
            y=1.d0+x*x



            tk=c*pot(x,l)*pot(4.d0,n)/(FACT(2*N)*a*a*a)**(.5d0)

            sfk=fact(n)/pot(y,n+1)

	do k=0,int((n-l)/2)
	sumind=sumind+tk*sfk
       sfk=-sfk*(n-l-2*k)*(n-l-2*k-1)*y
     &/((n-k)*(k+1.d0)*4.d0)
	enddo

 60      END DO
         
C***********FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ***************
         
c        sumind=sumind/pi
	   GammaF=GammaF+SUMind*SUMind


 52   END DO
 50   END DO
      
C***************FIN DE LA SUMA EN LOS ORBITALES *********************

         GAMMAF=1.d0*GammaF/(pi*pi*ne)

c     continue

      RETURN

 111  continue

	      bcf=0.d0

	      do 10 l=0,3
                 do 11 nn=l+1,l+norb(z,q,l)

		 nenl=mzi(z,q,nn,l)

		    dens=gam_nl(z,nn,l,pe)

		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    gammaf=bcf

	      return

 222          n=(-jx)/10
              l=-jx-10*n
              gammaf=gam_nl(Z,n,l,pe)

              return

              end



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function GAM_NL(jZ,n,l,p)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c densidad gammaf(p) para hidrogenoide (n,l) carga Z
c

      implicit real*8 (a-h,o-z)

      integer z

      z=abs(jz)

      x=(1.d0*z*z-n*n*p*p)/(1.d0*z*z+n*n*p*p)

      pi=4.d0*datan(1.d0)

      gamnl=0.d0

      do 10 i=0,(n-l-1)/2

         do 20 m=0,(n-l-1)/2

         term=pot(-1.d0,m+i)*pot(2.d0*x,2*(n-l-1-m-i))/(fact(m)*fact(i))
            term=term/(fact(n-l-1-2*m)*fact(n-l-1-2*i))
            term=term*(fact(n-m-1)*fact(n-i-1))

            gamnl=gamnl+term

 20         continue

 10         continue

            gamnl=gamnl*pot(2.d0,2*l-1)

            gamnl=gamnl*n*n*n*n

            gamnl=gamnl*fact(n-l-1)

            gamnl=gamnl/fact(n+l)

            gam_nl=gamnl*pot(1.d0-x,l)*pot(1.d0+x,l+4)/(z*z*z*pi*pi)

            return

            end


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function COMPTON(JZ,JX,q)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


c
c  Perfil Compton J(q) de atomo caracterizado por (JZ,JX)
c
C  (NORMA: J(0)=<p**-1>/2)
  
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

	pi=4.d0*datan(1.d0)

	   hnorm=0.d0

	DELTA=25.d0

           do 1 iv=0,31

	VMAX=q+iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=gammaf(JZ,JX,r(l))

		   hnorm=hnorm+R(L)*dens1*W*delta

88	CONTINUE
20	CONTINUE

 1      continue

	compton=hnorm*3.141592654d0

        return

	END











C#####################################################################
C
C
C
C                   ##########################################
C                   ##########################################
C                   ##                                      ##
C                   ##  LAS FUNCIONES DENSIDAD POR SUBCAPAS ##
C                   ##                                      ##
C                   ##########################################
C                   ##########################################
C
C
C
C#############################################################

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	DENSUB(JZ,JX,nlsub,R,ID)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


C     DENSIDAD SUBCAPAS

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)



      IF(ID.EQ.-2)DENSUB=BRsub(JZ,JX,nlsub,R)
      IF(ID.EQ.0)DENSUB=ROsub(JZ,JX,nlsub,R)
      IF(ID.EQ.1)DENSUB=GAMMAsub(JZ,JX,nlsub,R)
      IF(ID.EQ.2)DENSUB=COMPTONsub(JZ,JX,nlsub,R)
 
      RETURN

      END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	ROSUB(JZ,JX,nlsub,R)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


C     DENSIDAD DE CARGA ATOMICA SUBCAPA(S) 'nlsub'
c     NORMA = Nsub

c     nlsub>0 ===> Capas a considerar
c     nlsub<0 ===> Capas que no contribuyen

c     n=8 para todo n, l=4 para todo l

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)

      nsub=nlsub/10
      lsub=nlsub-10*nsub

      if((nsub.gt.0).and.(lsub.lt.4))then
         lmin=lsub
         lmax=lsub
         else
            lmin=0
            lmax=3
            endif

      z=abs(jz)

      if(jx.gt.1)then
         q=z-jx
         ne=jx
         else
            q=jx
            ne=z-jx
            endif

      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111

        RO=0.0D0
      
C***************SUMATORIA EN LOS ORBITALES **************************
      
      DO 50 l=lmin,lmax

         if((nsub.eq.-8).and.(l.eq.-lsub))goto 50


         if(nsub*(8-nsub).gt.0)then
            nmin=nsub
            nmax=nsub
            else
            nmin=l+1
            nmax=l+norb(z,q,l)
            endif
             

         do 51 nn=nmin,nmax
            if((nsub*(8+nsub).lt.0).and.(nn.eq.-nsub))then
               if(lsub.eq.-4)then
                  goto 51
                  else
                     if(l.eq.-lsub)goto 51
                     endif
                     endif
                  
         SUMind=0.0
         
C***************SUMATORIA EN LAS FUNCIONES DE SLATER ****************
         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0
            IF(R.EQ.0.0D0)THEN
               IF (N.EQ.1)THEN
                  ERRE=1.0D0
               ELSE
                  ERRE=0.0D0
               ENDIF
            ELSE
               ERRE=R**(N-1)
            ENDIF
            
            SUMind=SUMind+C*((2*A)**(N+.5))*ERRE*EXP(-A*R)/
     #           ((FACT(2*N))**(.5))
 60      END DO
         
C***********FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ****************
         
         RO=RO+SUMind*SUMind
 51   END DO
 50   END DO
      
C***********FIN DE LA SUMA EN LOS ORBITALES **************************

      ROsub=RO/12.5663706143591729539

      
      RETURN


 111  continue

	      bcf=0.d0

	      do 10 l=lmin,lmax

         if((nsub.eq.-8).and.(l.eq.-lsub))goto 10



         if(nsub*(8-nsub).gt.0)then
            nmin=nsub
            nmax=nsub
            else
            nmin=l+1
            nmax=l+norb(z,q,l)
            endif
             
                 do 11 nn=nmin,nmax
                    nenl=mzi(z,q,nn,l)

		    dens=rho_nl(z,nn,l,r)

		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    rosub=bcf

	      return


 222          n=(-jx)/10
              l=-jx-10*n
              rosub=rho_nl(Z,n,l,r)

              return

              end


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	GAMMASUB(JZ,JX,nlsub,PE)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      
C     DENSIDAD ATOMICA DE MOMENTOS SUBCAPA(S) 'nlsub' (NORMALIZADA A Nsub)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,s

      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3),d(7)
      DIMENSION NLVAL(103,-1:20)
      

      nsub=nlsub/10
      lsub=nlsub-10*nsub

      if((nsub.gt.0).and.(lsub.lt.4))then
         lmin=lsub
         lmax=lsub
         else
            lmin=0
            lmax=3
            endif

      z=abs(jz)

      if(jx.gt.1)then
         q=z-jx
         ne=jx
         else
            q=jx
            ne=z-jx
            endif

      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111

      GAMMAF=0.0D0
      
C**********SUMA EN LOS ORBITALES *********************************
      
      DO 50 l=lmin,lmax


         if((nsub.eq.-8).and.(l.eq.-lsub))goto 50



         if(nsub*(8-nsub).gt.0)then
            nmin=nsub
            nmax=nsub
            else
            nmin=l+1
            nmax=l+norb(z,q,l)
            endif
             

         do 52 nn=nmin,nmax

            if((nsub*(8+nsub).lt.0).and.(nn.eq.-nsub))then
               if(lsub.eq.-4)then
                  goto 52
                  else
                     if(l.eq.-lsub)goto 52
                     endif
                     endif

         SUMind=0.0D0
         
C**********SUMA EN LAS FUNCIONES DE SLATER ***********************
         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0
            D(1)=1.0D0
            D(3)=-(N-L)*(N-L-1)/(4*L+6.0D0)
            D(5)=-(N-L-2)*(N-L-3)*D(3)/(8*L+20.0D0)
            D(7)=-(N-L-4)*(N-L-5)*D(5)/(12*L+42.D0)
            SUMs=0.0D0
            
C*********CALCULO DE LA FUNCION HIPERGEOMETRICA ******************
            
            DO 51 S=0,3
               IF(D(2*S+1).EQ.0.0D0)GOTO 51
               
               IF (PE.NE.0.0D0)THEN
                  PEPE=PE**(2*S+L)
               ELSE
                  IF((2*S+L).EQ.0)THEN
                     PEPE=1.0D0
                  ELSE
                     PEPE=0.0D0
                  ENDIF
               ENDIF
               
               ax=a**(2*n-2*s-l+0.5)
               ay=(a*a+pe*pe)**(n+1)
               SUMs=SUMs+D(2*S+1)*ax*PEPE/ay
               
 51         CONTINUE
            
C**********FIN DEL CALCULO DE LA FUNCION HIPERGEOMETRICA ***********
            
            SUMind=SUMind+SUMs*(2**(N+L))*C*(FACT(N+L+1))/DSQRT(1.0D0*
     +           FACT(2*N))
 60      END DO
         
C***********FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ***************
         
         SUMind=SUMind*(FACT(L))/FACT(2*L+1)
         SUMind=SUMind/3.141592654
         GAMMAF=GAMMAF+SUMind*SUMind
 52   END DO
 50   END DO
      
C***************FIN DE LA SUMA EN LOS ORBITALES *********************

      gammasub=gammaF

      RETURN





 111  continue


	      bcf=0.d0

	      do 10 l=lmin,lmax
         if((nsub.eq.-8).and.(l.eq.-lsub))goto 10



         if(nsub*(8-nsub).gt.0)then
            nmin=nsub
            nmax=nsub
            else
            nmin=l+1
            nmax=l+norb(z,q,l)
            endif
             

                 do 11 nn=nmin,nmax

		 nenl=mzi(z,q,nn,l)

		    dens=gam_nl(z,nn,l,pe)

		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    gammasub=bcf

	      return

 222          n=(-jx)/10
              l=-jx-10*n
              gammasub=gam_nl(Z,n,l,pe)

              return

              end


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function COMPTONSUB(JZ,JX,nlsub,q)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


c
c  Perfil Compton J(q) de subcapa(s) 'nlsub' en 
c  atomo caracterizado por (JZ,JX)
c
c  Norma: J(0)=(<p**-1>_nl)/2

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      

	pi=4.d0*datan(1.d0)

	   hnorm=0.d0

	DELTA=25.d0

           do 1 iv=0,31

	VMAX=q+iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=gammasub(JZ,JX,nlsub,r(l))

		   hnorm=hnorm+R(L)*dens1*W*delta

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

 1      continue


	comptonsub=hnorm*3.141592654d0

        return

	END




CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function BRSUB(JZ,JX,nlsub,q)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


c
c  Factor de forma reciproco (numerico) B(q) de subcapa(s) 'nlsub' en
c  atomo caracterizado por (JZ,JX)
C
C  Norma: BRSUB(0)=Nsub

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      

	pi=4.d0*datan(1.d0)


        if(q.eq.0.d0)then
           brsub=1.d0*nesub(jz,jx,nlsub)
           return
           endif


	   hnorm=0.d0

	DELTA=1.d0

        ivmax=60

           do 1 iv=0,ivmax

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=gammasub(JZ,JX,nlsub,r(l))

		   hnorm=hnorm+R(L)*dens1*W*delta*dsin(q*R(L))

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

 1      continue

	brsub=hnorm*3.141592654d0*2.d0/q

        return

	END







C#############################################################
C
C
C                   ######################################
C                   ######################################
C                   ##                                  ##
C                   ##  DERIVADAS DE LA DENSIDAD TOTAL  ##
C                   ##                                  ##
C                   ######################################
C                   ######################################
C
C
C
C#############################################################

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	DERDEN(JZ,JX,R,ORD,ID)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


C     DERIVADA ORD-ESIMA DE LA DENSIDAD 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,ord
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      
C  LS = CONFIGURACION ESTADO FUNDAMENTAL
C  NORB = NUMERO DE ORBITALES CON UN 'L' DADO
C  NFUN = NUMERO DE FUNCIONES DE LA BASE DE SLATER PARA 'L' DADO
C  NZ,AZ = PARAMETROS DE LAS FUNCIONES DE SLATER
C  CZ = COEFICIENTES DE LOS DESARROLLOS EN FUNCIONES DE SLATER
C  MZI = NUMERO DE ELECTRONES EN UN ORBITAL DADO
C  NL = DESIGNACION DE LA CAPA (N,L)


      IF(ID.EQ.0)DERDEN=RON(JZ,JX,R,ORD)
      IF(ID.EQ.1)DERDEN=GAMMAD(JZ,JX,R,ORD)
     
      RETURN

      END

      
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	RON(JZ,JX,R,ORD)
      
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DERIVADA ORD-ESIMA DE LA DENSIDAD DE CARGA ATOMICA 

c     NORMA DE RO(R) = N

      
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      IMPLICIT	REAL*8	(A-H,O-Z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,ORD
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      
      z=abs(jz)


      if(jx.gt.1)then
         q=z-jx
         ne=jx
         else
            q=jx
            ne=z-jx
            endif

c      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111



      RON=0.0D0
      
C***************SUMATORIA EN LOS ORBITALES **************************
      
      DO 50 l=0,3
         do 51 nn=l+1,l+norb(z,q,l)
         SUMind=0.0d0
         
C***************SUMATORIA EN LAS FUNCIONES DE SLATER ****************
         
         DO 60 ind1=1,nfun(z,q,l)
            DO 61 ind2=1,nfun(z,q,l)
               N1=NZ(z,q,l,ind1)
               N2=NZ(z,q,l,ind2)
               A1=AZ(z,q,l,ind1)
               A2=AZ(z,q,l,ind2)
               C1=CZ(z,q,nn,l,ind1)
               C2=CZ(z,q,nn,l,ind2)
               N=N1+N2-2
               A=A1+A2
               C=C1*C2
               
               IF(ORD.Lt.N)THEN
                  LIM=ORD
                  rr=r**(n-lim)
               ELSE
                  LIM=N
                  rr=1.d0
               ENDIF
               
               TK=((-A)**(ORD-LIM))*rr*FACT(ORD)*FACT(N)/
     #              (FACT(LIM)*FACT(ORD-LIM)*FACT(N-LIM))
               ERRE=TK
               
               DO 65 K=LIM-1,0,-1
                  TK=-TK*A*(K+1.d0)*R/((ORD-K)*(N-K))
                  ERRE=ERRE+TK
 65            END DO
               
               SUMind=SUMind+C*((2*A1)**(N1+.5))*
     #              ((2*A2)**(N2+.5))*ERRE*DEXP(-A*R)/
     #              (DSQRT(FACT(2*N1)*FACT(2*N2)))
 61         END DO
 60      END DO
         
C***********FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ****************
         
         RON=RON+SUMind
 51   END DO
 50   END DO
      
C***********FIN DE LA SUMA EN LOS ORBITALES **************************
      
      RON=RON/(4.d0*3.141592654d0)


      RETURN
      

 111  continue


	
              return

              end






CXXXXXXXXXXXXXXXXXXXXXXXXXXXXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      REAL*8	FUNCTION	GAMMAD(JZ,JX,PE,ord)
      
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     DENSIDAD DE CARGA ATOMICA (NORMALIZADA A N)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,s,ord
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3),d(7)

      z=abs(jz)

      if(jx.gt.1)then
         ne=jx
         q=z-jx
         else
            ne=z-jx
            q=jx
            endif


      if(jx.lt.-9)goto 222

      if(jz.lt.0)goto 111


      GAMMAd=0.0D0
      
C**********SUMA EN LOS ORBITALES *********************************
      
             DO 50 l=0,3
             do 52 nn=l+1,l+norb(z,q,l) 


         sind=0d0
	 dersind=0d0
	 der2sind=0d0
C**********SUMA EN LAS FUNCIONES DE SLATER ***********************
         
         DO 60 ind=1,nfun(z,q,l)
            N=NZ(z,q,l,ind)*1.0D0
            A=AZ(z,q,l,ind)*1.0D0
            C=CZ(z,q,nn,l,ind)*1.0D0



	do k=0,int((n-l)/2)
	sind=sind+c*(2*A)**(N+.5)/FACT(2*N)**(.5d0)
     &	*fact(n-l)*(2*a)**n
     &	*(pe/a)**l*fact(n-k)/(fact(k)*
     &  fact(n-l-2*k)*(pe**2+a**2)**(n+1-k))
     &  *(-1/(4*a**2))**k*(2.d0/3.141592654D0)**0.5  

	dersind=dersind+c*(2*A)**(N+.5)/FACT(2*N)**(.5d0)
     &	*(-1/(4*a**2))**k*(2.d0/3.141592654D0)**0.5*fact(n-l)*(2*a)**n
     &	*fact(n-k)/(fact(k)*fact(n-l-2*k)*a**l)
     &  *(l*pe**(l-1)*(pe**2+a**2)-2*(n+1-k)*pe**(l+1))
     &  /((pe**2+a**2)**(n+2-k))

	der2sind=der2sind+c*(2*A)**(N+.5)/FACT(2*N)**(.5d0)
     &	*(-1/(4*a**2))**k*(2.d0/3.141592654D0)**0.5*fact(n-l)*(2*a)**n
     &	*fact(n-k)/(fact(k)*fact(n-l-2*k)*a**l)
     &  *((pe**2+a**2)*(l*(l-1)*(pe**2+a**2)*pe**(l-2)
     &  +2*l*pe**l-2*(l+1)*(n+1-k)*pe**l)-2*pe*(n-k+2)*
     &  (l*pe**(l-1)*(pe**2+a**2)-2*(n+1-k)*pe**(l+1)))
     &  /((pe**2+a**2)**(n+3-k))
	enddo

 60      END DO
         
C***********FIN DE LA SUMA EN LAS FUNCIONES DE SLATER ***************
         

	 if(ord.eq.1)GAMMAd=GAMMAd+derSind*Sind/(2.d0*3.141592654d0*jz)
	 if(ord.eq.2) GAMMAd=GAMMAd+2*(derSind**2+der2sind*Sind)
     &   /(4.d0*3.141592654d0*jz)
 52   END DO
 50   END DO
      
C***************FIN DE LA SUMA EN LOS ORBITALES *********************


      RETURN

 111  continue

	      bcf=0.d0

	      do 10 l=0,3
                 do 11 nn=l+1,l+norb(z,q,l)

		 nenl=mzi(z,q,nn,l)

		    dens=gam_nl(z,nn,l,pe)

		    bcf=bcf+nenl*dens

 11              continue

 10		    continue

                    gammad=bcf

	      return

 222          n=(-jx)/10
              l=-jx-10*n
              gammad=gam_nl(Z,n,l,pe)

              return

              end



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function GAMN_NL(JZ,n,l,p,nesima)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c Derivada nesima de gamma(p) de Hidrogeno (JZ,nl,l) 
c (nesima=1)   

      implicit real*8 (a-h,o-z)

      real*8 knl


      pi=4.d0*datan(1.d0)

      x=(jz*jz-n*n*p*p)/(jz*jz+n*n*p*p)
      knl=pot(2.d0,2*l-1)*n*n*n*n*fact(n-l-1)
      knl=knl/(jz*jz*jz*pi*pi*fact(n+l))

      sum=0.d0

      do 10 m1=0,(n-l-1)/2

      do 20 m2=0,(n-l-1)/2

      term=(l+4.d0)*pot(x,2*(n-l-1-m1-m2))*pot(1.d0-x,l)*pot(1.d0+x,l+3)
      if(l.gt.0)term=term-l*pot(x,2*(n-l-1-m1-m2))*pot(1.d0-x,l-1)
     #   *pot(1.d0+x,l+4)
      if((n-l-1-m1-m2).gt.0)term=term+2.d0*(n-l-1-m1-m2)
     #   *pot(x,2*(n-l-1-m1-m2)-1)*pot(1.d0-x,l)*pot(1.d0+x,l+4)

      c=pot(-1.d0,m1+m2)*pot(4.d0,n-l-1-m1-m2)*fact(n-m1-1)*fact(n-m2-1)
      c=c/(fact(m1)*fact(m2)*fact(n-l-1-2*m1)*fact(n-l-1-2*m2))

      sum=sum+c*term

 20    continue

 10    continue


       gamn_nl=-knl*sum*n/(1.d0*jz)*dsqrt((1.d0-x)*pot(1.d0+x,3))

       return

       end

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function DLNGAM_NL(Z,nn,ll,p)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c derivada logaritmica de densidad gammaf(p) para hidrogenoide (n,l) carga Z
c

      implicit real*8 (a-h,o-z)

      integer z

      x=(1.d0*z*z-nn*nn*p*p)/(1.d0*z*z+nn*nn*p*p)

      pi=4.d0*datan(1.d0)

      dlngamnl=0.d0


      snum=0.d0

      if(nn.eq.(ll+1))goto 11

      do 10 m=0,(nn-ll-2)/2

      term=pot(-1.d0,m)*pot(2.d0*x,nn-ll-2-2*m)*fact(nn-m-1)/fact(m)

         term=term/fact(nn-ll-2-2*m)

         snum=snum+term

 10         continue


 11         snum=snum/fact(ll+1)


      sden=0.d0

      do 20 i=0,(nn-ll-1)/2

      term=pot(-1.d0,i)*pot(2.d0*x,nn-ll-1-2*i)*fact(nn-i-1)/fact(i)

         term=term/fact(nn-ll-1-2*i)

         sden=sden+term

 20         continue


            sden=sden/fact(ll)


            dlngamnl=(4.d0-(2.d0*ll+4.d0)*x)/(1-x*x)

            dlngamnl=dlngamnl+4.d0*(ll+1.d0)*snum/sden

            dlngam_nl=-nn*dlngamnl*(1.d0+x)*dsqrt(1.d0-x*x)/z


            return

            end




CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function DFFHNL(JZ,n,l,vk)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c Derivada del factor de forma F(vk) del atomo de Hidrogeno (JZ,n,l)
c

      implicit real*8 (a-h,o-z)

      vk=vk/JZ

      sum=0.d0

      do 10 i=0,n-l-1

         do 20 m=0,n-l-1

            do 30 j=0,(m+i+1)/2+l

               term=pot(-1.d0,m+i+j)*fact(m+i+2*l+1)*fact(m+i+2*l+2)
               term=term/(fact(n-l-1-m)*fact(m+2*l+1)*fact(n-l-1-i))
         term=term/(fact(i+2*l+1)*pot(1.d0+0.25d0*vk*vk*n*n,m+i+2*l+2))
         term=term*pot(0.25d0*vk*vk*n*n,j)/fact(2*j+1)
         term=term/(fact(m+i+2*l-2*j+1)*fact(m)*fact(i))

         sum=sum+term

 30      continue

 20      continue

 10      continue

         dffhnl=sum*0.5d0*fact(n-l-1)*fact(n+l)/n

      return

      end








C#############################################################
C
C
C                        #################################
C                        #################################
C                        ##                             ##
C                        ##  VALORES ESPERADOS TOTALES  ##
C                        ##                             ##
C                        #################################
C                        #################################
C
C
C#############################################################

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function RADMOM(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Valor esperado radial de orden alfa
c de la densidad 'idens' caracterizada por (JZ,JX)
C con norma 1 para rho(r), gammaf(p) y h(u)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)


	   hnorm=0.d0

           ivmax=10.d0+90.d0*idens
           delta=5.d0+195.d0*idens

           do 1 iv=0,ivmax

              term=0.d0

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=den(JZ,JX,r(l),IDENS)
                      term=term+(R(L)**(alfa+2.d0))*dens1*W*delta

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

		   hnorm=hnorm+term

 1      enddo

	radmom=hnorm*2.d0*pi

      if((idens.eq.0).or.(idens.eq.1))radmom=radmom/
     #  (1.d0*nesub(jz,jx,84))

c        if(idens.eq.-2)radmom=radmom/(8.d0*pi*pi*pi*gammaf(jz,jx,0.d0))
           
c        if(idens.eq.2)radmom=3.d0*radmom/(radmom(jz,jx,2.d0,1)*(alfa+1.d0))

c        if(idens.eq.3)radmom=radmom/(8.d0*pi*pi*pi*ro(jz,jx,0.d0))
           
c        if(idens.eq.4)radmom=radmom/(8.d0*pi*pi*pi*hu(jz,0.d0))
           
        return

	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function EVLOG(JZ,JX,kk,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Valor esperado logaritmico de orden kk 
c de la densidad caracterizada por (JZ,JX)
c con norma 1
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=10.d0+90.d0*idens
           delta=5.d0+195.d0*idens

           do 1 iv=0,ivmax

              term=0.d0

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=den(JZ,JX,r(l),IDENS)

                      pot=1.d0

                      do 50 ii=1,kk
                         pot=pot*dlog(r(l))
 50                      continue

                         term=term+R(L)*R(L)*dens1*W*delta*pot

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

		   hnorm=hnorm+term

 1      continue



	evlog=hnorm*2.d0*pi

       if((idens.eq.0).or.(idens.eq.1))evlog=evlog/
     #  (1.d0*nesub(jz,jx,84))

c        if(idens.eq.-2)evlog=evlog/(8.d0*pi*pi*pi*gammaf(jz,jx,0.d0))

c        if(idens.eq.2)evlog=3.d0*evlog/(2.d0*pi*radmom(jz,jx,2.d0,1))

c        if(idens.eq.3)evlog=evlog/(8.d0*pi*pi*pi*ro(jz,jx,0.d0))

c        if(idens.eq.4)evlog=evlog/(8.d0*pi*pi*pi*hu(jz,0.d0))

	END







CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	REAL*8	FUNCTION	CUMUL(JZ,JX,DELTA)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C  CARGA ACUMULADA HASTA EL RADIO DELTA 
c  DEL SISTEMA CARACTERIZADO POR (JZ,JX)
C  CUMUL(infinito)=N

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

	pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           delta=delta/8.d0

           do 1 iv=0,7

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=ro(JZ,JX,r(l))


       if(dens1.gt.0.d0)hnorm=hnorm+R(L)*R(L)*dens1*W*delta

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

 1      continue


	cumul=hnorm*2.d0*pi


	END


C#####################################################################
C
C
C
C                        ######################################
C                        ######################################
C                        ##                                  ##
C                        ##  VALORES ESPERADOS POR SUBCAPAS  ##
C                        ##                                  ##
C                        ######################################
C                        ######################################
C
C
C#############################################################

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function RADMOMSUB(JZ,JX,nlsub,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Valor esperado radial de orden alfa
c de la densidad caracterizada por (JZ,JX)
C
   
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      
      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=10.d0+90.d0*idens
           delta=5.d0+195.d0*idens

           do 1 iv=0,ivmax

              term=0.d0

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=densub(JZ,JX,nlsub,r(l),IDENS)
                      term=term+(R(L)**(alfa+2.d0))*dens1*W*delta

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

		   hnorm=hnorm+term
 1      enddo

	radmomsub=hnorm*2.d0*pi

       if((idens.eq.0).or.(idens.eq.1))radmomsub=radmomsub/
     #(1.d0*nesub(jz,jx,nlsub))

c       if(idens.eq.-2)radmomsub=radmomsub/
c     #  (8.d0*pi*pi*pi*gammasub(jz,jx,nlsub,0.d0))
           
c        if(idens.eq.2)radmom=3.d0*radmom/(radmom(jz,jx,2.d0,1)*(alfa+1.d0))
c	   hnorm=0.d0

c        if(idens.eq.3)radmomsub=radmomsub/
c     #  (8.d0*pi*pi*pi*rosub(jz,jx,nlsub,0.d0))
           
	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function EVLOGSUB(JZ,JX,nlsub,kk,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Valor esperado logaritmico de orden kk 
c de la densidad caracterizada por (JZ,JX)
C 


C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      
      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=10.d0+90.d0*idens
           delta=5.d0+195.d0*idens

           do 1 iv=0,ivmax

              term=0.d0

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=densub(JZ,JX,nlsub,r(l),IDENS)

                      pot=1.d0

                      do 50 ii=1,kk
                         pot=pot*dlog(r(l))
 50                      continue

                         term=term+R(L)*R(L)*dens1*W*delta*pot

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

		   hnorm=hnorm+term

 1      enddo



	evlogsub=hnorm*2.d0*pi

       if((idens.eq.0).or.(idens.eq.1))evlogsub=evlogsub/
     #(1.d0*nesub(jz,jx,nlsub))

c        if(idens.eq.-2)evlogsub=evlogsub/
c     # (8.d0*pi*pi*pi*gammasub(jz,jx,nlsub,0.d0))

c        if(idens.eq.2)evlogsub=3.d0*evlogsub/
c     #  (2.d0*pi*radmomsub(jz,jx,nlsub,2.d0,1))
	   
c        if(idens.eq.3)evlogsub=evlogsub/
c     # (8.d0*pi*pi*pi*rosub(jz,jx,nlsub,0.d0))


	END






CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	REAL*8	FUNCTION	CUMULSUB(JZ,JX,nlsub,DELTA)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C  CARGA ACUMULADA HASTA EL RADIO DELTA 
c  DEL SISTEMA CARACTERIZADO POR (JZ,JX)
C
c  CUMULSUB(infinito)=Nsub

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)
      

	pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           delta=delta/8.d0

           do 1 iv=0,7

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=rosub(JZ,JX,nlsub,r(l))


       if(dens1.gt.0.d0)hnorm=hnorm+R(L)*R(L)*dens1*W*delta

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

 1      continue


	cumulsub=hnorm*2.d0*pi


	END


C#####################################################################
C
C
C
C                   ########################################
C                   ########################################
C                   ##                                    ##
C                   ##  FUNCIONALES DE LA DENSIDAD TOTAL  ##
C                   ##                                    ##
C                   ########################################
C                   ########################################
C
C
C#############################################################


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         real*8 function FRECMOM(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Momento de frecuencia de orden alfa 
c de la densidad 'idens' caracterizada por (JZ,JX)


C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

      a=0.0D0

      zrp=1.d0*JZ
      delta=1.d0/(zrp*alfa)
      if(idens.eq.1)then
        zrp=1.d0/zrp
        delta=1.d-2/(zrp*alfa)
        endif



       hnorm=0.d0

       test=1.d0

           
        iv=1

        do 1 while(test.gt.1.d-12)

c        write(*,*)'test=',test,iv,idens,JZ

    

        term=0.d0

         DO 20 K=1,48
         X=ABCIS(K)
         W=WEIGHT(K)

            R(1)=0.5d0*(DELTA*X+a)
            R(2)=0.5d0*(DELTA*(1.d0-X)+a)
        DO 88 L=1,2


        dens1=den(JZ,JX,R(L),IDENS)
                      term=term+R(L)*R(L)*(dens1**alfa)*W

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88       CONTINUE
20       CONTINUE

                   hnorm=hnorm+term*delta

                   test=term/hnorm

                   

                   iv=iv+1

                   a=a+delta

                   delta=delta*1.1d0

1       end do



         frecmom=hnorm*pi


c        write(*,*)jz,idens,iv,frecmom,test

         return

         END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         real*8 function FRECMOMDOB(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Momento de frecuencia de orden alfa 
c de la densidad DE PARES 'idens' caracterizada por (JZ,JX)
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R1(2),R2(2)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION test(0:1)

      pi=4.d0*datan(1.d0)

      iesp=idens-2*(idens/2)


      test(0)=1.d-6
      test(1)=1.d-6

      NI=10


        zrp=pot(1.d0*JZ,1-2*iesp)

      delta1=1.0d0/(NI*dsqrt(zrp)*alfa)
      if(iesp.eq.1)then
            if(alfa.lt.0.7d0)then
                delta1=20.0d0/(NI*dsqrt(zrp)*pot(alfa-3.d0/8.d0,3))
                else
        delta1=20.0d0/(NI*zrp*pot(4.d0*alfa-1.d0,2))
            endif
        endif
        delta10=delta1
        



      a1=0.0D0

       sumi1=0.d0

       test1=1.d0

       iv1=1

        do 1 while(test1.gt.test(iesp))


            sumk1=0.d0


       do 20 K1=1,48
        X1=ABCIS(K1)
        W1=WEIGHT(K1)

        suml1=0.d0

         do 88 L1=1,2

            R1(L1)=a1+0.5d0*DELTA1*(1.d0+(2.d0*L1-3.d0)*X1)

           

c empieza suma en segunda variable (con sufijo '2')    


      delta2=1.0d0/(NI*dsqrt(zrp)*alfa)
      if(iesp.eq.1)then
            if(alfa.lt.0.7d0)then
                delta2=20.0d0/(NI*dsqrt(zrp)*pot(alfa-3.d0/8.d0,3))
                else
        delta2=20.0d0/(NI*zrp*pot(4.d0*alfa-1.d0,2))
            endif
        endif
        delta20=delta2

    
       
    


        a2=0.0D0

       sumi2=0.d0

       test2=1.d0

       iv2=1

        do 2 while(test2.gt.test(iesp))


    

        sumk2=0.d0

         DO 22 K2=1,48
         X2=ABCIS(K2)
         W2=WEIGHT(K2)


         suml2=0.d0

        DO 82 L2=1,2

             R2(L2)=a2+0.5d0*DELTA2*(1.d0+(2.d0*L2-3.d0)*X2)


        densry=dendob(JZ,JX,R1(L1),R2(L2),IDENS)
C         if(idens.eq.0)then
C             densry=pot(zrp,6)*dexp(-2.d0*JZ*(R1(L1)+R2(L2)))/(pi*pi)
C             else
C        densry=64.d0*pot(zrp,10)/
C      &pot(pi*(zrp*zrp+R1(L1)*R1(L1))*(zrp*zrp+R2(L2)*R2(L2)),4)
C             endif

                      suml2=suml2+R2(L2)*R2(L2)*(densry**alfa)


82       CONTINUE

            sumk2=sumk2+W2*suml2

22       CONTINUE

                   sumi2=sumi2+sumk2*delta2

                   test2=sumk2*delta2/sumi2



                   iv2=iv2+1

            

                   a2=a2+delta2


      if(iv2.eq.(NI+1))delta2=3*pot(zrp,1-iesp)*delta2


2       end do


c                     
            

            suml1=suml1+R1(L1)*R1(L1)*sumi2

88          continue

        sumk1=sumk1+W1*suml1

20          continue            


                   sumi1=sumi1+delta1*sumk1

                    test1=delta1*sumk1/sumi1


                                     
C                                     v=(a1+delta1)*jz
    
C        exac=((pi/pot(1.d0*JZ,3))**(1.d0-alfa))/pot(alfa,3)
C        exac=exac*exac
C        exac=exac*(1.d0-(.5d0*v*v+v+1.d0)*dexp(-v))

C         exac=pi**(2.d0-2.d0*alfa)
C         exac=exac/pot(alfa*(jz**(1.d0-alfa)),6)
C                    exac=exac*(1.d0-(.5d0*v*v+v+1.d0)*dexp(-v))

C                    coc=sumi1*pi*pi/exac

       write(*,900)'idens=',idens,'Z=',JZ,'q=',alfa,'iv1=',iv1,'a1=',a1,
     &  'test1=',test1,'Norma= ',sumi1*4.d0*pi*pi
c     ,'Exacta= ',exac,'Coc= ',coc

900   format(1x,a6,i1,1x,a3,i3,2x,a2,f5.3,1x,a4,i4,1x,a3,f8.2,4x,a6,
     &  1pd8.2,4x,a7,1pd13.7,2x,a8,1pd13.7,2x,a5,1pd13.7)    



                    iv1=iv1+1

                    a1=a1+DELTA1


      if(iv1.eq.(NI+1))delta1=3*pot(zrp,1-iesp)*delta1

1       end do

       




         frecmomdob=sumi1*4.d0*pi*pi




         return

         END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         real*8 function QSMq1y2(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Medida QSM de orden alfa entre el producto de densidades a 1 cuerpo  
c y la densidad DE PARES en 'idens' caracterizadas por (JZ,JX)
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R1(2),R2(2)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION test(0:1)

      pi=4.d0*datan(1.d0)



      test(0)=1.d-6
      test(1)=1.d-6

      NI=10


        zrp=pot(1.d0*JZ,1-2*idens)

      delta1=1.0d0/(NI*alfa*2.d0*dsqrt(zrp))
      if(idens.eq.1)then
            if(alfa.lt.0.7d0)then
                delta1=20.0d0/(NI*dsqrt(zrp)*pot(alfa-3.d0/8.d0,3))
                else
        delta1=20.0d0/(NI*zrp*pot(4.d0*alfa-1.d0,2))
            endif
        endif
        delta10=delta1
        



      a1=0.0D0

       sumi1=0.d0

       test1=1.d0

       iv1=1

        do 1 while(test1.gt.test(idens))

c            print*, 'iv1=',iv1,'a1=',a1,'test1=',test1

            sumk1=0.d0


       do 20 K1=1,48
        X1=ABCIS(K1)
        W1=WEIGHT(K1)

        suml1=0.d0

         do 88 L1=1,2

            R1(L1)=a1+0.5d0*DELTA1*(1.d0+(2.d0*L1-3.d0)*X1)

           

c empieza suma en segunda variable (con sufijo '2')    


      delta2=1.0d0/(NI*alfa*2.d0*dsqrt(zrp))
      if(idens.eq.1)then
            if(alfa.lt.0.7d0)then
                delta2=20.0d0/(NI*dsqrt(zrp)*pot(alfa-3.d0/8.d0,3))
                else
        delta2=20.0d0/(NI*zrp*pot(4.d0*alfa-1.d0,2))
            endif
        endif
        delta20=delta2

    
       
    


        a2=0.0D0

       sumi2=0.d0

       test2=1.d0

       iv2=1

        do 2 while(test2.gt.test(idens))

c        print*, 'iv1=',iv1,'a1=',a1,'iv2=',iv2,'a2=',a2,'test2=',test2
    

        sumk2=0.d0

         DO 22 K2=1,48
         X2=ABCIS(K2)
         W2=WEIGHT(K2)


         suml2=0.d0

        DO 82 L2=1,2

             R2(L2)=a2+0.5d0*DELTA2*(1.d0+(2.d0*L2-3.d0)*X2)


        densprod=dendob(JZ,JX,R1(L1),R2(L2),IDENS+2)
        dens12=dendob(JZ,JX,R1(L1),R2(L2),IDENS)
        densint=densprod*dens12

C         if(idens.eq.0)then
C             densry=pot(zrp,6)*dexp(-2.d0*JZ*(R1(L1)+R2(L2)))/(pi*pi)
C             else
C        densry=64.d0*pot(zrp,10)/
C      &pot(pi*(zrp*zrp+R1(L1)*R1(L1))*(zrp*zrp+R2(L2)*R2(L2)),4)
C             endif

                      suml2=suml2+R2(L2)*R2(L2)*(densint**(0.5d0*alfa))


82       CONTINUE

            sumk2=sumk2+W2*suml2

22       CONTINUE

                   sumi2=sumi2+sumk2*delta2

                   test2=sumk2*delta2/sumi2



                   iv2=iv2+1



            

                   a2=a2+delta2


      if(iv2.eq.(NI+1))delta2=3.d0*pot(zrp,1-idens)*delta2


2       end do


c                     
            

            suml1=suml1+R1(L1)*R1(L1)*sumi2

88          continue

        sumk1=sumk1+W1*suml1

20          continue            


                   sumi1=sumi1+delta1*sumk1

                    test1=delta1*sumk1/sumi1


                                     
C                                     v=(a1+delta1)*jz
    
C        exac=((pi/pot(1.d0*JZ,3))**(1.d0-alfa))/pot(alfa,3)
C        exac=exac*exac
C        exac=exac*(1.d0-(.5d0*v*v+v+1.d0)*dexp(-v))

C         exac=pi**(2.d0-2.d0*alfa)
C         exac=exac/pot(alfa*(jz**(1.d0-alfa)),6)
C                    exac=exac*(1.d0-(.5d0*v*v+v+1.d0)*dexp(-v))

C                    coc=sumi1*pi*pi/exac

       write(*,900)'idens=',idens,'Z=',JZ,'q=',alfa,'iv1=',iv1,'a1=',a1,
     &  'test1=',test1,'Norma= ',sumi1*4.d0*pi*pi




c     ,'Exacta= ',exac,'Coc= ',coc

900   format(1x,a6,i1,1x,a3,i3,2x,a2,f5.3,1x,a4,i4,1x,a3,f8.2,4x,a6,
     &  1pd8.2,4x,a7,1pd13.7,2x,a8,1pd13.7,2x,a5,1pd13.7)    



                    iv1=iv1+1

                    

                    a1=a1+DELTA1


      if(iv1.eq.(NI+1))delta1=3.d0*pot(zrp,1-idens)*delta1

1       end do

       




         QSMq1y2=sumi1*4.d0*pi*pi




         return

         END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         real*8 function SHANNONDOB(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Momento de frecuencia de orden alfa 
c de la densidad DE PARES 'idens' caracterizada por (JZ,JX)
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R1(2),R2(2)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION test(0:1)

      pi=4.d0*datan(1.d0)

      iesp=idens-2*(idens/2)


      test(0)=1.d-6
      test(1)=1.d-6

      NI=10


        zrp=pot(1.d0*JZ,1-2*iesp)

      delta1=1.0d0/(NI*dsqrt(zrp))
      if(iesp.eq.1)delta1=20.0d0/(NI*zrp*9.d0)
        delta10=delta1
        



      a1=0.0D0

       sumi1=0.d0

       test1=1.d0

       iv1=1

        do 1 while((test1.gt.test(iesp)).or.(iv1.lt.15))


            sumk1=0.d0


       do 20 K1=1,48
        X1=ABCIS(K1)
        W1=WEIGHT(K1)

        suml1=0.d0

         do 88 L1=1,2

            R1(L1)=a1+0.5d0*DELTA1*(1.d0+(2.d0*L1-3.d0)*X1)

           

c empieza suma en segunda variable (con sufijo '2')    


      delta2=1.0d0/(NI*dsqrt(zrp))
      if(iesp.eq.1)delta2=20.0d0/(NI*zrp*9.d0)
    
        delta20=delta2

    
       
    


        a2=0.0D0

       sumi2=0.d0

       test2=1.d0

       iv2=1

        do 2 while((test2.gt.test(iesp)).or.(iv2.lt.15))


    

        sumk2=0.d0

         DO 22 K2=1,48
         X2=ABCIS(K2)
         W2=WEIGHT(K2)


         suml2=0.d0

        DO 82 L2=1,2

             R2(L2)=a2+0.5d0*DELTA2*(1.d0+(2.d0*L2-3.d0)*X2)


        densry=dendob(JZ,JX,R1(L1),R2(L2),IDENS)
        if(densry.lt.1.d-14)densry=1.d0
C         if(idens.eq.0)then
C             densry=pot(zrp,6)*dexp(-2.d0*JZ*(R1(L1)+R2(L2)))/(pi*pi)
C             else
C        densry=64.d0*pot(zrp,10)/
C      &pot(pi*(zrp*zrp+R1(L1)*R1(L1))*(zrp*zrp+R2(L2)*R2(L2)),4)
C             endif

                      suml2=suml2+R2(L2)*R2(L2)*densry*dlog(densry)


82       CONTINUE

            sumk2=sumk2+W2*suml2

22       CONTINUE

                   sumi2=sumi2+sumk2*delta2

                   test2=dabs(sumk2*delta2/sumi2)



                   iv2=iv2+1

            

                   a2=a2+delta2


      if(iv2.eq.(NI+1))delta2=3*pot(zrp,1-iesp)*delta2


2       end do


c                     
            

            suml1=suml1+R1(L1)*R1(L1)*sumi2

88          continue

        sumk1=sumk1+W1*suml1

20          continue            


                   sumi1=sumi1+delta1*sumk1

                    test1=dabs(delta1*sumk1/sumi1)


                                     
C                                     v=(a1+delta1)*jz
    
C        exac=((pi/pot(1.d0*JZ,3))**(1.d0-alfa))/pot(alfa,3)
C        exac=exac*exac
C        exac=exac*(1.d0-(.5d0*v*v+v+1.d0)*dexp(-v))

C         exac=pi**(2.d0-2.d0*alfa)
C         exac=exac/pot(alfa*(jz**(1.d0-alfa)),6)
C                    exac=exac*(1.d0-(.5d0*v*v+v+1.d0)*dexp(-v))

C                    coc=sumi1*pi*pi/exac

       write(*,900)'idens=',idens,'Z=',JZ,'iv1=',iv1,'a1=',a1,
     &  'test1=',test1,'Shannon= ',-sumi1*4.d0*pi*pi
c     ,'Exacta= ',exac,'Coc= ',coc

900   format(1x,a6,i1,1x,a3,i3,2x,a4,i4,1x,a3,f8.2,4x,a6,
     &  1pd8.2,4x,a9,1pd14.7,2x,a8,1pd13.7,2x,a5,1pd13.7)    



                    iv1=iv1+1

                    a1=a1+DELTA1


      if(iv1.eq.(NI+1))delta1=3*pot(zrp,1-iesp)*delta1

1       end do

       




         shannondob=-sumi1*4.d0*pi*pi




         return

         END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         real*8 function DJS1y2(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Divergencia JSD entre la densidad 'idens' a 2 cuerpos y el producto 
c de densidades a 1 cuerpo para el sistema (JZ,JX)
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R1(2),R2(2)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION test(0:1)

      pi=4.d0*datan(1.d0)


      test(0)=1.d-6
      test(1)=1.d-6

      NI=10


        zrp=pot(1.d0*JZ,1-2*idens)

      delta1=1.0d0/(NI*dsqrt(zrp))
      if(idens.eq.1)delta1=20.0d0/(NI*zrp*9.d0)
        delta10=delta1
        



      a1=0.0D0

       sumi1=0.d0

       test1=1.d0

       iv1=1

        do 1 while((test1.gt.test(idens)).or.(iv1.lt.15))


            sumk1=0.d0


       do 20 K1=1,48
        X1=ABCIS(K1)
        W1=WEIGHT(K1)

        suml1=0.d0

         do 88 L1=1,2

            R1(L1)=a1+0.5d0*DELTA1*(1.d0+(2.d0*L1-3.d0)*X1)

           

c empieza suma en segunda variable (con sufijo '2')    


      delta2=1.0d0/(NI*dsqrt(zrp))
      if(idens.eq.1)delta2=20.0d0/(NI*zrp*9.d0)
    
        delta20=delta2

    
       
    


        a2=0.0D0

       sumi2=0.d0

       test2=1.d0

       iv2=1

        do 2 while((test2.gt.test(idens)).or.(iv2.lt.15))


    

        sumk2=0.d0

         DO 22 K2=1,48
         X2=ABCIS(K2)
         W2=WEIGHT(K2)


         suml2=0.d0

        DO 82 L2=1,2

             R2(L2)=a2+0.5d0*DELTA2*(1.d0+(2.d0*L2-3.d0)*X2)


        dens2=dendob(JZ,JX,R1(L1),R2(L2),IDENS)
        dens12=den(JZ,JX,R1(L1),IDENS)*den(JZ,JX,R2(L2),IDENS)
        densmed=0.5d0*(dens2+dens12)

        densint=0.d0
        if(dens2.gt.1.d-14)densint=densint+0.5d0*dens2*dlog(dens2)
        if(dens12.gt.1.d-14)densint=densint+0.5d0*dens12*dlog(dens12)
        if(densmed.gt.1.d-14)densint=densint-densmed*dlog(densmed)



                      suml2=suml2+R2(L2)*R2(L2)*densint


82       CONTINUE

            sumk2=sumk2+W2*suml2

22       CONTINUE

                   sumi2=sumi2+sumk2*delta2

                   test2=dabs(sumk2*delta2/sumi2)



                   iv2=iv2+1

            

                   a2=a2+delta2


      if(iv2.eq.(NI+1))delta2=3*pot(zrp,1-idens)*delta2


2       end do


c                     
            

            suml1=suml1+R1(L1)*R1(L1)*sumi2

88          continue

        sumk1=sumk1+W1*suml1

20          continue            


                   sumi1=sumi1+delta1*sumk1

                    test1=dabs(delta1*sumk1/sumi1)


                                     

       write(*,900)'idens=',idens,'Z=',JZ,'iv1=',iv1,'a1=',a1,
     &  'test1=',test1,'JSD1y2= ',sumi1*4.d0*pi*pi


900   format(1x,a6,i1,1x,a3,i3,2x,a4,i4,1x,a3,f8.2,4x,a6,
     &  1pd8.2,4x,a8,1pd14.7,2x,a8,1pd13.7,2x,a5,1pd13.7)    



                    iv1=iv1+1

                    a1=a1+DELTA1


      if(iv1.eq.(NI+1))delta1=3*pot(zrp,1-idens)*delta1

1       end do

       




         DJS1y2=sumi1*4.d0*pi*pi




         return

         END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

         real*8 function DJT1y2(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Divergencia JTD_alfa entre la densidad 'idens' a 2 cuerpos y producto 
c de densidades a 1 cuerpo para el sistema (JZ,JX)
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R1(2),R2(2)
      DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION test(0:1)


      if(alfa.eq.1.d0)then
        DJT1y2= DJS1y2(JZ,JX,idens)
        goto 333
      endif

      pi=4.d0*datan(1.d0)


      test(0)=1.d-6
      test(1)=1.d-6

      NI=10


        zrp=pot(1.d0*JZ,1-2*idens)

      delta1=1.0d0/(NI*dsqrt(zrp))
      if(idens.eq.1)delta1=20.0d0/(NI*zrp*9.d0)
        delta10=delta1
        



      a1=0.0D0

       sumi1=0.d0

       test1=1.d0

       iv1=1

        do 1 while((test1.gt.test(idens)).or.(iv1.lt.15))


            sumk1=0.d0


       do 20 K1=1,48
        X1=ABCIS(K1)
        W1=WEIGHT(K1)

        suml1=0.d0

         do 88 L1=1,2

            R1(L1)=a1+0.5d0*DELTA1*(1.d0+(2.d0*L1-3.d0)*X1)

           

c empieza suma en segunda variable (con sufijo '2')    


      delta2=1.0d0/(NI*dsqrt(zrp))
      if(idens.eq.1)delta2=20.0d0/(NI*zrp*9.d0)
    
        delta20=delta2

    
       
    


        a2=0.0D0

       sumi2=0.d0

       test2=1.d0

       iv2=1

        do 2 while((test2.gt.test(idens)).or.(iv2.lt.15))


    

        sumk2=0.d0

         DO 22 K2=1,48
         X2=ABCIS(K2)
         W2=WEIGHT(K2)


         suml2=0.d0

        DO 82 L2=1,2

             R2(L2)=a2+0.5d0*DELTA2*(1.d0+(2.d0*L2-3.d0)*X2)


        dens2=dendob(JZ,JX,R1(L1),R2(L2),IDENS)
        dens12=den(JZ,JX,R1(L1),IDENS)*den(JZ,JX,R2(L2),IDENS)
        densmed=0.5d0*(dens2+dens12)

        densint=0.d0
        densint=densint+0.5d0*dens2**alfa
        densint=densint+0.5d0*dens12**alfa
        densint=densint-densmed**alfa



                      suml2=suml2+R2(L2)*R2(L2)*densint


82       CONTINUE

            sumk2=sumk2+W2*suml2

22       CONTINUE

                   sumi2=sumi2+sumk2*delta2

                   test2=dabs(sumk2*delta2/sumi2)



                   iv2=iv2+1

            

                   a2=a2+delta2


      if(iv2.eq.(NI+1))delta2=3*pot(zrp,1-idens)*delta2


2       end do


c                     
            

            suml1=suml1+R1(L1)*R1(L1)*sumi2

88          continue

        sumk1=sumk1+W1*suml1

20          continue            


                   sumi1=sumi1+delta1*sumk1

                    test1=dabs(delta1*sumk1/sumi1)


                                     

       write(*,900)'idens=',idens,'Z=',JZ,'q=',alfa,'iv1=',iv1,'a1=',a1,
     &  'test1=',test1,'JTD1y2= ',sumi1*4.d0*pi*pi


900   format(1x,a6,i1,1x,a3,i3,2x,a3,f7.5,2x,a4,i4,1x,a3,f8.2,4x,a6,
     &  1pd8.2,4x,a8,1pd14.7,2x,a8,1pd13.7,2x,a5,1pd13.7)    



                    iv1=iv1+1

                    a1=a1+DELTA1


      if(iv1.eq.(NI+1))delta1=3*pot(zrp,1-idens)*delta1

1       end do

       




         DJT1y2=sumi1*4.d0*pi*pi/(alfa-1.d0)


333      continue

         return

         END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function SHANNON(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropia de Shannon de la densidad 'idens'
c caracterizada por (JZ,JX) con norma 1
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=350d0!10.d0+90.d0*idens
           delta=10d0!5.d0+195.d0*idens

           do 1 iv=0,ivmax

	VMAX=iv*delta

        term=0.d0

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=den(JZ,JX,r(l),IDENS)
        if(dens1.gt.1.d-48)term=term+R(L)*R(L)*dens1*W*delta*dlog(dens1)

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

        		   hnorm=hnorm+term
 1      enddo

	shannon=-hnorm*2.d0*pi

c       if((idens.eq.0).or.(idens.eq.1))then
c           el=1.d0*nesub(jz,jx,84)
c           shannon=shannon/el+dlog(el)
c           endif

c        if(idens.eq.-2)then
c           vnorm=8.d0*pi*pi*pi*gammaf(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.2)then
c           vnorm=2.d0*pi*radmom(jz,jx,2.d0,1)/3.d0
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.3)then
c           vnorm=8.d0*pi*pi*pi*ro(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.4)then
c           vnorm=8.d0*pi*pi*pi*hu(jz,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

        return

	END




CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function HNORMADOB(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Norma de la densidad 'idens'
c caracterizada por (JZ,JX) (Supuestamente norma 1)
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

      z=1.d0*jz
      rz=z**0.5d0

	   hnorm=0.d0

	  ivmax=8
	  deltaw=6.d0*rz
      delta0=10.d0*rz
          


        rmin1=0.d0
        delta1=delta0

           do 1 iv1=0,ivmax


        rmin=0.d0
        delta=delta0

	  do 2 iv=0,ivmax
		term=0.d0



        
	DO J=1,48
	X1=ABCIS(J)
	W1=WEIGHT(J)

	R1(1)=rmin1+(0.5d0*(1.d0+X1))*DELTA1
	R1(2)=rmin1+(0.5d0*(1.d0-X1))*DELTA1

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=rmin+(0.5d0*(1.d0+X))*DELTA
	R(2)=rmin+(0.5d0*(1.d0-X))*DELTA
	DO L1=1,2
		DO 88 L=1,2


 	DENS1DOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

        term=term+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DENS1DOB





88	CONTINUE
	ENDDO
20	CONTINUE

	ENDDO

                hnorm=hnorm+term

                delta=deltaw
            rmin=delta0+iv*delta
            
              
2	ENDDO 

            delta1=deltaw
            rmin1=delta0+iv1*deltaw

1       ENDDO

    	hnormadob=hnorm*4.d0*pi*pi




        return
 
	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        real*8 function HMUTP(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Información mútua de (JZ,JX) en espacio p
c 
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

        CHARACTER         LS*2,LL*1,NL*2
         INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
        COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
         COMMON /GSS/ ABCIS,WEIGHT
        DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
         DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
          DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
         DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
         DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

         pi=4.d0*datan(1.d0)

         z=1.d0*jz
         rz=z**0.5d0

          hnorm=0.d0

         ivmax=8
         deltaw=6.d0*rz
         delta0=10.d0*rz
          


           rmin1=0.d0
           delta1=delta0

           do 1 iv1=0,ivmax


          rmin=0.d0
          delta=delta0

         do 2 iv=0,ivmax
          term=0.d0



        
        DO J=1,48
         X1=ABCIS(J)
         W1=WEIGHT(J)

        R1(1)=rmin1+(0.5d0*(1.d0+X1))*DELTA1
         R1(2)=rmin1+(0.5d0*(1.d0-X1))*DELTA1

        DO 20 K=1,48
         X=ABCIS(K)
          W=WEIGHT(K)

        R(1)=rmin+(0.5d0*(1.d0+X))*DELTA
        R(2)=rmin+(0.5d0*(1.d0-X))*DELTA

        DO L1=1,2

        den1=gammaf(JZ,JX,R1(L1))

        DO 88 L=1,2

        den0=gammaf(JZ,JX,R(L))

         DENS1DOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

         coc=dlog(DENS1DOB/(den0*den1))

c         print*,dens1dob,den0,den1,coc

        term=term+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DENS1DOB*coc





88      CONTINUE
        ENDDO
20      CONTINUE

         ENDDO

                hnorm=hnorm+term

                delta=deltaw
            rmin=delta0+iv*delta
            
              
2       ENDDO 

            delta1=deltaw
            rmin1=delta0+iv1*deltaw

1       ENDDO



          hmutp=hnorm*4.d0*pi*pi




        return
 
        END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        real*8 function HMUTRP(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Información mútua de (JZ,JX) en espacios r (idens=0) ó p (idens=1)
c 
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

        CHARACTER         LS*2,LL*1,NL*2
         INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

         COMMON /MUTUA/ HMUTUA,PI
        COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
         COMMON /GSS/ ABCIS,WEIGHT
        DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
         DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
          DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
         DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
         DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

         DIMENSION HMUTUA(0:1)

         pi=4.d0*datan(1.d0)

         z=1.d0*jz
         rz=z**0.5d0
         rz0=z**0.45d0
         rz1=z**1.1d0

c          hnorm0=0.d0
          hnorm1=0.d0

         ivmax=4
         deltaw=6.d0*rz*idens+(20.d0/rz0)*(1-idens)
         delta0=10.d0*rz*idens+(2.d0/rz1)*(1-idens)
          


           rmin1=0.d0
           delta1=delta0

           do 1 iv1=0,ivmax


          rmin=0.d0
          delta=delta0

         do 2 iv=0,ivmax

         print*,jz,iv,iv1,hnorm1*4.d0*pi*pi
    
c          term0=0.d0
          TERM1=0.D0



        
        DO J=1,48
         X1=ABCIS(J)
         W1=WEIGHT(J)

        R1(1)=rmin1+(0.5d0*(1.d0+X1))*DELTA1
         R1(2)=rmin1+(0.5d0*(1.d0-X1))*DELTA1

        DO 20 K=1,48
         X=ABCIS(K)
          W=WEIGHT(K)

        R(1)=rmin+(0.5d0*(1.d0+X))*DELTA
        R(2)=rmin+(0.5d0*(1.d0-X))*DELTA

        DO L1=1,2

        if(idens.eq.1)then
            den1=gammaf(JZ,JX,R1(L1))
            else
            den1=ro(JZ,JX,R1(L1))
            endif

        DO 88 L=1,2

        if(idens.eq.1)then
            den0=gammaf(JZ,JX,R(L))
            else
            den0=ro(JZ,JX,R(L))
            endif


         DENS1DOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

         coc=1.d0!dlog(DENS1DOB/(den0*den1))

c            term0=term0+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DENS1DOB

            term1=term1+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DENS1DOB





88      CONTINUE
        ENDDO
20      CONTINUE

         ENDDO

c                hnorm0=hnorm0+term0
                hnorm1=hnorm1+term1

                delta=deltaw
            rmin=delta0+iv*delta
            
              
2       ENDDO 

            delta1=deltaw
            rmin1=delta0+iv1*deltaw

1       ENDDO



          hmutrp=hnorm1*4.d0*pi*pi

          HMUTUA(0)=hnorm1*4.d0*pi*pi
c          HMUTUA(0)=HNORM0*4.D0*PI*PI




        return
 
        END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        real*8 function H2NORM(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Norma-doble de (JZ,JX) en espacios r (idens=0) ó p (idens=1)
c 
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

        CHARACTER         LS*2,LL*1,NL*2
         INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
        COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
         COMMON /GSS/ ABCIS,WEIGHT
        DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
         DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
          DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
         DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
         DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

         pi=4.d0*datan(1.d0)

         z=1.d0*jz
         rz=z**0.5d0
         rz0=z**0.45d0
         rz1=z**1.1d0

          hnorm=0.d0

         ivmax=8
         deltaw=6.d0*rz*idens+(10.d0/z**0.45d0)*(1-idens)
         delta0=10.d0*rz*idens+(2.d0/z**1.1d0)*(1-idens)
          


           rmin1=0.d0
           delta1=delta0

           do 1 iv1=0,ivmax


          rmin=0.d0
          delta=delta0

         do 2 iv=0,ivmax

         print*,jz,iv,iv1,hnorm*4.d0*pi*pi

          term=0.d0



        
        DO J=1,48
         X1=ABCIS(J)
         W1=WEIGHT(J)

        R1(1)=rmin1+(0.5d0*(1.d0+X1))*DELTA1
         R1(2)=rmin1+(0.5d0*(1.d0-X1))*DELTA1

        DO 20 K=1,48
         X=ABCIS(K)
          W=WEIGHT(K)

        R(1)=rmin+(0.5d0*(1.d0+X))*DELTA
        R(2)=rmin+(0.5d0*(1.d0-X))*DELTA

        DO L1=1,2

        

        DO 88 L=1,2

        


         DENS1DOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

c        DENS1DOB=(Z**6)/(pi*pi)*dexp(-2.d0*z*R(L))*dexp(-2.d0*z*R1(L1))

       

        term=term+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DENS1DOB





88      CONTINUE
        ENDDO
20      CONTINUE

         ENDDO

                hnorm=hnorm+term

                delta=deltaw
            rmin=delta0+iv*delta
            
              
2       ENDDO 

            delta1=deltaw
            rmin1=delta0+iv1*deltaw

1       ENDDO



          h2norm=hnorm*4.d0*pi*pi




        return
 
        END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        subroutine infmutua(JZ,JX,IDENS)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Información mútua de (JZ,JX) en espacios r (idens=0) ó p (idens=1)
c 
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

        CHARACTER         LS*2,LL*1,NL*2
         INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1
         integer*4      jz,jx,idens
         real*8        hmutua

         
        COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
         COMMON /GSS/ ABCIS,WEIGHT
         COMMON /HMUT/    HMUTUA
        DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
         DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
          DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
         DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
         DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

         DIMENSION HMUTUA(0:1)
        
        open(10,file='converg.dat',status='unknown')

         pi=4.d0*datan(1.d0)

         z=1.d0*(jz-jx)
         rz=z**0.5d0
         

          hnorm0=0.d0
          hnorm1=0.d0

         ivmax=5
         deltaw=6.d0*rz*idens+(15.d0/(1.d0*z**0.5d0))*(1-idens)
         delta0=10.d0*rz*idens+(2.d-1/(3.d0*z**1.5d0))*(1-idens)
          


           rmin1=0.d0
           delta1=delta0

           do 1 iv1=0,ivmax


          rmin=0.d0
          delta=delta0

         do 2 iv=0,ivmax

            write(*,900)jz,iv,iv1,hnorm0*4.d0*pi*pi,hnorm1*4.d0*pi*pi
         write(10,900)jz,iv,iv1,hnorm0*4.d0*pi*pi,hnorm1*4.d0*pi*pi
    
900      format(3(2x,i3),2(3x,e17.11))

          term0=0.d0
          TERM1=0.D0



        
        DO J=1,48
         X1=ABCIS(J)
         W1=WEIGHT(J)

        R1(1)=rmin1+(0.5d0*(1.d0+X1))*DELTA1
         R1(2)=rmin1+(0.5d0*(1.d0-X1))*DELTA1

        DO 20 K=1,48
         X=ABCIS(K)
          W=WEIGHT(K)

        R(1)=rmin+(0.5d0*(1.d0+X))*DELTA
        R(2)=rmin+(0.5d0*(1.d0-X))*DELTA

        DO L1=1,2

       if(idens.eq.1)then
           den1=gammaf(JZ,JX,R1(L1))
           else
           den1=ro(JZ,JX,R1(L1))
           endif

        DO 88 L=1,2

       if(idens.eq.1)then
            den0=gammaf(JZ,JX,R(L))
           else
            den0=ro(JZ,JX,R(L))
            endif


         DDOB=dendob(JZ,JX,R(L),R1(L1),IDENS)
c         DDOB=z*z*z*z*z*z*dexp(-2.d0*z*(R(L)+R1(L1)))/(pi*pi)
         coc=dlog(DDOB/(den0*den1))

        term0=term0+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DDOB


        term1=term1+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*DDOB*coc



88      CONTINUE
        ENDDO
20      CONTINUE

         ENDDO

                hnorm0=hnorm0+term0
                hnorm1=hnorm1+term1

                delta=deltaw
            rmin=delta0+iv*delta
            
              
2       ENDDO 

            delta1=deltaw
            rmin1=delta0+iv1*deltaw

1       ENDDO





          HMUTUA(1)=hnorm1*4.d0*pi*pi!*z*z*z*z*z*z
          HMUTUA(0)=HNORM0*4.D0*PI*PI!*z*z*z*z*z*z

c         print*,HMUTUA(0),HMUTUA(1)

c         if(hnorm0.lt.0.d0)read(*,*)

        close(10)

        return
 
        END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        real*8 function sdob(JZ,JX,IDENS)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropía de Shannon de (JZ,JX) en espacios r (idens=0) ó p (idens=1)
c 
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

        CHARACTER         LS*2,LL*1,NL*2
         INTEGER       Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1
         integer*4      jz,jx,idens
         

         
        COMMON /koga/ AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
         COMMON /GSS/ ABCIS,WEIGHT

        DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
         DIMENSION LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
          DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
         DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
         DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

         
        
        open(10,file='converg.dat',status='unknown')

         pi=4.d0*datan(1.d0)

         z=1.d0*jz
         rz=z**0.5d0
         

          hnorm1=0.d0

         ivmax=5
         deltaw=6.d0*rz*idens+(15.d0/(1.d0*z**0.5d0))*(1-idens)
         delta0=10.d0*rz*idens+(2.d-1/(3.d0*z**1.5d0))*(1-idens)
          


           rmin1=0.d0
           delta1=delta0

           do 1 iv1=0,ivmax


          rmin=0.d0
          delta=delta0

         do 2 iv=0,ivmax

            write(*,900)jz,iv,iv1,hnorm1*4.d0*pi*pi
         write(10,900)jz,iv,iv1,hnorm1*4.d0*pi*pi
    
900      format(3(2x,i3),2(3x,e17.11))

          
          TERM1=0.D0



        
        DO J=1,48
         X1=ABCIS(J)
         W1=WEIGHT(J)

        R1(1)=rmin1+(0.5d0*(1.d0+X1))*DELTA1
         R1(2)=rmin1+(0.5d0*(1.d0-X1))*DELTA1

        DO 20 K=1,48
         X=ABCIS(K)
          W=WEIGHT(K)

        R(1)=rmin+(0.5d0*(1.d0+X))*DELTA
        R(2)=rmin+(0.5d0*(1.d0-X))*DELTA

        DO L1=1,2

       

        DO 88 L=1,2

       


         DDOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

         vint=DDOB*dlog(DDOB)
        


        term1=term1+W*W1*delta*delta1*R(L)*R(L)*R1(L1)*R1(L1)*vint



88      CONTINUE
        ENDDO
20      CONTINUE

         ENDDO

                
                hnorm1=hnorm1-term1

                delta=deltaw
            rmin=delta0+iv*delta
            
              
2       ENDDO 

            delta1=deltaw
            rmin1=delta0+iv1*deltaw

1       ENDDO





          sdob=hnorm1*4.d0*pi*pi
          
            write(*,900)jz,iv,iv1,sdob
         write(10,900)jz,iv,iv1,sdob

        close(10)

        return
 
        END





CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function reservaSHANNONDOB(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropia de Shannon de la densidad 'idens'
c caracterizada por (JZ,JX) con norma 1
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

c	  ivmax=110d0!10.d0+90.d0*MOD(idens,2) 
c	  delta=20d0!5d0+195.d0*MOD(idens,2)

	ivmax=5*(4+idens*(JZ/25))
	delta=4.d0                      

c	   ivmax1=110d0!10.d0+90.d0*MOD(idens,2)
c           delta1=20d0!5d0+195.d0*MOD(idens,2)

	ivmax1=5*(4+idens*(JZ/25))
	delta1=4.d0  


           do 1 iv1=0,ivmax1

	  do 2 iv=0,ivmax
		term=0.d0



        
	DO J=1,48
	X1=ABCIS(J)
	W1=WEIGHT(J)

	R1(1)=iv1*delta1+0.5d0*DELTA1*(1.d0+X1)
	R1(2)=iv1*delta1+0.5d0*DELTA1*(1.d0-X1)

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=iv*delta+0.5d0*DELTA*(1.d0+X)
	R(2)=iv*delta+0.5d0*DELTA*(1.d0-X)
	DO L1=1,2
		DO 88 L=1,2


	DENS1DOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

		      
        if(dens1dob.gt.1.d-48)term=term+W*W1*delta**2*(R(L)*R1(L1))**2
     & 	*DENS1DOB*dlog(DENS1DOB)

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
	enddo
20	CONTINUE

	ENDDO

                hnorm=hnorm+term
2	ENDDO 
1      enddo

	shannondob=-hnorm*(2.d0*pi)**2

C        if((idens.eq.0).or.(idens.eq.1))then
C           el=1.d0*nesub(jz,jx,84)
C           shannon=shannon/el+dlog(el)
C           endif

c        if(idens.eq.-2)then
c           vnorm=8.d0*pi*pi*pi*gammaf(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.2)then
c           vnorm=2.d0*pi*radmom(jz,jx,2.d0,1)/3.d0
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.3)then
c           vnorm=8.d0*pi*pi*pi*ro(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.4)then
c           vnorm=8.d0*pi*pi*pi*hu(jz,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

        return

	END





CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function reservaFRECMOMDOB(JZ,JX,a,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Momento entropico de orden 'a' de la densidad 'idens'
c caracterizada por (JZ,JX) con norma 1
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

c	  ivmax=110d0!10.d0+90.d0*MOD(idens,2) 
c	  delta=20d0!5d0+195.d0*MOD(idens,2)

	ivmax=5*(4+idens*(JZ/25))
	delta=4.d0                      

c	   ivmax1=110d0!10.d0+90.d0*MOD(idens,2)
c           delta1=20d0!5d0+195.d0*MOD(idens,2)

	ivmax1=5*(4+idens*(JZ/25))
	delta1=4.d0  


           do 1 iv1=0,ivmax1

	  do 2 iv=0,ivmax
		term=0.d0



        
	DO J=1,48
	X1=ABCIS(J)
	W1=WEIGHT(J)

	R1(1)=iv1*delta1+0.5d0*DELTA1*(1.d0+X1)
	R1(2)=iv1*delta1+0.5d0*DELTA1*(1.d0-X1)

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=iv*delta+0.5d0*DELTA*(1.d0+X)
	R(2)=iv*delta+0.5d0*DELTA*(1.d0-X)
	DO L1=1,2
		DO 88 L=1,2


	DENS1DOB=dendob(JZ,JX,R(L),R1(L1),IDENS)

		      
        if(dens1dob.gt.1.d-48)term=term+W*W1*delta**2*(R(L)*R1(L1))**2
     & 	*(DENS1DOB**a)

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
	enddo
20	CONTINUE

	ENDDO

                hnorm=hnorm+term
2	ENDDO 
1      enddo

	frecmomdob=hnorm*(2.d0*pi)**2

        return

	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function DJSDOB(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c 
c Divergencia JSD (en espacio idens) entre la densidad a dos cuerpos y el
c producto de las densidades a un cuerpo del sistema (JZ,JX)
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

c	  ivmax=110d0!10.d0+90.d0*MOD(idens,2) 
c	  delta=20d0!5d0+195.d0*MOD(idens,2)

	ivmax=5*(4+idens*(JZ/25))
	delta=4.d0                      

c	   ivmax1=110d0!10.d0+90.d0*MOD(idens,2)
c           delta1=20d0!5d0+195.d0*MOD(idens,2)

	ivmax1=5*(4+idens*(JZ/25))
	delta1=4.d0  


           do 1 iv1=0,ivmax1

	  do 2 iv=0,ivmax
		term=0.d0



        
	DO J=1,48
	X1=ABCIS(J)
	W1=WEIGHT(J)

	R1(1)=iv1*delta1+0.5d0*DELTA1*(1.d0+X1)
	R1(2)=iv1*delta1+0.5d0*DELTA1*(1.d0-X1)

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=iv*delta+0.5d0*DELTA*(1.d0+X)
	R(2)=iv*delta+0.5d0*DELTA*(1.d0-X)
	DO L1=1,2
		DO 88 L=1,2

                   
	fx=den(JZ,JX,R(L),IDENS)
	fy=den(JZ,JX,R1(L1),IDENS)
	gxy=dendob(JZ,JX,R(L),R1(L1),IDENS)
        hxy=0.5d0*(fx*fy+gxy)
		      
        if(gxy.gt.1.d-48)term=term+W*W1*delta**2*(R(L)*R1(L1))**2
     & 	*(0.5d0*(gxy*dlog(gxy)+fx*fy*dlog(fx*fy))-hxy*dlog(hxy))

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88  	CONTINUE
	enddo
20  	CONTINUE

	ENDDO

                hnorm=hnorm+term
2   	ENDDO 
1      enddo

	djsdob=hnorm*(2.d0*pi)**2

C        if((idens.eq.0).or.(idens.eq.1))then
C           el=1.d0*nesub(jz,jx,84)
C           shannon=shannon/el+dlog(el)
C           endif

c        if(idens.eq.-2)then
c           vnorm=8.d0*pi*pi*pi*gammaf(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.2)then
c           vnorm=2.d0*pi*radmom(jz,jx,2.d0,1)/3.d0
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.3)then
c           vnorm=8.d0*pi*pi*pi*ro(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.4)then
c           vnorm=8.d0*pi*pi*pi*hu(jz,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

        return

	END





CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function HMUTINF(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Informacion mutua del sistema (JZ,JX) en espacio idens.
c
c O sea, KL de la densidad a dos cuerpos respecto al producto de 
c las densidades a un cuerpo.
c
 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

      implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,Li,Lj,J,ivi,ivj
c      real*16           fxy,dob,densi,densj
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),Ri(2),Rj(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)



	   sumj=0.d0
           termj=1.d0
              ivj=0
              vmaxj=0.d0
 
              if(idens.eq.1) then
               deltaj=1.d0
              else 
               deltaj=5.d-1
              endif


      do 1 while(((termj/(sumj+termj)).gt.1.d-12).or.(ivj.lt.10))

 
c	VMAXj=ivj*deltaj

         termj=0.d0

	      DO 20 Kj=1,48
	         Xj=ABCIS(Kj)
	         Wj=WEIGHT(Kj)

	         Rj(1)=VMAXj+0.5d0*DELTAj*(1.d0+Xj)
	         Rj(2)=VMAXj+0.5d0*DELTAj*(1.d0-Xj)

		      DO 88 Lj=1,2

	            DENSj=den(JZ,JX,rj(lj),IDENS)


	           sumi=0.d0
               termi=1.d0
               ivi=0

               if(idens.eq.1) then
                deltai=1.d0
               else 
                deltai=5.d-1
               endif

                vmaxi=0.d0

               do 2 while(((termi/(sumi+termi)).gt.1.d-12)
     +             .or.(ivi.lt.10))

c	VMAXi=ivi*deltai

                  termi=0.d0

	               DO 21 Ki=1,48
                  	  Xi=ABCIS(Ki)
	                  Wi=WEIGHT(Ki)

	                  Ri(1)=VMAXi+0.5d0*DELTAi*(1.d0+Xi)
	                  Ri(2)=VMAXi+0.5d0*DELTAi*(1.d0-Xi)

		               DO 89 Li=1,2

c                       fxy=0.d0

	                     DOB=dendob(JZ,JX,ri(li),Rj(Lj),IDENS)
	                     DENSi=den(JZ,JX,ri(li),IDENS)

c            fxy=dob*densi
c            fxy=-dob*dlog(dob)
c            print*,idens,dob,densi,densj,dob/(densi*densj)
c            if(dob.gt.1.d-48)fxy=dob*dlog(dob/(densi*densj))
                        fxy=dob*dlog(dob/(densi*densj))
c            fxy=dob

                        termi=termi+Ri(Li)*Ri(Li)*fxy*Wi*deltai
 89                  CONTINUE

 21               CONTINUE

                  ivi=ivi+1
                  sumi=sumi+termi
                  gx=sumi
                  vmaxi=vmaxi+deltai

c            deltai=deltai+5.d-1


c            print*,idens,vmaxi,vmaxj

 2             end do


               termj=termj+Rj(Lj)*Rj(Lj)*gx*Wj*deltaj
 88         CONTINUE

 20      CONTINUE

            ivj=ivj+1
            sumj=sumj+termj
            vmaxj=vmaxj+deltaj
c           deltaj=deltaj+5.d-1
 1    enddo

	   hmutinf=sumj*4.d0*pi*pi

c        print*,vmaxi,vmaxj

      return
	   END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function HNORM(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Informacion mutua del sistema (JZ,JX) en espacio idens.
c
c O sea, KL de la densidad a dos cuerpos respecto al producto de 
c las densidades a un cuerpo.
c
 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,Li,Lj,J,ivi,ivj

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),Ri(2),Rj(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)



	   sumj=0.d0
           termj=1.d0
              ivj=0

              vmaxj=0.d0

      deltaj=5.d-1

           do 1 while((termj/(sumj+termj)).gt.1.d-7)

 
c	VMAXj=ivj*deltaj

        termj=0.d0

	DO 20 Kj=1,48
	Xj=ABCIS(Kj)
	Wj=WEIGHT(Kj)

	Rj(1)=VMAXj+0.5d0*DELTAj*(1.d0+Xj)
	Rj(2)=VMAXj+0.5d0*DELTAj*(1.d0-Xj)

		DO 88 Lj=1,2

	DENSj=den(JZ,JX,rj(lj),IDENS)


	   sumi=0.d0
           termi=1.d0
              ivi=0
                 deltai=5.d-1


                 vmaxi=0.d0

           do 2 while((termi/(sumi+termi)).gt.1.d-7)

c	VMAXi=ivi*deltai

        termi=0.d0

	DO 21 Ki=1,48
	Xi=ABCIS(Ki)
	Wi=WEIGHT(Ki)

	Ri(1)=VMAXi+0.5d0*DELTAi*(1.d0+Xi)
	Ri(2)=VMAXi+0.5d0*DELTAi*(1.d0-Xi)

		DO 89 Li=1,2

                   fxy=0.d0

	DOB=dendob(JZ,JX,ri(li),Rj(Lj),IDENS)
	DENSi=den(JZ,JX,ri(li),IDENS)

c            if(dob.gt.1.d-48)fxy=dob*dlog(dob/(densi*densj))
            fxy=dob

        termi=termi+Ri(Li)*Ri(Li)*fxy*Wi*deltai
 89   CONTINUE

 21   CONTINUE

            ivi=ivi+1
            sumi=sumi+termi
            gx=sumi
            vmaxi=vmaxi+deltai

            deltai=deltai+5.d-1


c            print*,idens,vmaxi,vmaxj

 2          enddo


        termj=termj+Rj(Lj)*Rj(Lj)*gx*Wj*deltaj
88	CONTINUE

20	CONTINUE

            ivj=ivj+1
            sumj=sumj+termj
            vmaxj=vmaxj+deltaj
            deltaj=deltaj+5.d-1
 1      enddo

	hnorm=sumj*4.d0*pi*pi

        return

	END


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function HJSDDOB(JZ1,JX1,JZ2,JX2,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Divergencia jensen shannon de las densidades, en espacio 'idens', 
c caracterizadas por (JZ1,JX1) y (JZ2,JX2), con norma 1
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,L1,J,iv1

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2),R1(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=110d0!10.d0+90.d0*mod(idens,2)
           delta=20d0!5d0+195.d0*mod(idens,2)

           do 1 iv1=0,ivmax

	  do 2 iv=0,ivmax
		term=0.d0



        
	DO J=1,48
	X1=ABCIS(J)
	W1=WEIGHT(J)

	R1(1)=iv1*delta+0.5d0*DELTA*(1.d0+X1)
	R1(2)=iv1*delta+0.5d0*DELTA*(1.d0-X1)

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=iv*delta+0.5d0*DELTA*(1.d0+X)
	R(2)=iv*delta+0.5d0*DELTA*(1.d0-X)
	DO L1=1,2
		DO 88 L=1,2


	DENS1DOB=(dendob(JZ1,JX1,r(l),R1(L1),IDENS)+
     &  dendob(JZ2,JX2,r(l),R1(L1),IDENS))*0.5D0
		      
        if(dens1dob.gt.1.d-48)term=term+W*W1*delta**2*(R(L)*R1(L1))**2
     & 	*DENS1DOB*dlog(DENS1DOB)

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
	enddo
20	CONTINUE

	ENDDO

                hnorm=hnorm+term
2	ENDDO 
1      enddo

	HJSDdob=-hnorm*(2.d0*pi)**2-0.5D0*(SHANNONDOB(JZ1,JX1,IDENS)+
     &  SHANNONDOB(JZ2,JX2,IDENS))

C        if((idens.eq.0).or.(idens.eq.1))then
C           el=1.d0*nesub(jz,jx,84)
C           shannon=shannon/el+dlog(el)
C           endif

c        if(idens.eq.-2)then
c           vnorm=8.d0*pi*pi*pi*gammaf(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.2)then
c           vnorm=2.d0*pi*radmom(jz,jx,2.d0,1)/3.d0
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.3)then
c           vnorm=8.d0*pi*pi*pi*ro(jz,jx,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

c        if(idens.eq.4)then
c           vnorm=8.d0*pi*pi*pi*hu(jz,0.d0)
c           shannon=shannon/vnorm+dlog(vnorm)
c           endif

        return

	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function FISHER(JZ,JX,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropia de Fisher de la densidad 'idens'
c caracterizada por (JZ,JX)
C con norma 1

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=10.d0+90.d0*idens
           delta=5.d0+195.d0*idens

           do 1 iv=0,ivmax

	VMAX=iv*delta

        term=0.d0

	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens0=den(JZ,JX,r(l),IDENS)
		      dens1=derden(JZ,JX,r(l),1,IDENS)

        if(dens0.gt.(1.d-15))then
           term=term+R(L)*R(L)*(dens1*dens1/dens0)*W*delta
           endif

c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

        		   hnorm=hnorm+term

 1      enddo



	fisher=hnorm*2.d0*pi

       if((idens.eq.0).or.(idens.eq.1))fisher=fisher/
     #  (1.d0*nesub(jz,jx,84))

c        if(idens.eq.-2)then
c           vnorm=8.d0*pi*pi*pi*gammaf(jz,jx,0.d0)
c           fisher=fisher/vnorm
c           endif

c        if(idens.eq.2)then
c           vnorm=2.d0*pi*radmom(jz,jx,2.d0,1)/3.d0
c           fisher=fisher/vnorm
c           endif

c        if(idens.eq.3)then
c           vnorm=8.d0*pi*pi*pi*ro(jz,jx,0.d0)
c           fisher=fisher/vnorm
c           endif

c        if(idens.eq.4)then
c           vnorm=8.d0*pi*pi*pi*hu(jz,0.d0)
c           fisher=fisher/vnorm
c           endif


	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function RENYI(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropia de Renyi de orden alfa 
c de la densidad 'idens' caracterizada por (JZ,JX)
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      if (alfa.eq.1.d0)then
         renyi=shannon(JZ,JX,idens)
         else
         renyi=(dlog(frecmom(JZ,JX,alfa,idens)))/(1.d0-alfa)
            endif


	END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function TSALLIS(JZ,JX,alfa,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropia de Tsallis de orden alfa 
c de la densidad 'idens' caracterizada por (JZ,JX)
C 
C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ
      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)

      if (alfa.eq.1.d0)then
         tsallis=shannon(JZ,JX,idens)
         else
         tsallis=(1.d0-(frecmom(JZ,JX,alfa,idens)))/(alfa-1.d0)
            endif


	END





CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	real*8 function ENTREL(JZ1,JX1,JZ2,JX2,idens)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Entropia relativa de Kullback-Leibler de las densidades 
c caracterizadas por (JZ1,JX1) y (JZ2,JX2) 
C 

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      

      pi=4.d0*datan(1.d0)

	   hnorm=0.d0

           ivmax=10.d0+90.d0*idens
           delta=5.d0+195.d0*idens

           do 1 iv=0,ivmax

              term=0.d0

	VMAX=iv*delta
	DO 20 K=1,48
	X=ABCIS(K)
	W=WEIGHT(K)

	R(1)=VMAX+0.5d0*DELTA*(1.d0+X)
	R(2)=VMAX+0.5d0*DELTA*(1.d0-X)
		DO 88 L=1,2



		      dens1=den(JZ1,JX1,r(l),IDENS)
		      dens2=den(JZ2,JX2,r(l),IDENS)

                      
              term=term+R(L)*R(L)*dens1*W*delta*dlog(dens1/dens2)


c  EN GENERAL, PARA INTEGRAR ENTRE 0 E INFINITO F(R),
C  EL SEGUNDO SUMANDO ES  F(R)/(2*PI)


88	CONTINUE
20	CONTINUE

        hnorm=hnorm+term

 1      enddo



	entrel=hnorm*2.d0*3.141592654

              if((idens.eq.0).or.(idens.eq.1))then
                 el1=1.d0*nesub(jz1,jx1,84)
                 el2=1.d0*nesub(jz2,jx2,84)
                 entrel=entrel/el1+dlog(el2/el1)
                 endif


c        if(idens.eq.-2)then
c           entrel=entrel/(8.d0*pi*pi*pi*gammaf(jz1,jx1,0.d0))
c           entrel=entrel+dlog(gamma(jz2,jx2,0.d0)/gammaf(jz1,jx1,0.d0))
c           endif

c        if(idens.eq.2)then
c           entrel=3.d0*entrel/(2.d0*pi*radmom(jz1,jx1,2.d0,1))
c      entrel=entrel+dlog(radmom(jz2,jx2,2.d0,1)/radmom(jz1,jx1,2.d0,1))
c           endif

c        if(idens.eq.3)then
c           entrel=entrel/(8.d0*pi*pi*pi*ro(jz1,jx1,0.d0))
c           entrel=entrel+dlog(ro(jz2,jx2,0.d0)/ro(jz1,jx1,0.d0))
c           endif

c        if(idens.eq.4)then
c           entrel=entrel/(8.d0*pi*pi*pi*hu(jz1,0.d0))
c           entrel=entrel+dlog(hu(jz2,0.d0)/hu(jz1,0.d0))
c           endif

        return

	END






C#############################################################
C
C
C
C                   #############################
C                   #############################
C                   ##                         ##
C                   ##  SUBRUTINAS AUXILIARES  ##
C                   ##                         ##
C                   #############################
C                   #############################
C
C
C#############################################################


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      FUNCTION EULER(XX)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      implicit real*8 (a-h,o-z)
      dimension cof(6)
      parameter (pi=3.141592653)
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *     -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      if (xx.lt.0.0d0) then
         euler=pi/(rgam(1-xx)*sin(pi*xx))
      else
         if (xx.eq.1.0d0) then
            euler=1.0d0
         else
            X=XX-ONE
            TMP=X+FPF
            TMP=(X+HALF)*LOG(TMP)-TMP
            SER=ONE
            DO 11 J=1,6
               X=X+ONE
               SER=SER+COF(J)/X
 11         CONTINUE
            euler=EXP(TMP+LOG(STP*SER))
         end if
      end if
      RETURN
      END
      


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function FACT(n)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c factorial
c

      fact=1.d0

      if(n.lt.2)return

      do 10 k=n,2,-1

         fact=fact*k

 10      continue

         return

         end




CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function DFACT(n)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c doble factorial
c

      dfact=1.d0

      if(n.lt.2)return

      do 10 k=n,2,-2

         dfact=dfact*k

 10      continue

         return

         end




CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function POT(a,n)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c
c calcula  a**n  con a real y n entero
c

      implicit real*8 (a-h,o-z)
      pot=1.d0
      

      if(n.eq.0)return

      if(a.eq.0.d0)then
         pot=0.d0
         return
         endif

         m=abs(n)

         do 10 k=1,m
            pot=pot*a
 10         continue

            if(n.lt.0)pot=1.d0/pot

            return

            end




CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	SUBROUTINE GAUSS(N,ABCIS,WEIGHT)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c******** SUBRUTINA QUE DA LAS ABSCISAS Y PESOS DE LA INTEGRACION
c******** GAUSSIANA

C
        DOUBLE PRECISION ABCIS(48),WEIGHT(48),ABCP(48),WEIP(48)
C
        DO 10 I=1,48
        ABCIS(I)=0.D0
        WEIGHT(I)=0.D0
   10   CONTINUE
C  
C
C
	IF (N .EQ. 48) THEN
        ABCIS(1)=.016276744849602969579D0
        ABCIS(2)=.048812985136049731112D0
        ABCIS(3)=.081297495464425558994D0
        ABCIS(4)=.113695850110665920911D0
        ABCIS(5)=.145973714654896941989D0
        ABCIS(6)=.178096882367618602759D0
        ABCIS(7)=.210031310460567203603D0
        ABCIS(8)=.241743156163840012328D0
        ABCIS(9)=.273198812591049141487D0
        ABCIS(10)=.304364944354496353024D0
        ABCIS(11)=.335208522892625422616D0
        ABCIS(12)=.365696861472313635031D0
        ABCIS(13)=.395797649828908603285D0
        ABCIS(14)=.425478988407300545365D0
        ABCIS(15)=.454709422167743008636D0
        ABCIS(16)=.483457973920596359768D0
        ABCIS(17)=.511694177154667673586D0
        ABCIS(18)=.539388108324357436227D0
        ABCIS(19)=.566510418561397168404D0
        ABCIS(20)=.593032364777572080684D0
        ABCIS(21)=.618925840125468570386D0
        ABCIS(22)=.644163403784967106798D0
        ABCIS(23)=.668718310043916153953D0
        ABCIS(24)=.692564536642171561344D0
        ABCIS(25)=.715676812348967626225D0
        ABCIS(26)=.738030643744400132851D0
        ABCIS(27)=.759602341176647498703D0
        ABCIS(28)=.780369043867433217604D0
        ABCIS(29)=.800308744139140817229D0
        ABCIS(30)=.819400310737931675539D0
        ABCIS(31)=.837623511228187121494D0
        ABCIS(32)=.854959033434601455463D0
        ABCIS(33)=.871388505909296502874D0
        ABCIS(34)=.886894517402420416057D0
        ABCIS(35)=.901460635315852341319D0
        ABCIS(36)=.915071423120898074206D0
        ABCIS(37)=.927712456722308690965D0
        ABCIS(38)=.939370339752755216932D0
        ABCIS(39)=.950032717784437635756D0
        ABCIS(40)=.959688291448742539300D0
        ABCIS(41)=.968326828463264212174D0
        ABCIS(42)=.975939174585136466453D0
        ABCIS(43)=.982517263563014677447D0
        ABCIS(44)=.988054126329623799481D0
        ABCIS(45)=.992543900323762624572D0
        ABCIS(46)=.995981842987209290650D0
        ABCIS(47)=.998364375863181677724D0
        ABCIS(48)=.999689503883230766828D0
C  
        WEIGHT(1)=.032550614492363166242D0
        WEIGHT(2)=.032516118713868835987D0
        WEIGHT(3)=.032447163714064269364D0
        WEIGHT(4)=.032343822568575928429D0
        WEIGHT(5)=.032206204794030250669D0
        WEIGHT(6)=.032034456231992663218D0
        WEIGHT(7)=.031828758894411006535D0
        WEIGHT(8)=.031589330770727168558D0
        WEIGHT(9)=.031316425596861355813D0
        WEIGHT(10)=.031010332586313837423D0
        WEIGHT(11)=.030671376123669149014D0
        WEIGHT(12)=.030299915420827593794D0
        WEIGHT(13)=.029896344136328385984D0
        WEIGHT(14)=.029461089958167905970D0
        WEIGHT(15)=.028994614150555236543D0
        WEIGHT(16)=.028497411065085385646D0
        WEIGHT(17)=.027970007616848334440D0
        WEIGHT(18)=.027412962726029242823D0
        WEIGHT(19)=.026826866725591762198D0
        WEIGHT(20)=.026212340735672413913D0
        WEIGHT(21)=.025570036005349361499D0
        WEIGHT(22)=.024900633222483610288D0
        WEIGHT(23)=.024204841792364691282D0
        WEIGHT(24)=.023483399085926219842D0
        WEIGHT(25)=.022737069658329374001D0
        WEIGHT(26)=.021966644438744349195D0 
        WEIGHT(27)=.021172939892191298988D0
        WEIGHT(28)=.020356797154333324595D0
        WEIGHT(29)=.019519081140145022410D0
        WEIGHT(30)=.018660679627411467385D0
        WEIGHT(31)=.017782502316045260838D0
        WEIGHT(32)=.016885479864245172450D0
        WEIGHT(33)=.015970562902562291381D0
        WEIGHT(34)=.015038721026994938006D0
        WEIGHT(35)=.014090941772314860916D0
        WEIGHT(36)=.013128229566961572637D0
        WEIGHT(37)=.012151604671088319635D0
        WEIGHT(38)=.011162102099838498591D0
        WEIGHT(39)=.010160770535008415758D0
        WEIGHT(40)=.009148671230783386633D0
        WEIGHT(41)=.008126876925698759217D0
        WEIGHT(42)=.007096470791153865269D0
        WEIGHT(43)=.006058545504235961683D0
        WEIGHT(44)=.005014202742927517693D0
        WEIGHT(45)=.003964554338444686674D0
        WEIGHT(46)=.002910731817934946408D0
        WEIGHT(47)=.001853960788946921732D0
        WEIGHT(48)=.000796792065552012429D0
	ENDIF
C
C
	RETURN
        END



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      integer function KARTOT(JZ,JX)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c Carga total Q del sistema caracterizado por (JZ,JX)

      kartot=jx

      if(jx.lt.-9)kartot=abs(jz)-1
      if(jx.gt.1)kartot=abs(jz)-jx

      return

      end


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      SUBROUTINE TRANSFER
      
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C     TRANSFIERE LOS DATOS DE LA FUNCION DE ONDA DE LOS ATOMOS
C     A MATRICES


C     LS(Z,Q) = ESTADO FUNDAMENTAL DEL SISTEMA (Z,Q) EN HARTREE-FOCK
C     NORB(Z,Q,L) = NUMERO DE ORBITALES 'L' DEL SISTEMA (Z,Q)
C     NFUN(Z,Q,L) = NUMERO DE FUNCIONES SLATER EN BASE 'L' DEL SISTEMA (Z,Q)
C     NZ(Z,Q,L,IND) = EXPONENTE IND-ESIM0 DE R, MAS UNO, EN BASE 'L'
C     AZ(Z,Q,L,IND) = PARAMETRO IND-ESIMO DE LA EXPONENCIAL EN BASE 'L'
C     CZ(Z,Q,N,L,IND) = IND-ESIMO COEFICIENTE DE LA COMBINACION LINEAL 
C                       DEL DESARROLLO DEL ORBITAL (N.L)
C     MZI(Z,Q,N,L) = NUMERO DE ELECTRONES EN ORBITAL (N,L)
C     NEL(Z,Q,L) = NUMERO DE ELECTRONES TIPO 'L'
C     EN(Z,Q) = ENERGIA TOTAL
C     PION(Z,Q) = POTENCIAL DE IONIZACION
C     NLVAL(Z,Q) = NUMEROS CUANTICOS 'N,L' DE LA CAPA DE VALENCIA

      implicit real*8 (a-h,o-z)
      
      CHARACTER         LS*2,LL*1,NL*2,ruta*35
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ,zmin,zmax

      COMMON /energy/   EN,PION
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION EN(103,-1:20),PION(103,-1:20),NLVAL(103,-1:20)
      DIMENSION ruta(5)

C 888      write(*,*)'INDICAR ORDENADOR DE COMPILACION/EJECUCION'
C       WRITE(*,*)
C       WRITE(*,*)' 1: ORDENADOR NUEVO EN DESPACHO'
C       WRITE(*,*)' 2: ORDENADOR ANTIGUO EN DESPACHO'
C       WRITE(*,*)' 3: ORDENADOR DE TORRE EN CASA'
C       WRITE(*,*)' 4: PORTATIL'
C       WRITE(*,*)
C       iord=-20
C       READ(*,*)iord
C       if(iord.eq.0)stop
C       if((iord.lt.1).or.(iord.gt.4))then
C         write(*,*)'ERROR: INDICAR DE NUEVO o 0 PARA TERMINAR'
C         GOTO 888
C         endif


        ruta(1) = './'
        write(*,*)'ruta leida'

        ruta(2) = 'C:/Users/Usuario/Dropbox'        !ANTIGUO EN DESPACHO 
        ruta(3) = 'F:/Programas/Dropbox'       !ORDENADOR TORRE EN CASA
        ruta(4) = 'C:/Users/angul/Dropbox'     !PORTÁTIL 
        ruta(5) = '/home/angulo'                !PROTEUS 

        iord=1


666      OPEN(79,
     +  FILE=
     +trim(ruta(iord))//'COEFS/atoms.ani',
     +     STATUS='OLD',ACCESS='SEQUENTIAL',IOSTAT=ierror)

      if(ierror.ne.0)then
        iord=iord+1
        goto 666
         endif


      OPEN(80,
     +  FILE=
     +trim(ruta(iord))//'COEFS/control.ani',
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=12)

      OPEN(81,
     +  FILE=
     +trim(ruta(iord))//'COEFS/nalfa.ani',
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=13)
      OPEN(82,
     +  FILE=
     +trim(ruta(iord))//'COEFS/c.ani',
     +     STATUS='OLD',ACCESS='DIRECT',RECL=12,FORM='FORMATTED')
      OPEN(83,
     +  FILE=
     +trim(ruta(iord))//'COEFS/energy.ani',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(84,
     +  FILE=
     +trim(ruta(iord))//'COEFS/potioniz.ani',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(85,
     +  FILE=
     +trim(ruta(iord))//'COEFS/valen.ani',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(86,
     +  FILE=
     +trim(ruta(iord))//'COEFS/atoms.neu',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(87,
     +  FILE=
     +trim(ruta(iord))//'COEFS/control.neu',
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=12)
      OPEN(88,
     +  FILE=
     +trim(ruta(iord))//'COEFS/nalfa.neu',
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=13)
      OPEN(89,
     +  FILE=
     +trim(ruta(iord))//'COEFS/c.neu',
     +     STATUS='OLD',ACCESS='DIRECT',RECL=12,FORM='FORMATTED')
      OPEN(90,
     +  FILE=
     +trim(ruta(iord))//'COEFS/energy.neu',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
       OPEN(91,
     +  FILE=
     +trim(ruta(iord))//'COEFS/potioniz.neu',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(92,
     +  FILE=
     +trim(ruta(iord))//'COEFS/valen.neu',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(93,
     +  FILE=
     +trim(ruta(iord))//'COEFS/atoms.cat',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(94,
     +  FILE=
     +trim(ruta(iord))//'COEFS/control.cat',
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=12)
      OPEN(95,
     +  FILE=
     +trim(ruta(iord))//'COEFS/nalfa.cat',
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=13)
      OPEN(96,
     +  FILE=
     +trim(ruta(iord))//'COEFS/c.cat',
     +     STATUS='OLD',ACCESS='DIRECT',RECL=12,FORM='FORMATTED')
      OPEN(97,
     +  FILE=
     +trim(ruta(iord))//'COEFS/energy.cat',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(98,
     +  FILE=
     +trim(ruta(iord))//'COEFS/potioniz.cat',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(99,
     +  FILE=
     +trim(ruta(iord))//'COEFS/valen.cat',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      write(*,*)'abriendo archivos de la tabla periodica'

      do 12 i=2,10

         j=i-10*(i/10)

         iu1=13+6*i
         iu2=14+6*i
         iu3=15+6*i
         iu4=16+6*i
         iu5=17+6*i
         iu6=18+6*i

      OPEN(iu1,
     +  FILE=
     +trim(ruta(iord))//'COEFS/atoms.'
     +     //char(48+j),
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      OPEN(iu2,
     +  FILE=
     +trim(ruta(iord))//'COEFS/control.'
     +  //char(48+j),
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=12)
      OPEN(iu3,
     +   FILE=
     +trim(ruta(iord))//'COEFS/nalfa.'
     +//char(48+j),
     +     STATUS='OLD',ACCESS='DIRECT',FORM='FORMATTED',RECL=13)
      OPEN(iu4,
     +   FILE=
     +trim(ruta(iord))//'COEFS/c.'
     +//char(48+j),
     +     STATUS='OLD',ACCESS='DIRECT',RECL=12,FORM='FORMATTED')
      OPEN(iu5,
     +  FILE=
     +trim(ruta(iord))//'COEFS/energy.'
     +     //char(48+j),
     +     STATUS='OLD',ACCESS='SEQUENTIAL')



 12   continue

      do 6 z=1,103
      do 5 q=-1,20
      do 4 l=0,3
      norb(z,q,l)=0
      nfun(z,q,l)=0
      nel(z,q,l)=0
      do 7 n=1,7
         nl(z,q,n,l)='**'
         mzi(z,q,n,l)=0
 7       continue
 4    continue
 5    continue
 6    continue


      do 10 q=-1,20

         if(q.eq.-1)then
         iun1=79
         iun2=80
         iun3=81
         iun4=82
         iun5=83
         iun6=84
         iun7=85
         zmin=1
         zmax=53
         zz=0
         endif

         if(q.eq.0)then
         iun1=86
         iun2=87
         iun3=88
         iun4=89
         iun5=90
         iun6=91
         iun7=92
         zmin=1
         zmax=103
         zz=0
         endif

         if(q.eq.1)then
         iun1=93
         iun2=94
         iun3=95
         iun4=96
         iun5=97
         iun6=98
         iun7=99
         zmin=3
         zmax=55
         zz=2
         endif

         if(q.gt.1)then
         zmin=q+2
         zmax=q+10
         endif

      do 11 z=zmin,zmax

         if(q.eq.-1)then
         if(z.eq.2)goto 11
         if(z.eq.4)goto 11
         if(z.eq.10)goto 11
         if(z.eq.12)goto 11
         if(z.eq.18)goto 11
         if(z.eq.20)goto 11
         if(z.eq.30)goto 11
         if(z.eq.36)goto 11
         if(z.eq.38)goto 11
         if(z.eq.48)goto 11
         endif

         if(q.gt.1)then
            iun1=13+6*(z-q)
            iun2=14+6*(z-q)
            iun3=15+6*(z-q)
            iun4=16+6*(z-q)
            iun5=17+6*(z-q)
            iun6=18+6*(z-q)
            iun7=92
            endif

 1       REWIND(iun1)
         REWIND(iun5)
         REWIND(iun6)
         REWIND(iun7)

         if(q.lt.2)then

      DO WHILE(ZZ.LT.Z)
         READ(iun1,*)ZZ,LS(z,q),PZ,QZ,SZ
         READ(iun5,*)ZZ,EN(z,q)
         READ(iun6,*)ZZ,PION(z,q)
         READ(iun7,*)ZZ,nval,lval
      END DO

      nlval(z,q)=10*nval+lval

      else

         do 55 ii=1,q+1
         READ(iun1,*)ZZ,LS(z,q),PZ,QZ,SZ
         READ(iun5,*)ZZ,EN(z,q)
 55   continue

      nval=2
      lval=1
      if((z-q).eq.2)nval=1
      if((z-q).lt.5)lval=0

      nlval(z,q)=10*nval+lval

         endif

      do 3 j=0,3
      norb(z,q,j)=0
      nel(z,q,j)=0
 3    continue

      RZ=2

 222  DO WHILE(RZ.NE.0)
         READ(iun2,907,REC=PZ)KZ,RZ,NN,LL,MZ

         pz=pz+1
         IF(LL.EQ.'S')LZ=0
         IF(LL.EQ.'P')LZ=1
         IF(LL.EQ.'D')LZ=2
         IF(LL.EQ.'F')LZ=3
         norb(z,q,lz)=norb(z,q,lz)+1
         
         NL(z,q,nn,lz)=CHAR(48+NN)//LL
         mzi(z,q,nn,lz)=mz
         nel(z,q,lz)=nel(z,q,lz)+mz


         nfun(z,q,lz)=kz

         DO 20 I=1,KZ
            if(rz.gt.-1)then
               READ(iun3,908,REC=QZ)NZ(z,q,lz,i),AZ(z,q,lz,i)
            QZ=QZ+1
            endif
            READ(iun4,909,REC=SZ)CZ(z,q,nn,lz,i)
            CZ(z,q,nn,lz,i)=CZ(z,q,nn,lz,i)*SQRT(MZ*1.0)
            SZ=SZ+1
 20      CONTINUE

      END DO

 11   continue

      
 10   continue




      do 80 k=25,99
         close(k)
 80      enddo

      RETURN
      
 905  FORMAT(I3,1X,A2,1X,I4,1X,I4,1X,I5)
 907  FORMAT(1X,I2,1X,I2,1X,I1,A1,1X,I2)
 908  FORMAT(1X,I1,1X,F10.6)
 909  FORMAT(1X,F11.7)
      END






CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      function RGAM(x)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


      implicit real*8 (a-h,o-z)
      rgam=euler(x)
      return
      end
      

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	integer function NESUB(JZ,JX,nlsub)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


c Numero de electrones de la densidad 
c caracterizada por (JZ,JX,nlsub)

C |JZ| = CARGA NUCLEAR
C SI JZ>0  ROOTHAN-HARTREE-FOCK
C SI JZ<0  BARE COULOMB FIELD
C
C JX = -1  ANIONES (SIMPLES)
C JX =  0  NEUTROS
C JX = +1  CATIONES (SIMPLES)
C JX > +1  SERIE ISOELECTRONICA DE JX ELECTRONES
C JX < -9  HIDROGENO EN ESTADO (N,L) CON -JX = 10*N+L 

        implicit real*8 (a-h,o-z)

      CHARACTER         LS*2,LL*1,NL*2
      INTEGER		Z,Q,RZ,ZZ,PZ,QZ,SZ

      
      COMMON /capas/    NLVAL
      COMMON /koga/	AZ,CZ,NEL,NORB,NFUN,NZ,MZI,NL,LS
      COMMON /GSS/ ABCIS,WEIGHT
      DIMENSION ABCIS(48),WEIGHT(48),R(2)
      DIMENSION	LS(103,-1:20),NORB(103,-1:20,0:3),NFUN(103,-1:20,0:3)
      DIMENSION NZ(103,-1:20,0:3,15),AZ(103,-1:20,0:3,15)
      DIMENSION CZ(103,-1:20,7,0:3,15),NEL(103,-1:20,0:3)
      DIMENSION MZI(103,-1:20,7,0:3),NL(103,-1:20,7,0:3)
      DIMENSION NLVAL(103,-1:20)

      nesub=0

      if(jx.lt.-9)then
      njx=-jx/10
      ljx=-jx-10*njx

      if(nlsub.gt.0)then

         nesub=0

         if(nlsub.eq.84)nesub=1
         if(nlsub.eq.(80+ljx))nesub=1
         if(nlsub.eq.(10*njx+4))nesub=1
         if(nlsub.eq.(10*njx+ljx))nesub=1

         else

            nesub=1

         if(nlsub.eq.-84)nesub=0
         if(nlsub.eq.(-80-ljx))nesub=0
         if(nlsub.eq.(-10*njx-4))nesub=0
         if(nlsub.eq.(-10*njx-ljx))nesub=0

         endif

         return
         endif



         jza=abs(jz)

      if(jx.gt.1)then
         q=jza-jx
         else
            q=jx
         endif

         if(nlsub.gt.0)then
            nesub=0

            else

                  nesub=jza-q

                     endif

                     nsub=nlsub/10
                     lsub=abs(nlsub-10*nsub)
                     nsub=abs(nsub)

                     if(lsub.lt.4)then
                        lmin=lsub
                        lmax=lsub
                        else
                           lmin=0
                           lmax=3
                           endif

                     if(nsub.lt.8)then
                        nmin=nsub
                        nmax=nsub
                        else
                           nmin=1
                           nmax=7
                           endif

                           do 10 l=lmin,lmax

                              do 20 n=nmin,nmax

                                 if(nlsub.gt.0)then
                                    nesub=nesub+mzi(jza,q,n,l)
                                    else
                                    nesub=nesub-mzi(jza,q,n,l)
                                    endif

 20                                 continue

 10                                 continue

                                    return

                                    end



CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function YLM2(l,m,theta)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


c
c modulo cuadrado |Y|**2 del armonico esferico Y_lm(theta)
c

      implicit real*8(a-h,o-z)


      ylm2=-1.d0

      m=abs(m)

      pi=4.d0*datan(1.d0)

      s=dsin(theta)
      c=dcos(theta)

      if(l.eq.0)then
         ylm2=0.25d0/pi
         return
         endif

         if(l.eq.1)then
            if(m.eq.0)then
               ylm2=0.75d0*c*c/pi
               else
                  ylm2=0.375d0*s*s/pi
                  endif
                  return
                  endif

                  if(l.eq.2)then

                     if(m.eq.0)then
                     ylm2=3.d0*c*c-1.d0
                     ylm2=ylm2*ylm2*5.d0/(16.d0*pi)
                     endif

                     if(m.eq.1)ylm2=15.d0*s*s*c*c/(8.d0*pi)
                       
                     if(m.eq.2)ylm2=15.d0*s*s*s*s/(32.d0*pi)

                     return

                     endif

                     if(l.eq.3)then

                        if(m.eq.0)then
                           ylm2=5.d0*c*c-3.d0
                           ylm2=ylm2*ylm2*c*c*7.d0/(16.d0*pi)
                           endif

                           if(m.eq.1)then
                              ylm2=1.d0-5.d0*c*c
                              ylm2=ylm2*ylm2*s*s*21.d0/(64.d0*pi)
                              endif

                              return

                              endif

                              return

                              end


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      real*8 function DYLM2(l,m,theta)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


c
c derivada del modulo cuadrado |Y|**2 del armonico esferico Y_lm(theta)
c
      implicit real*8(a-h,o-z)


      dylm2=-1.d0

      m=abs(m)

      pi=4.d0*datan(1.d0)

      s=dsin(theta)
      c=dcos(theta)

      if(l.eq.0)then
         dylm2=0.d0
         return
         endif

         if(l.eq.1)then
            if(m.eq.0)then
               dylm2=-1.5d0*s*c/pi
               else
                  dylm2=-0.75d0*s*c/pi
                  endif
                  return
                  endif

                  if(l.eq.2)then

                     if(m.eq.0)then
                     dylm2=15.d0*c*(1.d0-3.d0*c*c)/(4.d0*pi)
                     endif

                     if(m.eq.1)dylm2=15.d0*s*c*(c*c-s*s)/(4.d0*pi)
                       
                     if(m.eq.2)dylm2=15.d0*s*s*s*c/(8.d0*pi)

                     return

                     endif

                     if(l.eq.3)then

      if(m.eq.0)dylm2=21.d0*s*c*(20.d0*c*c-25.d0*c*c*c*c-3.d0)/(8.d0*pi)
           
                           if(m.eq.1)then
           dylm2=1.d0-10.d0*c*c+10.d0*s*s+25.d0*c*c*c*c-50.d0*s*s*c*c
           dylm2=dylm2*s*c*21.d0/(32.d0*pi)
                              endif

                              return

                              endif

                              return

                              end






