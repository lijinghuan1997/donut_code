MODULE global
    DOUBLE PRECISION Emodel(500000),Bmodel(500000),Y1model(500000)
    DOUBLE PRECISION Y2model(500000),Y3model(500000),Rmodel(500000),Echannel(23)
    DOUBLE PRECISION pi,Re,lamda,tao_,tao,be,be_,te_,delta,kmodel,w,yibusen,ae,ai,n0,qe,me,v0,D,v00,omiga,omega,paramsm
    DOUBLE PRECISION B00,mmm,dir
    DOUBLE PRECISION eq,ttotal,eq2
    END MODULE
    
    
PROGRAM leakage
      USE global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION vstart(6),vvector(3),notrace(4)
      DOUBLE PRECISION A15(15),PS,AA11(11),BB8(8),XSTART,YSTART,ZSTART
	  DOUBLE PRECISION XSTARTGSM,YSTARTGSM,ZSTARTGSM,PSI
	  REAL PARMOD(10)
	  INTEGER J,IYR,IDAY,IHOUR,MIN,ISEC,IOPT,pwhere
	  INTEGER outputonce,fb,vol,GSEorGSM,vorp,qmsign
      INTEGER m,n,i,num,Tnum
      INTEGER alpha,beta
      CHARACTER *3 str_int(181),title
	
      
	COMMON /B_field/ B1,Bzimf,Vsw
      COMMON /basics/ qm,rad,E0,B0,vc,qmsign
	COMMON /red/ rp,ep,dpp
      COMMON /gyrora/ omega0,vd0
      COMMON /pitchgy/ energy,pitch,phi
	COMMON /ele/ Einv
	COMMON /calcu/ fb,vol,vorp
	COMMON /vvv/ vvector
	COMMON /oput/ outputonce
!	COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
!     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
!     * K,IY,CGST,SGST,BA
	
	
	
       data str_int/ &
     &'000','001','002','003','004','005','006','007','008','009',&
     &'010','011','012','013','014','015','016','017','018','019',&
     &'020','021','022','023','024','025','026','027','028','029',&
     &'030','031','032','033','034','035','036','037','038','039',&
     &'040','041','042','043','044','045','046','047','048','049',&
     &'050','051','052','053','054','055','056','057','058','059',&
     &'060','061','062','063','064','065','066','067','068','069',&
     &'070','071','072','073','074','075','076','077','078','079',&
     &'080','081','082','083','084','085','086','087','088','089',&
     &'090','091','092','093','094','095','096','097','098','099',&
     &'100','101','102','103','104','105','106','107','108','109',&
     &'110','111','112','113','114','115','116','117','118','119',&
     &'120','121','122','123','124','125','126','127','128','129',&
     &'130','131','132','133','134','135','136','137','138','139',&
     &'140','141','142','143','144','145','146','147','148','149',&
     &'150','151','152','153','154','155','156','157','158','159',&
     &'160','161','162','163','164','165','166','167','168','169',&
     &'170','171','172','173','174','175','176','177','178','179','180'/


!  qm - charge mass ratio for e (or qm=5986936.5d0 for O+ or qm=95790984 for H+)
!  qmsign - -1 for e, 1 for O+ or H+ 
      qmin=-175880284000.9d0					   
!      qmin=5986936.5d0
!       qmin=95790984.0d0
      
      IF (qmin.GT.0) THEN
	 qmsign=1
	ELSE
	 qmsign=-1
	ENDIF
	qm=DABS(qmin)

!  the defination of the speed of light:C
	vc=3.0d8

!  use internal function to define pi and radian conversion
      pi=4.0D0*datan(1.0D0)
      rad=pi/180.d0


	  
	  
!  scale parameters
      Re=6.3712d6  !Re - Earth raius in meter
      v0=1.d5	     !(m/s)
      E0=1.d-3     !(V/m)
      B0=1.d-9     !(nT)
!
      rp=qm*E0/v0
      ep=v0*B0/E0
      dpp=v0/Re

      omega0=qm*B0
      vd0=E0/B0

	print *,'Input B1(nT)' 
!	read *,B1
	B1=40
!      Byimf=-6.7
	print *,'Input Bzimf(nT)'
!	read *,Bzimf
!	Bzimf=7
      Bzimf=-10
	Print *,'Input the total time of Magnetic Storm' !磁暴一般持续1-3天
!	read *,tt
      tt=100000  
102	print *,'Electric field valid? '
	print *,'1 - no, 2 - yes '
!	read *,vol
      vol=2
      
	print *, 'Solar Wind Speed'  
!     read *, Vsw
      Vsw=300	
      print *, 'initial energy'
!      read *, energy
!      energy=5
      energy=0.1
!     print *, 'Radial distance (Re),azimu angle(rad), Z'

	print *, 'X,Y,Z(Re)'
!     read *, Req,FAI,vstart(6)
!	read *, Xstart,Ystart,Zstart

!      Xstart=6.34
!	Ystart=-1.4
!	Zstart=-7.32 

105   print *,'How to give the direction?'
	print *,'1 - by vector,such as: (vx/n,vy/n,vz/n)'
!      如果选择1，之后将要输入一个向量，三个分量同速度方向三分量成正比//zhou
	print *,'2 - by pitch angle and gyrophase'
!      如果选择2，之后将要输入投掷角和相位角//zhou
!	read *, vorp
	vorp=1

      if(vorp.NE.1) then
	 if(vorp.NE.2) then
	  print *,'wrong input!'
	  goto 105
	 endif
	endif

	if (vorp.EQ.2) then
      print *, 'Pitch angle (deg):'	 
!      read *, pitch
	read *,pitchh
      print *, 'Enter gyrophase (deg):'  !初始相位角
!      print *, 'Relative to the intersection line of GSM XY plane to'
!      print *, 'the plane that is perependicular to the magnetic field'
!      print *, 'vector and containes the intersection point of the'
!      print *, 'magnetic field line with GSM-XY plane if B_azimuth!=0'
!      print *, 'or relative to X axis if only Bz!=0)'
!      read *, phi
       phi=0
!       phi=180
      elseif(vorp.EQ.1) then
	print *,'Then input the vector please'
!	read *,vvector(1),vvector(2),vvector(3)
	endif

101   print *,'You can choose:'
	print *,'1: forwards, -1: backwards'
!	read *,fb
    fb=-1
	if (fb.NE.1.AND.fb.NE.-1) then
	 print *,'WRONG INPUT with fb'
       goto 101
      endif
      print *, 'Running time (in second):'
!      read *, T_end
!2       T_end=1
	PRINT *, 'Start time'
!	read *, t0
       t0=0

!     可能会输出过多的数据，所以可以考虑每计算多少步输出一次
!     推荐对于北向磁场取10，对于南向磁场取1
	print *, 'How many times calculation will you output once?'
!	read *, outputonce
!	outputonce=1
	outputonce=100

	pitch=pitchh	

      J=-1      

	vstart(4)=Xstart
	vstart(5)=Ystart
	vstart(6)=Zstart
	 
      qe=1.60218d-19
      me= 9.109389699999999d-31
      !  model parameters
	  lamda=1.5d0
      te_=20d0
	  tao_=205.0d0/te_
      be=0.084919730844829
      !be_=0.d0
      be_=0.d0
	  delta=0.4d0
	  kmodel=20d0
	  w=0d0
!	  tao=(290.0d0-240.0d0*delta)/(1-delta)/45.0d0/lamda
        tao=(235-delta*205)/(1-delta)/te_/lamda
        
        paramsm=5.4462d-4
	  yibusen=(w**2+paramsm)/kmodel**2/(w-1)**2
!        yibusen=5.838251281384d-3
	  ae=paramsm/(paramsm+w**2)
!        ae=0.9829
	  ai=w**2/(paramsm+w**2)
!        ai=0.017
	  n0=87d0
      D=4.88927539165936d7
      v00=sqrt(te_*lamda*1.6d-19*2/me)
!      v00=3.9404d6
      omega=0.016873812496050
      omiga=0
      B00=70.654957809081736
!  model results read data from data.txt

	  OPEN(1,FILE='data_miu_diedai_be=6_bz=23_lambda=1.5_k=20.txt')
      DO i=1,500000
        READ(1,*)Rmodel(i),Bmodel(i),Emodel(i),Y1model(i),Y2model(i),Y3model(i)
      ENDDO
      close(1)
      OPEN(50,FILE='jianyan.txt')
      DO i=1,500000
          WRITE(50,*)Bmodel(i)
      ENDDO 
      close(50)
! open files and record initial values of some parameters(结果存入处）
!      OPEN(10,FILE='H_L_E.txt')
!      OPEN(20,FILE='H_log.txt')
      if(vorp.EQ.2) then
	  title=str_int(pitchh+1)
	  OPEN(10,FILE='PA' // title //'H_L_E.txt')	  
        OPEN(20,FILE='PA' // title //'H_log.txt')	
	else
        OPEN(10,FILE='vecH_L_E.txt')
	  OPEN(20,FILE='vecH_log.txt')
      OPEN(40,FILE='FSD.txt')
      open(120,FILE='E+P.txt')
      open(180,FILE='tiduB.txt')
      open(200,FILE='FPSDdelta.txt')
      open(220,FILE='single_particlev.txt')
	endif

      write(20,*) 'Energy in keV, distance in Re, angle in degree, veloc&
     &ity in  m/s, time in sec'
!      write(20,'(A)') 'axe= 0, X0=-8 Re'

! start RK routine to calculate particle trajectory
	       
!       energy=Emin+float(i-1)*DE
!        vstart(4)=-Req*dsin(FAI*rad)   !x
!        vstart(5)=Req*dcos(FAI*rad)	 !y

       WRITE(*,'(/,I3,A)') j,' Start to trace the trajectory of the part&
     &icle'
       WRITE(20,'(/,I2,A)') j,' Start to trace the trajectory of the par&
     &ticle'

       energy=0.025
       eq=1.5d5
       eq2=2*eq
        mmm=75

        dir=1

       DO num=1,32,1

   
            T_end=2.0d0+0.1*(num-1)  !0.00001d0*(num-1)+0.15
            ttotal=T_end

    DO beta_=10,190,10
    DO alpha=45,360,45
            beta=beta_-10
            vvector(1)=sin(real(beta)/180.0*pi)*cos(real(alpha)/180.0*pi)
            vvector(2)=sin(real(beta)/180.0*pi)*sin(real(alpha)/180.0*pi)
            vvector(3)=cos(real(beta)/180.0*pi)
            Xstart=(-15+(num-1))*27600.d0/Re
            Ystart=-1.d4/Re
            Zstart=0.d0/Re;
            vstart(4)=Xstart
	        vstart(5)=Ystart
            vstart(6)=Zstart
        CALL rkdumb(t0,T_end,vstart)	!计算轨道的程序
      END DO
  !      END DO
        END DO
        END DO
      close(10)
      close(20)
      close(40)
      close(120)
      close(180)
      close(200)
      close(220)
!      stop 'Calculation completed!'
 	print*,'Calculation completed!'
	t=t0
199   continue
      end

!  -------- SUBROUTINE RKDUMB -----------
      SUBROUTINE rkdumb(t0,T_end,v)
!  Changed by zhou. I want to use 8 equations instead of 6,which can keep 
!  the conservation of energy. Independent variable is s instead of t,while
!  s=ct/gamma.                            ! //zhou 2002.10.12
  use global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	INTEGER pwhere,fb,vol,vorp
	REAL PARMOD(10)
	INTEGER IYR,IDAY,IHOUR,MIN,ISEC,IOPT,qmsign
!	REAL FAI1,FAI2,FAI3,FAI4
      DIMENSION dv(6),v(6),posi(3),Ef(3),Ec(3),Ecr(3),Ei(3),Bf(3)
      DIMENSION vm(6),notrace(4)
	DIMENSION p(8),dp(8) !//zhou 2002.10.12
      INTEGER*4 k,m

!      COMMON /efd/ Ec
      COMMON /basics/ qm,rad,E0,B0,vc,qmsign
      COMMON /B_field/ B1,Bzimf,Vsw
      COMMON /pitchgy/ energy,pitch,phi
      COMMON /outp/ h_mlt,RLm,enem,tm,vm
	COMMON /calcu/ fb,vol,vorp
      COMMON /FAIn/ FAI1,FAI2,FAI3,FAI4,FAI5,FAI6,FAI7,FAI8


      pi=180.*rad

      Xc=-7.6
      Yc=0.0
      Dx=15.0
      Dy=15.0

      ftf=1.0
      frf=1.0     !15 ??什么意思
      Af=0.0      !98	??

      AA=Af
      ftail=ftf
      frc=frf 

	m=0
      t=t0

      DO 10 i=1,3
	 posi(i)=v(i+3)
10    CONTINUE

! calculate initial magnetic field (计算开始的磁场）
      CALL Bfield(posi,t,Bf)

!  Calculate initial velocity（计算初始速度）
      CALL init_v(energy,pitch,Bf,phi,v)

! Calculate initial s !//zhou 2002.10.12 
      vv=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
	gam=dsqrt(1-(vv*v0*v0/(vc*vc)))
	s=t0*gam*vc

! Calculate initial p by the value of v !//zhou 2002.10.12
      p(1)=(1.d0/gam)*vc  ! p(1)=gamma*c
      DO 3 i=1,3
	 p(i+1)=v(i)
3     continue
      p(5)=vc*t0
	DO 4 i=6,8
	 p(i)=v(i-2)
4     continue
	  
!  Tc - gyroperiod （回旋周期）
      Tc=2.0*pi/(qm*B0*dsqrt(Bf(1)*Bf(1)+Bf(2)*Bf(2)+Bf(3)*Bf(3)))

!  Initial gyroraius (in Re)	（回旋半径）
      rgi=dsqrt(2000.*qm*energy)*Tc/(2.*pi*Re)

      WRITE(*, '(4x,3(A,f7.2),2(A,f6.2,A))') 'Bx=',Bf(1)*B0*1.d9,&
     &    '  By=',Bf(2)*B0*1.d9,'  Bz=',Bf(3)*B0*1.d9,'  Tc=',Tc,& 
     &    ' s','  rc=', rgi,' Re'
 
!  display GSM components of initial velocity and initial GSM coordinate
      WRITE(*,'(4x,A,f6.2,A,3(A,f6.2))') 'Eini=',energy,' keV',&
     &   '  X=',v(4),'  Y=',v(5),'  Z=',v(6) 
      WRITE(*,'(4x,3(A,e11.4))') 'Vx=',v(1)*v0/1.d3,'  Vy=',&
     &     v(2)*v0/1.d3,'  Vz=',v(3)*v0/1.d3

! initial MLT
!      CALL hmlt(v(4),v(5),v(6),h_mlt)
!      WRITE(20,'(A,e12.6,A)') ' (MLT=',h_mlt,')'

      WRITE(20,'(3x,A)') 'Initial position and velocity:'
      WRITE(20,'(5x,3(A,f6.2))') 'X=',v(4),'  Y=',v(5),'  Z=',v(6) 
      WRITE(20,'(5x,3(A,e11.4))') 'Vx=',v(1)*v0/1.d3,'  Vy=', &
     &     v(2)*v0/1.d3,'  Vz=',v(3)*v0/1.d3 
		write(20,18) 'x               ','y             ','z          ' &
     &,'  t       ','  energy 	   ','Ex            ','Ey            ' &
     &,' Ez            ','Bx            ','By            ','Bz         ' &
     &,'pitch angle            ','  potential          ','  latitude of  &
     &landing place''       where?'
18	format(10x,13A)
		
      k=0
	m=0

	DO while(dabs(t).LE.T_end)  !	 设定程序运算时间（可以是时间或运算步数）

! Calculate electric field at each position !//zhou 2002.10.10
       CALL Efield(posi,Bf,Ef,t)

! Save to output files
	 	 
       CALL outf(k,m,t,Ef,Bf,v,ks)
	 if (m.GT.30) then	! 用来控制每3000步输出一次
	    m=0
	 end if

       if(ks.eq.1) goto 2

! Calculate derivatives of v at time t	! 计算t时刻的速度
       CALL derivs(s, p, dp)

! 1/120 or 1/240 of cyclotron period for H+ or O+. Tc=2*pi/(q/m * B_size)
! 1/12 or 1/24 of cyclotron period for e-
       B_size=dsqrt(Bf(1)*Bf(1)+Bf(2)*Bf(2)+Bf(3)*Bf(3))
!      h=0.0261799/(qm*B0*B_size)
       vv=p(2)*p(2)+p(3)*p(3)+p(4)*p(4)
	 gam=dsqrt(1-(vv*v0*v0/(vc*vc)))
!      fb用于正推和反推的控制 //zhou
! 	 h=(fb*0.261799*2/(2*qm*B0*B_size))/50.d0 	 ! h=Tc/240      
       h=fb*0.002d0*(qmsign+1.1d0)/5
       hs=h*vc*gam  ! written by zhou. The step of s. //2002.10.12
! calculate position and velocity at time t+h	 ! 计算t+h时间后的速度和位置
       CALL rk4(s, hs, dp, p)

       s=s+hs ! //zhou 2002.10.12
       t=t+h
       k=k+1
	 m=m+1
!  change dp and p into dv and v !//zhou 2002.10.12
       DO 6 i=1,3
	  dv(i)=dp(i+1)
	  dv(i+3)=dp(i+5)
        v(i)=p(i+1)
	  v(i+3)=p(i+5)
6      continue

! position and magnetic field at time t+h
       DO 5 i=1,3
        posi(i)=v(i+3)
5      CONTINUE

	CALL Bfield(posi,t,Bf)

!       CALL Bfield(posi,t,Bf)
	 
! end of while loop
       end do
!c11	 Continue
      
	write(*,14) ' Trajectory tracing completed. Time over'
      write(10,14) 'Trajectory tracing completed. Time over'

2     write(10,14) 'Position and velocity at Lmin (Lmin=',RLm,'t=',tm,&
     &             ' s):'     
      write(10,15) 'X=',v(4),'   Y=',v(5),'   Z=',v(6)
      write(10,16) 'Vx=',v(1),'  Vy=',v(2),'  Vz=',v(3)
	write(10,19) 'Energy=',energy
      WRITE(10,17) h_mlt,energy,RLm,enem,tm
14    format(3x,A,f5.2,1x,A,f12.6,A)
15    format(5x,3(A,f10.6))
16    format(5x,3(A,e13.6))
17    format(1x,5e13.6)
19    format(5x,A,e13.6)
      return
      end

!  *********** SUBROUTINE RK4 ********************
!     这里应该用相对论公式吧
      SUBROUTINE rk4(s, hs, dp, y)
      use global
!  The subroutine is a fourth order Runge-Kutta integrator with the 
!  adaptive stepsize control

!  Changed by zhou. I want to use 8 equations instead of 6,which can keep 
!  the conservation of 
! Independent variable is s instead of t,while
!  s=ct/gamma.                            ! //zhou 2002.10.12

!  Changed by zhou for the second time. The equation to keep the conservation
!  of energy was changed from c**2*(p(2)**2+p(3)**2+p(4)**2)=W**2-m0**2 into
!  v=c*sqrt(1-c**2/p(1)**2)/v0            ! //zhou 2002.10.14

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION dp(8), y(8), dym(8), dyt(8), yt(8), dydx(8)
      INTEGER qmsign
	COMMON /basics/ qm,rad,E0,B0,vc,qmsign
      
	vc=3.0d8
      hh=hs*0.5d0
      h6=hs/6.0d0
      xh=s+hh

!  first step
      DO 10 i=1,8
        dydx(i)=dp(i)
10    CONTINUE
  	DO 20 i=1,8
        yt(i)=y(i)+hh*dydx(i)
20    CONTINUE

!  second step    
      CALL derivs(xh, yt, dp)
      DO 30 i=1,8
        dyt(i)=dp(i)
30    CONTINUE
      DO 40 i=1,8
        yt(i)=y(i)+hh*dyt(i)
40    CONTINUE

!  Third step
      CALL derivs(xh, yt, dp)
      DO 50 i=1,8
        dym(i)=dp(i)
50    CONTINUE
      DO 60 i=1,8
       yt(i)=y(i)+hs*dym(i)
       dym(i)=dym(i)+dyt(i)
60    CONTINUE
!  Fourth step
      CALL derivs(s+hs, yt, dp)
      DO 70 i=1,8
       dyt(i)=dp(i)
70    CONTINUE									   
!  Accumulate increments
      DO 80 i=1,8
       y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0d0*dym(i))
80    CONTINUE

	spe1=vc*dsqrt(1.d0-vc*vc/(y(1)*y(1)))/v0
	spe2=dsqrt(y(2)*y(2)+y(3)*y(3)+y(4)*y(4))
      EE=spe1/spe2

      DO 90 i=1,3
	 y(i+1)=y(i+1)*EE
90    continue

      return
      end

!  ********** SUBROUTINE DERIVS ********************
      SUBROUTINE derivs(s, p, dp)
    use global
!  The subroutine calculates the  derivatives of position and velocity
!  Note the Lonertz equation has been normalized to be dimentionaless

!  Changed by zhou. I want to use 8 equations instead of 6,which can keep 
!  the conservation of energy. Independent variable is s instead of t,while
!  s=ct/gamma.                            ! //zhou 2002.10.12

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION v(6),dv(6),posi(3),Bf(3),Ef(3),Ec(3),Ecr(3),Ei(3)
	INTEGER pwhere,vol,fb,vorp,notrace(4),qmsign
      DIMENSION p(8),dp(8) !//zhou 2002.10.12

!     COMMON /efd/Ec
      COMMON /basics/ qm,rad,E0,B0,vc,qmsign
	COMMON /B_field/ B1,Bzimf,Vsw
      COMMON /red/ rp,ep,dpp
	COMMON /ele/ Einv
	COMMON /calcu/ fb,vol,vorp
      DO 10 i=1,3
        posi(i)=p(i+5)
10    CONTINUE
      vv=p(2)*p(2)+p(3)*p(3)+p(4)*p(4) !//zhou 2002.10.12
	gam=dsqrt(1-vv*v0*v0/(vc*vc))
      
	t=s/(vc*gam)
!  Calculate mfd and efd at each time and position
      CALL Bfield(posi,t,Bf)
!  Calculate electric field at each position !//zhou 2002.10.10
      CALL Efield(posi,Bf,Ef,t)

!  Dv/Dt rp=q/m E0 /v0 ep=v0 B0/E0
!  相对论条件lorenz力的形式应该不变吧
      dp(1)=qmsign*rp*v0*v0*(p(2)*Ef(1)+p(3)*Ef(2)+p(4)*Ef(3))/(vc*vc*&
     &gam)
      dp(2)=qmsign*rp*(Ef(1)+ep*(p(3)*Bf(3)-p(4)*Bf(2)))/(vc*gam)					       
      dp(3)=qmsign*rp*(Ef(2)+ep*(p(4)*Bf(1)-p(2)*Bf(3)))/(vc*gam)
      dp(4)=qmsign*rp*(Ef(3)+ep*(p(2)*Bf(2)-p(3)*Bf(1)))/(vc*gam)

!  Dr/Dt  dpp=v0/Re	   !用来改变正推和反推的条件
!      dp(5)=dpp/gam
      dp(6)=dpp*p(2)/(vc*gam)
      dp(7)=dpp*p(3)/(vc*gam)
      dp(8)=dpp*p(4)/(vc*gam)

	continue
      return
      end
!  ------------------------------------------------------------------
!  Written by zhou in order to calculate the total electric field in 
!  GSM coordinates !! //zhou 2002.10.10
      SUBROUTINE Efield(posi,Bf,Ef,t)
	USE global
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DOUBLE PRECISION posi(3),Bf(3),Ef(3)
	INTEGER vol,fb,vorp,qmsign
	DOUBLE PRECISION R,angle,X,Y,Z,interp
    COMMON /basics/ qm,rad,E0,B0,vc,qmsign
    COMMON /B_field/ B1,Bzimf,Vsw
	COMMON /ele/ Einv
	COMMON /calcu/ fb,vol,vorp
	  X=posi(1)
      Y=posi(2)
      Z=posi(3)
      if (X==0.d0) then
          if (Y==0.d0) then
              R=0.d0
          else
               R=sqrt(X*X+Y*Y)
          endif
          else
            R=sqrt(X*X+Y*Y)
      endif
      if (R==0.d0) then
          Ef(1)=0.d0
          Ef(2)=0.d0
          Ef(3)=0.d0
          return 
          endif
	  if (X>=0.d0) then
	  angle=asin(Y/R)
	  else 
	  angle=pi-asin(Y/R)
      ENDIF
      R=R*Re
      ZRe=Z*Re
      if (R<1) then
          rr=R/eq
          zz=ZRe/eq2
          Ef(1)=0.d0-rr*exp(-rr*rr)*exp(-zz*zz)*sin(angle)/1000*mmm*dir
          Ef(2)=0.d0+rr*exp(-rr*rr)*exp(-zz*zz)*cos(angle)/1000*mmm*dir
          Ef(3)=0.d0
      else  !20km是存的数据文件的平衡态场的分布范围
        Emo=interp(R,'E ') !径向向里的电场
        rr=R/eq
        zz=ZRe/eq2
	    Ef(1)=Emo*cos(angle)*1000-rr*exp(-rr*rr)*exp(-zz*zz)*sin(angle)/1000*mmm*dir
	    Ef(2)=Emo*sin(angle)*1000+rr*exp(-rr*rr)*exp(-zz*zz)*cos(angle)/1000*mmm*dir
	    Ef(3)=0.d0 
	  endif
	  
	return
      end 
!  ------------------------------------------------------------------
!  Interface for rutines that calculate the magnetic field components
!  in GSM coordinates
      SUBROUTINE Bfield(posi,t,Bf)
      USE global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Bf(3),A15(15),AA11(11),BB8(8)
      REAL PARMOD(10)
	DOUBLE PRECISION kz,L,interp
      COMMON /dateut/ IYR,IDAY,IHOUR,MIN,ISEC,IOPT,PARMOD
      COMMON /GEOPACK/ A15,PS,AA11,K,IY,BB8
      COMMON /timepa/ tt,tp,tf
	COMMON /B_field/ B1,Bzimf,Vsw

      tt=ttotal-abs(t)
      X=posi(1)
	  Y=posi(2)
      Z=posi(3)
      if (X==0.d0) then
          if (Y==0.d0) then
              R=0.d0
          else 
              R=sqrt(X*X+Y*Y)
          endif
          else
      R=sqrt(X*X+Y*Y)
      endif
            if (R==0.d0) then
                rr=R/eq
          zz=ZRe/eq2
          Bf(1)=0.d0
          Bf(2)=0.d0
          Bf(3)=23.d0+(-2.0/eq*exp(-zz*zz)*tt*mmm)*1e3*dir
          return 
          endif
      if (X>=0.d0) then
	  angle=asin(Y/R)
	  else 
	  angle=pi-asin(Y/R)
      ENDIF
      ZRe=Z*Re
	    
      R=R*Re
      if (R<1) then 
          rr=R/eq
          zz=ZRe/eq2
          Bf(1)=0.d0-2*zz/eq2*rr*exp(-rr*rr)*exp(-zz*zz)*cos(angle)*tt*mmm*1e3*dir
          Bf(2)=0.d0-2*zz/eq2*rr*exp(-rr*rr)*exp(-zz*zz)*sin(angle)*tt*mmm*1e3*dir
          Bf(3)=23.d0+(2.0d0/eq*rr*rr*exp(-rr*rr)*exp(-zz*zz)*tt*mmm-2.0/eq*exp(-rr*rr)*exp(-zz*zz)*tt*mmm)*1e3*dir
      else 
        rr=R/eq
        zz=ZRe/eq2
	    Bmo=interp(R,'B ')
        Bf(1)=0.d0-2*zz/eq2*rr*exp(-rr*rr)*exp(-zz*zz)*cos(angle)*tt*mmm*1e3*dir
	    Bf(2)=0.d0-2*zz/eq2*rr*exp(-rr*rr)*exp(-zz*zz)*sin(angle)*tt*mmm*1e3*dir
	    Bf(3)=Bmo+(2.0d0/eq*rr*rr*exp(-rr*rr)*exp(-zz*zz)*tt*mmm-2.0/eq*exp(-rr*rr)*exp(-zz*zz)*tt*mmm)*1e3*dir
      endif
      return
      end

! --------------------------------------
      SUBROUTINE outf(k,m,t,Ef,Bf,v,ks)
    use global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	INTEGER pwhere,outputonce,ppwhere,qmsign
	DOUBLE PRECISION magne(4),LATIT,mu
    DOUBLE PRECISION FPSD,Fzhongjian,y1,y3,R,interp
      DIMENSION Ef(3),Bf(3),v(6),rg(3),pg(3),vm(6),speed(4),vx(6),vy(6),vz(6),xx(6),yy(6),zz(6)
      INTEGER*4 k,m
 
      COMMON /basics/ qm,rad,E0,B0,vc,qmsign
      COMMON /outp/ h_mlt,RLm,enem,tm,vm
	COMMON /Electricfield/ pwhere
	COMMON /ele/ Einv
	COMMON /oput/ outputonce
	COMMON /faii/ fai,LATIT
	COMMON /pp/ ppwhere
	COMMON /SIG/ AM01,X001,ASQ01,SPS,SIGMA,SIGMA2

!  determine the gyroradius vector
      CALL guid_c(Ef, Bf, v, rg)
      DO 8 i=1,3
        pg(i)=v(i+3)-rg(i)  !旋转中心的位置
8     CONTINUE

      if(k.eq.0) then	   !这里的k,kc,ks代表什么呀?kc代表转一圈的参数,k代表次数
       kc=0
       ks=0
!       CALL hmlt(v(4),v(5),v(6),h_mlt)
       RLm=dsqrt(pg(1)*pg(1)+pg(2)*pg(2)+pg(3)*pg(3))
	
      vv=(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v0*v0
	gam=dsqrt(1-vv/(vc*vc))
      enem=(vc*vc)/1000.d0/qm*(1/gam-1) 
!	 这里的能量是否应用相对论?
!      enem=(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v0*v0/2000.d0/qm
      
	 tm=t
        do 9 i=1,6
          vm(i)=v(i)
          if(i.le.3) vm(i)=vm(i)*v0/1.d3
9      continue
      endif
    r=sqrt(v(4)*v(4)+v(5)*v(5))
	 if(r>2000000.0/Re) then
	   write(*,'(4x,A,f12.6,A)')'Trajectory arrive at y=-5Re (t=',t,&
     &'s)'
         write(20,'(3x,A,f12.6,A)') 'Trajectory arrive at R=500km(t=', &
     &t, 's)'
         goto 2
       endif

!	  if(v(5).gt.4.d0) then
!	  write(*,'(4x,A,f12.6,A)')'Trajectory arrives at y=4 (t=',t,' 
!     !s)'
!         write(20,'(3x,A,f12.6,A)') 'Trajectory arrives at y=4 (t=',t
!     !	   ,'s)'
!         goto 2
!       endif

!  energy in keV
	  vv=(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v0*v0
	  gam=dsqrt(1-vv/(vc*vc))
      enem=(vc*vc)/1000.d0/qm*(1/gam-1) 

       energy=enem 
!      energy=(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))*v0*v0/2000.d0/qm

       RL=dsqrt(pg(1)*pg(1)+pg(2)*pg(2)+pg(3)*pg(3))
       if(RL.LT.RLm) then
	  tm=t
        RLm=RL
        enem=energy
        do 90 i=1,6
         vm(i)=v(i)
         if(i.le.3) then
          vm(i)=vm(i)*v0/1.d3
         endif
90      continue
       endif
!!      else
!  ensures program terminates
!!       write(*,'(4x,A,f12.6,A)') 'Outside magnetosphere (t=',t,' s)'
!!       write(20,'(3x,A,f12.6,A)') 'Outside magnetosphere (t=',t,' s)'
!!       goto 2
!!      endif
    
	if (m.GT.3000) then
	  write(*,15) 'X=',v(4),'   Y=',v(5),'   Z=',v(6)
	  write(*,17) '  t=',t
	  write(*,17) 'Energy=',enem
	end if
	if (dabs(v(6)).LT.1.0d-1) then
	     vv=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
	     Pangel=Acosd(v(3)/vv)
	     write(10,16) v(4),v(5),v(6),Pangel,t,enem
	  end if

!      计算pitch angle以便输出//zhou
	speed(4)=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
	speed(1)=v(1)/speed(4)
	speed(2)=v(2)/speed(4)
	speed(3)=v(3)/speed(4)
 
      magne(4)=dsqrt(Bf(1)*Bf(1)+Bf(2)*Bf(2)+Bf(3)*Bf(3))
	magne(1)=Bf(1)/magne(4)
	magne(2)=Bf(2)/magne(4)
	magne(3)=Bf(3)/magne(4)

    vdist=(speed(1)-magne(1))**2d0+(speed(2)-magne(2))**2d0
	vdist=vdist+(speed(3)-magne(3))**2d0

	pitchangle=(dacos(1-vdist/2))*180d0/3.1415927

	if(abs(mod(m-1,outputonce)).LE.0.001) then
        
!      每outputonce步输出一次 //zhou
   !     if (t<-1.97) then
	  write(20,18) v(1),v(2),v(3),v(4),v(5),v(6),t,&
     & Ef(1),Ef(2),enem,Bf(1),Bf(2)&
     &,Bf(3) ,pitchangle  !,fai,LATIT,SIGMA2,ppwhere
    Fzhongjian=FPSD(v)
      
      R=sqrt(v(4)*v(4)+v(5)*v(5))*Re
      ZRe=v(6)*Re
      bb=interp(R,'B ')
      write(40,19) t,Fzhongjian,Bf(3),bb,Ef(2),Ef(1)
      
      vx=v
      vx(1)=vx(1)+0.05/v0
      vy=v
      vy(2)=vy(2)+0.05/v0
      vz=v
      vz(3)=vz(3)+0.05/v0
      xx=v
      xx(4)=xx(4)+10/Re
      yy=v
      yy(5)=yy(5)+10/Re
      zz=v
      zz(6)=zz(6)+10/Re
	 write(200,20) FPSD(vx)*1d12,FPSD(vy)*1d12,FPSD(vz)*1d12,FPSD(xx)*1d12,FPSD(yy)*1d12,FPSD(zz)*1d12
     
      mu=1.0/2*me*(v(1)**2+v(2)**2)*v0**2/Bf(3)*1.0d9
      write(180,19) t,mu,tiduBperp(R,ZRe,t),tiduBpara(R,ZRe,t),tiduBt(R,ZRe,t)
   !   end if
    end if
15      format(5x,3(A,f10.6))
16	  format(5x,6  (e14.6))
17	  format(5x,A,e18.9)
18	  format(5x,15(e15.7),I2)
19    format(5x,6(e20.12))
20    format(5x,6(e20.12))
      RETURN

2     ks=1
      RETURN

      END

! ------------------------------------
      SUBROUTINE guid_c(Ef, Bf, v, rg)
    use global
!  calculate the gyroradius vector (in Re) at each time step

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	integer qmsign
      DIMENSION Ef(3), Bf(3), v(6), rg(3)
      COMMON /basics/ qm,rad,E0,B0,vc,qmsign
      COMMON /gyrora/ omega0,vd0

!  gyroradius vector in Re (dimensionless)
      B2=Bf(1)*Bf(1)+Bf(2)*Bf(2)+Bf(3)*Bf(3)
      tmp=1.d0/(omega0*B2*Re)
      rg(1)=tmp*((Bf(2)*v(3)-Bf(3)*v(2))*v0-Ef(1)*vd0)
      rg(2)=tmp*((Bf(3)*v(1)-Bf(1)*v(3))*v0-Ef(2)*vd0)
      rg(3)=tmp*((Bf(1)*v(2)-Bf(2)*v(1))*v0-Ef(3)*vd0)

      return
      end

!  -------------------------------------

! ******** calculate GSM comments of initial velocity	*********
      SUBROUTINE init_v(energy,pitch,Bf,phi,vstart)
    use global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION vstart(6), Bf(3),vvector(3)
	integer fb,vol,vorp,qmsign
      COMMON /basics/ qm,rad,E0,B0,vc,qmsign
	COMMON /calcu/ fb,vol,vorp
	COMMON /vvv/ vvector
	 
! determine the velocity components parallel and perpendicular to the 
! magnetic field 
!      v=dsqrt(energy*2000.d0*qm)/v0
      gam=1.d0/((energy*1000.d0*qm/(vc*vc))+1.0d0)
      v=vc*dsqrt(1.d0-gam*gam)/v0

	if(vorp.EQ.1) then
	 vtotal=dsqrt(vvector(1)**2.d0+vvector(2)**2.d0+vvector(3)**2.d0)
	 vstart(1)=v*vvector(1)/vtotal
	 vstart(2)=v*vvector(2)/vtotal
	 vstart(3)=v*vvector(3)/vtotal
	else !以下为原有的程序//zhou

      if(pitch.ne.90.0D0) then 
        v_para=v*dcos(pitch*rad)
      else 
       v_para=0.0D0
      endif
      if(pitch.ne.180.0) then
       v_perp=v*dsin(pitch*rad)
      else 
       v_perp=0.0D0
      endif

      if(phi.ne.90.) then
       vx=v_perp*dcos(phi*rad)
      else
       vx=0.0
      endif
      if(phi.ne.180.) then
   
       vy=v_perp*dsin(phi*rad)
      else
       vy=0.0
      endif

! Convert velocity components to GSM coordinate. v_perp=vz
      B_azimuth=dsqrt(Bf(1)*Bf(1)+Bf(2)*Bf(2))
      Bt=dsqrt(Bf(3)*Bf(3)+B_azimuth*B_azimuth)
      theta=dacos(Bf(3)/Bt)
      if(B_azimuth.ne.0.0) then
!       fai=dasin(Bf(2)/B_azimuth)     !改为下一行，似乎这里有严重错误！！//zhou
        fai1=dacos(Bf(1)/B_azimuth)
	  fai2=dasin(Bf(2)/B_azimuth)
      else
!      fai=phi*rad                     !//zhou 2001.5.18
       fai1=phi*rad
	 fai2=phi*rad
      endif

!     faideg=fai/rad                   !//zhou  
      faideg1=fai1/rad
	faideg2=fai2/rad
!     if(faideg.ne.90.0.and.faideg.ne.270.0) then 
      if(faideg1.ne.90.0.and.faideg2.ne.-90.0) then
!      cosfai=dcos(fai)
       cosfai=dcos(fai1)
      else
       cosfai=0.0
      endif
!     if(faideg.ne.180.0) then 
      if(faideg1.ne.180.0) then
!      sinfai=dsin(fai)
       sinfai=dsin(fai2)
      else
       sinfai=0.0
      endif

      if(B_azimuth.ne.0.0) then
!       vstart(1)=vx*sinfai+(vy*dcos(theta)+v_para*dsin(theta))*cosfai
      vstart(1)=vx*sinfai+(vy*dcos(theta)+v_para*dsin(theta))*cosfai
       vstart(2)=-vx*cosfai+(vy*dcos(theta)+v_para*dsin(theta))*sinfai
!	 vstart(2)=vx*cosfai+(vy*dcos(theta)+v_para*dsin(theta))*sinfai
      else
       vstart(1)=v_perp*cosfai
!	 vstart(1)=-v_perp*cosfai
       vstart(2)=v_perp*sinfai
!	 vstart(2)=-v_perp*sinfai
      endif
      vstart(3)=v_para*dcos(theta)-vy*dsin(theta)		

	endif


      return
      end

! -------------------------------------------------
      SUBROUTINE D_MODI(X,Y,F,DFX,DFY)
! Modify the half thickness determined by T96 model. See Lu, G., et.al, 
! J. Geophys. Res., 104, 12,327, 1999. 

! Also calculate the derivatives of modifting function F(x) with regard 
! to X. 
    use global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /thickne/ AA,Xc,Yc,Dx,Dy

      X1=(X-Xc)/Dx
      Y1=(Y-Yc)/Dy
      fx=dexp(-X1*X1)
      fy=dexp(-Y1*Y1)
      F=1.-AA*fx*fy
      DFX=2.0*X1/Dx*AA*fx*fy
      DFY=2.0*Y1/Dy*AA*fx*fy

      RETURN
      END      
	  
	  FUNCTION interp(R,biaoshi)
      use global
	  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	  INTEGER  weizhi
	  CHARACTER(len=2)  biaoshi
        DOUBLE PRECISION  interp,R

	  
! Rmodel precision equal to 1m	
	    weizhi=floor(R)+1
        if (weizhi==0) then
            interp=0

            return
        endif

! select
      select case (biaoshi)
	  
      case ('E ')
        if (R<499999) then
		interp=Emodel(weizhi)*(Rmodel(weizhi+1)-R) &
		& +Emodel(weizhi+1)*(R-Rmodel(weizhi))
        else
            interp=Emodel(500000)
        endif
        
        
		
      case ('B ')
          if (R<499999) then
	    interp=Bmodel(weizhi)*(Rmodel(weizhi+1)-R) &
		& +Bmodel(weizhi+1)*(R-Rmodel(weizhi))
          else
              interp=Bmodel(500000)
          endif
      
        case ('Y1')
         if (R<499999) then 
		interp=Y1model(weizhi)*(Rmodel(weizhi+1)-R) &
		& +Y1model(weizhi+1)*(R-Rmodel(weizhi))
         else
             interp=Y1model(500000)
2        endif
         
        case ('Y2')
          if (R<499999) then 
		  interp=Y2model(weizhi)*(Rmodel(weizhi+1)-R) &
		  & +Y2model(weizhi+1)*(R-Rmodel(weizhi))
		  else
          interp=Y2model(500000)
          endif
          
        case ('Y3')
              if (R<499999) then 
		    interp=Y3model(weizhi)*(Rmodel(weizhi+1)-R) &
		    & +Y3model(weizhi+1)*(R-Rmodel(weizhi))
              else
                  interp=Y3model(500000)
            endif

		
	  end select
	  end FUNCTION interp
	  
    
    FUNCTION tiduB(R)
        use global
	  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	  INTEGER  weizhi
	  DOUBLE PRECISION  tiduB,R
        weizhi=floor(R)+1
        if (weizhi==0) then
            tiduB=0
        return
        endif
        if (weizhi==1) then
            tiduB=(Bmodel(weizhi+1)-Bmodel(weizhi))
        return
        endif
        if(weizhi<499999) then
        tiduB=(Bmodel(weizhi+1)-Bmodel(weizhi-1))/2
        return
        else 
            tiduB=0
            return
            endif
    end function
       
   
    FUNCTION tiduBperp(R,ZRe,t)
    USE global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Bf(3),A15(15),AA11(11),BB8(8)
      REAL PARMOD(10)
	DOUBLE PRECISION kz,L,interp
      COMMON /dateut/ IYR,IDAY,IHOUR,MIN,ISEC,IOPT,PARMOD
      COMMON /GEOPACK/ A15,PS,AA11,K,IY,BB8
      COMMON /timepa/ tt,tp,tf
	COMMON /B_field/ B1,Bzimf,Vsw

      tt=ttotal-abs(t)
      rr=R/eq
      zz=ZRe/eq2
      Bmo=interp(R,'B ')
      Bnow=Bmo+(2.0d0/eq*rr*rr*exp(-rr*rr)*exp(-zz*zz)*tt*mmm-2.0/eq*exp(-rr*rr)*exp(-zz*zz)*tt*mmm)*1e3 
      rr2=(R+10)/eq
      zz2=ZRe/eq2
      Blater=interp(R+10,'B ')+(2.0d0/eq*rr2*rr2*exp(-rr2*rr2)*exp(-zz2*zz2)*tt*mmm-2.0/eq*exp(-rr2*rr2)*exp(-zz2*zz2)*tt*mmm)*1e3
      tiduBperp=(Blater-Bnow)/10

      return
    end
    
    FUNCTION tiduBpara(R,ZRe,t)
    USE global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Bf(3),A15(15),AA11(11),BB8(8)
      REAL PARMOD(10)
	DOUBLE PRECISION kz,L,interp
      COMMON /dateut/ IYR,IDAY,IHOUR,MIN,ISEC,IOPT,PARMOD
      COMMON /GEOPACK/ A15,PS,AA11,K,IY,BB8
      COMMON /timepa/ tt,tp,tf
	COMMON /B_field/ B1,Bzimf,Vsw

      tt=ttotal-abs(t)
      rr=R/eq
      zz=ZRe/eq2
      Bmo=interp(R,'B ')
      Bnow=Bmo+(2.0d0/eq*rr*rr*exp(-rr*rr)*exp(-zz*zz)*tt*mmm-2.0/eq*exp(-rr*rr)*exp(-zz*zz)*tt*mmm)*1e3 
      rr2=R/eq
      zz2=(ZRe+10)/eq2
      Blater=interp(R+10,'B ')+(2.0d0/eq*rr2*rr2*exp(-rr2*rr2)*exp(-zz2*zz2)*tt*mmm-2.0/eq*exp(-rr2*rr2)*exp(-zz2*zz2)*tt*mmm)*1e3
      tiduBpara=(Blater-Bnow)/10
    
      return
    end
    
    FUNCTION tiduBt(R,ZRe,t)
    USE global
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Bf(3),A15(15),AA11(11),BB8(8)
      REAL PARMOD(10)
	DOUBLE PRECISION kz,L,interp
      COMMON /dateut/ IYR,IDAY,IHOUR,MIN,ISEC,IOPT,PARMOD
      COMMON /GEOPACK/ A15,PS,AA11,K,IY,BB8
      COMMON /timepa/ tt,tp,tf
	COMMON /B_field/ B1,Bzimf,Vsw

      tt=ttotal-abs(t)
      rr=R/eq
      zz=ZRe/eq2
  
      tiduBt=(2.0d0/eq*rr*rr*exp(-rr*rr)*exp(-zz*zz)*mmm-2.0/eq*exp(-rr*rr)*exp(-zz*zz)*mmm)*1e3 
      return
      end
    
    FUNCTION FPSD(v)
      use global
        IMPLICIT DOUBLE PRECISION (A-H, O-Z)
        DOUBLE PRECISION interp,R,vz,vf,vr,y1,y3,E,B,rc_x,rc_y,rc,vE,vD,miu_guiyi1,miu_guiyi2
        DOUBLE PRECISION FPSD1,FPSD2,FPSD
        DIMENSION  v(6)
        if (v(4)==0.d0) then
          if (v(5)==0.d0) then
              R=0.d0
          else 
              R=sqrt(v(4)*v(4)+v(5)*v(5))*Re
          endif
        else
           R=sqrt(v(4)*v(4)+v(5)*v(5))*Re  !单位变成m
        endif
        
        vz=v(3)*v0/v00
        vx=v(1)*v0/v00
        vy=v(2)*v0/v00
        
        E=interp(R,'E ')
        B=interp(R,'B ')
        if (v(1)==0.d0) then
          if (v(2)==0.d0) then
              rc=0.d0
                rc_x=v(4)*Re/D*100
                rc_y=v(5)*Re/D*100
          else 
              rc=me*v00/B*1.d9/qe/D*100*sqrt(vx**2+vy**2)
        rc_x=v(4)*Re/D*100+rc*(-vy)/sqrt(vx**2+vy**2)
        rc_y=v(5)*Re/D*100+rc*vx/sqrt(vx**2+vy**2)
          endif
        else
        rc=me*v00/B*1.d9/qe/D*100*sqrt(vx**2+vy**2)
        rc_x=v(4)*Re/D*100+rc*(-vy)/sqrt(vx**2+vy**2)
        rc_y=v(5)*Re/D*100+rc*vx/sqrt(vx**2+vy**2)
           !单位变成m
        endif
        r_=sqrt(rc_x**2+rc_y**2)
        r_m=r_*D/100
        vE=-interp(r_m,'E ')/interp(r_m,'B ')*1.d9
        vf=(-v(1)*v(5)+v(2)*v(4))/sqrt(v(4)*v(4)+v(5)*v(5))*v0/v00
        tiduB_=tiduB(r_m)
        B_2=interp(r_m,'B ')**2
        v2=(vx**2+vy**2)*v00**2
        vD=1.0/2*me*v2/qe/B_2*1.d9
        vD=vD*tiduB_
        vdrift=vE-vD
        vdrift_x=vdrift*(-rc_y)/r_
        vdrift_y=vdrift*(rc_x)/r_
        miu_guiyi1=((vx*v00-vdrift_x)**2+(vy*v00-vdrift_y)**2)*be*B00
        miu_guiyi1=miu_guiyi1/interp(R,'B ')/v00**2
        miu_guiyi2=((vx*v00-vdrift_x)**2+(vy*v00-vdrift_y)**2)*be_*B00
        miu_guiyi2=miu_guiyi2/interp(R,'B ')/v00**2
        y1=interp(R,'Y1')
        y3=interp(R,'Y3')

        xx=R/D*100  !D单位cm 
        a0=n0*(1-delta)*(1/pi)**1.5*10**6*v00**(-3)
        a1=exp(-vx**2-vy**2+miu_guiyi1+vf*xx*sqrt(2*ae*yibusen*(1+tao)) &
        & -vz**2+ (1+tao)*(y3-1/(w-1)*y1))
        FPSD1=a0*a1
        a3=n0*delta*(lamda/pi)**1.5*10**6*v00**(-3)
        a4=exp(-lamda*(vx**2+vy**2+vz**2-miu_guiyi2) + lamda*(1+tao)*y3)
        FPSD2=a3*a4
        FPSD=FPSD1+FPSD2
        
    END FUNCTION FPSD