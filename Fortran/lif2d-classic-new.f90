program lif2d
	implicit real*8 (a-h,o-z)
	real*4 U1,V1
	integer*4 seed
	character (len=40) :: filename

!-------QUESTIONS-------!
!Connectivity Pattern
!(-) in sum calculation

!-------FIXES-------!
!dt and alpha common coefficient in sum
!calculate sum only for neighbours


!-------Parameter Declarations-------!
!uth: Potential threshold
!b: integrator floor (mi in paper)
!alphai: sigma in paper
	parameter (b=1.0,pi=3.14159)
	parameter (dt=0.001,itstep=1.0/dt)
	parameter (uth=0.980)
	parameter (N=100,ITIME=itstep*N*N)		!??????? WHY N*N?
	parameter (Nclassic=30)      !size of classic connectivity    
	parameter (alphai=0.7)
	parameter (nsteps=0,ITRANS=100*itstep)

!-------Array Declarations-------!
!u: Potential at time t
!unext: Potential at time t+dt
!f: Sum of differences from neighbors' potential at time t
!fnext: Sum of differences from neighbors' potential at time t+dt
!neigh: Adjacency matrix
!nk: Number of resets to resting potential	

!mk: ONLY INITIALIZED. NOT USED.
!itc: ?????????????????????????
	dimension unext(N,N),u(N,N),f(N,N),fnext(N,N),itc(N,N)
	dimension mk(N,N),nk(N,N),neigh(N,N)

!-------Files-------!
	open(unit=99,file='profile0.7-classic2D.dat3R30-minusp00')
	open(unit=98,file='wmegaa0.7-classic2D.dat3R30-minusp00')
	open(unit=97,file='liandf0.7-one-classic2DR30.dat3-minusp00')
	open(unit=96,file='liandf0.7-spt-classic2DR30.dat3-minusp00')

!-------Initialization-------!
!	seed=89732753      !seed1
!	seed=239354385     !seed2
	seed=765385665     !seed3
!	seed=930574656     !seed4
!	seed=573028468     !seed5
	ro=rand(seed)
	do i=1,N
		do j=1,N
			u(i,j)=rand()
		enddo
	enddo
	t=0.0
	do i=1,N
		do j=1,N
			mk(i,j)=0
			nk(i,j)=0
			itc(i,j)=0
		enddo
	enddo
	kk=0.
	K=0.
	nnn=0

!-------Random & Classic Connectivity Chimera-------!
	nnb=0 !2R+1
	do i=1,N
		do j=1,N
			neigh(i,j)=0
		enddo
	enddo
	do i=1,Nclassic
		do j=1,Nclassic     
			neigh(i,j)=1
			neigh(N+1-i,j)=1
			neigh(i,N+1-j)=1
			neigh(N+1-i,N+1-j)=1
			nnb=nnb+neigh(i,j)+neigh(N+1-i,j)
			nnb=nnb+neigh(i,N+1-j)+neigh(N+1-i,N+1-j)
		enddo
	enddo                                 

	do i=1,N
		do j=1,n
			nnn=nnn+neigh(i,j)
		enddo
	enddo
	if(nnn.ne.nnb) write(*,*) "PROBLEM!!!!!!",nnn,nnb
	alpha=alphai/dble(nnb*nnb-1) 					!Potential sum coefficient
	tnext=0.0
	write(*,*) 'classic-chimera with nnb=',nnb		!classicchimera

!-------Time Integration-------!
	do it=1, ITIME             !do on time
		Tc=0.
		Ti=it

		!-------DEBUGGING???-------!		
		if(mod(it,50000).eq.0) write(*,*) it,unext(2,3)
		!------------------------------!		

		do i=1,N                   !do on oscillators for refractory period
			do j=1,N

				!-------Refractory???-------!
				if(unext(i,j).eq.0.0.and.itc(i,j).lt.nsteps) then
					itc(i,j)=itc(i,j)+1
					exit !go to 33
				else
					itc(i,j)=0
				endif                     !end of refractory period
				
		                                  !oscillators

				!-------Sum calculation-------!
				f(i,j)=0.0				
				do kki=1,N-1
					kkpi=i+kki
					if(kkpi.gt.N) kkpi=kkpi-N
					do kkj=1,N-1
						kkpj=j+kkj
						if(kkpj.gt.N) kkpj=kkpj-N
						fnext(i,j)=f(i,j)-dt*alpha*neigh(kki,kkj)*(u(kkpi,kkpj)-u(i,j)) !??????? WHY (-)?
						f(i,j)=fnext(i,j)
					enddo           
				enddo
				
				!-------Euler Method Step-------!
				unext(i,j)=u(i,j)+dt*(b-u(i,j))+f(i,j)

				!-------Potential threshold is crossed-------!
				if(unext(i,j).gt.uth) then
					unext(i,j)=0.0
					nk(i,j)=nk(i,j)+1
				endif

				u(i,j)=unext(i,j)
		      
			enddo                         !oscillators
		enddo                         !oscillators

		tnext=tnext+dt
           
		!-------DEBUGGING???-------!		
		if(mod(it,itstep).eq.0) then
			do il=1,N
				write(96,*) il,tnext,unext(il,8)
			enddo
			write(96,*) 
		endif
		if(mod(it,itstep).eq.0) then
			write(97,*) tnext,unext(3,3),unext(1,5),unext(2,8),unext(3,1),unext(4,4)
		endif
        !------------------------------!

		!-------Write potentials to file-------!   
		if(mod(it,(ITRANS)).eq.0) then
			iit=dble(it)/dble(ITRANS)
			write(filename, '("lif-spt-class-dat3R30-",I0,".001")')  iit 
			open (file=filename,unit=90)
			do il=1,N
				do jl=1,N
					write(90,*) il,jl,u(il,jl)
				enddo
			enddo
			close(90)
		endif

		!-------Write mean phase velocities to file-------!
		if(mod(it,(ITRANS)).eq.0) then
			iit=dble(it)/dble(ITRANS)
			write(filename, '("wmegaa0.7-class-dat3R30-",I0,".001")')  iit 
			open (file=filename,unit=90)          
			do il=1,N
				do jl=1,N
					ww=2.0*pi*dble(nk(il,jl))/(dt*dble(it-ITRANS))
					write(90,*) il,jl,ww
				enddo
			enddo
			close(90)
		endif
        
	enddo                !time

	stop

end program lif2d
