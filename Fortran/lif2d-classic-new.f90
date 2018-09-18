program lif2d
	real*8 ro, t, tnext, sumCoeff
	integer*4 seed, R, totalIter, secIter, writeIter, debugCounter, testCounter
	character (len=40) :: filename

!-------FIXES-------!
!opening and closing files at writing

!-------Parameter Declarations-------!
!uth: Potential threshold
!mi: Integrator floor (mi in paper)
!sigma: Coupling strength (sigma in paper)
!secIter: Iterations per sec
!N: Grid dimension
!totalIter: Total iterations
!writeIter: Iterations per writing
!refracIter: Steps at refractory period
!R: R

	parameter (mi=1.0,pi=3.14159)
	parameter (dt=0.01,secIter=1.0/dt)
	parameter (uth=0.980)
	parameter (N=100,totalIter=50000)	!int(secIter*N*N)?
	parameter (R=10)      						!size of classic connectivity    
	parameter (sigma=0.1)
	parameter (refracIter=0,writeIter=100*secIter)

!-------Array Declarations-------!
!u: Potential at time t
!unext: Potential at time t+dt
!f: Sum of differences from neighbors' potential at time t
!fnext: Sum of differences from neighbors' potential at time t+dt
!neigh: Adjacency matrix
!nk: Number of resets to resting potential	
!itc: Current iterations already in refractory period

	dimension unext(N,N),u(N,N),f(N,N),fnext(N,N),itc(N,N)
	dimension nk(N,N),neigh(N,N)

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
			nk(i,j)=0
			itc(i,j)=0
		enddo
	enddo

	sumCoeff=sigma/dble((2*R+1)*(2*R+1)-1) 					!Potential sum coefficient
	tnext=0.0

!-------Time Integration-------!
	do it=1, totalIter             !do on time
		
		write(*,*) 'Iteration ',it,'of ',totalIter		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!-------DEBUGGING???-------!		
		if(mod(it,50000).eq.0) write(*,*) it,unext(2,3)
		!------------------------------!		
		
		testCounter=0
		do i=1,N                   !do on oscillators for refractory period
			do j=1,N
				
				!-------Refractory-------!
				if(unext(i,j).eq.0.0.and.itc(i,j).lt.refracIter) then
					itc(i,j)=itc(i,j)+1
					exit !go to 33
				else
					itc(i,j)=0
				endif                     !end of refractory period
				
		                                  !oscillators

				!-------Sum calculation-------!
				sumVar=0.0	
				iLeftCorner=mod(N+i-R,N)
				jLeftCorner=mod(N+j-R,N)

				!debugCounter=0
				do k=iLeftCorner,iLeftCorner+2*R
					do l=jLeftCorner, jLeftCorner+2*R
						debugCounter=debugCounter+1
						sumVar=sumVar+u(mod(k,N),mod(l,N))-u(i,j)
					enddo
				enddo
				
				!write(*,*) debugCounter
				!read(*,*)				

				!-------Euler Method Step-------!
				unext(i,j)=u(i,j)+dt*(mi-u(i,j)+sumCoeff*sumVar)

				!-------Potential threshold is crossed-------!
				if(unext(i,j).gt.uth) then
					unext(i,j)=0.0
					nk(i,j)=nk(i,j)+1
				endif
				
				u(i,j)=unext(i,j)
		      	
			enddo                         !oscillators
		enddo                         !oscillators
		tnext=tnext+dt
        
		!-------DEBUGGING-------!
		if (mod(it,secIter).eq.0) then		
			open (file='testing.dat',unit=90) !,position="append"
			do i=1,N
				do j=1,N
					ww=2.0*pi*dble(nk(i,j))/(dt*dble(it))
					write(90,fmt="(f10.8,a)",advance="no") ww,','
				enddo
				write(90,*)
			enddo						
			close(90)
		endif
   
		!-------RESULTS???-------!		
		if(mod(it,secIter).eq.0) then
			do il=1,N
				write(96,*) il,tnext,unext(il,8)
			enddo
			write(96,*) 
		endif
		if(mod(it,secIter).eq.0) then
			write(97,*) tnext,unext(3,3),unext(1,5),unext(2,8),unext(3,1),unext(4,4)
		endif
        !------------------------------!

		!-------Write potentials to file-------!   
		if(mod(it,(writeIter)).eq.0) then
			iit=dble(it)/dble(writeIter)
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
		if(mod(it,(writeIter)).eq.0) then
			iit=dble(it)/dble(writeIter)
			write(filename, '("wmegaa0.7-class-dat3R30-",I0,".001")')  iit 
			open (file=filename,unit=90)          
			do il=1,N
				do jl=1,N
					ww=2.0*pi*dble(nk(il,jl))/(dt*dble(it-writeIter))
					write(90,*) il,jl,ww
				enddo
			enddo
			close(90)
		endif
        
	enddo                !time

	stop

end program lif2d
