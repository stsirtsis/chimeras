program lif2d
	real*8 ro, t, sumCoeff
	integer*4 seed, R, totalIter, secIter, writeIter, currMPVIter
	character (len=80) :: filename

	!-------FIXES-------!
	!opening and closing files at writing
	!Euler step

	!-------Parameter Declarations-------!
	!uth: Potential threshold
	!mi: Integrator floor
	!sigma: Coupling strength
	!secIter: Iterations per sec
	!N: Grid dimension
	!totalIter: Total iterations
	!writeIter: Iterations per writing
	!maxRefracIter: Steps at refractory period
	!R: Square radius

	parameter (mi=1.0,pi=3.14159)
	parameter (dt=0.01,secIter=1.0/dt) !dt=0.001
	parameter (uth=0.98)
	parameter (N=100,totalIter=4000000)	!int(secIter*N*N)?
	parameter (R=22)      						!size of classic connectivity
	parameter (sigma=0.1)
	parameter (maxRefracIter=0,writeIter=100*secIter)
	parameter (maxMPVIter=20000, minMPVIter=50000)
	!-------Array Declarations-------!
	!u: Potential at time t
	!unext: Potential at time t+dt
	!currThresCrossings: Number of resets to resting potential
	!currRefracIter: Current iterations already in refractory period

	dimension unext(N,N),u(N,N),currRefracIter(N,N)
	dimension currThresCrossings(N,N)


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
			do while (u(i,j).ge.uth)
				u(i,j)=rand()
			enddo
			u(i,j)=rand()
		enddo
	enddo

	do i=1,N
		do j=1,N
			currThresCrossings(i,j)=0
			currRefracIter(i,j)=0
		enddo
	enddo

	currMPVIter=0
	sumCoeff=sigma/dble((2*R+1)*(2*R+1)-1) 					!Potential sum coefficient
	t=0.0

	!-------Time Integration-------!
	do it=1, totalIter             !do on time

		if (mod(it,1000).eq.0) then
			write(*,*) 'Iteration ',it,'of ',totalIter
		endif

		do i=1,N                   !do on oscillators for refractory period
			do j=1,N

				!-------Refractory-------!
				if(u(i,j).eq.0.0.and.currRefracIter(i,j).lt.maxRefracIter) then
					currRefracIter(i,j)=currRefracIter(i,j)+1
					cycle
				else
					currRefracIter(i,j)=0
				endif                     !end of refractory period

				!oscillators

				!-------Sum calculation-------!
				sumVar=0.0
				iLeftCorner=N+i-R
				jLeftCorner=N+j-R
				iLeftCorner=mod(N+i-R,N)
				jLeftCorner=mod(N+j-R,N)

				do k=iLeftCorner,iLeftCorner+2*R
					do l=jLeftCorner, jLeftCorner+2*R
						sumVar=sumVar+u(i,j)-u(mod(k,N),mod(l,N))
					enddo
				enddo

				!-------Euler Method Step-------!
				unext(i,j)=u(i,j)+dt*(mi-u(i,j)+sumCoeff*sumVar)

				!-------Potential threshold is crossed-------!
				if(unext(i,j).gt.uth) then
					unext(i,j)=0.0
					currThresCrossings(i,j)=currThresCrossings(i,j)+1
				endif

				u(i,j)=unext(i,j)

			enddo                         !oscillators
		enddo                         !oscillators


		if (it.ge.minMPVIter) then
			if (currMPVIter.eq.maxMPVIter) then
				write(filename, '("Results_MPV_LIF_2D_Classic_sigma_",f2.1,"_R_",I2,"_time_",f10.5,"_.dat")')  sigma, R, t
				open(file=filename, unit=90)
				do i=1,N
					do j=1,N
						ww=2.0*pi*dble(currThresCrossings(i,j))/(dt*(it-minMPVIter))
						write(90,fmt="(f10.8,a)",advance="no") ww,','
					enddo
					write(90,*)
				enddo
				close(90)
				currMPVIter=0
			else
				currMPVIter=currMPVIter+1
			endif
		endif

		if(mod(it,1000).eq.0) then
			write(filename, '("Results_POT_LIF_2D_Classic_sigma_",F8.8,"_R_",I2,"_time_",F8.8,"_.dat")')  sigma, R, t
			open(file=filename, unit=90)
			do i=1,N
				do j=1,N
					write(90,fmt="(f10.8,a)",advance="no") u(i,j),','
				enddo
				write(90,*)
			enddo
			close(90)
		endif

		t=t+dt
	enddo

	stop

end program lif2d
