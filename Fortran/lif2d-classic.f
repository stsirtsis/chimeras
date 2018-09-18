           program lif2d
           implicit real*8 (a-h,o-z)
           real*4 U1,V1
           integer*4 seed
           character (len=40) :: filename
           parameter (b=1.0,pi=3.14159)
           parameter (dt=0.001,itstep=1.0/dt)
           parameter (uth=0.980)
           parameter (N=100,ITIME=itstep*N*N)
           parameter (Nclassic=30)      !size of classic connectivity    
           parameter (alphai=0.7)
c           parameter (alphai=0.0)
           parameter (nsteps=0,ITRANS=100*itstep)

           dimension unext(N,N),u(N,N),f(N,N),fnext(N,N),itc(N,N)
c           dimension phi(N)
c           dimension du(N),ddu(N),neigh(N)
           dimension mk(N,N),nk(N,N),neigh(N,N)
           
c           open(unit=11,file='cant110011-det')
c           open(unit=12,file='profile-pattern.input')
         open(unit=99,file='profile0.7-classic2D.dat3R30-minusp00')
         open(unit=98,file='wmegaa0.7-classic2D.dat3R30-minusp00')
         open(unit=97,file='liandf0.7-one-classic2DR30.dat3-minusp00')
         open(unit=96,file='liandf0.7-spt-classic2DR30.dat3-minusp00')
c----------INITIALISE--------------------------------------------------
c          seed=89732753      !seed1
c           seed=239354385     !seed2
          seed=765385665     !seed3
c          seed=930574656     !seed4
c          seed=573028468     !seed5
           ro=rand(seed)
           do i=1,N
           do j=1,N
           u(i,j)=rand()
c!!           u(i)=0.0
c!!           read(12,*) ki,u(ki)
c!!           if(ki.ne.i) write(*,*) 'READING PROBLEMS'
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
ccccccccccFracatal Geometry ***************************
c           nnn=0
c           do i=1,N
c           read(11,*) j,neigh(j)
c           nnn=nnn+neigh(j)
c           enddo
c           alpha=alphai/dble(nnn)
c           write(*,*) 'nnn=',nnn
c           tnext=0.0
cccccccccccccccc-Random & Classic Connectivity Chimera-cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
                    nnb=0
                   do i=1,N
                   do j=1,N
                    neigh(i,j)=0
                   enddo
                   enddo
                   do i=1,Nclassic
                   do j=1,Nclassic     
                   neigh(i,j)=1                          !cclassicchimera
                   neigh(N+1-i,j)=1                      !classicchimera
                   neigh(i,N+1-j)=1                      !classicchimera
                   neigh(N+1-i,N+1-j)=1                  !cclassicchimera
                   nnb=nnb+neigh(i,j)+neigh(N+1-i,j)         !cclassicchimera
                   nnb=nnb+neigh(i,N+1-j)+neigh(N+1-i,N+1-j)
                   enddo
                   enddo                                 !cclassicchimera
c
                  do i=1,N
                  do j=1,n
                  nnn=nnn+neigh(i,j)
                  enddo
                  enddo
                  if(nnn.ne.nnb) write(*,*) "PROBLEM!!!!!!",nnn,nnb
                   alpha=alphai/dble(nnb*nnb-1)
                   tnext=0.0
          write(*,*) 'classic-chimera with nnb=',nnb     !cclassicchimera
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----------TIME INTEGRATION--------------------------------------------------

           do it=1, ITIME             !do on time
           Tc=0.
           Ti=it
           if(mod(it,50000).eq.0) write(*,*) it,unext(2,3)
           do i=1,N                   !do on oscillators for refractory period
           do j=1,N
           if(unext(i,j).eq.0.0.and.itc(i,j).lt.nsteps) then
            itc(i,j)=itc(i,j)+1
            go to 33
            else
            itc(i,j)=0
            endif                     !end of refractory period
           f(i,j)=0.0
                                      !oscillators
           do kki=1,N-1
            kkpi=i+kki
            if(kkpi.gt.N) kkpi=kkpi-N
           do kkj=1,N-1
            kkpj=j+kkj
            if(kkpj.gt.N) kkpj=kkpj-N
         fnext(i,j)=f(i,j)-dt*alpha*neigh(kki,kkj)*(u(kkpi,kkpj)-u(i,j))
            f(i,j)=fnext(i,j)
           enddo           
           enddo 
           unext(i,j)=u(i,j)+dt*(b-u(i,j))+f(i,j)

           if(unext(i,j).gt.uth) then
              unext(i,j)=0.0
              nk(i,j)=nk(i,j)+1
              endif

cc           ddu(i)=du(i)
cc           du(i)=(unext(i)-u(i))/dt
cc               if(u(i).eq.0.0) then
cc                   phi(i)=pi/2.0
cc               else
cc                   phi(i)=datan(unext(i)/u(i))
cc               endif
           u(i,j)=unext(i,j)
          
33           enddo                         !oscillators
             enddo                         !oscillators
cc            do k=1,N
cc            ddsign=du(k)*ddu(k)
cc           if (it.gt.ITRANS.and.ddsign.lt.0.and.du(k).gt.0) then
cc           mk(k)=mk(k)+1
cc           endif
cc           enddo

           tnext=tnext+dt
           
c           if(it.gt.150000.and.mod(it,1000).eq.0) then
           if(mod(it,itstep).eq.0) then
           do il=1,N
           write(96,*) il,tnext,unext(il,8)
           enddo
           write(96,*) 
           endif
           if(mod(it,itstep).eq.0) then
          write(97,*) tnext,unext(3,3),unext(1,5),
     &             unext(2,8),unext(3,1),unext(4,4)
            endif
           
           if(mod(it,(ITRANS)).eq.0) then
           iit=dble(it)/dble(ITRANS)
           write(filename, '("lif-spt-class-dat3R30-",I0,".001")')  iit 
c            write(filename,fmt=('a,I0,a')) 'lif-spt-class-dat3',it,'p001'
           open (file=filename,unit=90)
           do il=1,N
           do jl=1,N
            write(90,*) il,jl,u(il,jl)
c            write(*,*)  l,u(l)
           enddo
           enddo
cc           write(90,*)
           close(90)
           endif
c---------correlation functions-----------------------------
c---------A: Spatial correlations---------------------------
c          if(mod(it,100000).eq.0) then
c           cc=0.0
c           do i=1,N-1
c           cc=cc+unext(i)*unext(i+1)
c           enddo
c           cc=cc+unext(N)*unext(1)
c           cc=cc/dble(N)
c           write(95,*) tnext,cc
c           write(*,*) tnext,cc
c          endif
c-------END of Spatial Correlations--------------------------
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

           end
