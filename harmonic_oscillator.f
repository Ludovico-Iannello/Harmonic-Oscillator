	program harmonic_oscillator

	parameter (measures=2000000)
	real t1,t2
	real var_y2,var_Dy2,var_corr,lattice,y2,Dy2
	real O_mean_1(measures)
	real, dimension(:), allocatable :: field, O_arr
	real, dimension(:, :), allocatable :: corr
	common/various/eta,d_metro
	common/arr/y2(measures),Dy2(measures)

	call second(t1)
	open(1, file='parameters',status="unknown")
	open(2,file='meas_out',status='unknown')
	open(3,file='observables_all_paths',status='unknown')
	

CC====================================
CC LETTURA PARAMETRI DELLA SIMULAZIONE
CC====================================

	iflag=1              !!partenza caldo/freddo/precedente
	i_decorrel=10        !!updating fra una misura e l'altra
	i_term=500000        !!passi di termalizzazione
	beta_omega=10
	
	do i=1,20
		Npassi=i*10            
		eta=beta_omega/Npassi	     !!valore del parametro eta = omega * a
		write(*,*) 'eta', eta
		d_metro=2.0*sqrt(eta)        !!parametro del metropolis
		
CC===========================================	
CC TEST best delta
CC=========================================== 		
C do itd=1,10
C d_metro=0.0001+(2*sqrt(eta)-0.0001)/10*itd
C write(*,*) 'delta metropolis', d_metro
!! scrivo nel file best_delta il valore di delta del metropolis 
!! e il valore della dev calcolata su <y2> media su tutti i cammini.	
C call best_delta()
c enddo	
!!===========================================

CC=====================================
CC OPERAZIONI PRELIMINARI
CC=====================================
 		call ranstart                     !! initialize random number generator

		allocate(field(Npassi))
		allocate(O_arr(int(Npassi/2)))
		allocate(corr(measures,int(Npassi/2)))
		
		CALL initialize_lattice(field,Npassi,iflag)    
	
CC=============================
CC TERMALIZZAZIONE
CC=============================
		do it = 1,i_term
			call update_metropolis(field,Npassi)
		enddo
		
CC=============================================================
CC SESSIONE ALL'EQUILIBRIO CON MISURE
CC=============================================================
		do iter = 1,measures

CC AGGIORNAMENTO CONFIGURAZIONE
			do idec = 1,i_decorrel
				call update_metropolis(field,Npassi)	     
		     	enddo

		     	
	   		call measure(field,Npassi,obs1,obs2)
	   		
	   		call correlatore(field,Npassi,1,O_arr)
	   		
	   		call mean_O_corr(field,Npassi,1,mean_1_path)
	   		
	   		corr(iter,:)=O_arr
	   		
	   		O_mean_1(iter)=mean_1_path
	   		y2(iter)=obs1 !! array per fare le medie su tutti i cammini
	   		Dy2(iter)=obs2
			
		enddo !! measures
		
CC CAlcolo media ed errore su tutti i cammini
		y2_all_path=sum(y2)/float(measures)
		Dy2_all_path=sum(Dy2)/float(measures)
		
		
		open(i*10,status='unknown')
		
		
		do iter = 1, int(Npassi/2)
			media_corr=sum(corr(:,iter))
			call bootstrap(corr(:,iter),var_corr)
			write(i*10,*) iter, media_corr/float(measures), var_corr
		enddo
		
		
		call bootstrap(y2,var_y2)
		call bootstrap(Dy2,var_Dy2)
		
		write(1,*) 'eta,d,N,b_o:',eta,d_metro,Npassi,beta_omega
		write(1,*) '<O(q)>^2',(sum(O_mean_1)/float(measures))**2
	
		write(3,*) y2_all_path,var_y2,Dy2_all_path,var_Dy2
		
		
		!! creo file fort.i*10 dove i*10 rappresenta N passi
		!! riporto il valore della x del cammino

		open(i*11,status='unknown')
		write(i*11,*) field		
		close(i*11)
		deallocate(field)
		deallocate(O_arr)
		deallocate(corr)
	enddo

CC=========TERMINE SIMULAZIONE MONTE-CARLO===========


CC==============================================
CC PRENDO L'ULTIMO STATO DEL GENERATORE RANDOM
CC==============================================
	call ranfinish

CC==============================================
CC SALVO CONFIGURAZIONE E STATO GEN. RANDOM PER POTER RIPARTIRE
CC==============================================
	
	close(3)
	close(2)
	close(1)
	call second(t2)
	write(*,*)t2-t1
CC==============================================
	stop
	end
CC===========F   I   N   E======================




CC      INIZIO DEFINIZIONE SUBROUTINES


c*****************************************************************
	subroutine initialize_lattice(field,Npassi,iflag)
c*****************************************************************
c ASSEGNO LA CONFIGURAZIONE DI PARTENZA DELLA CATENA DI MARKOV
c=================================================================
      	real, dimension(Npassi) :: field
      	common/various/eta,d_metro

  	!!PARTENZA A FREDDO ...
      	if (iflag.eq.0) then
        	do i = 1,Npassi                  !! loop su tutti i siti
                	field(i) = 0.0
         	enddo
	!!A CALDO ...
      	elseif (iflag.eq.1) then
        	do i = 1,Npassi                  !! loop su tutti i siti
               		x = 1.0 - 2.*ran2()       !! frand() random fra -1 e 1
                  	field(i) = x
         	enddo
	!!O DA DOVE ERO RIMASTO L'ULTIMA VOLTA
      	else
        	open(9,file='lattice',status='old')
        	read(9,*) field

         	close(9)
      	endif

      	return
      	end
c=================================================================


c*****************************************************************
	subroutine update_metropolis(field,Npassi)
c*****************************************************************

      	real, dimension(Npassi) :: field
      	common/various/eta,d_metro

      	c1 = 1./eta
      	c2 = (1./eta + eta/2.)

      	do i = 1,Npassi         !! loop su tutti i siti, qui il sito
      	                        !! non e` scelto a caso ma faccio una spazzata
                                !! iterativa su tutti i siti, si puo` dimostrare
                		!! che va bene lo stesso per il bilancio, ma meno banale da provare
        	npp=i+1
         	nmm=i-1
         	if (i.eq.1) nmm=Npassi
         	if (i.eq.Npassi) npp=1
         	ip = npp            !! calcolo le coordinate
         	im = nmm            !! dei due primi vicini

         	force = field(ip) + field(im) !! costruisco la forza

         	phi =  field(i)        !! phi = valore attuale del campo.
         	phi_prova = phi + 2.*d_metro*(0.5-ran2())

         	p_rat = c1 * phi_prova * force - c2 * phi_prova**2
         	p_rat = p_rat - c1 * phi * force + c2 * phi**2


         	x = log(ran2())                      !! METRO-TEST! x = random (0,1)
                                              !! x < p_rat verifica anche caso
         	if (x.lt.p_rat) field(i) = phi_prova !! p_rat > 1 -> se si accetto

      	enddo                     !!  chiudo il loop

      	return
      	end
c=================================================================


c*****************************************************************
	subroutine measure(field,Npassi,obs1,obs2)
c*****************************************************************

      	real, dimension(Npassi) :: field
      	obs1 = 0.0
      	obs2 = 0.0
      	do i = 1,Npassi
		npp=i+1
		if (i.eq.Npassi) npp=1
		obs1 = obs1 + field(i)**2
		obs2 = obs2 + (field(i)-field(npp))**2
      	enddo

      	obs1 = obs1/float(Npassi) !! media sul singolo path di y^2
      	obs2 = obs2/float(Npassi) !! media sul singolo path di Delta y^2

      	!write(2,*) obs1,obs2
	return
      	end

c*****************************************************************
	subroutine correlatore(field,Npassi,power,O_arr)
c*****************************************************************
      	real, dimension(Npassi) :: field
      	real, dimension(int(Npassi/2)) :: O_arr
      	
      	integer tau_max,power
      	
      	tau_max = int(Npassi/2) !fino a dove calcolare la correlazione
    	!calcolo della correlazione
    	
	do k = 1, tau_max
		yc = 0.0
		do i = 1, Npassi-k
			yc = yc + (field(i)*field(i+k))**power
		enddo

		O_arr(k) = yc/float(Npassi-k) !!correlazione su singolo path
	enddo
	return
      	end      	
      
c*****************************************************************
	subroutine mean_O_corr(field,Npassi,power,mean_1path)
c*****************************************************************
      	integer power
      	real, dimension(Npassi) :: field
      	obs1 = 0.0
      	do i = 1,Npassi
		npp=i+1
		if (i.eq.Npassi) npp=1
		obs1 = obs1 + field(i)**power
      	enddo

      	mean_1path = obs1/float(Npassi) !! media sul singolo path di y^power

	return
      	end

C=================================================================
	subroutine bootstrap(x,dx)
        integer j,q,i,N,l,k,p
        parameter (k=5000,N=2000000,p=100)
        real dx,media_x
        real sim_arr(k)      		    !!media  simulazioni
        real a(N),x(N)        !!a array provvisorio

        dx=0.0
        do j=1,k
          	do i=1,N,p
            		l=int((N-p)*rand()+1)
            		do q=1,p !punti successivi della catena che taglio
              			a(i+q-1)=x(l+q)
            		end do
          	end do
          	sim_arr(j)=sum(a)/float(N)
        end do

        media_x=sum(sim_arr)/float(k)

        do i = 1,k
          	dx=(sim_arr(i)-media_x)**2+dx
        end do

        dx=sqrt(dx/((k-1)))
        return
      	end    
      
c=================================================================
      	subroutine best_delta()
        integer i,j,l,k,q,p,d
        parameter (measures=2000000,k=1000,p=100)
        real x,m,dm,u,dev,y2,Dy2
        real a(k)
        common/arr/y2(measures),Dy2(measures)
        common/various/eta,d_metro
        x=0
        u=0
        

        m=0
        dm=0
        do j=1,k
       		y=0
          	do i=1,measures,p
            		l=int((measures-p)*rand()+1)
            		do q=1,p !punti successivi della catena che taglio
              			y=y2(l+q)+y
            		end do
          	end do
         	a(j)=y/measures
         	m=a(j)+m
        end do
        
        m=m/float(k)
        do i = 1,k
        	dm=(a(i)-m)**2+dm
        end do
        
        dev=sqrt(dm/((k-1)))
        write(1,*) d_metro, dev


       return
       end


c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================
