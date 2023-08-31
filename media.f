      program media

      parameter(N=10000000)
      integer i,ios,k,j
      real*8 dm,m,tau,sum,C,m1,dm1
      real*8 x(N),Y(N)
      real t1,t2
      call second(t1)

      open(unit=1, file='meas_out', status="unknown")

      m=0
      m1=0
      dm=0
      dm1=0
      tau=0
      sum=0

      do i=1,N
        read (1,*) x(i),y(i)
        m=x(i)+m
        m1=y(i)+m1
      enddo



      m=m/float(N)
      m1=m1/float(N)



      do i=1,N
        dm=(x(i)-M)**2+dm
        dm1=(y(i)-m1)**2+dm1
      enddo

      open(unit=2, file='correlazione.txt', status="unknown")

      do k=1,5000
        sum=0
        C=0
        do i=1,N-k
          sum=(x(i)-M)*(x(i+k)-M)+sum
        enddo
        C=sum/float(N-k)
        write (unit=2, fmt=*) k, C
        tau=tau+C
      enddo


      close(unit=1, status="keep")
      write(*,*) 'tau ', tau
      write(*,*) 'y^2 medio ', m
c      write(*,*) 'il valor di tau è ', tau
      write(*,*) 'err correlato: ', sqrt(dm/(N-1)/N)*sqrt(1+2*tau)
      write(*,*) 'Delta y^2 ', m1
      write(*,*) 'err  ', sqrt(dm1/(N-1)/N)
      call second(t2)
      write(*,*) 'il tempo impiegato è ', t2-t1
      end
