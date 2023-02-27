C================================================
C     LBM2D - from python code by Jonas Latt
C================================================
      MODULE allocated
        INTEGER,  ALLOCATABLE :: obs(:,:)
        REAL,     ALLOCATABLE :: ang(:,:)
        REAL,     ALLOCATABLE :: vin(:,:)
        REAL,     ALLOCATABLE :: vel(:,:,:)
        REAL,     ALLOCATABLE :: fin(:,:,:)
        REAL,     ALLOCATABLE :: fout(:,:,:)
        REAL,     ALLOCATABLE :: feq(:,:,:)
        REAL,     ALLOCATABLE :: rho(:,:)
      END MODULE allocated
C===========================================

      program main
      use allocated
      implicit none

      integer, parameter :: ndim=2, nvec=9
      integer i,j,k,ni,nj,writeits,writevtk,maxits
      integer vec(nvec,ndim),col1(3),col2(3),col3(3)
      integer its,nexti,nextj,iseq

      real ls,re,uLB,nuLB,omega
      real sum2,sum3,cu,usqr
      real cx,cy,cr,xx,yy
      real wt(nvec)
      real vmag_sum,vmag

      CHARACTER(30) outfile
      CHARACTER(6)  string
      
C=====constants
      write(6,*) "constants"
      maxits   = 10000 ! max iterations
      writevtk = 100
      writeits = 100
      iseq = 1
C-----
      Re = 20.0         ! reynolds number
      ni = 420          ! lattice nodes
      nj = 200          ! lattice nodes
      cx = real(ni)/4.0 ! cyl x-coord
      cy = real(nj)/2.0 ! cyl y-coord
      cr = real(nj)/9.0 ! cyl radius
      ls = 20           ! length scale (prev cr = cylinder radius)
      uLB   = 0.03                   ! velocity in lattice units
      nuLB  = uLB*ls/Re              ! viscocity in lattice units
      omega = 1.0 / (3.0*nuLB + 0.5) ! relaxation parameter

C=====allocate
      ALLOCATE( obs(ni,nj) )
      ALLOCATE( ang(ni,nj) )
      ALLOCATE( vin(nj,ndim) )
      ALLOCATE( vel(ni,nj,ndim) )
      ALLOCATE( fin(ni,nj,nvec) )
      ALLOCATE( fout(ni,nj,nvec) )
      ALLOCATE( feq(ni,nj,nvec) )
      ALLOCATE( rho(ni,nj) )

C=====initialise
      vec(1,:) = (/ 1 ,  1 /)
      vec(2,:) = (/ 1 ,  0 /)
      vec(3,:) = (/ 1 , -1 /)
      vec(4,:) = (/ 0 ,  1 /) 
      vec(5,:) = (/ 0 ,  0 /)
      vec(6,:) = (/ 0 , -1 /) 
      vec(7,:) = (/-1 ,  1 /) 
      vec(8,:) = (/-1 ,  0 /)  
      vec(9,:) = (/-1 , -1 /)

      wt(1) = (1.0 / 36.0)
      wt(2) = (1.0 / 9.0)
      wt(3) = (1.0 / 36.0)
      wt(4) = (1.0 / 9.0)
      wt(5) = (4.0 / 9.0)
      wt(6) = (1.0 / 9.0)
      wt(7) = (1.0 / 36.0)
      wt(8) = (1.0 / 9.0)
      wt(9) = (1.0 / 36.0)

      col1 = (/ 1, 2, 3 /)
      col2 = (/ 4, 5, 6 /)
      col3 = (/ 7, 8, 9 /)

      obs = 0    ! int array
      ang = 0.0  ! real array
      vin = 0.0  ! real array
      fin = 0.0  ! real array
      fout= 0.0  ! real array
      feq = 0.0  ! real array
      rho = 0.0  ! real array
      vel = 0.0  ! real array

C=====obstacle

C------original cyclinder
       do i = 1,ni
       do j = 1,nj
         xx = real(i-1) + 0.5
         yy = real(j-1) + 0.5
         if (((xx-cx)**2 + (yy-cy)**2).LT.(cr**2)) obs(i,j) = 1
       enddo
       enddo

C=====inlet profile (vin)
      write(6,*)"inlet profile"
      do j = 1,nj
       vin(j,1) = 1.0 * uLB
       vin(j,2) = 1.0 * uLB
      enddo

C=====init velocity field
      write(6,*)"initial velocity"
      do i = 1,ni
      do j = 1,nj
       vel(i,j,1) = vin(j,1)
       vel(i,j,2) = vin(j,2)
      enddo
      enddo

C=====equilibrium distribution function
C-----fin = equilibrium(1,u)
      do i = 1,ni
      do j = 1,nj
        usqr = 1.5*(vel(i,j,1)**2 + vel(i,j,2)**2)
        do k = 1,nvec
          cu = 3.0*(vec(k,1)*vel(i,j,1) + vec(k,2)*vel(i,j,2))
          fin(i,j,k) = 1.0*wt(k)*(1.0 + cu + 0.5*cu**2 - usqr)
        enddo
      enddo
      enddo

      write(6,*)"start iterations"
C=====iteration loop
      do its = 1,maxits

C-------right wall outflow condition
        do j = 1,nj
        do k = 1,3
          fin(ni,j,col3(k)) = fin(ni-1,j,col3(k))
        enddo
        enddo

C-------compute macroscopic variables rho and u
        rho = 0.0
        vel = 0.0
        do i = 1,ni
        do j = 1,nj
          do k = 1,9
            rho(i,j) = rho(i,j) + fin(i,j,k)
            vel(i,j,1) = vel(i,j,1) + vec(k,1) * fin(i,j,k)
            vel(i,j,2) = vel(i,j,2) + vec(k,2) * fin(i,j,k)
          enddo
          vel(i,j,1) = vel(i,j,1) / rho(i,j)
          vel(i,j,2) = vel(i,j,2) / rho(i,j)
        enddo
        enddo

C-------left wall inflow condition
        do j = 1,nj
          vel(1,j,1) = vin(j,1)
          vel(1,j,2) = vin(j,2)
          sum2 = 0.0
          sum3 = 0.0
          do k=1,3
            sum2 = sum2 + fin(1,j,col2(k))
            sum3 = sum3 + fin(1,j,col3(k))
            rho(1,j) = (sum2 + 2.0*sum3) / (1.0-vel(1,j,1))
          enddo
        enddo

C-------compute equilibrium (rho, vel)
        do i = 1,ni
        do j = 1,nj
          usqr = 1.5*(vel(i,j,1)**2 + vel(i,j,2)**2)
          do k = 1,nvec
            cu = 3.0*(vec(k,1)*vel(i,j,1) + vec(k,2)*vel(i,j,2))
            feq(i,j,k) = rho(i,j)*wt(k)*(1.0 + cu + 0.5*cu**2 - usqr)
          enddo
        enddo
        enddo

C-------calculate populations (at inlet)       
        do j = 1,nj
        do k = 1,3
          fin(1,j,col1(k)) = feq(1,j,col1(k)) + fin(1,j,col3(4-k))
     &                                        - feq(1,j,col3(4-k))
        enddo
        enddo

C-------collision step
        do i = 1,ni
        do j = 1,nj
        do k = 1,nvec
          fout(i,j,k) = fin(i,j,k) - omega*(fin(i,j,k) - feq(i,j,k))
        enddo
        enddo
        enddo

C-------obstacle
        do i = 1,ni
        do j = 1,nj

C---------bounce-back condition
          if (obs(i,j).eq.1) then
            do k = 1,nvec
              fout(i,j,k) = fin(i,j,10-k)
            enddo
          endif

        enddo
        enddo

C-------streaming step
        do i = 1,ni
        do j = 1,nj
        do k = 1,nvec
          nexti = i + vec(k,1)
          if (nexti.lt.1)  nexti = ni
          if (nexti.gt.ni) nexti = 1
          nextj = j + vec(k,2)
          if (nextj.lt.1)  nextj = nj
          if (nextj.gt.nj) nextj = 1
          fin(nexti,nextj,k) = fout(i,j,k)
        enddo
        enddo
        enddo

C-------residiuals
        vmag_sum = 0.0
        do i = 1,ni
        do j = 1,nj
          vmag = sqrt(vel(i,j,1)**2 + vel(i,j,2)**2)
          vmag_sum = vmag_sum + vmag
        enddo
        enddo
        
C-------write iterations
        if (mod(its,writeits).eq.0) write(6,*)its,vmag_sum

C-------write ascii vtk
        if (mod(its,writevtk).eq.0) then

C---------set obstacle v=0
          if (.true.) then
          do i = 1,ni
          do j = 1,nj
            if (obs(i,j).eq.1) then
              vel(i,j,1) = 0.0
              vel(i,j,2) = 0.0
            endif
          enddo
          enddo
          endif

C---------vtk file
          if (iseq.eq.0) then
            outfile = './vtkfiles/lbm2d.vtk'
          elseif (iseq.eq.1) then
            write(unit=string, fmt='(I6.6)') its
            outfile = './vtkfiles/lbm2d_'//string//'.vtk'
          endif
          write(6,*)"outfile=",outfile
          open(UNIT=20,FILE=outfile)
          write(20,10)'# vtk DataFile Version 3.0'
          write(20,10)'vtk output'
          write(20,10)'ASCII'
          write(20,10)'DATASET RECTILINEAR_GRID'
          write(20,20)'DIMENSIONS ',ni+1,nj+1,2
          write(20,30)'X_COORDINATES ',ni+1,' float'
          write(20,*)   (real(i-1),i=1,ni+1)
          write(20,30)'Y_COORDINATES ',nj+1,' float'
          write(20,*)   (real(j-1),j=1,nj+1)
          WRITE(20,30)'Z_COORDINATES ', 2,' float'
          WRITE(20,*)   50.0, 51.0
          write(20,40)'CELL_DATA ',ni*nj
C---------density
          write(20,10)'SCALARS density float'
          write(20,10)'LOOKUP_TABLE default'
          write(20,*)((rho(i,j),i=1,ni),j=1,nj)
C---------velocity
          write(20,10)'VECTORS velocity float'
          write(20,*)((vel(i,j,1),vel(i,j,2),0.0,i=1,ni),j=1,nj)
          close(20)
        endif

      enddo ! maxits
C=====iteration loop
      write(6,*) "finished"

   10 format(A)
   20 format(A,3I4)
   30 format(A,I3,A)
   40 format(A,I9)

      end ! main







