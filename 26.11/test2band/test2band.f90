!global variable
module globalvar
    double precision, parameter :: pi=3.141592d0
    real, parameter :: varepsilon = 4.2
    real, parameter :: epsilon0 = 55.26349406e-3 
    ! [e^2 eV^-1 nm^-1]
    real, parameter :: tmin = -200.0
    ! time in fs
    real, parameter :: tmax = 400.0
    real, parameter :: deltat = 0.02  ![fs]
    integer, parameter :: ntmax=int((tmax-tmin)/deltat)
    real, parameter :: hbar=0.6582 
    ! fs
    real, parameter :: E0 = 0.22
    real, parameter :: qe = 1.
    ! V/nm
    real, parameter :: phi = 0.0 
    ! rad
    real, parameter :: TauL = 50
    ! fs (the TauL is not have the form x^2/2 \tauL**2)
    real, parameter :: epsilonL = 1.6
    ! eV
    real, parameter :: omegaL = epsilonL / hbar 
    ! 1/fs
    real, parameter :: me = 9.1094 / 1.6022 
    ! eV.fs^2/nm^2; E=mc^2;kg=J(s/m)^2=1/qe eV (s/m)^2
    ! Assign values to variables
    integer, parameter :: coulomb = 0 
    !coulomb=1: Yes, 0: No.
    integer, parameter :: Lightinteraction = 1
    !Lightinterac=1: Yes, 0: No.
    integer, parameter :: colid = 0
    save
    end module
    module defud
    integer, parameter :: Norbital=1
    integer, parameter :: nkamax=999
    double precision, parameter :: deltak=(3.+3.)/nkamax
    double precision, dimension(0:nkamax) :: grid
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax) :: Xi
    integer, dimension(0:Norbital,0:Norbital) :: Delta
    double precision, dimension(0:Norbital,0:Norbital,0:nkamax) :: deltaE
    double precision, dimension(0:Norbital,0:nkamax) :: E
    save
    contains
    subroutine grad(rho,gradrho)
    implicit none
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: rho
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(out) :: gradrho
    integer :: i
    do i=0,nkamax
    if (i==0) then
    gradrho(:,:,i)=(rho(:,:,i+1)-rho(:,:,i))/deltak
    elseif (i==nkamax) then
    gradrho(:,:,i)=(rho(:,:,i)-rho(:,:,i-1))/deltak
    else 
    gradrho(:,:,i)=(rho(:,:,i+1)-rho(:,:,i-1))/(2.*deltak)
    end if
    end do
    end subroutine
    subroutine variable(grids, Xis, Deltas, deltaEs)
    use globalvar
    implicit none
    integer :: i, nu, mu
    double precision, dimension(0:nkamax), intent(inout) :: grids
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(inout) :: Xis
    integer, dimension(0:Norbital,0:Norbital), intent(inout) ::  Deltas
    double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(inout) :: deltaEs
    double precision, dimension(0:Norbital,0:nkamax) :: Es
    real :: Eg=1.5
    Deltas(0,0)=1
    Deltas(0,1)=0
    Deltas(1,0)=0
    Deltas(1,1)=1
    open(unit=1,FILE="OUTPUT/grid.txt")
    open(unit=2,FILE="OUTPUT/Energy.txt")
    do i = 0, nkamax
        grids(i)=-3. + deltak*i
        write(1,*) grid(i)
        Es(0,i)=-(hbar*grid(i))**2./(2.*0.46*me)
        Es(1,i)=(hbar*grid(i))**2./(2.*0.067*me)+Eg
        write(2,'(E12.5, x,E12.5, x,E12.5)') grid(i), Es(0,i), Es(1,i)
        do mu = 0,Norbital
            do nu =0,Norbital
                deltaEs(nu,mu,i)=Es(nu,i)-Es(mu,i)
            end do
        end do
    end do
    close(1)
    close(2)
    do nu = 0,Norbital
        do mu = 0, Norbital
            !if (nu/=mu) then
            Xis(nu,mu,:)=0.3
            !end if
        end do
    end do
    
    end subroutine
    subroutine Elecrealt(ts,Et)
        use globalvar
        real, intent(in) :: ts
        double precision, intent(out) :: Et
        Et=E0*exp(-(ts)**2/TauL**2)*cos(omegaL*ts)
    end subroutine
    subroutine commu(A, B, C)
        use omp_lib
        implicit none
        double complex, dimension(:,:,:), intent(in) :: A, B
        double complex, dimension(:,:,:), intent(out) :: C
        integer :: n
            !$OMP PARALLEL DO
            do n = lbound(A,3), ubound(A,3)
                ! Calculate the commutator using matmul
                C(:,:,n) = matmul(A(:,:,n), B(:,:,n)) - matmul(B(:,:,n), A(:,:,n))
            end do
            !$OMP END PARALLEL DO
    end subroutine commu
    subroutine RHSsbeLG(RHSs, rhos , nt1s ,xis ,deltaEs)
        use globalvar
        use omp_lib
        implicit none
        double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: rhos
        double complex, allocatable, dimension(:,:,:) :: gradrhos
        double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: deltaEs
        double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(out) :: RHSs
        double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xis
        double complex, allocatable,dimension(:,:,:) :: lightinteractionparts,dummy1s
        double complex, allocatable, dimension(:,:,:)  :: deltaEparts
        integer :: nu,mu,nk1
        real, intent(in) :: nt1s
        double precision :: Ets
        !real :: T2 = 5.2
        allocate(lightinteractionparts(0:Norbital,0:Norbital,0:nkamax),dummy1s(0:Norbital,0:Norbital,0:nkamax)&
                ,deltaEparts(0:Norbital,0:Norbital,0:nkamax), gradrhos(0:Norbital,0:Norbital,0:nkamax))
        call Elecrealt(nt1s,Ets)
        write(3,'(F14.5, x,F14.5)') nt1s, Ets
            !$OMP PARALLEL DO
            do nk1 = 0, nkamax
                do mu = 0, Norbital
                    do nu = 0, Norbital
                        deltaEparts(nu, mu, nk1)= deltaEs(mu,nu, nk1)*rhos(nu, mu, nk1)
                        dummy1s(nu,mu,nk1)= xis(nu,mu,nk1)* Ets
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        call grad(rhos,gradrhos)
        call commu(dummy1s,rhos,lightinteractionparts)
        RHSs=-cmplx(0,1)/hbar*deltaEparts - cmplx(0,1)/hbar*qe*lightinteractionparts + gradrhos*qe/hbar*Ets
        deallocate(lightinteractionparts,dummy1s,deltaEparts,gradrhos)
    end subroutine

    subroutine RK4(rho, nt, xis, deltaEs)
        use globalvar
        implicit none
        double complex, dimension(:,:,:), intent(inout) :: rho
        double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xis
        double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: deltaEs
        double complex,allocatable, dimension(:,:,:):: RHS
        !double complex, dimension(:,:,:), intent(in) :: p
        double complex,allocatable, dimension(:,:,:,:) :: Rdummy
        integer, intent(in) :: nt
        real :: time
        allocate(RHS(0:Norbital,0:Norbital,0:nkamax),Rdummy(0:Norbital,0:Norbital,0:nkamax,0:4))
        ! Calculate k1
        time=tmin+deltat*nt
        Rdummy(:,:,:,0) = rho
        call RHSsbeLG(Rdummy(:,:,:,1), rho, time, xis, deltaEs)
    
        ! Calculate k2
        call RHSsbeLG(Rdummy(:,:,:,2), rho + deltat * 0.5* Rdummy(:,:,:,1), time+ 0.5*deltat,  xis, deltaEs)
    
        ! Calculate k3
        call RHSsbeLG(Rdummy(:,:,:,3), rho + deltat * 0.5* Rdummy(:,:,:,2), time+ 0.5*deltat,  xis, deltaEs)
    
        ! Calculate k4
        call RHSsbeLG(Rdummy(:,:,:,4), rho + deltat * Rdummy(:,:,:,3), time+ deltat,  xis, deltaEs)
    
        ! Update rho
        rho = Rdummy(:,:,:,0) + deltat *(Rdummy(:,:,:,1) + 2.0e0*Rdummy(:,:,:,2)&
                             + 2.0e0*Rdummy(:,:,:,3) + Rdummy(:,:,:,4)) / 6.0e0
        deallocate(RHS,Rdummy)
    end subroutine RK4
    subroutine calculatev( rhos, nts, xis)
        use globalvar
        implicit none
        double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: rhos
        double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xis
        double complex, allocatable, dimension(:):: N
        double complex, allocatable:: P
        integer :: nu,nk1
        double precision :: dS
        !double precision :: At
        integer, intent(in) :: nts
        allocate(N(0:Norbital),P)
        P=0.
        dS=deltak
        !call Elecrealt(nts*deltat+tmin,At)
            do nk1=0,nkamax
                P=P+1/(2*pi)*(xis(0,1,nk1)*rhos(1,0,nk1)+xis(1,0,nk1)*rhos(0,1,nk1))
            end do
            do nu=0,Norbital
                N(nu)=sum(rhos(nu,nu,:))*dS/(2*pi)
            end do
        write(1,'(E12.4, x,E12.5, x,E12.5)') tmin+nts*deltat, P
        write(2,'(E12.4, x,E12.5, x,E12.5, x,E12.5, x,E12.5)') tmin+nts*deltat, N(0), N(1)
        deallocate(N)
    end subroutine
    end module
    program main
        use globalvar
        use defud
        implicit none
        double complex,allocatable, dimension(:,:,:) :: rho
        integer :: count_start, count_end
        integer :: count_rate, count_max
        Real :: start, finish
        integer :: nt
        allocate(rho(0:Norbital,0:Norbital,0:nkamax))
        rho(0,0,:)=complex(1.,0.)
        rho(1,1,:)=complex(0.,0.)
        call cpu_time(start)
        print*,"begin calculate variable"
        call variable(grid, Xi, Delta, deltaE)
        print *, "done Calculate variable"
        !print*, p
        open(unit=1,FILE='OUTPUT/Polarize.txt')
        open(unit=2,FILE='OUTPUT/density.txt')
        open(unit=3,FILE='OUTPUT/Elec.txt')
        call cpu_time(finish)
        print *, finish-start
        do nt=1,ntmax
            call system_clock(count_start, count_rate, count_max)
            start = count_start*1.0/count_rate
            call RK4(rho, nt, xi, deltaE)
            call system_clock(count_end, count_rate, count_max)
            finish = count_end*1.0/count_rate
    
            if (mod(nt,int(10/deltat))==0) then
                print *, "-------------------------------------------------------------------------------------------------"
                write(*,'(A,F14.2,A,F14.2,A,F8.2)')"Real time:", real(tmin+nt*deltat)," on ",tmax, ".Time per step ", finish-start 
                write(*,'(F14.2,F5.2,F14.2,F5.2)') SUM(rho(0,0,:)),SUM(rho(1,1,:))
            end if
            if (mod(nt,int(0.2/deltat))==0) then
            call calculatev( rho, nt, xi)
            end if
        end do
        close(1)
        close(2)
        close(3) 
        deallocate(rho)
        print*, "deallocate rho"
    end program
    