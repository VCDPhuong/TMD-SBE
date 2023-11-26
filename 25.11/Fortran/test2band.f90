!global variable
module globalvar
!contains
!subroutine gbv()
!integer, parameter :: Norbital=5, nkamax=199, comid=0
!integer, parameter :: halfNorbital=2, qe=1
double precision, parameter :: pi=3.141592d0
real, parameter :: varepsilon = 4.2
real, parameter :: epsilon0 = 55.26349406e-3 
! [e^2 eV^-1 nm^-1]
real, parameter :: tmin = -400.0 
! time in fs
real, parameter :: tmax = 400.0
real, parameter :: deltat = 0.02
integer, parameter :: ntmax=int((tmax-tmin)/deltat)
real, parameter :: hbar=0.6582 
! fs
real, parameter :: E0 = 1.
real, parameter :: qe = 1.
! V/nm
real, parameter :: phi = 0.0 
! rad
real, parameter :: TauL = 60*sqrt(2.0) 
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
integer, parameter :: Lightinteraction = 0 
!Lightinterac=1: Yes, 0: No.
integer, parameter :: colid = 0
save
!end subroutine
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
subroutine variable(grid, Xi, Delta, deltaE)
use globalvar
implicit none
integer :: i, j, nu, mu
double precision, dimension(0:nkamax), intent(inout) :: grid
double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(inout) :: Xi
integer, dimension(0:Norbital,0:Norbital), intent(inout) ::  Delta
double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(inout) :: deltaE
double precision, dimension(0:Norbital,0:nkamax) :: E
real :: Eg=1.5
Delta(0,0)=1
Delta(0,1)=0
Delta(1,0)=0
Delta(1,1)=1
do i = 0, nkamax
    grid(i)=-3. + deltak*i
    E(0,i)=-(hbar*grid(i))**2./(2.*0.46*me)
    E(1,i)=(hbar*grid(i))**2./(2.*0.067*me)+Eg
    do nu = 0,Norbital
        do mu =0,Norbital
            deltaE(nu,mu,i)=E(nu,i)-E(mu,i)
        end do
    end do
end do
do nu = 0,Norbital
    do mu = 0, Norbital
        if (Delta(nu,mu)/=0) then
        Xi(nu,mu,:)=0.3
        end if
    end do
end do

end subroutine
subroutine EAlecntinteger(Elec, Alec)
    use globalvar
    implicit none
    integer:: nt
    double precision, dimension(0:ntmax), intent(inout) :: Elec, Alec
    double precision :: time
    Alec(0)=0.
    ! Calculate Elec
    do nt=0,ntmax
    time = nt*deltat + tmin
    Elec(nt) = E0*exp(-(time)**2/TauL**2)*cos(omegaL*time)*cos(phi)
!        Elec(nt,1) = E0*exp(-(time)**2/TauL**2)*cos(omegaL*time)*sin(phi)
    ! Calculate Alec
    if (nt /= 0) then
        Alec(nt) = Alec(nt-1) - Elec(nt-1)*deltat
    end if
    end do
end subroutine
subroutine Elect(nt1, Et,Elec)
    use globalvar
    implicit none
    real, intent(in) :: nt1
    double precision, intent(inout) :: Et
    double precision, dimension(0:ntmax):: Elec
    Et = Elec(int(nt1)) + (Elec(int(nt1)+1)-Elec(int(nt1)))*(nt1-int(nt1))
end subroutine Elect

subroutine Alect(nt1, At,Elec, Alec)
    use globalvar
    implicit none
    real, intent(in) :: nt1
    double precision, intent(inout) :: At
    double precision, dimension(0:ntmax), intent(in) :: Elec, Alec
    At = Alec(int(nt1)) - Elec(int(nt1)) * (nt1 - int(nt1)) * deltat
end subroutine Alect
subroutine RHSsbeLG(RHS, rho,nt1,xi,deltaE, Elec, Alec)
    use globalvar
use omp_lib
implicit none
double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: rho
double complex, dimension(0:Norbital,0:Norbital,0:nkamax) :: gradrho
double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: deltaE
double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(out) :: RHS
double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xi
double complex,allocatable,dimension(:,:,:) :: lightinteractionpart,dummy1
double precision,allocatable, dimension(:,:,:)  :: deltaEpart
double precision, dimension(0:ntmax) :: Elec, Alec
integer :: nu,mu,nk1,nk2
real, intent(in) :: nt1
double precision :: Et
real :: T2 = 5.2
allocate(lightinteractionpart(0:Norbital,0:Norbital,0:nkamax),dummy1(0:Norbital,0:Norbital,0:nkamax)&
        ,deltaEpart(0:Norbital,0:Norbital,0:nkamax))
call Elect(nt1, Et,Elec)
    !!$OMP PARALLEL DO
    !do nk2 = 0, nkamax
        !$OMP PARALLEL DO
        do nk1 = 0, nkamax
            do mu = 0, Norbital
                do nu = 0, Norbital
                    deltaEpart(nu, mu, nk1)= deltaE(nu, mu, nk1)*rho(nu, mu, nk1)
                    dummy1(nu, mu, nk1) = xi(nu, mu, nk1) * Et
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    !end do
    !!$OMP END PARALLEL DO
    call grad(rho,gradrho)
    call commu(dummy1,rho,lightinteractionpart)
    RHS=-cmplx(0,1)/hbar*deltaEpart - cmplx(0,1)/hbar*qe*lightinteractionpart+gradrho*qe/hbar*Et
    deallocate(lightinteractionpart,dummy1,deltaEpart)
end subroutine
subroutine RK4(rho1, nt1, xi1, deltaE1, Elec1, Alec1)
    use globalvar
    implicit none
    double complex, dimension(:,:,:), intent(inout) :: rho1
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xi1
    double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: deltaE1
    double complex,allocatable, dimension(:,:,:):: RHS
    !double complex, dimension(:,:,:), intent(in) :: p
    double complex,allocatable, dimension(:,:,:) :: k1, k2, k3, k4
    integer, intent(in) :: nt1
    double precision, dimension(0:ntmax) :: Elec1, Alec1
    allocate(RHS(0:Norbital,0:Norbital,0:nkamax),k1(0:Norbital,0:Norbital,0:nkamax)&
            ,k2(0:Norbital,0:Norbital,0:nkamax),k3(0:Norbital,0:nkamax,0:nkamax)&
            ,k4(0:Norbital,0:Norbital,0:nkamax))
    ! Calculate k1
    
    call RHSsbeLG(RHS, rho1, real(nt1), xi1, deltaE1, Elec1, Alec1)
    k1 = deltat * RHS

    ! Calculate k2
    call RHSsbeLG(RHS, rho1 + 0.5d0*k1, real(nt1) + 0.5e0,  xi1, deltaE1,  Elec1, Alec1)
    k2 = deltat * RHS

    ! Calculate k3
    call RHSsbeLG(RHS, rho1 + 0.5d0*k2, real(nt1) + 0.5e0,  xi1, deltaE1, Elec1, Alec1)
    k3 = deltat * RHS

    ! Calculate k4
    call RHSsbeLG(RHS, rho1 + k3, real(nt1) + 1.0e0,  xi, deltaE, Elec1, Alec1)
    k4 = deltat * RHS

    ! Update rho
    rho1 = rho1 + (k1 + 2.0e0*k2 + 2.0e0*k3 + k4) / 6.0e0
    deallocate(RHS,k1,k2,k3,k4)
end subroutine RK4
subroutine calculatev( rho, nt, xi, Elec, Alec)
    use globalvar
    use para_MX2
    implicit none
    double complex :: I, Polar
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: rho
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xi
    !double complex, allocatable, dimension(:,:,:,:):: velo
    double complex, allocatable, dimension(:):: N
    double complex, allocatable:: P
    double precision, dimension(0:ntmax) :: Elec, Alec
    integer :: j,nu,mu,nk1,nk2
    double precision :: dS
    double precision :: At
    integer, intent(in) :: nt
    allocate(N(0:Norbital),P)
    dS =(4.*pi/(sqrt(3.)*a*nkamax))**2*sqrt(3.)/2
    call Elect(real(nt),At,Elec)
    
    !do j=0,1
    !    velo(:,:,:,:,j)=(p(:,:,:,:,j))/me
    !end do
    
    !do j=0,1
    
    !do nk2=0,nkamax
        do nk1=0,nkamax
            do nu=0,Norbital
    !            do mu=0,Norbital
    !                I(j)=qe*velo(nu,mu,nk1,nk2,j)*rho(mu,nu,nk1,nk2)*dS
    !            end do
                N(nu)=rho(nu,nu,nk1)*dS/(2*pi)
            end do
            P=1/(2*pi)*(xi(0,1,nk1)*rho(1,0,nk1)+xi(1,0,nk1)*rho(0,1,nk1))
        end do
    !end do
    !end do
    !write(1,'(E12.5, x,E12.5,x,E12.5,x,E12.5,x, E12.5)') tmin+nt*deltat, I
    write(1,'(E12.5, x,E12.5, x,E12.5)') tmin+nt*deltat, P
    write(2,'(E12.5, x,E12.5, x,E12.5)') tmin+nt*deltat, real(N(0)), real(N(1))
    
    !open(unit=1,FILE='OUTPUT/polarization.txt')
    !do j=0,1
    !do nk2=0,nkamax
    !    do nk1=0,nkamax
    !        do nu=0,Norbital
    !            do mu=0,Norbital
    !                Polar(j)=xi(nu,mu,nk1,nk2)*rho(mu,nu,nk1,nk2)*dS
    !            end do
    !        end do
    !    end do
    !end do
    !end do
    !write(1,'(E12.5, E12.5, E12.5)') tmin+nt*deltat, Polar(0), Polar(1)
    !close(1)    
    deallocate(N)
end subroutine
end module
program main
    use globalvar
        !use varia
    use defud
    !use calculatevariable
    !use prepare
        implicit none
    double complex,allocatable, dimension(:,:,:) :: rho
    integer :: count_start, count_end
    integer :: count_rate, count_max
    Real :: start, finish
    double precision :: At
    double precision, dimension(0:ntmax):: Elec, Alec
    integer :: nt
    allocate(rho(0:Norbital,0:Norbital,0:nkamax))
    rho(0,0,:)=complex(1.,0.)
    rho(1,1,:)=complex(0.,0.)
    !rho(2,2,:,:)=complex(0.,0.)
    !rho(3,3,:,:)=complex(1.,0.)
    !rho(4,4,:,:)=complex(0.,0.)
    !rho(5,5,:,:)=complex(0.,0.)
    call cpu_time(start)
    print*,"begin calculate variable"
        call variable(grid, Xi, Delta, deltaE)
    print *, "done Calculate variable"
    !print*, p
    open(unit=1,FILE='OUTPUT/Polarize.txt')
    open(unit=2,FILE='OUTPUT/density.txt')
    open(unit=3,FILE='OUTPUT/Alec.txt')
    call cpu_time(finish)
    print *, finish-start
    call EAlecntinteger(Elec, Alec)
    do nt=1,ntmax
        !call Alect(real(nt), At,Elec, Alec)
        !write(3,'(F14.5, x,F14.5, x,F14.5)') real(nt*deltat +tmin), At(0), At(1)
        call system_clock(count_start, count_rate, count_max)
        start = count_start*1.0/count_rate
        call RK4(rho, nt, xi, deltaE, Elec, Alec)
        call system_clock(count_end, count_rate, count_max)
        finish = count_end*1.0/count_rate

        if (mod(nt,int(10/deltat))==0) then
            print *, "-------------------------------------------------------------------------------------------------"
            write(*,'(A,F14.2,A,F14.2,A,F8.2)')"Real time:", real(tmin+nt*deltat)," on ",tmax, ".Time per step ", finish-start 
            write(*,'(F14.2,F5.2,F5.2)') REAL(SUM(rho(0,0,:)+rho(1,1,:))),&
                sum(real(rho(1,1,:)/((rho(0,0,:)+rho(1,1,:)))))
            !write(*,'(F14.2,F5.2,F5.2)') REAL(SUM(rho(3,3,:,:)+rho(4,4,:,:)+rho(5,5,:,:))),&
            !    sum(real(rho(4,4,:,:)/((rho(3,3,:,:)+rho(4,4,:,:)+rho(5,5,:,:))))),&
            !    sum(real(rho(5,5,:,:)/((rho(3,3,:,:)+rho(4,4,:,:)+rho(5,5,:,:)))))
        end if
        if (mod(nt,int(0.2/deltat))==0) then
        call calculatev( rho, nt, xi, Elec, Alec)
        end if
    end do
    close(1)
    close(2)
    close(3) 
    deallocate(rho)
    print*, "deallocate rho"
end program
