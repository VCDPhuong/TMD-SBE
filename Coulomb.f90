!global variable
module globalvar
    integer, parameter :: halfNorbital = 2
    integer, parameter :: Norbital = 5
    integer, parameter :: nkamax = 249
    double precision, parameter :: pi=3.141592d0
    real, parameter :: varepsilon = 4.2
    real, parameter :: epsilon0 = 55.26349406e-3 
    ! [e^2 eV^-1 nm^-1]
    real, parameter :: tmin = -200.0
    ! time in fs
    real, parameter :: tmax = 400.0
    real, parameter :: deltat = 0.01  ![fs]
    integer, parameter :: ntmax=int((tmax-tmin)/deltat)
    real, parameter :: hbar=0.6582 
    ! fs
    real, parameter :: E0 = 1e-2 ![V/m]
    real, parameter :: qe = 1.
    ! V/nm
    real, parameter :: phi = 0e0
    ! rad
    real, parameter :: TauL = 50
    ! fs (the TauL is not have the form x^2/2 \tauL**2)
    real, parameter :: epsilonL = 2.0
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
    integer, parameter :: comid = 0 !0 MoS2
    save
end module globalvar
    module para_MX2
        double precision :: a,e1,e2,t0,t1,t2,t11,t12,t22,r0,r1,r2,r11,r12,&
                            u0,u1,u2,u11,u12,u22,lambda
        save
     contains
        subroutine tb_parameters(material)
           implicit none
           integer :: material
           double precision, dimension(0:5) :: &
                  atab=(/0.3190,0.3191,0.3326,0.3325,0.3557,0.3560/),&
                  e1tab=(/0.683,0.717,0.684,0.728,0.588,0.697/),&
                  e2tab=(/1.707,1.916,1.546,1.655,1.303,1.380/),&
                  t0tab=(/-0.146,-0.152,-0.146,-0.146,-0.226,-0.109/),&
                  t1tab=(/-0.114,-0.097,-0.130,-0.124,-0.234,-0.164/),&
                  t2tab=(/0.506,0.590,0.432,0.507,0.036,0.368/),&
                  t11tab=(/0.085,0.047,0.144,0.117,0.400,0.204/),&
                  t12tab=(/0.162,0.178,0.117,0.127,0.098,0.093/),&
                  t22tab=(/0.073,0.016,0.075,0.015,0.017,0.038/),&
                  r0tab=(/0.060,0.069,0.039,0.036,0.003,-0.015/),&
                  r1tab=(/-0.236,-0.261,-0.209,-0.234,-0.025,-0.209/),&
                  r2tab=(/0.067,0.107,0.069,0.107,-0.169,0.107/),&
                  r11tab=(/0.016,-0.003,0.052,0.044,0.082,0.115/),&
                  r12tab=(/0.087,0.109,0.060,0.075,0.051,0.009/),&
                  u0tab=(/-0.038,-0.054,-0.042,-0.061,0.057,-0.066/),&
                  u1tab=(/0.046,0.045,0.036,0.032,0.103,0.011/),&
                  u2tab=(/0.001,0.002,0.008,0.007,0.187,-0.013/),&
                  u11tab=(/0.266,0.325,0.272,0.329,-0.045,0.312/),&
                  u12tab=(/-0.176,-0.206,-0.172,-0.202,-0.141,-0.177/),&
                  u22tab=(/-0.150,-0.163,-0.150,-0.164,0.087,-0.132/),&
                  lambdatab=(/0.073,0.211,0.091,0.228,0.107,0.237/)
           a=atab(material)
           e1=e1tab(material)
           e2=e2tab(material)
           t0=t0tab(material)
           t1=t1tab(material)
           t2=t2tab(material)
           t11=t11tab(material)
           t12=t12tab(material)
           t22=t22tab(material)
           r0=r0tab(material)
           r1=r1tab(material)
           r2=r2tab(material)
           r11=r11tab(material)
           r12=r12tab(material)
           u0=u0tab(material)
           u1=u1tab(material)
           u2=u2tab(material)
           u11=u11tab(material)
           u12=u12tab(material)
           u22=u22tab(material)
           lambda=lambdatab(material)
        end subroutine
     end module
     !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     module prepare
     use para_MX2
     use globalvar
     contains
     subroutine tbham(kx,ky,hamu,hamd)
        use globalvar
        use para_MX2
        implicit none
        double precision, intent(IN) :: kx, ky
        complex*16, intent(out), dimension(0:halfNorbital,0:halfNorbital) :: hamu,hamd
        double precision :: alpha,beta
        double complex, ALLOCATABLE ,dimension(:,:) :: ham0,ham1
        call tb_parameters(comid)
        allocate(ham0(0:halfNorbital,0:halfNorbital),ham1(0:halfNorbital,0:halfNorbital))
        alpha = kx *a/2.
        beta = sqrt(3.) * ky * a/2.
        ham0(0,0)=e1 + 2.*t0*(2 *cos(alpha) *  cos(beta) + cos(2*alpha))&
                     + 2.*r0*(2 *cos(3*alpha)* cos(beta) + cos(2*beta))&
                     + 2.*u0*(2 *cos(2*alpha)* cos(2*beta) + cos(4*alpha))
        ham0(0,1)=cmplx(-2.*sqrt(3.)*t2*sin(alpha)*sin(beta)&
                        +2.*(r1+r2)*sin(3*alpha)*sin(beta)&
                        -2.*sqrt(3.)*u2*sin(2*alpha)*sin(2*beta),2.*t1*sin(alpha)*(2*cos(alpha)+cos(beta))&
                        +2.*(r1-r2)*sin(3*alpha)*cos(beta)&
                        +2.*u1*sin(2*alpha)*(2*cos(2*alpha)+cos(2*beta)), kind=8)
        ham0(0,2)=cmplx(2.*t2*(cos(2*alpha)-cos(alpha)*cos(beta))&
                        -2.*(r1+r2)*(cos(3*alpha)*cos(beta)-cos(2*beta))/sqrt(3.)&
                        +2.*u2*(cos(4*alpha)-cos(2*alpha)*cos(2*beta)),2.*sqrt(3.)*t1*cos(alpha)*sin(beta)&
                        +2.*(r1-r2)*sin(beta)*(cos(3*alpha)+2*cos(beta))/sqrt(3.)&
                        +2.*sqrt(3.)*u1*cos(2*alpha)*sin(2*beta), kind=8)
        ham0(1,0)=conjg(ham0(0,1))
        ham0(1,1)=e2+(t11+3.*t22)*cos(alpha)*cos(beta)+2.*t11*cos(2.*alpha)&
                    +4.*r11*cos(3.*alpha)*cos(beta)+2.*(r11+sqrt(3.)*r12)*cos(2.*beta)&
                    +(u11+3.*u22)*cos(2.*alpha)*cos(2.*beta)+2.*u11*cos(4.*alpha)
        ham0(1,2)=cmplx(sqrt(3.)*(t22-t11)*sin(alpha)*sin(beta)&
                        +4.*r12*sin(3.*alpha)*sin(beta)&
                        +sqrt(3.)*(u22-u11)*sin(2.*alpha)*sin(2.*beta), 4*t12*sin(alpha)*(cos(alpha)-cos(beta))&
                        +4.*u12*sin(2.*alpha)*(cos(2.*alpha)-cos(2.*beta)), kind=8)
        ham0(2,0)=conjg(ham0(0,2))
        ham0(2,1)=conjg(ham0(1,2))
        ham0(2,2)=e2+(3*t11+t22)*cos(alpha)*cos(beta)+2*t22*cos(2*alpha)&
                    +2*r11*(2*cos(3*alpha)*cos(beta)+cos(2*beta))&
                    +2*r12*(4*cos(3*alpha)*cos(beta)-cos(2*beta))/sqrt(3.)&
                    +(3*u11+u22)*cos(2*alpha)*cos(2*beta)+2*u22*cos(4*alpha)
        ham1(0,0)=0.
        ham1(0,1)=0.
        ham1(0,2)=0.
        ham1(1,0)=0.
        ham1(1,1)=0.
        ham1(1,2)=lambda*cmplx(0.,1.)
        ham1(2,0)=0.
        ham1(2,1)=lambda*cmplx(0.,-1.)
        ham1(2,2)=0.
        hamu=ham0+ham1
        hamd=ham0-ham1
        deallocate(ham0,ham1)
     end subroutine
     subroutine gradient_H(kx, ky, grad_H)
        use globalvar
        implicit none
        double precision, intent(In) :: kx, ky
        double complex, allocatable, dimension(:,:) :: hamu_plus, hamu_minus, hamd_plus, hamd_minus
        double complex, dimension(6,6,2), INTENT(OUT) :: grad_H
        real :: delta = 5e-2
        allocate(hamd_plus(3,3), hamd_minus(3,3), hamu_plus(3,3), hamu_minus(3,3))
           call tbham(kx + delta, ky, hamu_plus , hamd_plus )
           call tbham(kx - delta, ky, hamu_minus, hamd_minus)
           grad_H(1:3,1:3,1) = (hamu_plus - hamu_minus) / (2 * delta)
           grad_H(4:6,4:6,1) = (hamd_plus - hamd_minus) / (2 * delta)
           call tbham(kx, ky + delta, hamu_plus , hamd_plus )
           call tbham(kx, ky - delta, hamu_minus, hamd_minus)
           grad_H(1:3,1:3,2) = (hamu_plus - hamu_minus) / (2 * delta)
           grad_H(4:6,4:6,2) = (hamd_plus - hamd_minus) / (2 * delta)
        deallocate(hamd_plus, hamd_minus, hamu_plus, hamu_minus)
        end subroutine gradient_H
     end module prepare
     module varia
        use globalvar
        double precision,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax, 0:1) :: deltaE
        double precision,dimension(0:Norbital,0:nkamax,0:nkamax, 0:1) :: E
        complex*16 ,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax, 0:1) :: V
        complex*16 ,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1, 0:1) :: p,Xi
        double precision,dimension(0:nkamax,0:nkamax,0:1, 0:1) :: grid
        integer, dimension(0:Norbital,0:Norbital) :: Delta
        save
        contains
        subroutine variable(grids,ps,Xis,Deltas,deltaEs,Es,Vs)
        use prepare
        use omp_lib
        use para_MX2
        use globalvar
        implicit none
        double precision,dimension(0:Norbital,0:nkamax,0:nkamax,0:1),intent(out) :: Es
        double precision,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1),intent(out) :: deltaEs
        double complex ,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1,0:1),intent(out) :: ps, xis
        double precision,dimension(0:nkamax,0:nkamax,0:1,0:1),intent(out) :: grids
        integer,dimension(0:Norbital,0:Norbital),intent(out) :: Deltas
        double precision,dimension(0:1,0:1) :: B
        double complex,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1),intent(out) :: Vs
        double complex, dimension(:,:), Allocatable :: Vpkx, Vpky
        double complex, dimension(:,:,:), Allocatable  :: Vpkxy
        double complex, Allocatable ,dimension( :, :) :: hamu,hamd
        double precision, dimension(0:halfNorbital) :: E_sorted
        double complex :: WORK(5)
        double precision :: kx,ky,RWORK(7)
        integer :: nk1, nk2, nu, mu,INFO 
        double precision, dimension(0:1) :: array
        external :: zheev
        allocate(hamu(0:halfNorbital,0:halfNorbital),hamd(0:halfNorbital,0:halfNorbital), &
                 Vpkx(0:Norbital,0:Norbital), Vpky(0:Norbital,0:Norbital), &
                 Vpkxy(0:Norbital,0:Norbital,0:1))
        call tb_parameters(comid)
        ! Initialize nk1, nk2, and nkamax here
        B(0,0) =  2.*pi/a
        B(0,1) =  2.*pi/a
        B(1,0) =  2.*pi/(sqrt(3.)*a)
        B(1,1) = -2.*pi/(sqrt(3.)*a)
        print *, "begin variable"
        print *, "grid begin"
        do nk2 = 0,nkamax
            do nk1 =0,nkamax
                array(0) = real(nk1)/real(2*nkamax)
                array(1) = real(nk2)/real(2*nkamax)
                grid( nk1, nk2, :, 0) = matmul(B, array)
                grids( nk1, nk2, 0, 0) = grids( nk1, nk2, 0, 0) - 2.*pi/a
                grids( nk1, nk2, 1, 0) = grids( nk1, nk2, 1, 0) - 0.
                array(0) = real(nk1 + nkamax)/real(2*nkamax)
                array(1) = real(nk2 + nkamax)/real(2*nkamax)
                grid( nk1, nk2, :, 1) = matmul(B, array)
                grids( nk1, nk2, 0, 1) = grids( nk1, nk2, 0, 1) - 2.*pi/a
                grids( nk1, nk2, 1, 1) = grids( nk1, nk2, 1, 1) - 0.
            end do
        end do
            open(unit=1,FILE='OUTPUT/grid.txt')
            !$OMP PARALLEL DO
            do nk2 = 0,nkamax
                do nk1 =0,nkamax
                write(1,'(E12.5, 2x, E12.5)') grids( nk1, nk2, 0, 0), grids( nk1, nk2, 1, 0)
                write(1,'(E12.5, 2x, E12.5)') grids( nk1, nk2, 0, 1), grids( nk1, nk2, 1, 1)
                end do
            end do
            !$OMP END PARALLEL DO
            close(1)
            !open(unit=1,FILE='OUTPUT/Energy.txt')
            !open(unit=2,FILE='OUTPUT/p.txt')
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(kx, ky, hamu, hamd, E_sorted, WORK, RWORK, INFO, Vpkxy, Vpkx, Vpky)
            do nk2 = 0,nkamax
               do nk1 =0,nkamax
                    kx = grids( nk1, nk2, 0 ,0)
                    ky = grids( nk1, nk2, 1 ,0)
                    call tbham( kx, ky, hamu, hamd)
                    call zheev('V','U',3,hamu,3,E_sorted,WORK,5,RWORK,INFO)
                    !sort eigenvalue and eigenvector
                    Es(0:2,nk1,nk2, 0) = E_sorted
                    Vs(0:2,0:2,nk1,nk2, 0)=hamu
                    call zheev('V','U',3,hamd,3,E_sorted,WORK,5,RWORK,INFO)
                    !final index is 0 mean eigenvalue is ascending order
                    !sort eigenvalue and eigenvector
                    Es(3:5,nk1,nk2, 0) = E_sorted
                    Vs(3:5,3:5,nk1,nk2, 0)=hamd
                    write(1,'(F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6)')&
                    kx, ky, Es(0,nk1,nk2, 0),Es(1,nk1,nk2, 0),Es(2,nk1,nk2, 0),Es(3,nk1,nk2, 0),Es(4,nk1,nk2, 0),Es(5,nk1,nk2, 0)

                    kx = grids( nk1, nk2, 0 ,1)
                    ky = grids( nk1, nk2, 1 ,1)
                    call tbham( kx, ky, hamu, hamd)
                    call zheev('V','U',3,hamu,3,E_sorted,WORK,5,RWORK,INFO)
                    !sort eigenvalue and eigenvector
                    Es(0:2,nk1,nk2, 1) = E_sorted
                    Vs(0:2,0:2,nk1,nk2, 1)=hamu
                    call zheev('V','U',3,hamd,3,E_sorted,WORK,5,RWORK,INFO)
                    !final index is 0 mean eigenvalue is ascending order
                    !sort eigenvalue and eigenvector
                    Es(3:5,nk1,nk2, 1) = E_sorted
                    Vs(3:5,3:5,nk1,nk2, 1)=hamd

                    write(1,'(F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6, 2x, F12.6)')&
                    kx, ky, Es(0,nk1,nk2, 1),Es(1,nk1,nk2, 1),Es(2,nk1,nk2, 1),Es(3,nk1,nk2, 1),Es(4,nk1,nk2, 1),Es(5,nk1,nk2, 1)
                    
                    kx = grids( nk1, nk2, 0 ,0)
                    ky = grids( nk1, nk2, 1 ,0)
                    call gradient_H(kx ,ky ,Vpkxy)
                    Vpkx=Vpkxy(:,:,0)
                    Vpky=Vpkxy(:,:,1)
                    do mu = 0, Norbital
                        do nu = 0, Norbital
                        ps(nu,mu,nk1,nk2, 0, 0)=me /hbar *DOT_PRODUCT(Vs(:, nu, nk1, nk2, 0),matmul(Vpkx, Vs(:, mu, nk1, nk2, 0)))
                        ps(nu,mu,nk1,nk2, 1, 0)=me /hbar *DOT_PRODUCT(Vs(:, nu, nk1, nk2, 0),matmul(Vpky, Vs(:, mu, nk1, nk2, 0)))
                        if (nu>=mu) then
                        !write(2,'(I5, 2x, I5, 2x, F12.6, 2x, F12.6, 2x, F15.10, 2x, F15.10, 2x, F15.10, 2x, F15.10)')&
                        !    nu ,mu ,kx ,ky ,real(ps(nu,mu,nk1,nk2,0, 0)) ,imag(ps(nu,mu,nk1,nk2,0, 0)),&
                        !    real(ps(nu,mu,nk1,nk2,1, 0)) ,imag(ps(nu,mu,nk1,nk2,1, 0))
                        end if
                        deltaEs(nu,mu,nk1,nk2,0)=Es(nu,nk1,nk2,0)-Es(mu,nk1,nk2,0)
                            if (abs(deltaEs(nu, mu, nk1, nk2 ,0)) < 1e-2) then
                                Deltas(nu, mu) = Deltas(nu, mu) + 1
                            end if
                        end do
                    end do
                    kx = grids( nk1, nk2, 0 ,1)
                    ky = grids( nk1, nk2, 1 ,1)
                    call gradient_H(kx ,ky ,Vpkxy)
                    Vpkx=Vpkxy(:,:,0)
                    Vpky=Vpkxy(:,:,1)
                    do mu = 0, Norbital
                        do nu = 0, Norbital
                        ps(nu,mu,nk1,nk2,0,1)=me/hbar*DOT_PRODUCT(Vs(:,nu,nk1,nk2, 1),matmul(Vpkx, Vs(:, mu, nk1, nk2, 1)))
                        ps(nu,mu,nk1,nk2,1,1)=me/hbar*DOT_PRODUCT(Vs(:,nu,nk1,nk2, 1),matmul(Vpky, Vs(:, mu, nk1, nk2, 1)))
                        if (nu>=mu) then
                        !write(2,'(I5, 2x, I5, 2x, F12.6, 2x, F12.6, 2x, F15.10, 2x, F15.10, 2x, F15.10, 2x, F15.10)')&
                        !    nu ,mu ,kx ,ky ,real(ps(nu,mu,nk1,nk2,0, 1)) ,imag(ps(nu,mu,nk1,nk2,1, 1)),&
                        !    real(ps(nu,mu,nk1,nk2,1, 1)) ,imag(ps(nu,mu,nk1,nk2,1, 1))
                        end if
                        deltaEs(nu,mu,nk1,nk2, 1)=Es(nu,nk1,nk2, 1)-Es(mu,nk1,nk2, 1)
                            if (abs(deltaEs(nu, mu, nk1, nk2 , 1)) < 1e-2) then
                                Deltas(nu, mu) = Deltas(nu, mu) + 1
                            end if
                        end do
                    end do
                end do  !nk2 do
            end do !nk1 do
            !$OMP END PARALLEL DO
            !close(1)
            !close(2)
            do nu = 0, Norbital
                do mu = 0, Norbital
                    if (Deltas(nu, mu) /= 0) then
                        Deltas(nu, mu) = 1
                    end if
                end do
            end do
            !open(unit=1,FILE='OUTPUT/Xi.txt')
            !do nk2 = 0, nkamax
            !    do nk1 = 0, nkamax
            !        do nu = 0, Norbital
            !            do mu = 0, Norbital
            !                if (Deltas(nu, mu) == 0) then
            !                    xis(nu, mu, nk1, nk2, 0, :) = hbar / me * ps(nu, mu, nk1, nk2, 0, :) / deltaEs(nu, mu, nk1, nk2, :)
            !                    xis(nu, mu, nk1, nk2, 1, :) = hbar / me * ps(nu, mu, nk1, nk2, 1, :) / deltaEs(nu, mu, nk1, nk2, :)
            !                else
            !                    xis(nu, mu, nk1, nk2, :, :) = 0.
            !                end if
            !                if (nu > mu) then
            !                    write(1, '(I5, 2x, I5, 2x, F12.6, 2x, F12.6, 2x, F15.10, 2x, F15.10, 2x, F15.10, 2x, F15.10)') &
            !                        nu, mu, grids(nk1, nk2, 0, 0), grids(nk1, nk2, 1, 0), real(xis(nu, mu, nk1, nk2, 0, 0)), &
            !                        imag(xis(nu, mu, nk1, nk2, 0, 0)), real(xis(nu, mu, nk1, nk2, 1, 0)), &
            !                        imag(xis(nu, mu, nk1, nk2, 1, 0))
            !                    write(1, '(I5, 2x, I5, 2x, F12.6, 2x, F12.6, 2x, F15.10, 2x, F15.10, 2x, F15.10, 2x, F15.10)') &
            !                        nu, mu, grids(nk1, nk2, 0, 1), grids(nk1, nk2, 1, 1), real(xis(nu, mu, nk1, nk2, 0, 1)), &
            !                        imag(xis(nu, mu, nk1, nk2, 0, 1)), real(xis(nu, mu, nk1, nk2, 1, 1)), &
            !                        imag(xis(nu, mu, nk1, nk2, 1, 1))
            !                end if
            !            end do
            !        end do
            !    end do
            !end do
            !close(1)
        end subroutine variable
    end module varia
    module Elect
    use globalvar
    implicit none
    double precision, dimension(0:ntmax,0:1) :: Alec,Elec
    contains
    subroutine El(Elecs,Alecs)
        use globalvar
        implicit none
        double precision, dimension(0:ntmax,0:1), intent(inout) :: Alecs,Elecs
        integer, allocatable :: nt
        real, allocatable :: time
        allocate(nt,time)
        open(unit=1,FILE='OUTPUT/Elec.txt')
        !open(unit=2,FILE='OUTPUT/Alec.txt')
        !$OMP PARALLEL DO PRIVATE(nt,time)
        do nt = 0, ntmax
            time=tmin+nt*deltat
            Elecs(nt,0)=E0*exp(-(time)**2/TauL**2)*cos(omegaL*time)*cos(phi)
            Elecs(nt,1)=E0*exp(-(time)**2/TauL**2)*cos(omegaL*time)*sin(phi)
            write(1, '(F14.5, x, F14.5, x, F14.5)') time, Elecs(nt,0), Elecs(nt,1)
        end do
        !$OMP END PARALLEL DO
        do nt = 0,ntmax
            time=tmin+nt*deltat
            if (nt == 0) then
                Alecs(nt,0) = 0.
                Alecs(nt,1) = 0.
            else
                Alecs(nt,:)=Alecs(nt-1,:)-Elecs(nt,:)*deltat
            end if
            !write(2, '(F14.5, x, F14.5, x, F14.5)') time, Alecs(nt,0), Alecs(nt,1)
        end do
        close(1)
        !close(2)
        deallocate(nt,time)
    end subroutine
    end module Elect
    module def
        use globalvar
        !use varia
        contains
    subroutine commu(A, B, C)
        use omp_lib
        implicit none
        double complex, dimension(:,:,:,:), intent(in) :: A, B
        double complex, dimension(:,:,:,:), intent(out) :: C
        integer :: nk1,nk2
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nk1,nk2)
        do nk2 = lbound(A,3), ubound(A,3)
            do nk1 = lbound(A,3), ubound(A,3)
                ! Calculate the commutator using matmul
                C(:,:,nk1,nk2) = matmul(A(:,:,nk1,nk2),B(:,:,nk1,nk2))-matmul(B(:,:,nk1,nk2),A(:,:,nk1,nk2))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine commu
        subroutine grad(rhos,gradrhos)
            use omp_lib
            use para_MX2
            use globalvar
            implicit none
            double precision :: deltak
            double complex, dimension(:,:,:,:), intent(in) :: rhos
            double complex, dimension(:,:,:,:,:), intent(out) :: gradrhos
            integer :: i
            call tb_parameters(comid)
            deltak=4.*pi/(sqrt(3.)*a*nkamax)
            ! Calculate the gradient along the third dimension and save it in gradrhos(:,:,:,:,1)
            !$OMP PARALLEL DO
            do i=lbound(rhos,3),ubound(rhos,3)
                if (i==lbound(rhos,3)) then
                    gradrhos(:,:,i,:,1)=(rhos(:,:,i+1,:)-rhos(:,:,i,:))/deltak
                elseif (i==ubound(rhos,3)) then
                    gradrhos(:,:,i,:,1)=(rhos(:,:,i,:)-rhos(:,:,i-1,:))/deltak
                else 
                    gradrhos(:,:,i,:,1)=(rhos(:,:,i+1,:)-rhos(:,:,i-1,:))/(2.*deltak)
                end if
            end do
            !$OMP END PARALLEL DO
        
            ! Calculate the gradient along the fourth dimension and save it in gradrhos(:,:,:,:,2)
            !$OMP PARALLEL DO
            do i=lbound(rhos,4),ubound(rhos,4)
                if (i==lbound(rhos,4)) then
                    gradrhos(:,:,:,i,2)=(rhos(:,:,:,i+1)-rhos(:,:,:,i))/deltak
                elseif (i==ubound(rhos,4)) then
                    gradrhos(:,:,:,i,2)=(rhos(:,:,:,i)-rhos(:,:,:,i-1))/deltak
                else 
                    gradrhos(:,:,:,i,2)=(rhos(:,:,:,i+1)-rhos(:,:,:,i-1))/(2.*deltak)
                end if
            end do
            !$OMP END PARALLEL DO
        end subroutine
        subroutine Alecrealt(ts,Ats,As)
            use globalvar
            implicit none
            double precision, dimension(:),INTENT(OUT) :: As
            double precision, dimension(:,:), INTENT(IN) :: Ats
            real :: ts
            integer, allocatable :: nts
            real, allocatable :: ntsreal
            allocate(nts, ntsreal)
            ntsreal=(ts - tmin) /deltat
            nts = int(ntsreal)
            As=Ats(nts,:)*(ntsreal-nts)+Ats(nts+1,:)*(1-(ntsreal-nts))
            deallocate(nts, ntsreal)
        end subroutine
        !subroutine Elecrealt(ts,Et)
    !    use globalvar
    !    real, intent(in) :: ts
    !    double precision, intent(out) :: Et
    !    Et=E0*exp(-(ts)**2/TauL**2)*cos(omegaL*ts)
    !end subroutine
    
    !subroutine RHSsbeLG(RHSs, rhos , nt1s ,xis ,deltaEs)
    !    use globalvar
    !    use omp_lib
    !    implicit none
    !    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: rhos
    !    double complex, allocatable, dimension(:,:,:) :: gradrhos
    !    double precision, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in)  :: deltaEs
    !    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(out) :: RHSs
    !    double complex, dimension(0:Norbital,0:Norbital,0:nkamax), intent(in) :: xis
    !    double complex, allocatable,dimension(:,:,:) :: lightinteractionparts,dummy1s
    !    double complex, allocatable, dimension(:,:,:)  :: deltaEparts
    !    integer :: nu,mu,nk1
    !    real, intent(in) :: nt1s
    !    double precision :: Ets
    !    !real :: T2 = 5.2
    !    allocate(lightinteractionparts(0:Norbital,0:Norbital,0:nkamax),dummy1s(0:Norbital,0:Norbital,0:nkamax)&
    !            ,deltaEparts(0:Norbital,0:Norbital,0:nkamax), gradrhos(0:Norbital,0:Norbital,0:nkamax))
    !    call Elecrealt(nt1s,Ets)
    !    write(3,'(F14.5, x,F14.5)') nt1s, Ets
    !        !$OMP PARALLEL DO
    !        do nk1 = 0, nkamax
    !            do mu = 0, Norbital
    !                do nu = 0, Norbital
    !                    deltaEparts(nu, mu, nk1)= deltaEs(mu,nu, nk1)*rhos(nu, mu, nk1)
    !                    dummy1s(nu,mu,nk1)= xis(nu,mu,nk1)* Ets
    !                end do
    !            end do
    !        end do
    !        !$OMP END PARALLEL DO
    !    call grad(rhos,gradrhos)
    !    call commu(dummy1s,rhos,lightinteractionparts)
    !    RHSs=-cmplx(0,1)/hbar*deltaEparts - cmplx(0,1)/hbar*qe*lightinteractionparts + gradrhos*qe/hbar*Ets
    !    deallocate(lightinteractionparts,dummy1s,deltaEparts,gradrhos)
    !end subroutine
        subroutine RHSsbeVG(RHSs, rhos , nt1s ,ps ,deltaEs, Ats)
            use globalvar
            use omp_lib
            implicit none
            double complex,   dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: rhos
            double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: deltaEs
            double complex,   dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(out) :: RHSs
            double complex,   dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: ps
            double precision, dimension(:,:), intent(in) :: Ats
            double complex,   allocatable, dimension(:,:,:,:) :: lightinteractionparts,dummy1s
            double complex,   allocatable, dimension(:,:,:,:)  :: deltaEparts
            real, intent(in) :: nt1s
            double precision, dimension(0:1) :: As
            integer :: nu,mu,nk1,nk2
            !real :: T2 = 5.2
            allocate(lightinteractionparts(0:Norbital,0:Norbital,0:nkamax,0:nkamax)&
                    ,dummy1s(0:Norbital,0:Norbital,0:nkamax,0:nkamax)&
                    ,deltaEparts(0:Norbital,0:Norbital,0:nkamax,0:nkamax))
            call Alecrealt(nt1s, Ats, As)
            !deltaEparts= deltaEs*rhos
            !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(nk2, nk1, mu, nu) 
            do nk2 = 0, nkamax 
                do nk1 = 0, nkamax
                    do mu = 0, Norbital
                        do nu = 0, mu
                            if (nu == mu) then
                                dummy1s(nu, nu, nk1, nk2) = ps(nu, nu, nk1, nk2, 0) * As(0) + ps(nu, nu, nk1, nk2, 1) * As(1)
                            else
                                dummy1s(nu, mu, nk1, nk2) = ps(nu, mu, nk1, nk2, 0) * As(0) + ps(nu, mu, nk1, nk2, 1) * As(1)
                                deltaEparts(nu, mu, nk1, nk2) = deltaEs(nu, mu, nk1, nk2) * rhos(nu, mu, nk1, nk2)
                                dummy1s(mu, nu, nk1, nk2) = CONJG(dummy1s(nu, mu, nk1, nk2))
                            end if
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            call commu(dummy1s,rhos,lightinteractionparts)
            !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(nk2, nk1, mu, nu) 
            do nk2 = 0, nkamax 
                do nk1 = 0, nkamax
                    do mu = 0, Norbital
                        do nu = 0, mu
                            if (nu == mu) then
                                RHSs(nu, nu, nk1, nk2) = &
                                - cmplx(0,1)/hbar* deltaEparts(nu, nu, nk1, nk2)&
                                - cmplx(0,1)/(hbar*me)* qe* lightinteractionparts(nu, nu, nk1, nk2)
                            else
                                RHSs(nu, mu, nk1, nk2) = &
                                - cmplx(0,1)/hbar* deltaEparts(nu, mu, nk1, nk2)&
                                - cmplx(0,1)/(hbar*me)* qe* lightinteractionparts(nu, mu, nk1, nk2)
                                RHSs(mu,nu,nk1,nk2) = CONJG(RHSs(nu,mu,nk1,nk2))
                            end if
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            !RHSs=-cmplx(0,1)/hbar*deltaEparts - cmplx(0,1)/(hbar*me)*qe*lightinteractionparts
            deallocate(lightinteractionparts,dummy1s,deltaEparts)
        end subroutine
        subroutine RK4(rhos, nts, ps, deltaEs, Ats)
        use globalvar
        implicit none
        double complex, dimension(0:,0:,0:,0:,0:), intent(inout) :: rhos
        double complex, dimension(0:,0:,0:,0:,0:,0:), intent(in) :: ps
        double precision, dimension(0:,0:,0:,0:,0:), intent(in) :: deltaEs
        double complex,allocatable, dimension(:,:,:,:,:,:) :: Rdummy
        integer, intent(in) :: nts
        double precision, dimension(:,:), intent(in) :: Ats
        real :: time
        integer :: j
        allocate(Rdummy(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:4,0:1))
        ! Calculate k1
        time = tmin + deltat * nts
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
        do j = 0, 1
            Rdummy(:,:,:,:,0,j) = rhos(:,:,:,:,j)
            call RHSsbeVG(Rdummy(:,:,:,:,1,j), rhos(:,:,:,:,j), time, ps, deltaEs, Ats)

            ! Calculate k2
            call RHSsbeVG(Rdummy(:,:,:,:,2,j), rhos(:,:,:,:,j) + deltat * 0.5* Rdummy(:,:,:,:,1,j), time+ 0.5*deltat,&
                    ps(:,:,:,:,:,j), deltaEs(:,:,:,:,j), Ats)

            ! Calculate k3
            call RHSsbeVG(Rdummy(:,:,:,:,3,j), rhos(:,:,:,:,j) + deltat * 0.5* Rdummy(:,:,:,:,2,j), time+ 0.5*deltat,&
                    ps(:,:,:,:,:,j), deltaEs(:,:,:,:,j), Ats)
            ! Calculate k4
            call RHSsbeVG(Rdummy(:,:,:,:,4,j), rhos(:,:,:,:,j) + deltat * Rdummy(:,:,:,:,3,j), time+ deltat,&
                    ps(:,:,:,:,:,j), deltaEs(:,:,:,:,j), Ats)
        end do
        !$OMP END PARALLEL DO
        ! Update rho
        rhos = Rdummy(:,:,:,:,0,:) + deltat *(Rdummy(:,:,:,:,1,:) + 2.0e0*Rdummy(:,:,:,:,2,:)&
                             + 2.0e0*Rdummy(:,:,:,:,3,:) + Rdummy(:,:,:,:,4,:)) / 6.0e0
        
        deallocate(Rdummy)
    end subroutine RK4
    end module def
    module calva
        implicit none
        contains
        subroutine calculatev( rhos, nts, xis, ps, Ats)
            use omp_lib
            use globalvar
            use para_MX2
            use def
            implicit none
            double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in)  :: rhos
            double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1,0:1), intent(in) :: xis, ps
            double complex, allocatable, dimension(:,:,:,:,:,:) :: velo
            double complex, allocatable, dimension(:):: N,P,I
            double precision, dimension(:,:), intent(in) :: Ats
            integer :: nu, mu, nk1, nk2, j
            double precision :: dS, deltak
            double precision, dimension(0:1) :: At
            integer, intent(in) :: nts
            real :: time
            time = tmin + nts*deltat
            call Alecrealt(time,Ats,At)
            call tb_parameters(comid)
            deltak=4.*pi/(sqrt(3.)*a*2*nkamax)
            allocate(N(0:Norbital) ,P(0:1) ,I(0:1) , velo(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1,0:1))
            P = complex(1,0)
            dS = deltak**2*sin(pi/3)
            !$OMP PARALLEL DO PRIVATE(j, nk2, nk1, mu, nu) REDUCTION(+:P, I)
            do j=0,1
                velo(:,:,:,:,j,:)=(ps(:,:,:,:,j,:)-qe*At(j))/me
                do nk2=0, nkamax
                    do nk1=0, nkamax
                        do mu = 0 ,Norbital
                            do nu = 0, Norbital
                                if (nu/=mu) then
                                P(j)=P(j)+1/(2*pi)**2*(sum(xis(nu,mu,nk1,nk2,j,:)*rhos(mu,nu,nk1,nk2,:)))*dS
                                end if
                                I(j)=I(j)+1/(2*pi)**2*(sum(velo(nu,mu,nk1,nk2,j,:)*rhos(mu,nu,nk1,nk2,:)))*dS
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            do nu=0,Norbital
                N(nu)=sum(rhos(nu,nu,:,:,:))*dS/(2*pi)**2
            end do
            write(1,'(F12.4, 2x,F12.5, 2x,F12.5, 2x,F12.5, 2x,F12.5)') time, P(0), P(1)
            write(2,'(F12.4,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5,x,F12.5)')&
                      time, N(0), N(1),N(2), N(3), N(4), N(5)
            write(3,'(F12.4, x,F12.5, x,F12.5)') time, At(0), At(1)
            write(4,'(F12.4, 2x,F12.5, 2x,F12.5, 2x,F12.5, 2x,F12.5)') time, I(0), I(1)
            deallocate(N,P,I)
        end subroutine
    end module calva
    program main
        use globalvar
        use varia
        use def
        use Elect
        use calva
        implicit none
        double complex,allocatable, dimension(:,:,:,:,:) :: rho
        integer :: nt, start, finish
        real :: rate
        allocate(rho(0:Norbital,0:Norbital,0:nkamax,0:nkamax, 0:1))
        rho(0,0,:,:,:)=1.
        rho(1,1,:,:,:)=0.
        rho(2,2,:,:,:)=0.
        rho(3,3,:,:,:)=1.
        rho(4,4,:,:,:)=0.
        rho(5,5,:,:,:)=0.
        call system_clock(start, rate)
        call El(Elec,Alec)
        print*,"begin calculate variable"
        call variable(grid , p , Xi , Delta , deltaE , E , V)
        call system_clock(finish)
        print *, "done Calculate variable. Time:", real(finish-start)/rate
        open(unit=1,FILE='OUTPUT/Interbandpolarize.txt')
        open(unit=2,FILE='OUTPUT/density.txt')
        open(unit=3,FILE='OUTPUT/Alecruntime.txt')
        open(unit=4,FILE='OUTPUT/currentdensity.txt')
        !print *, finish-start
        print*, "begin calculate"
        write(*,'(F12.3,x,F12.3,x,F12.3,x,F12.3,x,F12.3,x,F12.3)') &
                    real(SUM(rho(0,0,:,:,:))),real(SUM(rho(1,1,:,:,:))),real(SUM(rho(2,2,:,:,:))),&
                    real(SUM(rho(3,3,:,:,:))),real(SUM(rho(4,4,:,:,:))),real(SUM(rho(5,5,:,:,:)))
        do nt=1,ntmax
            call system_clock(start, rate)
            call RK4(rho, nt, p, deltaE, Alec)
            call system_clock(finish)
            if (mod(nt,int(10/deltat))==0) then
                print *, "-------------------------------------------------------------------------------------------------"
                write(*,'(A,F8.2,A,F8.2,A,F4.2)')"Real time:", real(tmin+nt*deltat)," on ",tmax,&
                                                ". Time per step ", real(finish-start)/real(rate)
                write(*,'(F12.3,x,F12.3,x,F12.3,x,F12.3,x,F12.3,x,F12.3)') &
                            real(SUM(rho(0,0,:,:,:))),real(SUM(rho(1,1,:,:,:))),real(SUM(rho(2,2,:,:,:))),&
                            real(SUM(rho(3,3,:,:,:))),real(SUM(rho(4,4,:,:,:))),real(SUM(rho(5,5,:,:,:)))
            end if
            if (mod(nt,int(0.2/deltat))==0) then
                call calculatev(rho, nt, xi, p, Alec)
            end if
        end do
        close(1)
        close(2)
        close(3)
        close(4)
        deallocate(rho)
        print*, "end calculate"
    end program
    