
!global variable
module globalvar
!contains
!subroutine gbv()
integer, parameter :: Norbital=5, nkamax=199, comid=0
integer, parameter :: halfNorbital=2, qe=1
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
real, parameter :: E0 = 1e-2 
! V/nm
real, parameter :: phi = 0.0 
! rad
real, parameter :: TauL = 60*sqrt(2.0) 
! fs (the TauL is not have the form x^2/2\tauL**2)
real, parameter :: epsilonL = 1.81 
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
!eV.fs; unit charge; eV.fs^2/nm^2
!watch-out all array start from 0 !!!!!!!!
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
	double precision :: kx,ky,alpha,beta
	!call gbv()
	double complex, dimension(0:2,0:2) :: ham0,ham1,hamu,hamd
	call tb_parameters(comid)
	alpha=kx*a/2.
	beta=sqrt(3.)*ky*a/2.
	ham0(0,0)=e1+2.*t0*(2*cos(alpha)*cos(beta)+cos(2*alpha))&
				+2.*r0*(2*cos(3*alpha)*cos(beta)+cos(2*beta))&
				+2.*u0*(2*cos(2*alpha)*cos(2*beta)+cos(4*alpha))
	ham0(0,1)=cmplx(-2.*sqrt(3.)*t2*sin(alpha)*sin(beta)&
					+2.*(r1+r2)*sin(3*alpha)*sin(beta)&
					-2.*sqrt(3.)*u2*sin(2*alpha)*sin(2*beta),2.*t1*sin(alpha)*(2*cos(alpha)+cos(beta))&
					+2.*(r1-r2)*sin(3*alpha)*cos(beta)&
					+2.*u1*sin(2*alpha)*(2*cos(2*alpha)+cos(2*beta)))
	ham0(0,2)=cmplx(2.*t2*(cos(2*alpha)-cos(alpha)*cos(beta))&
					-2.*(r1+r2)*(cos(3*alpha)*cos(beta)-cos(2*beta))/sqrt(3.)&
					+2.*u2*(cos(4*alpha)-cos(2*alpha)*cos(2*beta)),2.*sqrt(3.)*t1*cos(alpha)*sin(beta)&
					+2.*(r1-r2)*sin(beta)*(cos(3*alpha)+2*cos(beta))/sqrt(3.)&
					+2.*sqrt(3.)*u1*cos(2*alpha)*sin(2*beta))
	ham0(1,0)=conjg(ham0(0,1))
	ham0(1,1)=e2+(t11+3.*t22)*cos(alpha)*cos(beta)+2.*t11*cos(2.*alpha)&
				+4.*r11*cos(3.*alpha)*cos(beta)+2.*(r11+sqrt(3.)*r12)*cos(2.*beta)&
				+(u11+3.*u22)*cos(2.*alpha)*cos(2.*beta)+2.*u11*cos(4.*alpha)
	ham0(1,2)=cmplx(sqrt(3.)*(t22-t11)*sin(alpha)*sin(beta)&
					+4.*r12*sin(3.*alpha)*sin(beta)&
					+sqrt(3.)*(u22-u11)*sin(2.*alpha)*sin(2.*beta), 4*t12*sin(alpha)*(cos(alpha)-cos(beta))&
					+4.*u12*sin(2.*alpha)*(cos(2.*alpha)-cos(2.*beta)))
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
end subroutine
subroutine gradient_H(kx, ky, grad_H)
use globalvar
    implicit none
    double precision :: kx, ky
	double complex, dimension(3,3) :: hamu_plus, hamu_minus, hamd_plus, hamd_minus
    double complex, dimension(6,6,2) :: grad_H
    real :: delta = 1.0e-2

    call tbham(kx + delta, ky, hamu_plus,hamd_plus)
    call tbham(kx - delta, ky, hamu_minus,hamd_minus)
    grad_H(1:3,1:3,1) = (hamu_plus - hamu_minus) / (2 * delta)
	grad_H(4:6,4:6,1) = (hamd_plus - hamd_minus) / (2 * delta)

    call tbham(kx, ky + delta, hamu_plus,hamd_plus)
    call tbham(kx, ky - delta, hamu_minus,hamd_minus)
    grad_H(1:3,1:3,2) = (hamu_plus - hamu_minus) / (2 * delta)
	grad_H(4:6,4:6,2) = (hamd_plus - hamd_minus) / (2 * delta)
end subroutine
end module prepare
module varia
	use omp_lib
	use globalvar
    use prepare
	use para_MX2
	double precision,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax) :: deltaE
	double precision,dimension(0:Norbital,0:nkamax,0:nkamax) :: E
    double complex,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax) :: V
    double complex,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1) :: p,Xi
    double precision,dimension(0:nkamax,0:nkamax,0:1) :: grid
    integer, dimension(0:Norbital,0:Norbital) :: Delta
    save
	contains
    subroutine variable(grid,p,Xi,Delta,deltaE,E,V)
    	use prepare
    	use omp_lib
	use para_MX2
	use globalvar
        implicit none
        double complex,dimension(0:Norbital,0:Norbital) :: Vpkx, Vpky ,ham, ham_plus, ham_minus
		double complex,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax),intent(inout) :: V
        double complex,dimension(0:Norbital,0:Norbital,0:1) :: Vpkxy
        double complex,dimension(0:2,0:2) :: hamu,hamd
        double precision,dimension(0:2) :: E_sorted
	double precision, dimension(0:1,0:1) :: B
	double precision,dimension(0:Norbital,0:nkamax,0:nkamax),intent(inout) :: E
	double precision,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax),intent(inout) :: deltaE
	double complex,dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1),intent(inout) :: p,Xi
	double precision,dimension(0:nkamax,0:nkamax,0:1),intent(inout) :: grid
	integer, dimension(0:Norbital,0:Norbital),intent(inout) :: Delta
        integer :: nk1, nk2, lwork, i, j, nu, mu,INFO 	
        double precision , dimension(0:1) :: array
	double complex :: WORK(5)
	real ::  da = 1.0e-5
	double precision :: kx,ky,RWORK(7)
        external :: zheev
	call tb_parameters(1)
        lwork = max(1, 5)
        ! Initialize nk1, nk2, and nkamax here
        B(0,0) = 2.*pi/a
		B(0,1) = 2.*pi/a
		B(1,0) = 2.*pi/(sqrt(3.)*a)
		B(1,1) = -2.*pi/(sqrt(3.)*a)
		print *, "begin variable"
		print *, "grid begin"
		do nk2 = 0,nkamax
			do nk1 =0,nkamax
				array(0) = real(nk1)/real(nkamax)
				array(1) = real(nk2)/real(nkamax)
				grid(nk1,nk2,:)=matmul(B,array)
				grid(nk1,nk2,0)=grid(nk1,nk2,0)-2*pi/a
				grid(nk1,nk2,1)=grid(nk1,nk2,1)-0
			end do
		end do
		open(unit=1,FILE='grid.txt')
		do nk2 = 0,nkamax
			do nk1 =0,nkamax
			write(1,'(E12.5, 2x, E12.5, 2x, E12.5, 2x, E12.5)') real(nk1), real(nk2), grid(nk1,nk2,0), grid(nk1,nk2,1)
			end do
		end do
		close(1)
		print *, "grid done"
		print *,"Calculate p,Xi,E,V,.."
        do nk2 = 0,nkamax
            do nk1 =0,nkamax
!                kx=grid(nk1,nk2,0)
!                ky=grid(nk1,nk2,1) 
!                call tbham(kx,ky,hamu,hamd)
!				!print*, "ham okay"
                call zheev('V','U',3,hamu,3,E_sorted,WORK,5,RWORK,INFO)
				!print*, "ZHEEV okay"
!                !sort eigenvalue and eigenvector
                E(0:2,nk1,nk2)=E_sorted
                V(0:2,0:2,nk1,nk2)=hamu
                call zheev('V','U',3,hamd,3,E_sorted,WORK,5,RWORK,INFO)
!				!final index is 0 mean eigenvalue is ascending order
!                !sort eigenvalue and eigenvector
                E(3:5,nk1,nk2)=E_sorted
                V(3:5,3:5,nk1,nk2)=hamd
                call gradient_H(kx,ky,Vpkxy)
                Vpkx=Vpkxy(:,:,0)
                Vpky=Vpkxy(:,:,1)
                do mu = 0, Norbital
                    do nu = 0, Norbital
                    p(nu,mu,nk1,nk2,0)=me/hbar*DOT_PRODUCT(V(:,nu,nk1,nk2),matmul(Vpkx, V(:, mu, nk1, nk2)))
                    p(nu,mu,nk1,nk2,1)=me/hbar*DOT_PRODUCT(V(:,nu,nk1,nk2),matmul(Vpky, V(:, mu, nk1, nk2)))
                    deltaE(nu,mu,nk1,nk2)=E(nu,nk1,nk2)-E(mu,nk1,nk2)
                    if (nu /= mu) then
                        xi(nu,mu,nk1,nk2,0)=hbar/me*p(nu,mu,nk1,nk2,0)/deltaE(nu,mu,nk1,nk2)
                        xi(nu,mu,nk1,nk2,1)=hbar/me*p(nu,mu,nk1,nk2,1)/deltaE(nu,mu,nk1,nk2)
                    else
                        xi(nu,mu,nk1,nk2,0)=0.
                        xi(nu,mu,nk1,nk2,1)=0.
                    end if
                        if (abs(deltaE(nu, mu, nk1, nk2)) < 1e-2) then
                            Delta(nu, mu) = Delta(nu, mu) + 1
                        end if
                    end do
                end do
            end do  !nk2 do
        end do !nk1 do
        do nu = 1, Norbital
            do mu = 1, Norbital
                if (Delta(nu, mu) /= 0) then
                    Delta(nu, mu) = 1
                end if
            end do
        end do
		open(unit=1,FILE='band1.txt')
			do i=1,nk1
				do j=1,nk2
					write(1,'(8E12.5)')grid(i,j,0),grid(i,j,1),&
									E(0,i,j),E(1,i,j),E(2,i,j),&
									E(3,i,j),E(4,i,j),E(5,i,j)
				end do
			end do
		close(1)
		
		open(unit=1,FILE='px.txt')
		do i=0,nkamax
		do j=0,nkamax
			if (i==j) then
			write(1,'(E12.5, x, E12.5, x, E12.5, x, E12.5, x, E12.5,x)')& 
			grid(i,j,0),grid(i,j,1),real(p(0,0,i,j,0)),real(p(0,1,i,j,0)),real(p(0,2,i,j,0))
			end if
			end do
		end do
		close(1)
    end subroutine
end module
module defud
use globalvar
contains
subroutine EAlecntinteger(Elec, Alec, nt)
    implicit none
    integer, intent(in) :: nt
    double precision, dimension(2), intent(inout) :: Elec, Alec
    double precision :: time

    ! Calculate time
    time = nt*deltat + tmin

    ! Calculate Elec
    Elec(1) = E0*exp(-(time)**2/TauL**2)*cos(omegaL*time)*cos(phi)
    Elec(2) = E0*exp(-(time)**2/TauL**2)*cos(omegaL*time)*sin(phi)

    ! Calculate Alec
    if (nt /= 0) then
        Alec = Alec - Elec*deltat
    end if
end subroutine
subroutine Elect(nt1, Et)
    implicit none
    real, intent(in) :: nt1
    double precision, dimension(2), intent(inout) :: Et
	double precision, dimension(2) :: Elec, Alec
    call EAlecntinteger(Elec, Alec, int(nt1))
    Et = Elec + (Elec - Elec) * (nt1 - int(nt1))
end subroutine Elect

subroutine Alect(nt1, At)
    implicit none
    real, intent(in) :: nt1
    double precision, dimension(2), intent(inout) :: At
	double precision, dimension(2) :: Elec, Alec
    call EAlecntinteger(Elec, Alec, int(nt1))
    At = Alec - Elec * (nt1 - int(nt1)) * deltat
end subroutine Alect
subroutine commu(A, B, C)
    use omp_lib
    implicit none
    double complex, dimension(:,:,:,:), intent(in) :: A, B
    double complex, dimension(:,:,:,:), intent(out) :: C
    integer :: m, n

    ! Initialize C to zero
    C = 0.0

    ! Loop over the last two indices
    !!$OMP PARALLEL DO PRIVATE(m)
    do m = lbound(A,3), ubound(A,3)
    	!!$OMP PARALLEL DO PRIVATE(m)
        do n = lbound(A,4), ubound(A,4)
            ! Calculate the commutator using matmul
            C(:,:,m,n) = matmul(A(:,:,m,n), B(:,:,m,n)) - matmul(B(:,:,m,n), A(:,:,m,n))
        end do
        !!$OMP END PARALLEL DO
    end do
    !!$OMP END PARALLEL DO
end subroutine commu
subroutine RHSsbeVG(RHS, rho,nt1,p,deltaE)
	use omp_lib
    	implicit none
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: rho
	double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: deltaE
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(out) :: RHS
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: p
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax) :: lightinteractionpart,dummy1
	double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax)  :: deltaEpart
	integer :: nu,mu,nk1,nk2
	real, intent(in) :: nt1
	double precision, dimension(0:1) :: At
	real :: T2 = 5.2
	call Alect(nt1, At)
	!!$OMP PARALLEL DO PRIVATE(nk1)
	do nk1 = 0, nkamax
		!!$OMP PARALLEL DO PRIVATE(nk2,mu,nu)
		do nk2 = 0, nkamax
			do mu = 0, Norbital
				do nu = 0, Norbital
					deltaEpart(nu, mu, nk1, nk2)= deltaE(nu, mu, nk1, nk2)*rho(nu, mu, nk1, nk2)
					dummy1(nu, mu, nk1, nk2) = p(nu, mu, nk1, nk2, 0) * At(0) + p(nu, mu, nk1, nk2, 1) * At(1)
				end do
			end do
		end do
		!!$OMP END PARALLEL DO
	end do
	!!$OMP END PARALLEL DO
	call commu(dummy1,rho,lightinteractionpart)
	RHS=-cmplx(0,1)/hbar*deltaEpart - cmplx(0,1)/hbar*qe/me*lightinteractionpart
end subroutine
subroutine RK4(rho, nt, p, deltaE)
    implicit none
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(inout) :: rho
    double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in) :: deltaE
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax):: RHS
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: p
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax) :: k1, k2, k3, k4
    integer, intent(in) :: nt

    ! Calculate k1
    call RHSsbeVG(RHS, rho, real(nt), p, deltaE)
    k1 = deltat * RHS

    ! Calculate k2
    call RHSsbeVG(RHS, rho + 0.5d0*k1, real(nt) + 0.5e0, p, deltaE)
    k2 = deltat * RHS

    ! Calculate k3
    call RHSsbeVG(RHS, rho + 0.5d0*k2, real(nt) + 0.5e0, p, deltaE)
    k3 = deltat * RHS

    ! Calculate k4
    call RHSsbeVG(RHS, rho + k3, real(nt) + 1.0e0, p, deltaE)
    k4 = deltat * RHS

    ! Update rho
    rho = rho + (k1 + 2.0e0*k2 + 2.0e0*k3 + k4) / 6.0e0
end subroutine RK4
subroutine calcuRK4(x,y)
	use varia
    	implicit none
    	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: x
    	double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in) :: y
	double complex, dimension(0:1) :: I
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax) :: rho
	integer :: nt
	Real :: start, finish
	!call variable(grid, p, Xi, Delta, deltaE, E, V)
	rho(0,0,:,:)=1
	rho(1,1,:,:)=0
	rho(2,2,:,:)=0
	rho(3,3,:,:)=1
	rho(4,4,:,:)=0
	rho(5,5,:,:)=0
	I=0
	!write(*,*)'.....'
	do nt=1,ntmax
		call cpu_time(start)
		call RK4(rho, nt, x, y)
		!write(*,*)'.....'
		call cpu_time(finish)
		if (mod(nt,int(10/deltat))==0) then
			print *, "-------------------------------------------------------------------------------------------------"
			print *,"Real time:", tmin+nt*deltat," on ",tmax, ".Time per step ", finish-start 
			print*, REAL(SUM(rho(0,0,:,:))), REAL(SUM(rho(1,1,:,:))), REAL(SUM(rho(2,2,:,:)))
			print*,	REAL(SUM(rho(3,3,:,:))), REAL(SUM(rho(4,4,:,:))), REAL(SUM(rho(5,5,:,:)))
		end if 
		!call calculatev(I,rho,nt,p)
	end do 

end subroutine
end module

module calculatevariable
	use omp_lib
	use para_MX2
	use globalvar
	use varia
	use defud
	implicit none
	contains
	subroutine Calculatev(I,rho,nt1,p)
	double complex, dimension(0:1), intent(inout) :: I
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: rho
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: p
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1):: velocity
	integer :: j,nu,mu,nk1,nk2
	double precision :: dS
	double precision, dimension(0:1) :: At
	real, intent(in) :: nt1
	dS =(4.*pi/(sqrt(3.)*a*nkamax))**2*sqrt(3.)/2
	call Alect(nt1,At)
	!!$OMP PARALLEL DO
	do j=0,1
		velocity(:,:,:,:,j)=(p(:,:,:,:,j)-qe*At(j))/me
	end do
	!!$OMP END PARALLEL DO
	open(unit=1,FILE='OUTPUT/Currentdensity.txt')
	!!$OMP PARALLEL DO PRIVATE(j,nu,mu,nk1,nk2)
	do j=0,1
	do nu=0,Norbital
		do mu=0,Norbital
			do nk1=0,nkamax
				do nk2=0,nkamax
					I(j)=I(j)+qe*velocity(nu,mu,nk1,nk2,j)*rho(mu,nu,nk1,nk2)*dS
				end do
			end do
		end do
	end do
	end do
	!!$OMP END PARALLEL DO
	
	close(1)	
	end subroutine
end module

program main
	use omp_lib
	use globalvar	
    	use varia
	use defud
	use calculatevariable
	use prepare
    	implicit none
    	double complex, dimension(0:1) :: I
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax) :: rho
	integer :: nt
	Real :: start, finish
	call cpu_time(start)
	call variable(grid, p, Xi, Delta, deltaE, E, V)
	call calcuRK4(p,deltaE)
	call cpu_time(finish)
	print *, finish-start
	
    ! Print or use the arrays here...
end program
