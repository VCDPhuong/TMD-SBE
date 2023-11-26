module main
integer :: Norbital=2
integer :: nkamax=10
contains
subroutine RHSsbeVG1(RHS, rho,nt1,p,deltaE, Elec, Alec)
    use omp_lib
    implicit none
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: rho
    double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: deltaE
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(out) :: RHS
    double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: p
    double complex,allocatable,dimension(:,:,:,:) :: lightinteractionpart,dummy1
    double precision,allocatable, dimension(:,:,:,:)  :: deltaEpart
    double precision, dimension(0:ntmax,0:1) :: Elec, Alec
    integer :: nu,mu,nk1,nk2
    real, intent(in) :: nt1
    double precision, dimension(0:1) :: At
    real :: T2 = 5.2
	allocate(lightinteractionpart(0:Norbital,0:Norbital,0:nkamax,0:nkamax),dummy1(0:Norbital,0:Norbital,0:nkamax,0:nkamax)&
			,deltaEpart(0:Norbital,0:Norbital,0:nkamax,0:nkamax))
	call Alect(nt1, At,Elec, Alec)
	!$OMP PARALLEL DO PRIVATE
	do nk2 = 0, nkamax
		!$OMP PARALLEL DO PRIVATE
		do nk1 = 0, nkamax
			do mu = 0, Norbital
				do nu = 0, Norbital
					deltaEpart(nu, mu, nk1, nk2)= deltaE(nu, mu, nk1, nk2)*rho(nu, mu, nk1, nk2)
					dummy1(nu, mu, nk1, nk2) = p(nu, mu, nk1, nk2, 0) * At(0) + p(nu, mu, nk1, nk2, 1) * At(1)
				end do
			end do
		end do
		!$OMP END PARALLEL DO
	end do
	!$OMP END PARALLEL DO
	call commu(dummy1,rho,lightinteractionpart)
	RHS=-cmplx(0,1)/hbar*deltaEpart - cmplx(0,1)/hbar*qe/me*lightinteractionpart
	deallocate(lightinteractionpart,dummy1,deltaEpart)
end subroutine
subroutine RHSsbeVG2(RHS, rho,nt1,p,deltaE, Elec, Alec)
	use omp_lib
    	implicit none
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: rho
	double precision, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(in)  :: deltaE
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax), intent(out) :: RHS
	double complex, dimension(0:Norbital,0:Norbital,0:nkamax,0:nkamax,0:1), intent(in) :: p
	double complex,allocatable,dimension(:,:,:,:) :: lightinteractionpart,dummy1
	double precision,allocatable, dimension(:,:,:,:)  :: deltaEpart
	double precision, dimension(0:ntmax,0:1) :: Elec, Alec
	integer :: nu,mu,nk1,nk2
	real, intent(in) :: nt1
	double precision, dimension(0:1) :: At
	real :: T2 = 5.2
	allocate(lightinteractionpart(0:Norbital,0:Norbital,0:nkamax,0:nkamax),dummy1(0:Norbital,0:Norbital,0:nkamax,0:nkamax)&
			,deltaEpart(0:Norbital,0:Norbital,0:nkamax,0:nkamax))
	call Alect(nt1, At,Elec, Alec)
	do nk2 = 0, nkamax
		do nk1 = 0, nkamax
			do mu = 0, Norbital
				do nu = 0, Norbital
					deltaEpart(nu, mu, nk1, nk2)= deltaE(nu, mu, nk1, nk2)*rho(nu, mu, nk1, nk2)
					dummy1(nu, mu, nk1, nk2) = p(nu, mu, nk1, nk2, 0) * At(0) + p(nu, mu, nk1, nk2, 1) * At(1)
				end do
			end do
		end do
	end do
	call commu(dummy1,rho,lightinteractionpart)
	RHS=-cmplx(0,1)/hbar*deltaEpart - cmplx(0,1)/hbar*qe/me*lightinteractionpart
	deallocate(lightinteractionpart,dummy1,deltaEpart)
end subroutine
program parralel_test
use main
double precision, dimension(5,5,10,20) :: A1, A2, C3,C4

    do i = 1, 5
        do j = 1, 5
            do k = 1, 10
                do l = 1, 20
                    call random_seed()
                    call random_number(A1(i,j,k,l))
                    call random_seed()
                    call random_number(A2(i,j,k,l))
                end do
            end do
        end do
    end do
call commu1(A1,A2,C3)
call commu2(A1,A2,C4)
print *, (C4)/C3


end program
