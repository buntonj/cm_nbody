module morse
!MORSE POTENTIAL MODULE
use constants
implicit none
contains

!FORCE SUBROUTINE
!At a given timestep, calculates the total force in each direction for all particles in system
subroutine morseforce(r, F)
integer :: nbodies
real*8, intent(in), dimension(:,:) :: r
real*8, intent(out), dimension(:,:) :: F ! first dimension is x, y, z, second is particle #
integer ::  ipart, jpart, i 
real*8 :: length, fmag

nbodies = size(r,2)

do ipart = 1, nbodies
	do jpart = ipart+1, nbodies
	!DISTANCE BETWEEN IPART and JPART
	length = norm2(r(:,ipart) - r(:,jpart))
	!MAGNITUDE of FORCE from POSITIONS
	fmag = -2.0*D*a*dexp(a*(re-length))*(1-dexp(a*(re-length)))
	!COMPONENTS of FORCE on IPART from JPART
	F(:,ipart) = F(:,ipart) + fmag*(r(:,ipart)-r(:,jpart))/length
	F(:,jpart) = F(:,jpart) - F(:,ipart)
	end do
end do
end subroutine

!POTENTIAL ENERGY SUBROUTINE
!At a given timestep, calculates all potential energy contributions for all particles in system
subroutine morsepotential(r, U)
real*8, dimension(:,:) :: r ! first dimension is x, y, z, second is particle #
integer :: ipart, jpart, nbodies
real*8 :: length, U
nbodies = size(r,2)
do ipart = 1, nbodies
	do jpart = ipart+1, nbodies
	length = norm2(r(:,ipart) - r(:,jpart))
	U = U + D*(1-dexp(-a*(length-re)))**2
	end do
end do
end subroutine

!MOMENTUM SUBROUTINE
!At a given timestep, sums momentum from all particles in system
subroutine morsemomentum(v, p)
real*8, dimension(:,:) :: v
real*8, dimension(3) :: p 
integer :: nbodies, ipart
nbodies = size(v,2)
do ipart = 1,nbodies
	p(:) = p(:) + m*v(:, ipart)
end do

end subroutine

!ENERGY SUBROUTINE
!At a given timestep, calculates the total energy in the system
subroutine morseenergy(r, v, E)
real*8, dimension(:,:) :: r, v
integer :: nbodies,  i
real*8 :: E
nbodies = size(v,2)

call morsepotential(r, E) !Stores potential of system in E

do i = 1,nbodies
	E = E + 0.5*m*norm2(v(:,i))**2 !Adds KE of each particle
end do


end subroutine

end module
