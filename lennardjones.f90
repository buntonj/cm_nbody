module lennardjones
!LENNARD JONES POTENTIAL MODULE
use constants
implicit none
contains

!FORCE SUBROUTINE
!At a given timestep, calculates all components of force between all particles in system
subroutine ljforce(r, F)
real*8, dimension(:,:) :: r, F ! first dimension is x, y, z, second is particle #
integer :: nbodies, ipart, jpart
real*8 :: length, fmag
nbodies = size(r,2)

do ipart = 1, nbodies
	do jpart = ipart+1, nbodies
	!DISTANCE BETWEEN IPART and JPART
	length = norm2(r(:,ipart) - r(:,jpart))
	!MAGNITUDE of FORCE from POSITIONS
	fmag = eps*(12/rm)*((rm/length)**13-(rm/length)**7)
	!COMPONENTS of FORCE on IPART from JPART
	F(:,ipart) = F(:,ipart) + fmag*(r(:,ipart)-r(:,jpart))/length
	F(:,jpart) = F(:,jpart) - F(:,ipart)
	end do
end do
end subroutine

!POTENTIAL ENERGY SUBROUTINE
!At a given timestep, calculates all potential energy contributions for all particles in system
subroutine ljpotential(r, U)
real*8, dimension(:,:) :: r ! first dimension is x, y, z, second is particle #
integer :: ipart, jpart, nbodies
real*8 :: length, U
nbodies = size(r,2)

do ipart = 1, nbodies
	do jpart = ipart+1, nbodies
	length = norm2(r(:,ipart) - r(:,jpart))
	U = U + eps*((rm/length)**12-2.0*(rm/length)**6)
	end do
end do
end subroutine


subroutine ljmomentum(v, p)
real*8, dimension(:,:) :: v
real*8, dimension(3) :: p 
integer :: nbodies, ipart
nbodies = size(v,2)

do ipart = 1,nbodies
	p(:) = p(:) + m*v(:, ipart)
end do

end subroutine


subroutine ljenergy(r, v, E)
real*8, dimension(:,:) :: r, v
integer :: nbodies, i
real*8 :: E
nbodies = size(r,2)

call ljpotential(r, E) !Stores potential of system in E

do i = 1,nbodies
	E = E + 0.5*m*norm2(v(:,i))**2 !Adds KE of each particle
end do

end subroutine

end module
