! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
module test_vibration
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io, only : structure_type, new
   use mctc_io_codata2018, only : electron_relative_atomic_mass, hartree_inverse_meter_relationship
   use tblite_xtb_vibration, only : get_vibrational_analysis, get_standard_atomic_mass
   implicit none
   private
   public :: collect_vibration
contains

subroutine collect_vibration(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      new_unittest("h2-stretch", test_h2_stretch), &
      new_unittest("h2-imaginary", test_h2_imaginary), &
      new_unittest("rigid-rank", test_rigid_rank), &
      new_unittest("custom-mass", test_custom_mass) &
      ]
end subroutine collect_vibration

subroutine test_h2_stretch(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp), allocatable :: freq(:), modes(:,:)
   real(wp) :: hess(6,6), expected, mass_au
   integer :: nimag, nrigid
   real(wp), parameter :: kforce = 0.2_wp
   real(wp), parameter :: xyz(3,2) = reshape([-0.7_wp,0.0_wp,0.0_wp, 0.7_wp,0.0_wp,0.0_wp],[3,2])

   call new(mol, [1,1], xyz)
   hess = 0.0_wp
   hess(1,1)=kforce; hess(1,4)=-kforce
   hess(4,1)=-kforce; hess(4,4)=kforce
   call get_vibrational_analysis(mol,hess,freq,modes,nimag,nrigid,error)
   if (allocated(error)) return
   if (nrigid /= 5 .or. size(freq) /= 1 .or. nimag /= 0) then
      call test_failed(error,"Unexpected H2 vibrational dimensions")
      return
   end if
   mass_au = get_standard_atomic_mass(1)/electron_relative_atomic_mass
   expected = sqrt(2.0_wp*kforce/mass_au)*hartree_inverse_meter_relationship*1.0e-2_wp
   call check(error, freq(1), expected, thr=1.0e-8_wp*max(1.0_wp,expected))
   if (allocated(error)) return
   call check(error, sqrt(sum(modes(:,1)**2)), 1.0_wp, thr=1.0e-12_wp)
end subroutine test_h2_stretch

subroutine test_h2_imaginary(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp), allocatable :: freq(:), modes(:,:)
   real(wp) :: hess(6,6)
   integer :: nimag, nrigid
   real(wp), parameter :: kforce = -0.2_wp
   real(wp), parameter :: xyz(3,2) = reshape([-0.7_wp,0.0_wp,0.0_wp, 0.7_wp,0.0_wp,0.0_wp],[3,2])

   call new(mol, [1,1], xyz)
   hess = 0.0_wp
   hess(1,1)=kforce; hess(1,4)=-kforce
   hess(4,1)=-kforce; hess(4,4)=kforce
   call get_vibrational_analysis(mol,hess,freq,modes,nimag,nrigid,error)
   if (allocated(error)) return
   if (nimag /= 1 .or. freq(1) >= 0.0_wp) call test_failed(error,"Imaginary H2 mode not detected")
end subroutine test_h2_imaginary

subroutine test_rigid_rank(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp), allocatable :: freq(:), modes(:,:)
   real(wp), allocatable :: hess(:,:)
   integer :: nimag, nrigid
   real(wp), parameter :: xyz(3,3) = reshape([ &
      0.0_wp,0.0_wp,0.0_wp, 1.4_wp,1.1_wp,0.0_wp, -1.4_wp,1.1_wp,0.0_wp],[3,3])

   call new(mol,[8,1,1],xyz)
   allocate(hess(9,9)); hess=0.0_wp
   call get_vibrational_analysis(mol,hess,freq,modes,nimag,nrigid,error)
   if (allocated(error)) return
   if (nrigid /= 6 .or. size(freq) /= 3 .or. size(modes,2) /= 3) then
      call test_failed(error,"Nonlinear molecule rigid-body rank is incorrect")
   end if
end subroutine test_rigid_rank

subroutine test_custom_mass(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp), allocatable :: f1(:), f2(:), m1(:,:), m2(:,:)
   real(wp) :: hess(6,6), masses(2)
   integer :: ni, nr
   real(wp), parameter :: kforce=0.2_wp
   real(wp), parameter :: xyz(3,2)=reshape([-0.7_wp,0.0_wp,0.0_wp,0.7_wp,0.0_wp,0.0_wp],[3,2])
   call new(mol,[1,1],xyz)
   hess=0.0_wp; hess(1,1)=kforce; hess(1,4)=-kforce; hess(4,1)=-kforce; hess(4,4)=kforce
   call get_vibrational_analysis(mol,hess,f1,m1,ni,nr,error)
   if (allocated(error)) return
   masses=[2.016_wp,2.016_wp]
   call get_vibrational_analysis(mol,hess,f2,m2,ni,nr,error,masses=masses)
   if (allocated(error)) return
   if (abs(f2(1)/f1(1)-sqrt(0.5_wp)) > 2.0e-3_wp) then
      call test_failed(error,"Custom isotope masses do not rescale frequency correctly")
   end if
end subroutine test_custom_mass

end module test_vibration
