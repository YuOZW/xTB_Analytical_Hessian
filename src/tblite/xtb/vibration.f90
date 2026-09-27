! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Vibrational analysis utilities for analytic GFN1-xTB Hessians.
module tblite_xtb_vibration
   use mctc_env, only : wp, dp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_codata2018, only : electron_relative_atomic_mass, &
      & hartree_inverse_meter_relationship
   implicit none
   private

   public :: get_vibrational_analysis, get_standard_atomic_mass

   ! Standard atomic weights (u) for Z=1..86. For elements without a standard
   ! atomic weight, the mass number of a representative long-lived isotope is
   ! used. User supplied masses can override these values.
   real(wp), parameter :: standard_mass(86) = [ &
      1.008_wp, 4.002602_wp, 6.94_wp, 9.0121831_wp, 10.81_wp, 12.011_wp, &
      14.007_wp, 15.999_wp, 18.998403163_wp, 20.1797_wp, 22.98976928_wp, &
      24.305_wp, 26.9815385_wp, 28.085_wp, 30.973761998_wp, 32.06_wp, &
      35.45_wp, 39.948_wp, 39.0983_wp, 40.078_wp, 44.955908_wp, 47.867_wp, &
      50.9415_wp, 51.9961_wp, 54.938044_wp, 55.845_wp, 58.933194_wp, &
      58.6934_wp, 63.546_wp, 65.38_wp, 69.723_wp, 72.630_wp, 74.921595_wp, &
      78.971_wp, 79.904_wp, 83.798_wp, 85.4678_wp, 87.62_wp, 88.90584_wp, &
      91.224_wp, 92.90637_wp, 95.95_wp, 98.0_wp, 101.07_wp, 102.90550_wp, &
      106.42_wp, 107.8682_wp, 112.414_wp, 114.818_wp, 118.710_wp, 121.760_wp, &
      127.60_wp, 126.90447_wp, 131.293_wp, 132.90545196_wp, 137.327_wp, &
      138.90547_wp, 140.116_wp, 140.90766_wp, 144.242_wp, 145.0_wp, &
      150.36_wp, 151.964_wp, 157.25_wp, 158.92535_wp, 162.500_wp, &
      164.93033_wp, 167.259_wp, 168.93422_wp, 173.045_wp, 174.9668_wp, &
      178.49_wp, 180.94788_wp, 183.84_wp, 186.207_wp, 190.23_wp, 192.217_wp, &
      195.084_wp, 196.966569_wp, 200.592_wp, 204.38_wp, 207.2_wp, &
      208.98040_wp, 209.0_wp, 210.0_wp, 222.0_wp ]

   interface
      pure subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsyev
   end interface

contains

pure function get_standard_atomic_mass(number) result(mass)
   integer, intent(in) :: number
   real(wp) :: mass
   if (number >= 1 .and. number <= size(standard_mass)) then
      mass = standard_mass(number)
   else
      mass = 0.0_wp
   end if
end function get_standard_atomic_mass


!> Convert a Cartesian Hessian to harmonic vibrational frequencies and normal
!> modes. Translational and rotational vectors are projected in mass-weighted
!> coordinates. Frequencies are returned in cm^-1; imaginary modes are
!> represented by negative frequencies. Normal modes are Cartesian displacement
!> vectors (3*nat,nvib), each normalized to unit Euclidean norm.
subroutine get_vibrational_analysis(mol, hessian, frequencies, modes, n_imag, &
      & n_rigid, error, masses, imag_cutoff)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: hessian(:, :)
   real(wp), allocatable, intent(out) :: frequencies(:)
   real(wp), allocatable, intent(out) :: modes(:, :)
   integer, intent(out) :: n_imag
   integer, intent(out) :: n_rigid
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(in), optional :: masses(:)
   real(wp), intent(in), optional :: imag_cutoff

   integer :: nat, ndim, iat, ic, j, k, nvib, info, lwork, nq
   real(wp) :: cutoff, scale, normv
   real(wp), allocatable :: mass_u(:), mass_au(:), mw_hess(:, :), rigid(:, :)
   real(wp), allocatable :: basis(:, :), vib_basis(:, :), red_hess(:, :)
   real(wp), allocatable :: eval(:), work(:), vec(:)

   nat = mol%nat
   ndim = 3*nat
   n_imag = 0
   n_rigid = 0
   cutoff = 20.0_wp
   if (present(imag_cutoff)) cutoff = max(0.0_wp, imag_cutoff)

   if (size(hessian,1) /= ndim .or. size(hessian,2) /= ndim) then
      call fatal_error(error, "Invalid Hessian dimensions for vibrational analysis")
      return
   end if
   if (any(mol%periodic)) then
      call fatal_error(error, "Vibrational analysis currently supports non-periodic systems only")
      return
   end if

   allocate(mass_u(nat), mass_au(nat))
   if (present(masses)) then
      if (size(masses) /= nat .or. any(masses <= 0.0_wp)) then
         call fatal_error(error, "Invalid atomic masses for vibrational analysis")
         return
      end if
      mass_u = masses
   else
      do iat = 1, nat
         mass_u(iat) = get_standard_atomic_mass(mol%num(mol%id(iat)))
         if (mass_u(iat) <= 0.0_wp) then
            call fatal_error(error, "No standard atomic mass available for vibrational analysis")
            return
         end if
      end do
   end if

   ! 1 u / m_e. Atomic-unit Hessians require nuclear masses in electron masses.
   mass_au = mass_u/electron_relative_atomic_mass

   allocate(rigid(ndim,6))
   call build_rigid_basis(mol, mass_au, rigid, n_rigid)
   nvib = ndim - n_rigid
   if (nvib < 0) then
      call fatal_error(error, "Invalid rigid-body rank in vibrational analysis")
      return
   end if
   if (nvib == 0) then
      allocate(frequencies(0), modes(ndim,0))
      return
   end if

   allocate(mw_hess(ndim,ndim))
   do iat = 1, nat
      do ic = 1, 3
         j = 3*(iat-1)+ic
         do k = 1, ndim
            mw_hess(j,k) = hessian(j,k)/sqrt(mass_au(iat)*mass_au((k-1)/3+1))
         end do
      end do
   end do
   mw_hess = 0.5_wp*(mw_hess + transpose(mw_hess))

   ! Construct an explicit orthonormal complement to the rigid-body space.
   ! This is robust even when genuine vibrational eigenvalues are exactly zero.
   allocate(basis(ndim,ndim), vec(ndim))
   basis = 0.0_wp
   if (n_rigid > 0) basis(:,:n_rigid) = rigid(:,:n_rigid)
   nq = n_rigid
   do j = 1, ndim
      if (nq == ndim) exit
      vec = 0.0_wp
      vec(j) = 1.0_wp
      call append_orthonormal(vec, basis, nq, 1.0e-10_wp)
   end do
   if (nq /= ndim) then
      call fatal_error(error, "Failed to construct vibrational subspace")
      return
   end if
   allocate(vib_basis(ndim,nvib), red_hess(nvib,nvib))
   vib_basis = basis(:,n_rigid+1:ndim)
   red_hess = matmul(transpose(vib_basis), matmul(mw_hess, vib_basis))
   red_hess = 0.5_wp*(red_hess + transpose(red_hess))

   allocate(eval(nvib))
   lwork = max(1, 3*nvib-1)
   allocate(work(lwork))
   call dsyev('V', 'U', nvib, red_hess, nvib, eval, work, lwork, info)
   if (info /= 0) then
      call fatal_error(error, "LAPACK diagonalization failed in vibrational analysis")
      return
   end if

   allocate(frequencies(nvib), modes(ndim,nvib))
   scale = hartree_inverse_meter_relationship*1.0e-2_wp
   do k = 1, nvib
      if (eval(k) < 0.0_wp) then
         frequencies(k) = -sqrt(abs(eval(k)))*scale
      else
         frequencies(k) = sqrt(eval(k))*scale
      end if

      vec = matmul(vib_basis, red_hess(:,k))
      do iat = 1, nat
         vec(3*(iat-1)+1:3*iat) = vec(3*(iat-1)+1:3*iat)/sqrt(mass_au(iat))
      end do
      normv = sqrt(sum(vec*vec))
      if (normv > 0.0_wp) vec = vec/normv
      modes(:,k) = vec
   end do

   n_imag = count(frequencies < -cutoff)
end subroutine get_vibrational_analysis


subroutine build_rigid_basis(mol, mass_au, rigid, n_rigid)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: mass_au(:)
   real(wp), intent(out) :: rigid(:, :)
   integer, intent(out) :: n_rigid

   integer :: nat, ndim, iat, iaxis
   real(wp) :: mtot, normv
   real(wp) :: com(3), r(3), v3(3)
   real(wp), allocatable :: vec(:)
   real(wp), parameter :: tol = 1.0e-10_wp

   nat = mol%nat
   ndim = 3*nat
   rigid = 0.0_wp
   n_rigid = 0
   allocate(vec(ndim))

   mtot = sum(mass_au)
   com = 0.0_wp
   do iat = 1, nat
      com = com + mass_au(iat)*mol%xyz(:,iat)
   end do
   com = com/mtot

   ! Translations.
   do iaxis = 1, 3
      vec = 0.0_wp
      do iat = 1, nat
         vec(3*(iat-1)+iaxis) = sqrt(mass_au(iat))
      end do
      call append_orthonormal(vec, rigid, n_rigid, tol)
   end do

   ! Rotations about Cartesian axes around the center of mass.
   do iaxis = 1, 3
      vec = 0.0_wp
      do iat = 1, nat
         r = mol%xyz(:,iat)-com
         select case(iaxis)
         case(1)
            v3 = [0.0_wp, -r(3), r(2)]
         case(2)
            v3 = [r(3), 0.0_wp, -r(1)]
         case default
            v3 = [-r(2), r(1), 0.0_wp]
         end select
         vec(3*(iat-1)+1:3*iat) = sqrt(mass_au(iat))*v3
      end do
      call append_orthonormal(vec, rigid, n_rigid, tol)
   end do
end subroutine build_rigid_basis


subroutine append_orthonormal(vec, basis, nvec, tol)
   real(wp), intent(inout) :: vec(:)
   real(wp), intent(inout) :: basis(:, :)
   integer, intent(inout) :: nvec
   real(wp), intent(in) :: tol

   integer :: j
   real(wp) :: normv

   ! Reorthogonalized modified Gram-Schmidt is sufficient for at most six
   ! rigid-body vectors and avoids a dependency on an additional QR wrapper.
   if (nvec > 0) then
      do j = 1, nvec
         vec = vec - dot_product(basis(:,j),vec)*basis(:,j)
      end do
      do j = 1, nvec
         vec = vec - dot_product(basis(:,j),vec)*basis(:,j)
      end do
   end if
   normv = sqrt(sum(vec*vec))
   if (normv > tol) then
      nvec = nvec + 1
      basis(:,nvec) = vec/normv
   end if
end subroutine append_orthonormal

end module tblite_xtb_vibration
