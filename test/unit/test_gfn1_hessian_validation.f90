! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Broad validation suite for the analytic GFN1-xTB Cartesian Hessian.
module test_gfn1_hessian_validation
   use mctc_env, only : wp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container, only : container_type
   use tblite_context_type, only : context_type
   use tblite_data_spin, only : get_spin_constant
   use tblite_lapack_solver, only : lapack_solver, lapack_algorithm
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_hessian, only : get_gfn1_hessian
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_gfn1_hessian_validation

   real(wp), parameter :: racc = 2.0e-4_wp
   real(wp), parameter :: step = 2.0e-4_wp

contains

subroutine collect_gfn1_hessian_validation(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("diatomics", test_diatomics), &
      new_unittest("main-group", test_main_group), &
      new_unittest("charged", test_charged), &
      new_unittest("open-shell", test_open_shell), &
      new_unittest("frontier-degeneracy-guard", test_frontier_degeneracy_guard), &
      new_unittest("diverse-mb01", test_diverse_mb01) &
      ]
end subroutine collect_gfn1_hessian_validation


subroutine add_spin_polarization(calc, mol)
   type(xtb_calculator), intent(inout) :: calc
   type(structure_type), intent(in) :: mol
   class(container_type), allocatable :: cont
   type(spin_polarization), allocatable :: spinpol
   real(wp), allocatable :: wll(:, :, :)
   integer :: izp, ish, jsh, il, jl

   allocate(wll(calc%bas%nsh, calc%bas%nsh, mol%nid), source=0.0_wp)
   do izp = 1, mol%nid
      do ish = 1, calc%bas%nsh_id(izp)
         il = calc%bas%cgto(ish,izp)%ang
         do jsh = 1, calc%bas%nsh_id(izp)
            jl = calc%bas%cgto(jsh,izp)%ang
            wll(jsh,ish,izp) = get_spin_constant(jl,il,mol%num(izp))
         end do
      end do
   end do
   allocate(spinpol)
   call new_spin_polarization(spinpol, mol, wll, calc%bas%nsh_id)
   call move_alloc(spinpol, cont)
   call calc%push_back(cont)
end subroutine add_spin_polarization


subroutine validate_hessian(mol, nspin, etemp, columns, hthr, error, with_spin)
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: nspin
   real(wp), intent(in) :: etemp, hthr
   integer, intent(in) :: columns(:)
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in), optional :: with_spin

   integer :: k, col, jat, jc, row, iat, ic, iaxis
   logical :: spin_on
   real(wp) :: energy, errmax, asym, tsum, trans_err, rot_err
   real(wp), allocatable :: hessian(:, :), g0(:, :), gp(:, :), gm(:, :), sigma(:, :), refcol(:), rvec(:), hrot(:), grot(:)
   type(context_type) :: ctx
   type(structure_type) :: molp, molm
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfnp, wfnm

   spin_on = .false.
   if (present(with_spin)) spin_on = with_spin

   ctx%solver = lapack_solver(lapack_algorithm%gvr)
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   if (spin_on) call add_spin_polarization(calc, mol)

   allocate(hessian(3*mol%nat, 3*mol%nat), g0(3,mol%nat), gp(3,mol%nat), gm(3,mol%nat), &
      & sigma(3,3), refcol(3*mol%nat))
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, etemp)
   g0 = 0.0_wp
   sigma = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, racc, energy, gradient=g0, sigma=sigma, verbosity=0)
   if (ctx%failed()) then
      call test_failed(error, "Central GFN1 validation SCF failed")
      return
   end if

   call get_gfn1_hessian(mol, calc, wfn, racc, hessian, error)
   if (allocated(error)) return

   if (.not.all(ieee_is_finite(hessian))) then
      call test_failed(error, "GFN1 validation Hessian contains non-finite values")
      return
   end if

   asym = maxval(abs(hessian-transpose(hessian)))
   if (asym > max(1.0e-8_wp, 0.1_wp*hthr)) then
      print '(a,es14.6)', 'max Hessian asymmetry = ', asym
      call test_failed(error, "GFN1 validation Hessian is not symmetric")
      return
   end if

   ! Translational invariance: the sum over a rigid Cartesian translation must vanish.
   trans_err = 0.0_wp
   do row = 1, 3*mol%nat
      do ic = 1, 3
         tsum = 0.0_wp
         do iat = 1, mol%nat
            tsum = tsum + hessian(row, 3*(iat-1)+ic)
         end do
         trans_err = max(trans_err, abs(tsum))
      end do
   end do
   if (trans_err > max(2.0e-7_wp, 0.2_wp*hthr)) then
      print '(a,es14.6)', 'max translational residual = ', trans_err
      call test_failed(error, "GFN1 validation Hessian violates translational invariance")
      return
   end if

   ! Rotational covariance. Away from a stationary point an infinitesimal
   ! rigid rotation also rotates the Cartesian gradient, therefore the exact
   ! identity is H*(omega x R) = omega x g rather than H*(omega x R) = 0.
   ! Removing translation from the rotation vector makes the check independent
   ! of the coordinate origin, given the acoustic sum rule above.
   allocate(rvec(3*mol%nat), hrot(3*mol%nat), grot(3*mol%nat))
   rot_err = 0.0_wp
   do iaxis = 1, 3
      rvec = 0.0_wp
      grot = 0.0_wp
      do iat = 1, mol%nat
         select case(iaxis)
         case(1)
            rvec(3*(iat-1)+2) = -mol%xyz(3,iat)
            rvec(3*(iat-1)+3) =  mol%xyz(2,iat)
            grot(3*(iat-1)+2) = -g0(3,iat)
            grot(3*(iat-1)+3) =  g0(2,iat)
         case(2)
            rvec(3*(iat-1)+1) =  mol%xyz(3,iat)
            rvec(3*(iat-1)+3) = -mol%xyz(1,iat)
            grot(3*(iat-1)+1) =  g0(3,iat)
            grot(3*(iat-1)+3) = -g0(1,iat)
         case default
            rvec(3*(iat-1)+1) = -mol%xyz(2,iat)
            rvec(3*(iat-1)+2) =  mol%xyz(1,iat)
            grot(3*(iat-1)+1) = -g0(2,iat)
            grot(3*(iat-1)+2) =  g0(1,iat)
         end select
      end do
      do ic = 1, 3
         tsum = 0.0_wp
         do iat = 1, mol%nat
            tsum = tsum + rvec(3*(iat-1)+ic)
         end do
         tsum = tsum/real(mol%nat,wp)
         do iat = 1, mol%nat
            rvec(3*(iat-1)+ic) = rvec(3*(iat-1)+ic)-tsum
         end do
      end do
      hrot = matmul(hessian, rvec)
      rot_err = max(rot_err, maxval(abs(hrot-grot)))
   end do
   if (rot_err > max(2.0e-6_wp, 0.5_wp*hthr)) then
      print '(a,es14.6)', 'max rotational covariance residual = ', rot_err
      call test_failed(error, "GFN1 validation Hessian violates rotational covariance")
      return
   end if

   errmax = 0.0_wp
   do k = 1, size(columns)
      col = columns(k)
      if (col < 1 .or. col > 3*mol%nat) then
         call test_failed(error, "Invalid validation Hessian column")
         return
      end if
      jat = (col-1)/3 + 1
      jc = mod(col-1, 3) + 1
      molp = mol
      molm = mol
      molp%xyz(jc,jat) = molp%xyz(jc,jat) + step
      molm%xyz(jc,jat) = molm%xyz(jc,jat) - step

      call new_wavefunction(wfnp, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, etemp)
      gp = 0.0_wp
      sigma = 0.0_wp
      call xtb_singlepoint(ctx, molp, calc, wfnp, racc, energy, gradient=gp, sigma=sigma, verbosity=0)
      if (ctx%failed()) then
         call test_failed(error, "Positive displaced GFN1 validation SCF failed")
         return
      end if

      call new_wavefunction(wfnm, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, etemp)
      gm = 0.0_wp
      sigma = 0.0_wp
      call xtb_singlepoint(ctx, molm, calc, wfnm, racc, energy, gradient=gm, sigma=sigma, verbosity=0)
      if (ctx%failed()) then
         call test_failed(error, "Negative displaced GFN1 validation SCF failed")
         return
      end if

      refcol = reshape((gp-gm)/(2.0_wp*step), [3*mol%nat])
      errmax = max(errmax, maxval(abs(hessian(:,col)-refcol)))
   end do

   if (errmax > hthr) then
      print '(a,es14.6)', 'max Hessian finite-difference error = ', errmax
      call test_failed(error, "GFN1 analytic Hessian failed broad finite-difference validation")
      return
   end if
end subroutine validate_hessian


subroutine test_diatomics(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call validate_hessian(mol, 1, 0.0_wp, [1,2,3,4,5,6], 2.0e-6_wp, error)
   if (allocated(error)) return

   call get_structure(mol, "MB16-43", "LiH")
   call validate_hessian(mol, 1, 0.0_wp, [1,3,4,6], 4.0e-6_wp, error)
   if (allocated(error)) return

   call get_structure(mol, "MB16-43", "F2")
   call validate_hessian(mol, 1, 0.0_wp, [1,3,4,6], 8.0e-6_wp, error)
   if (allocated(error)) return

   call get_structure(mol, "MB16-43", "N2")
   call validate_hessian(mol, 1, 0.0_wp, [1,2,3,4,5,6], 8.0e-6_wp, error)
   if (allocated(error)) return

   ! Closed-shell S2 has an exactly degenerate frontier at zero electronic
   ! temperature; use a finite-T state-averaged response for validation.
   call get_structure(mol, "MB16-43", "S2")
   call validate_hessian(mol, 1, 2.0e-2_wp, [1,3,4,6], 2.0e-5_wp, error)
   if (allocated(error)) return

   call get_structure(mol, "MB16-43", "Cl2")
   call validate_hessian(mol, 1, 0.0_wp, [1,2,3,4,5,6], 1.0e-5_wp, error)
   if (allocated(error)) return

   ! PCl has an exactly degenerate frontier pair at zero electronic temperature;
   ! finite-T occupations make the state-averaged response differentiable.
   call get_structure(mol, "MB16-43", "PCl")
   call validate_hessian(mol, 1, 2.0e-2_wp, [1,3,4,6], 2.0e-5_wp, error)
end subroutine test_diatomics


subroutine test_main_group(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "CH4")
   call validate_hessian(mol, 1, 0.0_wp, [1,2,4,8,12], 6.0e-6_wp, error)
   if (allocated(error)) return

   call get_structure(mol, "MB16-43", "SiH4")
   call validate_hessian(mol, 1, 0.0_wp, [1,3,5,9,13], 1.0e-5_wp, error)
end subroutine test_main_group


subroutine test_charged(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp), parameter :: xyz(3,3) = reshape([ &
      &  0.0000000000_wp,  1.6500000000_wp, 0.0000000000_wp, &
      & -1.4289419162_wp, -0.8250000000_wp, 0.0000000000_wp, &
      &  1.4289419162_wp, -0.8250000000_wp, 0.0000000000_wp], [3,3])

   ! Equilateral H3+ (two-electron closed-shell cation), deliberately not optimized.
   call new(mol, [1,1,1], xyz, charge=1.0_wp, uhf=0)
   call validate_hessian(mol, 1, 0.0_wp, [1,2,3,4,7], 8.0e-6_wp, error)
   if (allocated(error)) return

   ! Closed-shell hydroxide anion exercises the negative-charge path.
   call new(mol, [8,1], reshape([ &
      & 0.0000000000_wp, 0.0000000000_wp, 0.0000000000_wp, &
      & 0.0000000000_wp, 0.0000000000_wp, 1.8500000000_wp], [3,2]), &
      & charge=-1.0_wp, uhf=0)
   call validate_hessian(mol, 1, 0.0_wp, [1,2,3,4,5,6], 1.5e-5_wp, error)
end subroutine test_charged


subroutine test_open_shell(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   real(wp), parameter :: xyz_o2(3,2) = reshape([ &
      & 0.0000000000_wp, 0.0000000000_wp, -1.1521537400_wp, &
      & 0.0000000000_wp, 0.0000000000_wp,  1.1521537400_wp], [3,2])
   real(wp), parameter :: xyz_oh(3,2) = reshape([ &
      & 0.0000000000_wp, 0.0000000000_wp, 0.0000000000_wp, &
      & 0.0000000000_wp, 0.0000000000_wp, 1.8340000000_wp], [3,2])

   ! Triplet O2 exercises exactly degenerate equally occupied pi manifolds at T=0.
   call new(mol, [8,8], xyz_o2, charge=0.0_wp, uhf=2)
   call validate_hessian(mol, 2, 0.0_wp, [1,2,3,4,5,6], 2.0e-5_wp, error, with_spin=.true.)
   if (allocated(error)) return

   ! OH(2Pi) has a zero-T frontier degeneracy. Finite-T provides the smooth,
   ! state-averaged occupation response and validates the degenerate Fermi path.
   call new(mol, [8,1], xyz_oh, charge=0.0_wp, uhf=1)
   call validate_hessian(mol, 2, 2.0e-2_wp, [1,2,3,4,5,6], 3.0e-5_wp, error, with_spin=.true.)
end subroutine test_open_shell


subroutine test_frontier_degeneracy_guard(error)
   type(error_type), allocatable, intent(out) :: error
   type(error_type), allocatable :: expected_error
   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), allocatable :: hessian(:, :)

   call get_structure(mol, "MB16-43", "PCl")
   ctx%solver = lapack_solver(lapack_algorithm%gvr)
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   allocate(hessian(3*mol%nat,3*mol%nat))
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 0.0_wp)
   call xtb_singlepoint(ctx, mol, calc, wfn, racc, energy, verbosity=0)
   if (ctx%failed()) then
      call test_failed(error, "PCl central SCF failed in degeneracy guard")
      return
   end if
   call get_gfn1_hessian(mol, calc, wfn, racc, hessian, expected_error)
   if (.not.allocated(expected_error)) then
      call test_failed(error, "Zero-T frontier degeneracy was not detected")
      return
   end if
end subroutine test_frontier_degeneracy_guard


subroutine test_diverse_mb01(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol

   ! MB16-43/01 simultaneously exercises H, B, N, O, F, Na, Al and Cl.
   ! Only representative columns are finite-differenced to keep CI cost bounded;
   ! the complete analytic Hessian is still formed and checked for symmetry and ASR.
   call get_structure(mol, "MB16-43", "01")
   call validate_hessian(mol, 1, 0.0_wp, [1,8,17,26,35,44], 3.0e-5_wp, error)
end subroutine test_diverse_mb01

end module test_gfn1_hessian_validation
