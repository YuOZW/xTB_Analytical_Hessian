! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Analytic Cartesian Hessian assembly for the initial GFN1-xTB implementation.
module tblite_xtb_hessian
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_container, only : container_cache
   use tblite_coulomb_charge_effective, only : effective_coulomb
   use tblite_cutoff, only : get_lattice_points
   use tblite_disp_d3, only : d3_dispersion
   use tblite_integral_native_integrals, only : overlap_hess_cgto, msao
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_ncoord, only : add_coordination_number_hessian
   use tblite_partition, only : work_partition, owns_index, owns_pair
   use tblite_scf_potential, only : potential_type, new_potential, add_pot_to_h1
   use tblite_wavefunction_type, only : wavefunction_type, get_density_matrix
   use tblite_wavefunction_spin, only : updown_to_magnet
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_hamiltonian_gradient, &
      & get_hamiltonian_hessian, get_hamiltonian_cn_derivatives
   use tblite_xtb_response, only : add_gfn1_response_hessian
   implicit none
   private

   public :: get_gfn1_hessian

   real(wp), parameter :: gfn1_cn_cutoff = 25.0_wp

contains

!> Assemble the analytic Cartesian Hessian for the initial GFN1-xTB scope:
!> finite molecular systems with one or two spin channels, finite electronic
!> temperature, GFN1 halogen correction, and spin-polarization interactions.
subroutine get_gfn1_hessian(mol, calc, wfn, accuracy, hessian, error)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn
   real(wp), intent(in) :: accuracy
   real(wp), intent(out) :: hessian(:, :)
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, nao, nsh, nspin, iat, jat, ic, jc, kat, row, col, iao, ish, spin
   real(wp) :: cutoff
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: selfenergy(:), dsedcn(:), lattr(:, :), cntrans(:, :), rcov(:)
   real(wp), allocatable :: density(:, :, :), wdensity(:, :, :)
   real(wp), allocatable :: occ(:), dEdcn(:), mixed(:, :, :)
   real(wp), allocatable :: dEdcn_resp(:), gradtmp(:, :), sigma(:, :)
   real(wp), allocatable :: xzero(:, :, :)
   real(wp), allocatable :: zeroself(:), zerodsed(:), h1dummy(:, :, :)
   type(adjacency_list) :: list
   type(integral_type) :: ints
   type(container_cache) :: ccache, icache
   type(potential_type) :: pot, potg, potb

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   nspin = wfn%nspin

   hessian = 0.0_wp
   if (size(hessian,1) /= 3*nat .or. size(hessian,2) /= 3*nat) then
      call fatal_error(error, "Invalid GFN1 Hessian dimensions")
      return
   end if
   if (trim(calc%method) /= "gfn1") then
      call fatal_error(error, "Analytic Hessian is currently implemented for GFN1-xTB only")
      return
   end if
   if (any(mol%periodic)) then
      call fatal_error(error, "GFN1 analytic Hessian currently supports non-periodic systems only")
      return
   end if
   if (nspin < 1 .or. nspin > 2) then
      call fatal_error(error, "GFN1 analytic Hessian supports one or two spin channels")
      return
   end if
   if (calc%h0%do_diat_scale) then
      call fatal_error(error, "GFN1 analytic Hessian does not support diatomic-frame H0 scaling")
      return
   end if
   if (.not.allocated(calc%coulomb) .or. .not.allocated(calc%coulomb%es2)) then
      call fatal_error(error, "GFN1 analytic Hessian requires the SCC Coulomb model")
      return
   end if
   if (allocated(calc%coulomb%aes2)) then
      call fatal_error(error, "GFN1 analytic Hessian does not support anisotropic SCC terms")
      return
   end if

   ! Coordination-number data and its first derivatives.
   allocate(cn(nat), dcndr(3,nat,nat), dcndL(3,3,nat))
   cn = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp
   if (allocated(calc%ncoord)) call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)

   allocate(selfenergy(nsh), dsedcn(nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, &
      & cn=cn, selfenergy=selfenergy, dsedcn=dsedcn)

   cutoff = get_cutoff(calc%bas, accuracy)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   call new_integral(ints, nao)
   call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)

   ! Central density and energy-weighted density in charge/magnetization
   ! representation, matching get_hamiltonian_gradient.
   allocate(density(nao,nao,nspin), wdensity(nao,nao,nspin), occ(nao))
   density = wfn%density
   if (nspin == 1) then
      occ = wfn%focc(:,1) + wfn%focc(:,2)
      call get_density_matrix(occ*wfn%emo(:,1), wfn%coeff(:,:,1), wdensity(:,:,1))
   else
      do spin = 1, nspin
         call get_density_matrix(wfn%focc(:,spin)*wfn%emo(:,spin), &
            & wfn%coeff(:,:,spin), wdensity(:,:,spin))
      end do
      call updown_to_magnet(density)
      call updown_to_magnet(wdensity)
   end if

   ! Classical explicit second derivatives.  D3 writes its result rather than
   ! accumulating it, so evaluate it first directly into the final Hessian and
   ! avoid a second O((3N)^2) temporary matrix.
   if (allocated(calc%dispersion)) then
      select type(disp => calc%dispersion)
      type is(d3_dispersion)
         call disp%get_hessian(mol, hessian, error)
         if (allocated(error)) return
      class default
         call fatal_error(error, "GFN1 analytic Hessian currently supports D3 dispersion only")
         return
      end select
   end if
   if (allocated(calc%repulsion)) call calc%repulsion%get_hessian(mol, hessian)
   if (allocated(calc%halogen)) call calc%halogen%get_hessian(mol, hessian)

   ! Fixed-charge SCC2 Hessian.  Charge-response terms are added after the
   ! streamed CP solve has populated wfn%dqshdr.
   select type(es2 => calc%coulomb%es2)
   type is(effective_coulomb)
      call es2%get_hessian(mol, wfn%qat(:,1), wfn%qsh(:,1), hessian)
   class default
      call fatal_error(error, "GFN1 analytic Hessian requires effective Coulomb electrostatics")
      return
   end select

   ! Central SCC/solvation potential at the converged electronic state.  This is
   ! available before the CP solve and is used by the streamed P'/X' accumulator.
   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%interactions)) call calc%interactions%update(mol, icache)
   call new_potential(pot, mol, calc%bas, nspin)
   call pot%reset
   call calc%coulomb%get_potential(mol, ccache, wfn, pot)
   if (allocated(calc%interactions)) call calc%interactions%get_potential(mol, icache, wfn, pot)
   allocate(h1dummy(nao,nao,nspin))
   call add_pot_to_h1(calc%bas, ints, pot, h1dummy)

   ! Explicit H0/Pulay second derivative at fixed electronic state and fixed CN.
   call get_hamiltonian_hessian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & density, wdensity, hessian, calc%partition)

   ! The SCC-potential contribution to the Pulay term contains an additional
   ! fixed-potential overlap second derivative: -P(V_mu+V_nu) S''.
   call add_overlap_potential_hessian(mol, lattr, list, calc, pot, density, hessian)

   ! H0 coordination-number chain rule at fixed density.
   allocate(dEdcn(nat), mixed(3,nat,nat))
   dEdcn = 0.0_wp
   mixed = 0.0_wp
   call get_hamiltonian_cn_derivatives(mol, lattr, list, calc%bas, calc%h0, &
      & dsedcn, density, dEdcn, mixed, calc%partition)

   do iat = 1, nat
      do ic = 1, 3
         row = 3*(iat-1) + ic
         do jat = 1, nat
            do jc = 1, 3
               col = 3*(jat-1) + jc
               do kat = 1, nat
                  hessian(row,col) = hessian(row,col) &
                     & + mixed(ic,iat,kat)*dcndr(jc,jat,kat) &
                     & + dcndr(ic,iat,kat)*mixed(jc,jat,kat)
               end do
            end do
         end do
      end do
   end do

   allocate(rcov(mol%nid))
   rcov = get_covalent_rad(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, gfn1_cn_cutoff, cntrans)
   call add_coordination_number_hessian(mol, cntrans, gfn1_cn_cutoff, rcov, &
      & dEdcn, hessian, calc%partition)

   ! Electronic response columns.  Stream P'/X' directly into the Hessian
   ! instead of retaining two nao*nao*nspin*3*nat tensors.  The response solver
   ! still persists q'/m' on wfn for SCC and solvation response terms below.
   allocate(dEdcn_resp(nat), gradtmp(3,nat), sigma(3,3))
   call add_gfn1_response_hessian(mol, calc, wfn, accuracy, hessian, error)
   if (allocated(error)) return

   ! SCC2 charge response now that dq/dR is available.
   select type(es2 => calc%coulomb%es2)
   type is(effective_coulomb)
      call add_coulomb_charge_response(es2, mol, wfn%qsh(:,1), &
         & wfn%dqshdr(:,:,:,1), hessian, calc%partition)
   end select

   ! Solvation/container Hessians include their own q' response terms.
   if (allocated(calc%interactions)) then
      call calc%interactions%get_interaction_hessian(mol, icache, wfn, hessian)
   end if

   ! Full SCC/solvation potential derivative (explicit geometry + q' response).
   call new_potential(potg, mol, calc%bas, nspin, grad=.true.)
   call potg%reset
   call calc%coulomb%get_potential_gradient(mol, ccache, wfn, potg)
   if (allocated(calc%interactions)) call calc%interactions%get_potential_gradient(mol, icache, wfn, potg)

   allocate(xzero(nao,nao,nspin))
   allocate(zeroself(nsh), zerodsed(nsh))
   zeroself = 0.0_wp
   zerodsed = 0.0_wp
   xzero = 0.0_wp
   call new_potential(potb, mol, calc%bas, nspin)

   ! Response of the SCC potential in -P(V_mu+V_nu)S'.
   do jat = 1, nat
      do jc = 1, 3
         col = 3*(jat-1) + jc
         call potb%reset
         do spin = 1, nspin
            do iao = 1, nao
               ish = calc%bas%ao2sh(iao)
               iat = calc%bas%ao2at(iao)
               potb%vao(iao,spin) = potg%dvshdr(jc,jat,ish,spin) &
                  & + potg%dvatdr(jc,jat,iat,spin)
            end do
         end do
         dEdcn_resp = 0.0_wp
         gradtmp = 0.0_wp
         sigma = 0.0_wp
         call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, &
            & zeroself, zerodsed, potb, density, xzero, dEdcn_resp, &
            & gradtmp, sigma, calc%partition)
         hessian(:,col) = hessian(:,col) + reshape(gradtmp, [3*nat])
      end do
   end do

   ! The exact Cartesian Hessian is symmetric.  Reject a genuinely inconsistent
   ! assembly, but remove harmless floating-point skew before exposing the
   ! result to downstream normal-mode analysis.
   call stabilize_hessian(hessian, error)
   if (allocated(error)) return

end subroutine get_gfn1_hessian


!> Validate basic numerical properties and remove roundoff-level antisymmetry.
subroutine stabilize_hessian(hessian, error)
   real(wp), intent(inout) :: hessian(:, :)
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: scale, asym, tol

   if (.not.all(ieee_is_finite(hessian))) then
      call fatal_error(error, "GFN1 analytic Hessian contains non-finite values")
      return
   end if

   scale = max(1.0_wp, maxval(abs(hessian)))
   asym = maxval(abs(hessian-transpose(hessian)))
   tol = max(5.0e-8_wp, 5.0e-8_wp*scale)
   if (asym > tol) then
      call fatal_error(error, "GFN1 analytic Hessian violates permutation symmetry")
      return
   end if

   hessian = 0.5_wp*(hessian + transpose(hessian))
end subroutine stabilize_hessian


!> Add -P(V_mu+V_nu) S'' at fixed P and fixed SCC potential.
subroutine add_overlap_potential_hessian(mol, trans, list, calc, pot, pmat, hessian)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   type(adjacency_list), intent(in) :: list
   type(xtb_calculator), intent(in) :: calc
   type(potential_type), intent(in) :: pot
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, jat, izp, jzp, itr, img, inl, is, js, ish, jsh
   integer :: ii, jj, iao, jao, naoi, naoj, ic, jc, irow, jcol, mmax, spin
   real(wp) :: vec(3), r2, sval, block(3,3)
   real(wp), allocatable :: stmp(:, :), dstmp(:, :, :), hstmp(:, :, :, :)

   mmax = msao(calc%bas%maxl)
   allocate(stmp(mmax,mmax), dstmp(3,mmax,mmax), hstmp(3,3,mmax,mmax))

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = calc%bas%ish_at(iat)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         if (.not.owns_index(calc%partition, img+inl)) cycle
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = calc%bas%ish_at(jat)
         vec = mol%xyz(:,iat) - mol%xyz(:,jat) - trans(:,itr)
         r2 = dot_product(vec,vec)
         if (r2 <= epsilon(1.0_wp)) cycle

         block = 0.0_wp
         do ish = 1, calc%bas%nsh_id(izp)
            ii = calc%bas%iao_sh(is+ish)
            naoi = msao(calc%bas%cgto(ish,izp)%ang)
            do jsh = 1, calc%bas%nsh_id(jzp)
               jj = calc%bas%iao_sh(js+jsh)
               naoj = msao(calc%bas%cgto(jsh,jzp)%ang)
               call overlap_hess_cgto(calc%bas%cgto(jsh,jzp), calc%bas%cgto(ish,izp), &
                  & r2, vec, calc%bas%intcut, stmp(1:naoj,1:naoi), &
                  & dstmp(:,1:naoj,1:naoi), hstmp(:,:,1:naoj,1:naoi))
               do iao = 1, naoi
                  do jao = 1, naoj
                     do spin = 1, size(pmat,3)
                        sval = -pmat(jj+jao,ii+iao,spin) &
                           & *(pot%vao(jj+jao,spin)+pot%vao(ii+iao,spin))
                        block = block + sval*hstmp(:,:,jao,iao)
                     end do
                  end do
               end do
            end do
         end do

         do ic = 1, 3
            irow = 3*(iat-1)+ic
            do jc = 1, 3
               jcol = 3*(jat-1)+jc
               hessian(irow,3*(iat-1)+jc) = hessian(irow,3*(iat-1)+jc) + block(ic,jc)
               hessian(3*(jat-1)+ic,jcol) = hessian(3*(jat-1)+ic,jcol) + block(ic,jc)
               hessian(irow,jcol) = hessian(irow,jcol) - block(ic,jc)
               hessian(3*(jat-1)+ic,3*(iat-1)+jc) = &
                  & hessian(3*(jat-1)+ic,3*(iat-1)+jc) - block(ic,jc)
            end do
         end do
      end do
   end do
end subroutine add_overlap_potential_hessian


!> Add the SCC2 Hessian terms proportional to the first-order shell-charge
!> response, q'^T A_R q.  The fixed-charge A_RR term is provided by Phase 5.
subroutine add_coulomb_charge_response(es2, mol, qsh, dqshdr, hessian, partition)
   type(effective_coulomb), intent(in) :: es2
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(in) :: dqshdr(:, :, :)
   real(wp), intent(inout) :: hessian(:, :)
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh, ic, jc, kat, row, col
   real(wp) :: vec(3), r1, r2, gam, t, pref, dG(3), qresp

   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = es2%offset(iat)
      do jat = 1, iat-1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         jj = es2%offset(jat)
         vec = mol%xyz(:,iat) - mol%xyz(:,jat)
         r2 = dot_product(vec,vec)
         if (r2 < 1.0e-16_wp) cycle
         r1 = sqrt(r2)

         do ish = 1, es2%nshell(iat)
            do jsh = 1, es2%nshell(jat)
               gam = es2%hubbard(jsh,ish,jzp,izp)
               t = r1**es2%gexp + gam**(-es2%gexp)
               pref = -r1**(es2%gexp-2.0_wp) &
                  & *(1.0_wp/t) * (1.0_wp/t)**(1.0_wp/es2%gexp)
               dG = pref*vec
               do kat = 1, mol%nat
                  do jc = 1, 3
                     col = 3*(kat-1)+jc
                     qresp = dqshdr(jc,kat,ii+ish)*qsh(jj+jsh) &
                        & + qsh(ii+ish)*dqshdr(jc,kat,jj+jsh)
                     do ic = 1, 3
                        row = 3*(iat-1)+ic
                        hessian(row,col) = hessian(row,col) + dG(ic)*qresp
                        row = 3*(jat-1)+ic
                        hessian(row,col) = hessian(row,col) - dG(ic)*qresp
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine add_coulomb_charge_response

end module tblite_xtb_hessian
