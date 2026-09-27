! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/xtb/h0.f90
!> Provides the effective core Hamiltonian for xTB.

!> Implementation of the effective core Hamiltonian used in the extended tight binding.
module tblite_xtb_h0
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list
   use tblite_basis_type, only : basis_type
   use tblite_integral_diat_trafo, only : diat_trafo_cache, setup_diat_trafo, &
      & diat_trafo
   use tblite_integral_native_integrals, only : multipole_cgto, multipole_grad_cgto, &
      & overlap_grad_cgto, overlap_hess_cgto, msao, smap, sdim
   use tblite_partition, only : work_partition, owns_index
   use tblite_scf_potential, only : potential_type
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: new_hamiltonian
   public :: get_selfenergy, get_hamiltonian, get_occupation, get_hamiltonian_gradient
   public :: get_hamiltonian_derivative, get_hamiltonian_hessian
   public :: get_hamiltonian_cn_derivatives


   type, public :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Electronegativity scaled coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn_en(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq1(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq2(:, :)
      !> Enhancement factor to scale the Hamiltonian elements
      real(wp), allocatable :: hscale(:, :, :, :)
      !> Polynomial coefficients for distance dependent enhancement factor
      real(wp), allocatable :: shpoly(:, :)
      !> Atomic radius for polynomial enhancement
      real(wp), allocatable :: rad(:)
      !> Reference occupation numbers
      real(wp), allocatable :: refocc(:, :)
      !> Diatomic frame scaling of sigma bonding contribution
      real(wp), allocatable :: ksig(:, :)
      !> Diatomic frame scaling of pi bonding contribution
      real(wp), allocatable :: kpi(:, :)
      !> Diatomic frame scaling of delta bonding contribution
      real(wp), allocatable :: kdel(:, :)
      !> Perform diatomic frame scaling for the overlap in the Hamiltonian
      logical :: do_diat_scale = .false.
   end type tb_hamiltonian


contains


!> Constructor for a new Hamiltonian object, consumes a Hamiltonian specification
subroutine new_hamiltonian(self, mol, bas, spec)
   !> Hamiltonian object
   type(tb_hamiltonian), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   integer :: mshell

   mshell = maxval(bas%nsh_id)
   allocate(self%selfenergy(mshell, mol%nid), &
      & self%kcn(mshell, mol%nid), self%kcn_en(mshell, mol%nid), &
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
   call spec%get_selfenergy(mol, bas, self%selfenergy)
   call spec%get_cnshift(mol, bas, self%kcn)
   call spec%get_cnenshift(mol, bas, self%kcn_en)
   call spec%get_q1shift(mol, bas, self%kq1)
   call spec%get_q2shift(mol, bas, self%kq2)

   allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
   call spec%get_hscale(mol, bas, self%hscale)

   allocate(self%shpoly(mshell, mol%nid), self%rad(mol%nid))
   call spec%get_rad(mol, bas, self%rad)
   call spec%get_shpoly(mol, bas, self%shpoly)

   allocate(self%refocc(mshell, mol%nid))
   call spec%get_reference_occ(mol, bas, self%refocc)

   allocate(self%ksig(mol%nid, mol%nid), self%kpi(mol%nid, mol%nid), &
      & self%kdel(mol%nid, mol%nid))
   call spec%get_diat_scale(mol, bas, self%ksig, self%kpi, self%kdel, &
      & self%do_diat_scale)

end subroutine new_hamiltonian


subroutine get_selfenergy(h0, id, ish_at, nshell, cn, cn_en, qat, selfenergy, &
   & dsedcn, dsedcn_en, dsedq)
   !> Hamiltonian object
   type(tb_hamiltonian), intent(in) :: h0
   !> Species identifiers
   integer, intent(in) :: id(:)
   !> Index of first shell per atom
   integer, intent(in) :: ish_at(:)
   !> Number of shells per element
   integer, intent(in) :: nshell(:)
   !> Optional coordination number
   real(wp), intent(in), optional :: cn(:)
   !> Optional electronegativity scaled coordination number
   real(wp), intent(in), optional :: cn_en(:)
   !> Optional atomic partial charge
   real(wp), intent(in), optional :: qat(:)
   !> Reference energy levels
   real(wp), intent(out) :: selfenergy(:)
   !> Derivative of energy levels w.r.t. coordination number
   real(wp), intent(out), optional :: dsedcn(:)
   !> Derivative of energy levels w.r.t. EN-scaled coordination number
   real(wp), intent(out), optional :: dsedcn_en(:)
   !> Derivative of energy levels w.r.t. atomic charge
   real(wp), intent(out), optional :: dsedq(:)

   integer :: iat, izp, ish, ii

   selfenergy(:) = 0.0_wp
   if (present(dsedcn)) dsedcn(:) = 0.0_wp
   if (present(dsedq)) dsedq(:) = 0.0_wp
   do iat = 1, size(id)
      izp = id(iat)
      ii = ish_at(iat)
      do ish = 1, nshell(izp)
         selfenergy(ii+ish) = h0%selfenergy(ish, izp)
      end do
   end do
   if (present(cn)) then
      if (present(dsedcn)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn(ish, izp) * cn(iat)
               dsedcn(ii+ish) = -h0%kcn(ish, izp)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn(ish, izp) * cn(iat)
            end do
         end do
      end if
   end if
   if (present(cn_en)) then
      if (present(dsedcn_en)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn_en(ish, izp) * cn_en(iat)
               dsedcn_en(ii+ish) = -h0%kcn_en(ish, izp)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn_en(ish, izp) * cn_en(iat)
            end do
         end do
      end if
   end if
   if (present(qat)) then
      if (present(dsedq)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
               dsedq(ii+ish) = -h0%kq1(ish, izp) - h0%kq2(ish, izp)*2*qat(iat)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
            end do
         end do
      end if
   end if

end subroutine get_selfenergy


subroutine get_hamiltonian(mol, trans, list, bas, h0, selfenergy, overlap, &
   & dpint, qpint, hamiltonian, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Diagonal elememts of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Overlap integral matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)
   !> Effective Hamiltonian
   real(wp), intent(out) :: hamiltonian(:, :)
   !> Share of the neighbour list entries evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, itr, img, inl
   integer :: ish, jsh, is, js, nsi, nsj, ii, jj, iao, jao, nao, ij, iaosh, jaosh
   real(wp) :: rr, r2, vec(3), hij, shpolyi, shpoly, dtmpj(3), qtmpj(6)
   real(wp) :: mod_h0_fraction
   real(wp), allocatable :: stmp(:), dtmpi(:, :), qtmpi(:, :), block_overlap(:, :)
   type(diat_trafo_cache) :: dt_cache

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp
   hamiltonian(:, :) = 0.0_wp

   ! Select if we construct the Hamiltonian with modifications
   mod_h0_fraction = 0.0_wp
   if(h0%do_diat_scale) then
      mod_h0_fraction = 1.0_wp
   end if

   allocate(stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2), &
      & qtmpi(6, msao(bas%maxl)**2), block_overlap(sdim(bas%maxl), sdim(bas%maxl)))

   !$omp parallel do schedule(runtime) default(none) &
   !$omp firstprivate(mod_h0_fraction) &
   !$omp shared(mol, bas, trans, list, overlap, dpint, qpint, hamiltonian, h0, selfenergy) &
   !$omp shared(partition) &
   !$omp private(iat, jat, izp, jzp, itr, inl, img, is, js, ish, jsh, nsi, nsj, ii, jj) &
   !$omp private(iao, jao, iaosh, jaosh, nao, ij, r2, vec, hij, shpolyi, shpoly, rr) &
   !$omp private(stmp, dtmpi, qtmpi, dtmpj, qtmpj, block_overlap, dt_cache)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         ! the offset into the neighbour list enumerates the diatomic blocks
         if (.not.owns_index(partition, img+inl)) cycle
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
         ! Get the overlap and multipole integrals for the current diatomic pair
         if (h0%do_diat_scale) block_overlap = 0.0_wp
         do ish = 1, nsi
            ii = bas%iao_sh(is+ish)
            iaosh = smap(ish-1)

            ! Shell polynomial enhancement factor
            shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr
            do jsh = 1, nsj
               jj = bas%iao_sh(js+jsh)
               jaosh = smap(jsh-1)
               call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)

               shpoly = shpolyi * (1.0_wp + h0%shpoly(jsh, jzp)*rr)

               hij = 0.5_wp * (1.0_wp - mod_h0_fraction) &
                  & * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                  & * h0%hscale(jsh, ish, jzp, izp) * shpoly

               nao = msao(bas%cgto(jsh, jzp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     call shift_operator(vec, stmp(ij), dtmpi(:, ij), qtmpi(:, ij), &
                        & dtmpj, qtmpj)

                     ! Save overlap for possible diatomic frame trafo
                     if (h0%do_diat_scale) block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)

                     overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                        + stmp(ij)

                     dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                        + dtmpi(:, ij)

                     qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                        + qtmpi(:, ij)

                     hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                        + stmp(ij) * hij

                     if (iat /= jat) then
                        overlap(ii+iao, jj+jao) = overlap(ii+iao, jj+jao) &
                           + stmp(ij)

                        dpint(:, ii+iao, jj+jao) = dpint(:, ii+iao, jj+jao) &
                           + dtmpj

                        qpint(:, ii+iao, jj+jao) = qpint(:, ii+iao, jj+jao) &
                           + qtmpj

                        hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                           + stmp(ij) * hij
                     end if
                  end do
               end do

            end do
         end do

         ! Optional diatomic frame transformation and scaling of the overlap
         if (h0%do_diat_scale) then
            call setup_diat_trafo(dt_cache, vec, nsj-1, nsi-1)
            call diat_trafo(dt_cache, h0%ksig(izp, jzp), h0%kpi(izp, jzp), &
               & h0%kdel(izp, jzp), block_overlap)
         end if

         ! Optional repeated setup of the Hamiltonian after modifications
         if (h0%do_diat_scale) then
            do ish = 1, nsi
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)

               ! Recalculate shell polynomial enhancement factor
               shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr

               do jsh = 1, nsj
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)

                  ! Recalculate shell polynomial enhancement factor
                  shpoly = shpolyi * (1.0_wp + h0%shpoly(jsh, jzp)*rr)

                  ! Recalculate Hamiltonian elements for modified H0
                  hij = 0.5_wp * mod_h0_fraction &
                     & * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                     & * h0%hscale(jsh, ish, jzp, izp) * shpoly

                  ! Distribute shell block to Hamiltonian matrix
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        ! Add modified Hamiltonian contribution
                        hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao) * hij

                        if (iat /= jat) then
                           hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                              + block_overlap(jaosh+jao, iaosh+iao) * hij
                        end if
                     end do
                  end do

               end do
            end do
         end if
      end do

      ! Onsite contribution to the Hamiltonian
      if (.not.owns_index(partition, iat)) cycle
      vec(:) = 0.0_wp
      do ish = 1, nsi
         ii = bas%iao_sh(is+ish)
         do jsh = 1, nsi
            jj = bas%iao_sh(is+jsh)
            call multipole_cgto(bas%cgto(jsh, izp), bas%cgto(ish, izp), &
               & 0.0_wp, vec, bas%intcut, stmp, dtmpi, qtmpi)

            ! shpoly is always 1.0, because rr is always 0.0
            hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(is+jsh))

            nao = msao(bas%cgto(jsh, izp)%ang)
            do iao = 1, msao(bas%cgto(ish, izp)%ang)
               do jao = 1, nao
                  ij = jao + nao*(iao-1)
                  overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                     + stmp(ij)

                  dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                     + dtmpi(:, ij)

                  qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                     + qtmpi(:, ij)

                  hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                     + stmp(ij) * hij
               end do
            end do

         end do
      end do
   end do

end subroutine get_hamiltonian



!> First Cartesian derivatives of the GFN1 zero-order Hamiltonian and overlap.
!>
!> The returned Hamiltonian derivative contains the explicit two-centre
!> geometry dependence and the chain rule through the GFN1 coordination-number
!> dependent self energies.  Density-dependent SCC potentials are deliberately
!> excluded and are added by the response module.  This routine currently
!> targets the non-periodic GFN1 Hamiltonian, for which diatomic-frame scaling
!> is disabled.
subroutine get_hamiltonian_derivative(mol, trans, list, bas, h0, selfenergy, &
      & dsedcn, dcndr, doverlap, dhamiltonian, partition)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   type(adjacency_list), intent(in) :: list
   type(basis_type), intent(in) :: bas
   type(tb_hamiltonian), intent(in) :: h0
   real(wp), intent(in) :: selfenergy(:)
   real(wp), intent(in) :: dsedcn(:)
   !> CN derivatives, dcndr(component, displaced atom, CN atom)
   real(wp), intent(in) :: dcndr(:, :, :)
   !> dS/dR, shape [nao, nao, 3, nat]
   real(wp), intent(out) :: doverlap(:, :, :, :)
   !> dH0/dR, shape [nao, nao, 3, nat]
   real(wp), intent(out) :: dhamiltonian(:, :, :, :)
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, kat, izp, jzp, itr, img, inl
   integer :: ish, jsh, is, js, nsi, nsj, ii, jj, iao, jao, naoi, naoj
   integer :: ic
   real(wp) :: r2, vec(3), rr, r0, api, apj
   real(wp) :: shpolyi, shpolyj, shpoly, hscale, pref, sval
   real(wp) :: dsh(3), dhpair(3), dse
   real(wp), allocatable :: stmp(:, :), dstmp(:, :, :)
   integer :: mmax

   doverlap = 0.0_wp
   dhamiltonian = 0.0_wp

   ! GFN1 has no diatomic-frame overlap scaling.  Keeping this guard makes an
   ! accidental use with another Hamiltonian fail closed rather than silently
   ! returning an inconsistent derivative.
   if (h0%do_diat_scale) return
   if (any(mol%periodic)) return

   mmax = msao(bas%maxl)
   allocate(stmp(mmax, mmax), dstmp(3, mmax, mmax))

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         if (.not.owns_index(partition, img+inl)) cycle
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = dot_product(vec, vec)
         if (r2 <= epsilon(1.0_wp)) cycle

         r0 = h0%rad(jzp) + h0%rad(izp)
         rr = sqrt(sqrt(r2)/r0)

         do ish = 1, nsi
            ii = bas%iao_sh(is+ish)
            api = h0%shpoly(ish, izp)
            shpolyi = 1.0_wp + api*rr
            naoi = msao(bas%cgto(ish, izp)%ang)
            do jsh = 1, nsj
               jj = bas%iao_sh(js+jsh)
               apj = h0%shpoly(jsh, jzp)
               shpolyj = 1.0_wp + apj*rr
               shpoly = shpolyi*shpolyj
               naoj = msao(bas%cgto(jsh, jzp)%ang)

               call overlap_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp(:naoj, :naoi), &
                  & dstmp(:, :naoj, :naoi))

               ! d(shpoly)/dR_i.  vec = R_i - R_j.
               dsh = (shpolyi*apj + shpolyj*api) * 0.5_wp*rr/r2 * vec
               hscale = h0%hscale(jsh, ish, jzp, izp)
               pref = 0.5_wp*(selfenergy(is+ish) + selfenergy(js+jsh))*hscale

               do iao = 1, naoi
                  do jao = 1, naoj
                     sval = stmp(jao, iao)
                     dhpair = pref*(shpoly*dstmp(:, jao, iao) + sval*dsh)

                     ! Relative-coordinate derivatives.
                     do ic = 1, 3
                        doverlap(jj+jao, ii+iao, ic, iat) = &
                           doverlap(jj+jao, ii+iao, ic, iat) + dstmp(ic, jao, iao)
                        dhamiltonian(jj+jao, ii+iao, ic, iat) = &
                           dhamiltonian(jj+jao, ii+iao, ic, iat) + dhpair(ic)
                        if (iat /= jat) then
                           doverlap(jj+jao, ii+iao, ic, jat) = &
                              doverlap(jj+jao, ii+iao, ic, jat) - dstmp(ic, jao, iao)
                           dhamiltonian(jj+jao, ii+iao, ic, jat) = &
                              dhamiltonian(jj+jao, ii+iao, ic, jat) - dhpair(ic)
                        end if
                     end do

                     ! CN dependence of both shell self energies.  A CN can
                     ! respond to displacement of any atom, so this term is not
                     ! restricted to the two atoms forming this AO block.
                     do kat = 1, mol%nat
                        do ic = 1, 3
                           dse = dsedcn(is+ish)*dcndr(ic, kat, iat) &
                              & + dsedcn(js+jsh)*dcndr(ic, kat, jat)
                           dhamiltonian(jj+jao, ii+iao, ic, kat) = &
                              dhamiltonian(jj+jao, ii+iao, ic, kat) &
                              & + 0.5_wp*hscale*shpoly*sval*dse
                        end do
                     end do

                     if (iat /= jat) then
                        do ic = 1, 3
                           doverlap(ii+iao, jj+jao, ic, iat) = &
                              doverlap(ii+iao, jj+jao, ic, iat) + dstmp(ic, jao, iao)
                           dhamiltonian(ii+iao, jj+jao, ic, iat) = &
                              dhamiltonian(ii+iao, jj+jao, ic, iat) + dhpair(ic)
                           doverlap(ii+iao, jj+jao, ic, jat) = &
                              doverlap(ii+iao, jj+jao, ic, jat) - dstmp(ic, jao, iao)
                           dhamiltonian(ii+iao, jj+jao, ic, jat) = &
                              dhamiltonian(ii+iao, jj+jao, ic, jat) - dhpair(ic)
                        end do
                        do kat = 1, mol%nat
                           do ic = 1, 3
                              dse = dsedcn(is+ish)*dcndr(ic, kat, iat) &
                                 & + dsedcn(js+jsh)*dcndr(ic, kat, jat)
                              dhamiltonian(ii+iao, jj+jao, ic, kat) = &
                                 dhamiltonian(ii+iao, jj+jao, ic, kat) &
                                 & + 0.5_wp*hscale*shpoly*sval*dse
                           end do
                        end do
                     end if
                  end do
               end do
            end do
         end do
      end do

      ! Onsite overlap is translationally invariant.  Only the CN dependent
      ! self energies contribute to dH0 for onsite AO blocks.
      do ish = 1, nsi
         ii = bas%iao_sh(is+ish)
         naoi = msao(bas%cgto(ish, izp)%ang)
         do jsh = 1, nsi
            jj = bas%iao_sh(is+jsh)
            naoj = msao(bas%cgto(jsh, izp)%ang)
            vec = 0.0_wp
            call overlap_grad_cgto(bas%cgto(jsh, izp), bas%cgto(ish, izp), &
               & 0.0_wp, vec, bas%intcut, stmp(:naoj, :naoi), &
               & dstmp(:, :naoj, :naoi))
            do iao = 1, naoi
               do jao = 1, naoj
                  sval = stmp(jao, iao)
                  do kat = 1, mol%nat
                     do ic = 1, 3
                        dse = dsedcn(is+ish)*dcndr(ic, kat, iat) &
                           & + dsedcn(is+jsh)*dcndr(ic, kat, iat)
                        dhamiltonian(jj+jao, ii+iao, ic, kat) = &
                           dhamiltonian(jj+jao, ii+iao, ic, kat) &
                           & + 0.5_wp*sval*dse
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine get_hamiltonian_derivative

subroutine get_hamiltonian_gradient(mol, trans, list, bas, h0, selfenergy, dsedcn, &
      & pot, pmat, xmat, dEdcn, gradient, sigma, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Diagonal elememts of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Derivative of the diagonal elements of the Hamiltonian w.r.t. the coordination number
   real(wp), intent(in) :: dsedcn(:)
   !> Density dependent potential shifts on the Hamiltonian
   type(potential_type), intent(in) :: pot
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Energy weighted density matrix
   real(wp), intent(in) :: xmat(:, :, :)

   !> Derivative of the electronic energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   !> Share of the neighbour list entries evaluated here, absent selects the complete work
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, itr, img, inl, spin, nspin
   integer :: ish, jsh, is, js, nsi, nsj, ii, jj, iao, jao, iaosh, jaosh, nao, ij
   real(wp) :: rr, r2, vec(3), hij, dG(3), hscale, hs
   real(wp) :: shpolyi, shpolyj, shpoly, dshpoly, dsv(3)
   real(wp) :: sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij, mod_h0_fraction
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)
   real(wp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
   real(wp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)
   real(wp), allocatable :: block_overlap(:, :), block_doverlap(:, :, :)
   type(diat_trafo_cache) :: dt_cache

   nspin = size(pmat, 3)

   ! Select if we construct the Hamiltonian with modifications
   mod_h0_fraction = 0.0_wp
   if(h0%do_diat_scale) then
      mod_h0_fraction = 1.0_wp
   end if

   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmp(3, msao(bas%maxl)**2), ddtmpi(3, 3, msao(bas%maxl)**2), &
      & qtmp(6, msao(bas%maxl)**2), dqtmpi(3, 6, msao(bas%maxl)**2), &
      & ddtmpj(3, 3, msao(bas%maxl)**2), dqtmpj(3, 6, msao(bas%maxl)**2), &
      & block_overlap(sdim(bas%maxl), sdim(bas%maxl)), &
      & block_doverlap(sdim(bas%maxl), sdim(bas%maxl), 3))

   !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn, gradient, sigma) &
   !$omp shared(nspin, mol, bas, trans, h0, selfenergy, dsedcn, pot, pmat, xmat, list) &
   !$omp shared(partition) &
   !$omp firstprivate(mod_h0_fraction) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, nsi, nsj, ii, jj) &
   !$omp private(iaosh, jaosh, iao, jao, nao, ij, inl, img, spin, r2, vec) &
   !$omp private(stmp, dstmp, dtmp, ddtmpi, ddtmpj, qtmp, dqtmpi, dqtmpj) &
   !$omp private(dG, sval, hscale, hs, hij, hpij, pij, dcni, dcnj, dhdcni, dhdcnj) &
   !$omp private(rr, shpolyi, shpolyj, shpoly, dshpoly, dsv) &
   !$omp private(block_overlap, block_doverlap, dt_cache)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         ! the offset into the neighbour list enumerates the diatomic blocks
         if (.not.owns_index(partition, img+inl)) cycle
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         if (r2 <= epsilon(1.0_wp)) cycle
         rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
         if (h0%do_diat_scale) then
            block_overlap(:, :) = 0.0_wp
            block_doverlap(:, :, :) = 0.0_wp
         end if
         dG(:) = 0.0_wp
         dcni = 0.0_wp
         dcnj = 0.0_wp
         do ish = 1, nsi
            ii = bas%iao_sh(is+ish)
            iaosh = smap(ish-1)
            shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr
            do jsh = 1, nsj
               jj = bas%iao_sh(js+jsh)
               jaosh = smap(jsh-1)
               call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddtmpj, dqtmpj, &
                  & ddtmpi, dqtmpi)
               shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr
               shpoly = shpolyi * shpolyj
               dshpoly = (shpolyi * h0%shpoly(jsh, jzp) + shpolyj * &
                  & h0%shpoly(ish, izp)) * 0.5_wp * rr / r2
               dsv(:) = dshpoly / shpoly * vec

               hscale = h0%hscale(jsh, ish, jzp, izp)
               hs = hscale * shpoly * (1.0_wp - mod_h0_fraction)
               hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hs

               dhdcni = dsedcn(is+ish) * hs
               dhdcnj = dsedcn(js+jsh) * hs

               nao = msao(bas%cgto(jsh, jzp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)

                     ! Overlap and derivatives for possible diatomic frame trafo
                     if (h0%do_diat_scale) then
                        block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)
                        block_doverlap(jaosh+jao, iaosh+iao, :) = dstmp(:, ij)
                     end if
                     do spin = 1, nspin
                        pij = pmat(jj+jao, ii+iao, spin)
                        sval = - pij * (pot%vao(jj+jao, spin) + pot%vao(ii+iao, spin))

                        dG(:) = dG + sval * dstmp(:, ij)
                     end do
                     pij = pmat(jj+jao, ii+iao, 1)
                     hpij = pij * hij
                     sval = 2*hpij - 2*xmat(jj+jao, ii+iao, 1)

                     dG(:) = dG + sval * dstmp(:, ij) &
                        + 2*hpij*stmp(ij) * dsv &
                        - pij * matmul(ddtmpi(:, :, ij), pot%vdp(:, iat, 1)) &
                        - pij * matmul(ddtmpj(:, :, ij), pot%vdp(:, jat, 1)) &
                        - pij * matmul(dqtmpi(:, :, ij), pot%vqp(:, iat, 1)) &
                        - pij * matmul(dqtmpj(:, :, ij), pot%vqp(:, jat, 1))

                     dcni = dcni + dhdcni * pmat(jj+jao, ii+iao, 1) * stmp(ij)
                     dcnj = dcnj + dhdcnj * pmat(jj+jao, ii+iao, 1) * stmp(ij)
                  end do
               end do
            end do
         end do

         ! Optional diatomic frame scaling transformation for overlap and derivatives
         if (h0%do_diat_scale) then
            call setup_diat_trafo(dt_cache, vec, nsj-1, nsi-1, grad=.true.)
            call diat_trafo(dt_cache, h0%ksig(izp, jzp), h0%kpi(izp, jzp), &
               & h0%kdel(izp, jzp), block_overlap, block_doverlap)
         end if

         ! Optional repeated setup of the Hamiltonian after modifications
         if (h0%do_diat_scale) then
            do ish = 1, nsi
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)

               ! Recalculate shell polynomial enhancement factor and derivative
               shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr

               do jsh = 1, nsj
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)

                  ! Recalculate shell polynomial enhancement factor and derivative
                  shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr
                  shpoly = shpolyi * shpolyj
                  dshpoly = (shpolyi * h0%shpoly(jsh, jzp) + shpolyj * &
                  & h0%shpoly(ish, izp)) * 0.5_wp * rr / r2

                  dsv(:) = dshpoly / shpoly * vec

                  ! Recalculate Hamiltonian elements
                  hscale = h0%hscale(jsh, ish, jzp, izp)
                  hs = hscale * shpoly * mod_h0_fraction
                  hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hs

                  ! CN derivative of modified Hamiltonian elements
                  dhdcni = dsedcn(is+ish) * hs
                  dhdcnj = dsedcn(js+jsh) * hs

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        ! Add only H0 gradient with the modified overlap
                        pij = pmat(jj+jao, ii+iao, 1)
                        hpij = 2*pij * hij

                        ! Accumulate gradient from overlap and shell-polynomials
                        dG(:) = dG + hpij * (block_doverlap(jaosh+jao, iaosh+iao, :) &
                           + block_overlap(jaosh+jao, iaosh+iao) * dsv)

                        ! Accumulate CN derivatives
                        dcni = dcni + dhdcni * pij * block_overlap(jaosh+jao, iaosh+iao)
                        dcnj = dcnj + dhdcnj * pij * block_overlap(jaosh+jao, iaosh+iao)
                     end do
                  end do

               end do
            end do
         end if

         if (iat == jat) then
            dEdcn(iat) = dEdcn(iat) + 0.5_wp * (dcni + dcnj)
         else
            dEdcn(iat) = dEdcn(iat) + dcni
            dEdcn(jat) = dEdcn(jat) + dcnj
         end if
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
         if (iat == jat) then
            sigma(:, :) = sigma + 0.25_wp * (spread(vec, 1, 3) &
               & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
         else
            sigma(:, :) = sigma + 0.5_wp * (spread(vec, 1, 3) &
               & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
         end if

      end do

      ! Onsite contributions
      if (.not.owns_index(partition, iat)) cycle
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         dhdcni = dsedcn(is+ish)
         dcni = 0.0_wp
         do iao = 1, msao(bas%cgto(ish, izp)%ang)
            dcni = dcni + dhdcni * pmat(ii+iao, ii+iao, 1)
         end do
         dEdcn(iat) = dEdcn(iat) + dcni
      end do
   end do

end subroutine get_hamiltonian_gradient


!> Add the explicit Cartesian Hessian of the GFN1 zero-order Hamiltonian
!> contribution for a fixed electronic state and fixed coordination numbers.
!>
!> This routine differentiates the H0/Pulay part of get_hamiltonian_gradient
!> while keeping pmat, xmat and selfenergy fixed.  The CN chain-rule terms
!> (dE/dCN * d2CN/dR2 and the mixed d2E/dR dCN terms) and the SCF response
!> (dP/dR, dX/dR) are intentionally not included here and are assembled by
!> later Hessian stages.  GFN1 does not use diatomic-frame overlap scaling;
!> Hamiltonians that enable it are currently outside this routine's scope.
subroutine get_hamiltonian_hessian(mol, trans, list, bas, h0, selfenergy, &
      & pmat, xmat, hessian, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Fixed diagonal elements of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Fixed density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Fixed energy-weighted density matrix
   real(wp), intent(in) :: xmat(:, :, :)
   !> Cartesian Hessian, flattened as (3*nat, 3*nat)
   real(wp), intent(inout) :: hessian(:, :)
   !> Share of the neighbour list entries evaluated here
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, itr, img, inl
   integer :: ish, jsh, is, js, nsi, nsj, ii, jj, iao, jao, naoi, naoj
   integer :: ic, jc, irow, jcol
   real(wp) :: r2, vec(3), rr, r0, hscale, shpolyi, shpolyj, shpoly
   real(wp) :: api, apj, fp1, fp2, arr, hij, pij, xij, sval, pref
   real(wp) :: grr(3), hrr(3, 3), gsp(3), hsp(3, 3), gh(3), hh(3, 3)
   real(wp) :: block(3, 3)
   real(wp), allocatable :: stmp(:, :), dstmp(:, :, :), hstmp(:, :, :, :)
   integer :: mmax

   ! GFN1 uses the ordinary overlap Hamiltonian.  The diatomic-frame scaled
   ! branch (used by other Hamiltonians) requires derivatives of that
   ! transformation and is deliberately deferred.
   if (h0%do_diat_scale) return

   mmax = msao(bas%maxl)
   allocate(stmp(mmax, mmax), dstmp(3, mmax, mmax), &
      & hstmp(3, 3, mmax, mmax))

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         if (.not.owns_index(partition, img+inl)) cycle
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = dot_product(vec, vec)
         if (r2 <= epsilon(1.0_wp)) cycle

         r0 = h0%rad(jzp) + h0%rad(izp)
         rr = sqrt(sqrt(r2)/r0)
         arr = 0.5_wp*rr/r2
         grr(:) = arr*vec
         hrr(:, :) = -1.5_wp*arr/r2 * spread(vec, 2, 3)*spread(vec, 1, 3)
         do ic = 1, 3
            hrr(ic, ic) = hrr(ic, ic) + arr
         end do

         block(:, :) = 0.0_wp
         do ish = 1, nsi
            ii = bas%iao_sh(is+ish)
            api = h0%shpoly(ish, izp)
            shpolyi = 1.0_wp + api*rr
            naoi = msao(bas%cgto(ish, izp)%ang)
            do jsh = 1, nsj
               jj = bas%iao_sh(js+jsh)
               apj = h0%shpoly(jsh, jzp)
               shpolyj = 1.0_wp + apj*rr
               shpoly = shpolyi*shpolyj
               fp1 = api + apj + 2.0_wp*api*apj*rr
               fp2 = 2.0_wp*api*apj
               gsp(:) = fp1*grr
               hsp(:, :) = fp2*spread(grr, 2, 3)*spread(grr, 1, 3) + fp1*hrr

               hscale = h0%hscale(jsh, ish, jzp, izp)
               pref = 0.5_wp*(selfenergy(is+ish) + selfenergy(js+jsh))*hscale
               hij = pref*shpoly
               gh(:) = pref*gsp
               hh(:, :) = pref*hsp

               naoj = msao(bas%cgto(jsh, jzp)%ang)
               call overlap_hess_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp(1:naoj, 1:naoi), &
                  & dstmp(:, 1:naoj, 1:naoi), hstmp(:, :, 1:naoj, 1:naoi))

               do iao = 1, naoi
                  do jao = 1, naoj
                     pij = pmat(jj+jao, ii+iao, 1)
                     xij = xmat(jj+jao, ii+iao, 1)
                     sval = 2.0_wp*pij*hij - 2.0_wp*xij
                     block(:, :) = block + sval*hstmp(:, :, jao, iao) &
                        & + 2.0_wp*pij*( &
                        & spread(gh, 2, 3)*spread(dstmp(:, jao, iao), 1, 3) &
                        & + spread(dstmp(:, jao, iao), 2, 3)*spread(gh, 1, 3) &
                        & + stmp(jao, iao)*hh)
                  end do
               end do
            end do
         end do

         do ic = 1, 3
            irow = 3*(iat - 1) + ic
            do jc = 1, 3
               jcol = 3*(jat - 1) + jc
               hessian(irow, 3*(iat - 1) + jc) = &
                  & hessian(irow, 3*(iat - 1) + jc) + block(ic, jc)
               hessian(3*(jat - 1) + ic, jcol) = &
                  & hessian(3*(jat - 1) + ic, jcol) + block(ic, jc)
               hessian(irow, jcol) = hessian(irow, jcol) - block(ic, jc)
               hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) = &
                  & hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) - block(ic, jc)
            end do
         end do
      end do
   end do

end subroutine get_hamiltonian_hessian


!> Evaluate the CN derivative of the fixed-density H0 energy and its
!> mixed Cartesian/CN derivative for GFN1.
!>
!> dEdcndr(:, B, A) = d/dR_B [dE/dCN_A].  Together with dcndr this
!> supplies the two mixed chain-rule terms in the Cartesian Hessian.
subroutine get_hamiltonian_cn_derivatives(mol, trans, list, bas, h0, dsedcn, &
      & pmat, dEdcn, dEdcndr, partition)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   type(adjacency_list), intent(in) :: list
   type(basis_type), intent(in) :: bas
   type(tb_hamiltonian), intent(in) :: h0
   real(wp), intent(in) :: dsedcn(:)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(inout) :: dEdcn(:)
   real(wp), intent(inout) :: dEdcndr(:, :, :)
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, itr, img, inl
   integer :: ish, jsh, is, js, nsi, nsj, ii, jj, iao, jao, naoi, naoj, ic
   real(wp) :: r2, vec(3), rr, r0, api, apj, shpolyi, shpolyj, shpoly
   real(wp) :: fp1, arr, hscale, pij, dcni, dcnj
   real(wp) :: grr(3), gsp(3), gci(3), gcj(3)
   real(wp), allocatable :: stmp(:, :), dstmp(:, :, :)
   integer :: mmax

   if (h0%do_diat_scale) return

   mmax = msao(bas%maxl)
   allocate(stmp(mmax, mmax), dstmp(3, mmax, mmax))

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         if (.not.owns_index(partition, img+inl)) cycle
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = dot_product(vec, vec)
         if (r2 <= epsilon(1.0_wp)) cycle

         r0 = h0%rad(jzp) + h0%rad(izp)
         rr = sqrt(sqrt(r2)/r0)
         arr = 0.5_wp*rr/r2
         grr(:) = arr*vec
         dcni = 0.0_wp
         dcnj = 0.0_wp
         gci(:) = 0.0_wp
         gcj(:) = 0.0_wp

         do ish = 1, nsi
            ii = bas%iao_sh(is+ish)
            api = h0%shpoly(ish, izp)
            shpolyi = 1.0_wp + api*rr
            naoi = msao(bas%cgto(ish, izp)%ang)
            do jsh = 1, nsj
               jj = bas%iao_sh(js+jsh)
               apj = h0%shpoly(jsh, jzp)
               shpolyj = 1.0_wp + apj*rr
               shpoly = shpolyi*shpolyj
               fp1 = api + apj + 2.0_wp*api*apj*rr
               gsp(:) = fp1*grr
               hscale = h0%hscale(jsh, ish, jzp, izp)

               naoj = msao(bas%cgto(jsh, jzp)%ang)
               call overlap_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp(1:naoj, 1:naoi), &
                  & dstmp(:, 1:naoj, 1:naoi))

               do iao = 1, naoi
                  do jao = 1, naoj
                     pij = pmat(jj+jao, ii+iao, 1)
                     dcni = dcni + dsedcn(is+ish)*hscale*pij*shpoly*stmp(jao, iao)
                     dcnj = dcnj + dsedcn(js+jsh)*hscale*pij*shpoly*stmp(jao, iao)
                     gci(:) = gci + dsedcn(is+ish)*hscale*pij*( &
                        & shpoly*dstmp(:, jao, iao) + stmp(jao, iao)*gsp)
                     gcj(:) = gcj + dsedcn(js+jsh)*hscale*pij*( &
                        & shpoly*dstmp(:, jao, iao) + stmp(jao, iao)*gsp)
                  end do
               end do
            end do
         end do

         if (iat == jat) then
            dEdcn(iat) = dEdcn(iat) + 0.5_wp*(dcni+dcnj)
            dEdcndr(:, iat, iat) = dEdcndr(:, iat, iat) + 0.5_wp*(gci+gcj)
         else
            dEdcn(iat) = dEdcn(iat) + dcni
            dEdcn(jat) = dEdcn(jat) + dcnj
            dEdcndr(:, iat, iat) = dEdcndr(:, iat, iat) + gci
            dEdcndr(:, jat, iat) = dEdcndr(:, jat, iat) - gci
            dEdcndr(:, iat, jat) = dEdcndr(:, iat, jat) + gcj
            dEdcndr(:, jat, jat) = dEdcndr(:, jat, jat) - gcj
         end if
      end do

      ! Onsite CN derivative.  It has no explicit Cartesian derivative.
      if (.not.owns_index(partition, iat)) cycle
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         do iao = 1, msao(bas%cgto(ish, izp)%ang)
            dEdcn(iat) = dEdcn(iat) + dsedcn(is+ish)*pmat(ii+iao, ii+iao, 1)
         end do
      end do
   end do
end subroutine get_hamiltonian_cn_derivatives


subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Occupation number
   real(wp), intent(out) :: nocc
   !> Reference occupation for each atom
   real(wp), intent(out) :: n0at(:)
   !> Reference occupation for each shell
   real(wp), intent(out) :: n0sh(:)

   integer :: iat, ish, izp, ii

   nocc = -mol%charge
   n0at(:) = 0.0_wp
   n0sh(:) = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(izp)
         nocc = nocc + h0%refocc(ish, izp)
         n0at(iat) = n0at(iat) + h0%refocc(ish, izp)
         n0sh(ii+ish) = n0sh(ii+ish) + h0%refocc(ish, izp)
      end do
   end do

end subroutine get_occupation


!> Shift multipole operator from Ket function (center i) to Bra function (center j),
!> the multipole operator on the Bra function can be assembled from the lower moments
!> on the Ket function and the displacement vector using horizontal shift rules.
pure subroutine shift_operator(vec, s, di, qi, dj, qj)
   !> Displacement vector of center i and j
   real(wp),intent(in) :: vec(:)
   !> Overlap integral between basis functions
   real(wp),intent(in) :: s
   !> Dipole integral with operator on Ket function (center i)
   real(wp),intent(in) :: di(:)
   !> Quadrupole integral with operator on Ket function (center i)
   real(wp),intent(in) :: qi(:)
   !> Dipole integral with operator on Bra function (center j)
   real(wp),intent(out) :: dj(:)
   !> Quadrupole integral with operator on Bra function (center j)
   real(wp),intent(out) :: qj(:)

   real(wp) :: tr

   ! Create dipole operator on Bra function from Ket function and shift contribution
   ! due to monopol displacement
   dj(1) = di(1) + vec(1)*s
   dj(2) = di(2) + vec(2)*s
   dj(3) = di(3) + vec(3)*s

   ! For the quadrupole operator on the Bra function we first construct the shift
   ! contribution from the dipole and monopol displacement, since we have to remove
   ! the trace contribution from the shift and the moment integral on the Ket function
   ! is already traceless
   qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
   qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
   qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
   qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
   qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
   qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
   ! Now collect the trace of the shift contribution
   tr = 0.5_wp * (qj(1) + qj(3) + qj(6))

   ! Finally, assemble the quadrupole operator on the Bra function from the operator
   ! on the Ket function and the traceless shift contribution
   qj(1) = qi(1) + 1.5_wp * qj(1) - tr
   qj(2) = qi(2) + 1.5_wp * qj(2)
   qj(3) = qi(3) + 1.5_wp * qj(3) - tr
   qj(4) = qi(4) + 1.5_wp * qj(4)
   qj(5) = qi(5) + 1.5_wp * qj(5)
   qj(6) = qi(6) + 1.5_wp * qj(6) - tr
end subroutine shift_operator


end module tblite_xtb_h0
