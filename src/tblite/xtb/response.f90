! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Analytical first-order electronic response for GFN1-xTB.
module tblite_xtb_response
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_container, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_lapack, only : getrf, getrs
   use tblite_scf_potential, only : potential_type, new_potential, add_pot_to_h1
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, &
      & get_hamiltonian_derivative, get_hamiltonian_gradient
   implicit none
   private

   public :: get_gfn1_response, add_gfn1_response_hessian

   interface get_gfn1_response
      module procedure :: get_gfn1_response_closed
      module procedure :: get_gfn1_response_spin
   end interface get_gfn1_response

   !> Reusable workspace for one orbital spin channel.
   type :: response_workspace
      integer :: ispin = 1
      real(wp), allocatable :: c(:, :), hmo(:, :), smo(:, :), u(:, :)
      real(wp), allocatable :: tmp(:, :), cw(:, :), ce(:, :), cd(:, :)
      real(wp), allocatable :: p0(:, :), x0(:, :), h0(:, :), sinv(:, :), bmat(:, :)
      real(wp), allocatable :: focc(:), weight(:), eps(:), deps(:)
   end type response_workspace

contains

!> Backward-compatible closed-shell response interface.
subroutine get_gfn1_response_closed(mol, calc, wfn, accuracy, dpdr, dxdr, error)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn
   real(wp), intent(in) :: accuracy
   real(wp), intent(out) :: dpdr(:, :, :, :)
   real(wp), intent(out) :: dxdr(:, :, :, :)
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: dp5(:, :, :, :, :), dx5(:, :, :, :, :)

   if (wfn%nspin /= 1) then
      call fatal_error(error, "Rank-4 GFN1 response output is only valid for one spin channel")
      return
   end if
   allocate(dp5(size(dpdr,1),size(dpdr,2),1,size(dpdr,3),size(dpdr,4)))
   allocate(dx5(size(dxdr,1),size(dxdr,2),1,size(dxdr,3),size(dxdr,4)))
   call get_gfn1_response_spin(mol, calc, wfn, accuracy, dp5, dx5, error)
   if (allocated(error)) return
   dpdr = dp5(:,:,1,:,:)
   dxdr = dx5(:,:,1,:,:)
end subroutine get_gfn1_response_closed


!> Solve the GFN1 response and accumulate the electronic P'/X' contribution
!> directly into a Cartesian Hessian without retaining full response tensors.
subroutine add_gfn1_response_hessian(mol, calc, wfn, accuracy, hessian, error)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn
   real(wp), intent(in) :: accuracy
   real(wp), intent(inout) :: hessian(:, :)
   type(error_type), allocatable, intent(out) :: error

   call get_gfn1_response_spin(mol, calc, wfn, accuracy, error=error, &
      & hessian_response=hessian)
end subroutine add_gfn1_response_hessian


!> Compute the first-order GFN1-xTB response for closed- or open-shell systems.
!>
!> For unrestricted calculations the fixed-point variables are shell charge and
!> shell magnetization, giving a 2*nsh CP space.  Orbital quantities are treated
!> in alpha/beta representation, while qsh/dqsh and the returned density
!> responses are stored in charge/magnetization representation to match the
!> gradient implementation.
subroutine get_gfn1_response_spin(mol, calc, wfn, accuracy, dpdr, dxdr, error, hessian_response)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn
   real(wp), intent(in) :: accuracy
   !> Density derivative [nao, nao, spin(charge/magnet), 3, nat]
   real(wp), intent(out), optional :: dpdr(:, :, :, :, :)
   !> Energy-weighted density derivative [nao, nao, spin(charge/magnet), 3, nat]
   real(wp), intent(out), optional :: dxdr(:, :, :, :, :)
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(inout), optional :: hessian_response(:, :)

   integer :: nat, nsh, nao, nspin, nvar, ncoord
   integer :: iat, ic, k, ish, iao, jao, info, ivar, jspin, ispin
   integer, allocatable :: ipiv(:)
   real(wp) :: cutoff, cp_residual, cp_scale
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: selfenergy(:), dsedcn(:), lattr(:, :)
   real(wp), allocatable :: ds(:, :, :, :), dh0(:, :, :, :)
   real(wp), allocatable :: kernel(:, :), jac(:, :), amat(:, :)
   real(wp), allocatable :: rhs(:, :), rhs0(:, :), dvresp(:, :)
   real(wp), allocatable :: dhspin(:, :, :), dhtot(:, :, :)
   real(wp), allocatable :: dpab(:, :, :), dxab(:, :, :), dpcm(:, :, :), dxcm(:, :, :)
   real(wp), allocatable :: dqraw(:, :), dqcm(:, :), dqcheck(:, :)
   real(wp), allocatable :: dEdcn_resp(:), gradtmp(:, :), sigma_resp(:, :)
   real(wp), allocatable :: dvao(:, :, :, :), h1dummy(:, :, :)
   type(integral_type) :: ints
   type(adjacency_list) :: list
   type(container_cache) :: ccache, icache
   type(potential_type) :: pot0, potg
   type(wavefunction_type) :: rwfn
   type(response_workspace), allocatable :: work(:)

   nat = mol%nat
   nsh = calc%bas%nsh
   nao = calc%bas%nao
   nspin = wfn%nspin
   nvar = nsh*nspin
   ncoord = 3*nat

   if (any(mol%periodic)) then
      call fatal_error(error, "GFN1 analytic response currently supports non-periodic systems only")
      return
   end if
   if (nspin < 1 .or. nspin > 2) then
      call fatal_error(error, "GFN1 analytic response supports one or two spin channels")
      return
   end if
   if (.not.allocated(calc%coulomb) .or. .not.allocated(calc%coulomb%es2)) then
      call fatal_error(error, "GFN1 analytic response requires the SCC Coulomb model")
      return
   end if
   if (allocated(calc%coulomb%aes2)) then
      call fatal_error(error, "GFN1 analytic response does not support anisotropic SCC terms")
      return
   end if
   if (present(dpdr) .neqv. present(dxdr)) then
      call fatal_error(error, "GFN1 response requires both density response outputs or neither")
      return
   end if
   if (present(dpdr)) then
      if (size(dpdr,1) /= nao .or. size(dpdr,2) /= nao .or. &
          size(dpdr,3) /= nspin .or. size(dpdr,4) /= 3 .or. size(dpdr,5) /= nat .or. &
          any(shape(dxdr) /= shape(dpdr))) then
         call fatal_error(error, "Invalid GFN1 response output dimensions")
         return
      end if
   end if
   if (present(hessian_response)) then
      if (size(hessian_response,1) /= 3*nat .or. size(hessian_response,2) /= 3*nat) then
         call fatal_error(error, "Invalid GFN1 response Hessian dimensions")
         return
      end if
   end if

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

   allocate(work(nspin))
   do ispin = 1, nspin
      call init_response_workspace(calc, wfn, ispin, ints%overlap, work(ispin))
   end do

   allocate(ds(nao,nao,3,nat), dh0(nao,nao,3,nat))
   call get_hamiltonian_derivative(mol, lattr, list, calc%bas, calc%h0, &
      & selfenergy, dsedcn, dcndr, ds, dh0)

   ! Central potential in charge/magnetization representation.
   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%interactions)) call calc%interactions%update(mol, icache)
   call new_potential(pot0, mol, calc%bas, nspin)
   call pot0%reset
   call calc%coulomb%get_potential(mol, ccache, wfn, pot0)
   if (allocated(calc%interactions)) call calc%interactions%get_potential(mol, icache, wfn, pot0)
   allocate(h1dummy(nao,nao,nspin))
   call add_pot_to_h1(calc%bas, ints, pot0, h1dummy)

   ! Explicit potential derivative at fixed electronic populations.
   rwfn = wfn
   call ensure_charge_derivative_storage(rwfn, nat, nsh)
   rwfn%dqshdr = 0.0_wp
   rwfn%dqatdr = 0.0_wp
   rwfn%dqshdL = 0.0_wp
   rwfn%dqatdL = 0.0_wp
   call new_potential(potg, mol, calc%bas, nspin, grad=.true.)
   call potg%reset
   call calc%coulomb%get_potential_gradient(mol, ccache, rwfn, potg)
   if (allocated(calc%interactions)) call calc%interactions%get_potential_gradient(mol, icache, rwfn, potg)

   allocate(dvao(nao,3,nat,nspin), source=0.0_wp)
   do jspin = 1, nspin
      do iao = 1, nao
         ish = calc%bas%ao2sh(iao)
         iat = calc%bas%ao2at(iao)
         do k = 1, nat
            do ic = 1, 3
               dvao(iao,ic,k,jspin) = potg%dvshdr(ic,k,ish,jspin) &
                  & + potg%dvatdr(ic,k,iat,jspin)
            end do
         end do
      end do
   end do

   ! Full response kernel in charge/magnetization shell space.
   allocate(kernel(nvar,nvar), source=0.0_wp)
   call build_response_kernel(mol, calc, wfn, ccache, icache, kernel)

   allocate(jac(nvar,nvar), amat(nvar,nvar), rhs(nvar,ncoord))
   allocate(dhspin(nao,nao,nspin), dhtot(nao,nao,nspin))
   allocate(dpab(nao,nao,nspin), dxab(nao,nao,nspin), &
      & dpcm(nao,nao,nspin), dxcm(nao,nao,nspin))
   if (present(hessian_response)) then
      allocate(dEdcn_resp(nat), gradtmp(3,nat), sigma_resp(3,3))
   end if
   allocate(dqraw(nsh,nspin), dqcm(nsh,nspin), dqcheck(nsh,nspin))
   ! Fixed-point Jacobian d(q,m)_out / d(q,m)_in.  Only the Mulliken
   ! population response is needed here, so avoid constructing X' entirely and
   ! skip the overlap-derivative transformation for this zero-S' perturbation.
   do ivar = 1, nvar
      call build_response_hamiltonians(calc, ints%overlap, kernel(:,ivar), nspin, dhspin)
      do ispin = 1, nspin
         call orbital_charge_response_one(calc, wfn, ispin, ints%overlap, &
            & dhspin(:,:,ispin), work(ispin), dqraw(:,ispin), error)
         if (allocated(error)) return
      end do
      call combine_spin_response(nspin, dpab, dxab, dqraw, dqcm=dqcm)
      jac(:,ivar) = reshape(dqcm, [nvar])
   end do

   ! Explicit geometry RHS at fixed input populations.  Again only q' is
   ! required for the CP right-hand side; build one Hamiltonian-derivative
   ! column at a time rather than storing the full O(nao^2*3N) tensor.
   k = 0
   do iat = 1, nat
      do ic = 1, 3
         k = k + 1
         call build_explicit_hamiltonian_derivative_column(calc, ints%overlap, &
            & dh0, ds, pot0%vao, dvao, nspin, ic, iat, dhtot)
         do ispin = 1, nspin
            call orbital_charge_response_one(calc, wfn, ispin, ints%overlap, &
               & dhtot(:,:,ispin), work(ispin), dqraw(:,ispin), error, ds(:,:,ic,iat))
            if (allocated(error)) return
         end do
         call combine_spin_response(nspin, dpab, dxab, dqraw, dqcm=dqcm)
         rhs(:,k) = reshape(dqcm, [nvar])
      end do
   end do

   allocate(rhs0(nvar,ncoord))
   rhs0 = rhs
   amat = -jac
   do ivar = 1, nvar
      amat(ivar,ivar) = amat(ivar,ivar) + 1.0_wp
   end do
   allocate(ipiv(nvar))
   call getrf(amat, ipiv, info)
   if (info /= 0) then
      call fatal_error(error, "GFN1 CP-xTB response matrix is singular")
      return
   end if
   call getrs(amat, rhs, ipiv, info)
   if (info /= 0) then
      call fatal_error(error, "Failed to solve the GFN1 CP-xTB response equations")
      return
   end if

   cp_scale = max(1.0_wp, maxval(abs(rhs0)))
   cp_residual = maxval(abs(rhs - matmul(jac, rhs) - rhs0))/cp_scale
   if (cp_residual > 1.0e-9_wp) then
      call fatal_error(error, "GFN1 CP-xTB linear solve residual is too large")
      return
   end if

   ! Persist charge/magnetization derivatives.
   call ensure_charge_derivative_storage(wfn, nat, nsh)
   wfn%dqshdr = 0.0_wp
   wfn%dqatdr = 0.0_wp
   k = 0
   do iat = 1, nat
      do ic = 1, 3
         k = k + 1
         wfn%dqshdr(ic,iat,:,:) = reshape(rhs(:,k), [nsh,nspin])
      end do
   end do
   do ish = 1, nsh
      iat = calc%bas%sh2at(ish)
      wfn%dqatdr(:,:,iat,:) = wfn%dqatdr(:,:,iat,:) + wfn%dqshdr(:,:,ish,:)
   end do

   ! Reconstruct final self-consistent density responses.  They can either be
   ! stored for the public response API or streamed directly to a Hessian
   ! Hessian accumulator to keep peak memory independent of the number of coordinates.
   if (present(dpdr)) dpdr = 0.0_wp
   if (present(dxdr)) dxdr = 0.0_wp
   allocate(dvresp(nvar,ncoord))
   dvresp = matmul(kernel, rhs)
   k = 0
   do iat = 1, nat
      do ic = 1, 3
         k = k + 1
         call build_response_hamiltonians(calc, ints%overlap, dvresp(:,k), nspin, dhspin)
         call build_explicit_hamiltonian_derivative_column(calc, ints%overlap, &
            & dh0, ds, pot0%vao, dvao, nspin, ic, iat, dhtot)
         do ispin = 1, nspin
            dhtot(:,:,ispin) = dhtot(:,:,ispin) + dhspin(:,:,ispin)
            call orbital_response_one(calc, wfn, ispin, ints%overlap, dhtot(:,:,ispin), &
               & ds(:,:,ic,iat), work(ispin), dpab(:,:,ispin), dxab(:,:,ispin), &
               & dqraw(:,ispin), error)
            if (allocated(error)) return
         end do
         call combine_spin_response(nspin, dpab, dxab, dqraw, dpcm, dxcm, dqcheck)
         if (maxval(abs(reshape(dqcheck,[nvar])-rhs(:,k))) > 7.5e-7_wp) then
            call fatal_error(error, "GFN1 CP-xTB response residual is too large")
            return
         end if
         if (present(dpdr)) then
            dpdr(:,:,:,ic,iat) = dpcm
            dxdr(:,:,:,ic,iat) = dxcm
         end if
         if (present(hessian_response)) then
            dEdcn_resp = 0.0_wp
            gradtmp = 0.0_wp
            sigma_resp = 0.0_wp
            call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, &
               & selfenergy, dsedcn, pot0, dpcm, dxcm, dEdcn_resp, gradtmp, &
               & sigma_resp, calc%partition)
            hessian_response(:,k) = hessian_response(:,k) + reshape(gradtmp, [3*nat])
            do jao = 1, nat
               do jspin = 1, 3
                  hessian_response(3*(jao-1)+jspin,k) = &
                     hessian_response(3*(jao-1)+jspin,k) &
                     & + sum(dcndr(jspin,jao,:)*dEdcn_resp(:))
               end do
            end do
         end if
      end do
   end do

end subroutine get_gfn1_response_spin


!> Build the shell-potential response kernel in charge/magnetization space.
subroutine build_response_kernel(mol, calc, wfn, ccache, icache, kernel)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   type(container_cache), intent(inout) :: ccache, icache
   real(wp), intent(out) :: kernel(:, :)

   integer :: nsh, nspin, nvar, ish, jsh, iat, jat, izp, jin, sin, sout
   type(wavefunction_type) :: twfn
   type(potential_type) :: pot
   real(wp), allocatable :: charge_kernel(:, :)

   nsh = calc%bas%nsh
   nspin = wfn%nspin
   nvar = nsh*nspin
   kernel = 0.0_wp

   allocate(charge_kernel(nsh,nsh), source=0.0_wp)
   call build_charge_kernel(mol, calc, wfn, ccache, charge_kernel)
   kernel(1:nsh,1:nsh) = charge_kernel

   if (.not.allocated(calc%interactions)) return

   twfn = wfn
   call new_potential(pot, mol, calc%bas, nspin)
   do sin = 1, nspin
      do jsh = 1, nsh
         twfn%qsh = 0.0_wp
         twfn%qat = 0.0_wp
         twfn%qsh(jsh,sin) = 1.0_wp
         jat = calc%bas%sh2at(jsh)
         twfn%qat(jat,sin) = 1.0_wp
         call pot%reset
         call calc%interactions%get_potential(mol, icache, twfn, pot)
         jin = (sin-1)*nsh + jsh
         do sout = 1, nspin
            do ish = 1, nsh
               iat = calc%bas%sh2at(ish)
               kernel((sout-1)*nsh+ish,jin) = kernel((sout-1)*nsh+ish,jin) &
                  & + pot%vsh(ish,sout) + pot%vat(iat,sout)
            end do
         end do
      end do
   end do

end subroutine build_response_kernel


!> Charge-channel SCC response kernel dV_s/dq_t, including SCC3.
subroutine build_charge_kernel(mol, calc, wfn, ccache, kernel)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   type(container_cache), intent(inout) :: ccache
   real(wp), intent(out) :: kernel(:, :)

   integer :: nsh, ish, jsh, iat, jat, izp
   type(wavefunction_type) :: twfn
   type(potential_type) :: pot

   nsh = calc%bas%nsh
   kernel = 0.0_wp
   twfn = wfn
   call new_potential(pot, mol, calc%bas, wfn%nspin)

   do jsh = 1, nsh
      twfn%qsh = 0.0_wp
      twfn%qat = 0.0_wp
      twfn%qsh(jsh,1) = 1.0_wp
      jat = calc%bas%sh2at(jsh)
      twfn%qat(jat,1) = 1.0_wp
      call pot%reset
      call calc%coulomb%es2%get_potential(mol, ccache, twfn, pot)
      kernel(:,jsh) = pot%vsh(:,1)
   end do

   if (allocated(calc%coulomb%es3)) then
      if (calc%coulomb%es3%shell_resolved) then
         do jsh = 1, nsh
            jat = calc%bas%sh2at(jsh)
            izp = mol%id(jat)
            ish = jsh - calc%coulomb%es3%ish_at(jat)
            kernel(jsh,jsh) = kernel(jsh,jsh) &
               & + 2.0_wp*wfn%qsh(jsh,1) &
               & * calc%coulomb%es3%hubbard_derivs(ish,izp)
         end do
      else
         do jsh = 1, nsh
            jat = calc%bas%sh2at(jsh)
            izp = mol%id(jat)
            do ish = 1, nsh
               iat = calc%bas%sh2at(ish)
               if (iat == jat) then
                  kernel(ish,jsh) = kernel(ish,jsh) &
                     & + 2.0_wp*wfn%qat(iat,1) &
                     & * calc%coulomb%es3%hubbard_derivs(1,izp)
               end if
            end do
         end do
      end if
   end if
end subroutine build_charge_kernel


!> Convert one charge/magnetization shell-potential vector into the actual
!> alpha/beta Hamiltonian perturbations solved by the unrestricted SCF.
subroutine build_response_hamiltonians(calc, overlap, dvflat, nspin, dh)
   type(xtb_calculator), intent(in) :: calc
   real(wp), intent(in) :: overlap(:, :)
   real(wp), intent(in) :: dvflat(:)
   integer, intent(in) :: nspin
   real(wp), intent(out) :: dh(:, :, :)

   integer :: iao, jao, is, js, ispin, nsh
   real(wp) :: vi, vj

   nsh = calc%bas%nsh
   dh = 0.0_wp
   do ispin = 1, nspin
      do iao = 1, calc%bas%nao
         is = calc%bas%ao2sh(iao)
         do jao = 1, calc%bas%nao
            js = calc%bas%ao2sh(jao)
            if (nspin == 1) then
               vi = dvflat(is)
               vj = dvflat(js)
            else if (ispin == 1) then
               vi = dvflat(is) + dvflat(nsh+is)
               vj = dvflat(js) + dvflat(nsh+js)
            else
               vi = dvflat(is) - dvflat(nsh+is)
               vj = dvflat(js) - dvflat(nsh+js)
            end if
            dh(jao,iao,ispin) = -0.5_wp*overlap(jao,iao)*(vi+vj)
         end do
      end do
   end do
end subroutine build_response_hamiltonians


!> Build explicit alpha/beta H'(R) at fixed charge/magnetization populations.
subroutine build_explicit_hamiltonian_derivative_column(calc, overlap, dh0, ds, vao, dvao, &
      & nspin, ic, iat, dh)
   type(xtb_calculator), intent(in) :: calc
   real(wp), intent(in) :: overlap(:, :), dh0(:, :, :, :), ds(:, :, :, :)
   real(wp), intent(in) :: vao(:, :), dvao(:, :, :, :)
   integer, intent(in) :: nspin, ic, iat
   real(wp), intent(out) :: dh(:, :, :)

   integer :: iao, jao, ispin
   real(wp) :: vi, vj, dvi, dvj

   do ispin = 1, nspin
      do iao = 1, calc%bas%nao
         do jao = 1, calc%bas%nao
            if (nspin == 1) then
               vi = vao(iao,1); vj = vao(jao,1)
               dvi = dvao(iao,ic,iat,1); dvj = dvao(jao,ic,iat,1)
            else if (ispin == 1) then
               vi = vao(iao,1)+vao(iao,2); vj = vao(jao,1)+vao(jao,2)
               dvi = dvao(iao,ic,iat,1)+dvao(iao,ic,iat,2)
               dvj = dvao(jao,ic,iat,1)+dvao(jao,ic,iat,2)
            else
               vi = vao(iao,1)-vao(iao,2); vj = vao(jao,1)-vao(jao,2)
               dvi = dvao(iao,ic,iat,1)-dvao(iao,ic,iat,2)
               dvj = dvao(jao,ic,iat,1)-dvao(jao,ic,iat,2)
            end if
            dh(jao,iao,ispin) = dh0(jao,iao,ic,iat) &
               & - 0.5_wp*ds(jao,iao,ic,iat)*(vj+vi) &
               & - 0.5_wp*overlap(jao,iao)*(dvj+dvi)
         end do
      end do
   end do
end subroutine build_explicit_hamiltonian_derivative_column


subroutine init_response_workspace(calc, wfn, ispin, overlap, work)
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: ispin
   real(wp), intent(in) :: overlap(:, :)
   type(response_workspace), intent(out) :: work

   integer :: nao, n

   nao = calc%bas%nao
   work%ispin = ispin
   allocate(work%c(nao,nao), work%hmo(nao,nao), work%smo(nao,nao), &
      & work%u(nao,nao), work%tmp(nao,nao), &
      & work%cw(nao,nao), work%ce(nao,nao), work%cd(nao,nao), &
      & work%p0(nao,nao), work%x0(nao,nao), work%h0(nao,nao), &
      & work%sinv(nao,nao), work%bmat(nao,nao))
   allocate(work%focc(nao), work%weight(nao), work%eps(nao), work%deps(nao))

   work%c = wfn%coeff(:,:,ispin)
   work%eps = wfn%emo(:,ispin)
   if (wfn%nspin == 1) then
      work%focc = wfn%focc(:,1) + wfn%focc(:,2)
   else
      work%focc = wfn%focc(:,ispin)
   end if
   work%weight = 0.0_wp
   if (wfn%kt > 0.0_wp) then
      if (wfn%nspin == 1) then
         work%weight = wfn%focc(:,1)*(1.0_wp-wfn%focc(:,1))/wfn%kt &
            & + wfn%focc(:,2)*(1.0_wp-wfn%focc(:,2))/wfn%kt
      else
         work%weight = wfn%focc(:,ispin)*(1.0_wp-wfn%focc(:,ispin))/wfn%kt
      end if
   end if
   do n = 1, nao
      work%cw(:,n) = work%c(:,n)*work%focc(n)
      work%ce(:,n) = work%c(:,n)*work%focc(n)*work%eps(n)
   end do

   ! Central AO matrices used by the gauge-invariant X response.  The
   ! generalized eigenvectors satisfy C^T S C = I, hence S^{-1}=C C^T and
   ! H = S C eps C^T S.
   work%p0 = matmul(work%cw, transpose(work%c))
   work%x0 = matmul(work%ce, transpose(work%c))
   work%sinv = matmul(work%c, transpose(work%c))
   work%tmp = matmul(overlap, work%c)
   do n = 1, nao
      work%cd(:,n) = work%tmp(:,n)*work%eps(n)
   end do
   work%h0 = matmul(work%cd, transpose(work%tmp))
end subroutine init_response_workspace


!> Analytic response of one alpha/beta orbital channel.
!> Build the gauge-invariant MO density response matrix for one spin channel.
!> If ds is absent, the perturbation is assumed to have zero overlap derivative;
!> this fast path is used while constructing the CP fixed-point Jacobian.
subroutine build_mo_density_response(wfn, dh, work, error, ds)
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(in) :: dh(:, :)
   type(response_workspace), intent(inout) :: work
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(in), optional :: ds(:, :)

   integer :: nao, m, n
   real(wp) :: den, tol, fm, fn, dmu, denom, kval, ldiv, occ_tol

   nao = size(dh,1)
   work%tmp = matmul(dh, work%c)
   work%hmo = matmul(transpose(work%c), work%tmp)
   if (present(ds)) then
      work%tmp = matmul(ds, work%c)
      work%smo = matmul(transpose(work%c), work%tmp)
   else
      work%smo = 0.0_wp
   end if
   work%u = 0.0_wp
   tol = 1000.0_wp*epsilon(1.0_wp)
   occ_tol = 100.0_wp*epsilon(1.0_wp)

   do n = 1, nao
      work%deps(n) = work%hmo(n,n) - work%eps(n)*work%smo(n,n)
   end do

   dmu = 0.0_wp
   denom = sum(work%weight)
   if (wfn%kt > 0.0_wp .and. denom > sqrt(tiny(1.0_wp))) then
      dmu = dot_product(work%weight, work%deps)/denom
   end if

   do n = 1, nao
      fn = work%focc(n)
      do m = 1, nao
         fm = work%focc(m)
         kval = work%hmo(m,n) &
            & - 0.5_wp*(work%eps(m)+work%eps(n))*work%smo(m,n)
         if (m == n) then
            ldiv = -work%weight(n)
            work%u(n,n) = -work%smo(n,n)*fn + ldiv*kval &
               & + work%weight(n)*dmu
         else
            den = work%eps(m)-work%eps(n)
            if (abs(den) > tol*max(1.0_wp,abs(work%eps(m)),abs(work%eps(n)))) then
               ldiv = (fm-fn)/den
            else if (wfn%kt > 0.0_wp) then
               ldiv = -0.5_wp*(work%weight(m)+work%weight(n))
            else if (abs(fm-fn) <= occ_tol) then
               ldiv = 0.0_wp
            else
               call fatal_error(error, &
                  "Zero-temperature occupied/virtual degeneracy in GFN1 response; use finite electronic temperature")
               return
            end if
            work%u(m,n) = -0.5_wp*work%smo(m,n)*(fm+fn) + ldiv*kval
         end if
      end do
   end do
end subroutine build_mo_density_response


!> Compute only the Mulliken shell-population response for a perturbation.
!> This intentionally skips the energy-weighted density response, which is not
!> needed while forming the CP Jacobian or its explicit geometry right-hand side.
subroutine orbital_charge_response_one(calc, wfn, ispin, overlap, dh, work, dq, error, ds)
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: ispin
   real(wp), intent(in) :: overlap(:, :), dh(:, :)
   type(response_workspace), intent(inout) :: work
   real(wp), intent(out) :: dq(:)
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(in), optional :: ds(:, :)

   integer :: nao, ish, ii, iao, jao

   nao = calc%bas%nao
   if (present(ds)) then
      call build_mo_density_response(wfn, dh, work, error, ds)
   else
      call build_mo_density_response(wfn, dh, work, error)
   end if
   if (allocated(error)) return

   ! Reuse bmat as a temporary AO density response to avoid another workspace.
   work%tmp = matmul(work%c, work%u)
   work%bmat = matmul(work%tmp, transpose(work%c))
   work%bmat = 0.5_wp*(work%bmat + transpose(work%bmat))

   dq = 0.0_wp
   do ish = 1, calc%bas%nsh
      ii = calc%bas%iao_sh(ish)
      do iao = 1, calc%bas%nao_sh(ish)
         do jao = 1, nao
            dq(ish) = dq(ish) - work%bmat(jao,ii+iao)*overlap(jao,ii+iao)
            if (present(ds)) then
               dq(ish) = dq(ish) - wfn%density(jao,ii+iao,ispin)*ds(jao,ii+iao)
            end if
         end do
      end do
   end do
end subroutine orbital_charge_response_one


subroutine orbital_response_one(calc, wfn, ispin, overlap, dh, ds, work, dp, dx, dq, error)
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   integer, intent(in) :: ispin
   real(wp), intent(in) :: overlap(:, :), dh(:, :), ds(:, :)
   type(response_workspace), intent(inout) :: work
   real(wp), intent(out) :: dp(:, :), dx(:, :), dq(:)
   type(error_type), allocatable, intent(out) :: error

   integer :: nao, iao, jao, ish, ii

   nao = calc%bas%nao
   call build_mo_density_response(wfn, dh, work, error, ds)
   if (allocated(error)) return

   work%tmp = matmul(work%c, work%u)
   dp = matmul(work%tmp, transpose(work%c))
   dp = 0.5_wp*(dp + transpose(dp))

   ! Obtain the energy-weighted density response from the differentiated
   ! generalized-eigenvalue identity H P = S X.  This is likewise invariant to
   ! rotations within degenerate subspaces:
   !   X' = S^{-1} (H' P + H P' - S' X).
   work%bmat = matmul(dh, work%p0) + matmul(work%h0, dp) &
      & - matmul(ds, work%x0)
   dx = matmul(work%sinv, work%bmat)
   dx = 0.5_wp*(dx + transpose(dx))
   ! Raw alpha/beta Mulliken population derivative.  The caller converts these
   ! to charge/magnetization representation when nspin=2.
   dq = 0.0_wp
   do ish = 1, calc%bas%nsh
      ii = calc%bas%iao_sh(ish)
      do iao = 1, calc%bas%nao_sh(ish)
         do jao = 1, nao
            dq(ish) = dq(ish) &
               & - dp(jao,ii+iao)*overlap(jao,ii+iao) &
               & - wfn%density(jao,ii+iao,ispin)*ds(jao,ii+iao)
         end do
      end do
   end do
end subroutine orbital_response_one


!> Convert alpha/beta response matrices and populations to charge/magnetization.
subroutine combine_spin_response(nspin, dpab, dxab, dqraw, dpcm, dxcm, dqcm)
   integer, intent(in) :: nspin
   real(wp), intent(in) :: dpab(:, :, :), dxab(:, :, :), dqraw(:, :)
   real(wp), intent(out), optional :: dpcm(:, :, :), dxcm(:, :, :), dqcm(:, :)

   if (present(dpcm)) then
      if (nspin == 1) then
         dpcm(:,:,1) = dpab(:,:,1)
      else
         dpcm(:,:,1) = dpab(:,:,1) + dpab(:,:,2)
         dpcm(:,:,2) = dpab(:,:,1) - dpab(:,:,2)
      end if
   end if
   if (present(dxcm)) then
      if (nspin == 1) then
         dxcm(:,:,1) = dxab(:,:,1)
      else
         dxcm(:,:,1) = dxab(:,:,1) + dxab(:,:,2)
         dxcm(:,:,2) = dxab(:,:,1) - dxab(:,:,2)
      end if
   end if
   if (present(dqcm)) then
      if (nspin == 1) then
         dqcm(:,1) = dqraw(:,1)
      else
         dqcm(:,1) = dqraw(:,1) + dqraw(:,2)
         dqcm(:,2) = dqraw(:,1) - dqraw(:,2)
      end if
   end if
end subroutine combine_spin_response


subroutine ensure_charge_derivative_storage(wfn, nat, nsh)
   type(wavefunction_type), intent(inout) :: wfn
   integer, intent(in) :: nat, nsh

   if (.not.allocated(wfn%dqshdr)) allocate(wfn%dqshdr(3,nat,nsh,wfn%nspin))
   if (.not.allocated(wfn%dqatdr)) allocate(wfn%dqatdr(3,nat,nat,wfn%nspin))
   if (.not.allocated(wfn%dqshdL)) allocate(wfn%dqshdL(3,3,nsh,wfn%nspin))
   if (.not.allocated(wfn%dqatdL)) allocate(wfn%dqatdL(3,3,nat,wfn%nspin))
   wfn%dqshdL = 0.0_wp
   wfn%dqatdL = 0.0_wp
end subroutine ensure_charge_derivative_storage

end module tblite_xtb_response
