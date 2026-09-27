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

!> @file tblite/solvation/born.f90
!> Provides a Born radii integrator

!> Integrator for Born radii based on the Onufriev-Bashford-Case model
module tblite_solvation_born
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_solvation_hyperdual
   implicit none
   private

   public :: new_born_integrator

   !> Implementation of GBOBC integrator
   type, public :: born_integrator
      !> van der Waals radii of the particles
      real(wp), allocatable :: vdwr(:)
      !> pair descreening approximation radii
      real(wp), allocatable :: rho(:)
      !> offset van der Waals radii
      real(wp), allocatable :: svdw(:)
      !> cut-off radius for the Born radius NN list
      real(wp) :: lrcut
      !> Scaling factor for Born radii
      real(wp) :: born_scale
      !> Volume polynome correction, default parameters correspond to GBOBCII
      real(wp) :: obc(3)
   contains
      !> Calculate Born radii for a given geometry
      procedure :: get_rad
   end type born_integrator

   real(wp), parameter :: lrcut_default = 35.0_wp * aatoau
   real(wp), parameter :: born_scale_default = 1.0_wp
   real(wp), parameter :: born_offset_default = 0.0_wp
   real(wp), parameter :: descreening_default = 0.8_wp
   real(wp), parameter :: obc_default(3) = [1.0_wp, 0.8_wp, 4.85_wp]

contains

!> Create new Born radii integrator
subroutine new_born_integrator(self, mol, vdwrad, descreening, born_scale, born_offset, &
      & obc, rcutoff)
   !> Instance of the Born integrator
   type(born_integrator), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Van-der-Waals Radii
   real(wp), intent(in) :: vdwRad(:)
   !> Dielectric descreening parameter
   real(wp), intent(in), optional :: descreening(:)
   !> Scaling factor for Born radii
   real(wp), intent(in), optional :: born_scale
   !> Offset parameter for Born radii integration
   real(wp), intent(in), optional :: born_offset
   !> GBOBC integrator parameters
   real(wp), intent(in), optional :: obc(3)
   !> Real-space cutoff for Born radii integration
   real(wp), intent(in), optional :: rCutoff

   self%lrcut = lrcut_default
   if (present(rCutoff)) then
      self%lrcut = rCutoff
   end if

   self%born_scale = born_scale_default
   if (present(born_scale)) then
      self%born_scale = born_scale
   end if

   self%obc = obc_default
   if (present(obc)) then
      self%obc = obc
   end if

   self%vdwr = vdwRad(mol%id)

   if (present(descreening)) then
      self%rho = self%vdwr * descreening(mol%id)
   else
      self%rho = self%vdwr * descreening_default
   end if

   if (present(born_offset)) then
      self%svdw = self%vdwr - born_offset
   else
      self%svdw = self%vdwr - born_offset_default
   end if
end subroutine new_born_integrator

!> Calculate Born radii
subroutine get_rad(self, mol, rad, draddr, d2raddr2)
   !> Instance of the Born integrator
   class(born_integrator), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Born radii
   real(wp), intent(out) :: rad(:)
   !> Derivative of Born radii w.r.t. cartesian displacements
   real(wp), intent(out), optional :: draddr(:, :, :)
   !> Second derivatives of Born radii
   real(wp), intent(out), optional :: d2raddr2(:, :, :, :, :)

   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: brdr(:, :, :), brhess(:, :, :, :, :)

   call new_adjacency_list(list, mol, trans, self%lrcut)
   allocate(brdr(3, mol%nat, mol%nat))
   if (present(d2raddr2)) allocate(brhess(3,mol%nat,3,mol%nat,mol%nat))

   if (present(d2raddr2)) then
      call compute_bornr(mol%nat, mol%xyz, list, &
         & self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, rad, brdr, brhess)
   else
      call compute_bornr(mol%nat, mol%xyz, list, &
         & self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, rad, brdr)
   end if

   if (present(draddr)) draddr = brdr
   if (present(d2raddr2)) d2raddr2 = brhess
end subroutine get_rad

subroutine compute_bornr(nat, xyz, list, vdwr, rho, svdw, c1, obc, &
      & brad, brdr, brhess)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Neighbourlist
   type(adjacency_list), intent(in) :: list
   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwr(:)
   !> Descreened van-der-Waals radii
   real(wp), intent(in) :: rho(:)
   !> van-der-Waals radii with offset
   real(wp), intent(in) :: svdw(:)
   !> Scaling factor for the Born radii
   real(wp), intent(in) :: c1
   !> Volume polynome correction
   real(wp), intent(in) :: obc(3)
   !> Born radii
   real(wp), intent(out) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(out) :: brdr(:, :, :)
   !> Second derivative of Born radii
   real(wp), intent(out), optional :: brhess(:, :, :, :, :)

   integer :: iat, jat, ic, jc
   real(wp) :: br, dpsi, svdwi, vdwri, s1, v1, s2, arg, arg2
   real(wp) :: th, ch

   if (present(brhess)) then
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr, brhess)
   else
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr)
   end if

   do iat = 1, nat
      call born_transform_derivs(brad(iat), svdw(iat), vdwr(iat), c1, obc, br, dpsi, arg)

      if (present(brhess)) then
         ! arg carries F'' here; brdr still contains dPsi/dR.
         do jat = 1, nat
            do jc = 1, 3
               do ic = 1, 3
                  brhess(ic,:,jc,jat,iat) = dpsi*brhess(ic,:,jc,jat,iat) &
                     + arg*brdr(ic,:,iat)*brdr(jc,jat,iat)
               end do
            end do
         end do
      end if
      brdr(:,:,iat) = dpsi*brdr(:,:,iat)
      brad(iat) = br
   end do

end subroutine compute_bornr


pure subroutine compute_psi(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2)
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   type(adjacency_list), intent(in) :: list
   real(wp), intent(in) :: vdwr(:), rho(:)
   real(wp), intent(out) :: psi(:)
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), intent(out), optional :: d2psidr2(:, :, :, :, :)

   integer :: iat, jat, img, inl
   real(wp) :: vec(3), r, rhoi, rhoj, rvdwi, rvdwj
   real(wp) :: g, gp, gpp
   logical :: ijov, jiov

   psi = 0.0_wp
   dpsidr = 0.0_wp
   if (present(d2psidr2)) d2psidr2 = 0.0_wp

   do iat = 1, nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)
         vec = xyz(:,iat)-xyz(:,jat)
         r = norm2(vec)
         rhoi = rho(iat); rhoj = rho(jat)
         rvdwi = vdwr(iat); rvdwj = vdwr(jat)
         ijov = r < (rvdwi+rhoj)
         jiov = r < (rhoi+rvdwj)

         if (.not.(ijov .or. jiov)) then
            call psi_nonoverlap(r, rhoj, g, gp, gpp)
            call accumulate_psi(iat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            if (abs(rhoi-rhoj) < 1.0e-8_wp) then
               call accumulate_psi(jat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            else
               call psi_nonoverlap(r, rhoi, g, gp, gpp)
               call accumulate_psi(jat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            end if
         else if (.not.ijov .and. jiov) then
            call psi_nonoverlap(r, rhoj, g, gp, gpp)
            call accumulate_psi(iat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            if ((r+rhoi) > rvdwj) then
               call psi_overlap(r, rhoi, rvdwj, g, gp, gpp)
               call accumulate_psi(jat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            end if
         else if (ijov .and. .not.jiov) then
            if ((r+rhoj) > rvdwi) then
               call psi_overlap(r, rhoj, rvdwi, g, gp, gpp)
               call accumulate_psi(iat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            end if
            call psi_nonoverlap(r, rhoi, g, gp, gpp)
            call accumulate_psi(jat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
         else
            if ((r+rhoj) > rvdwi) then
               call psi_overlap(r, rhoj, rvdwi, g, gp, gpp)
               call accumulate_psi(iat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            end if
            if ((r+rhoi) > rvdwj) then
               call psi_overlap(r, rhoi, rvdwj, g, gp, gpp)
               call accumulate_psi(jat, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
            end if
         end if
      end do
   end do
end subroutine compute_psi

pure subroutine accumulate_psi(target, iat, jat, vec, r, g, gp, gpp, psi, dpsidr, d2psidr2)
   integer, intent(in) :: target, iat, jat
   real(wp), intent(in) :: vec(3), r, g, gp, gpp
   real(wp), intent(inout) :: psi(:), dpsidr(:, :, :)
   real(wp), intent(inout), optional :: d2psidr2(:, :, :, :, :)
   real(wp) :: a, b, h(3,3)
   integer :: ic,jc
   psi(target)=psi(target)+g
   a=gp/r
   dpsidr(:,iat,target)=dpsidr(:,iat,target)+a*vec
   dpsidr(:,jat,target)=dpsidr(:,jat,target)-a*vec
   if (present(d2psidr2)) then
      b=(gpp-a)/(r*r)
      h=0.0_wp
      do ic=1,3
         h(ic,ic)=a
         do jc=1,3
            h(ic,jc)=h(ic,jc)+b*vec(ic)*vec(jc)
         end do
      end do
      d2psidr2(:,iat,:,iat,target)=d2psidr2(:,iat,:,iat,target)+h
      d2psidr2(:,jat,:,jat,target)=d2psidr2(:,jat,:,jat,target)+h
      d2psidr2(:,iat,:,jat,target)=d2psidr2(:,iat,:,jat,target)-h
      d2psidr2(:,jat,:,iat,target)=d2psidr2(:,jat,:,iat,target)-h
   end if
end subroutine accumulate_psi

pure subroutine psi_nonoverlap(r, rho, g, gp, gpp)
   real(wp), intent(in) :: r, rho
   real(wp), intent(out) :: g, gp, gpp
   type(hyperdual) :: x, ap, am, val
   x=hd(r,1.0_wp,1.0_wp)
   ap=x+rho; am=x-rho
   val=rho/(ap*am) + 0.5_wp*hd_log(am/ap)/x
   g=val%v; gp=val%dx; gpp=val%dxy
end subroutine psi_nonoverlap

pure subroutine psi_overlap(r, rho, rvdw, g, gp, gpp)
   real(wp), intent(in) :: r, rho, rvdw
   real(wp), intent(out) :: g, gp, gpp
   type(hyperdual) :: x, ap, am, val
   x=hd(r,1.0_wp,1.0_wp)
   ap=x+rho; am=x-rho
   val=1.0_wp/rvdw - 1.0_wp/ap + 0.5_wp/x * &
      (0.5_wp*am*(1.0_wp/ap-ap/(rvdw*rvdw)) - hd_log(ap/rvdw))
   g=val%v; gp=val%dx; gpp=val%dxy
end subroutine psi_overlap

pure subroutine born_transform_derivs(psi, svdw, vdwr, scale, obc, value, d1, d2)
   real(wp), intent(in) :: psi, svdw, vdwr, scale, obc(3)
   real(wp), intent(out) :: value, d1, d2
   type(hyperdual) :: x, b, arg2, arg, out
   real(wp) :: s1,v1,s2
   s1=1.0_wp/svdw; v1=1.0_wp/vdwr; s2=0.5_wp*svdw
   x=hd(psi,1.0_wp,1.0_wp)
   b=x*s2
   arg2=b*(obc(3)*b-obc(2))
   arg=b*(obc(1)+arg2)
   out=scale/(s1-v1*hd_tanh(arg))
   value=out%v; d1=out%dx; d2=out%dxy
end subroutine born_transform_derivs

end module tblite_solvation_born
