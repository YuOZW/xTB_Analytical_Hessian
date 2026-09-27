! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Analytic second derivatives for the exponential coordination number used by GFN1-xTB.
module tblite_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_partition, only : work_partition, owns_pair
   implicit none
   private

   public :: add_coordination_number_hessian

   real(wp), parameter :: default_kcn = 16.0_wp

contains

!> Add the Cartesian Hessian of the exponential coordination number contracted
!> with a fixed derivative dE/dCN. This is the second derivative contribution
!> sum_A dE/dCN_A * d2CN_A/dR dR.
subroutine add_coordination_number_hessian(mol, trans, cutoff, rcov, dEdcn, hessian, partition)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: rcov(:)
   real(wp), intent(in) :: dEdcn(:)
   real(wp), intent(inout) :: hessian(:, :)
   type(work_partition), intent(in), optional :: partition

   integer :: iat, jat, izp, jzp, itr, ic, jc, ii, jj
   real(wp) :: vec(3), r2, r1, cutoff2, rc, expterm, cf, dcf, d2cf
   real(wp) :: dEdcnij, block(3, 3), tmp

   cutoff2 = cutoff*cutoff

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         rc = rcov(izp) + rcov(jzp)
         dEdcnij = dEdcn(iat) + dEdcn(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = dot_product(vec, vec)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            expterm = exp(-default_kcn*(rc/r1 - 1.0_wp))
            cf = 1.0_wp/(1.0_wp + expterm)
            tmp = cf*(1.0_wp - cf)
            dcf = -tmp*default_kcn*rc/r2
            d2cf = tmp*(1.0_wp - 2.0_wp*cf)*(default_kcn*rc)**2/(r2*r2) &
               & + 2.0_wp*tmp*default_kcn*rc/(r2*r1)

            do ic = 1, 3
               do jc = 1, 3
                  block(ic, jc) = dEdcnij*( &
                     & d2cf*vec(ic)*vec(jc)/r2 &
                     & - dcf*vec(ic)*vec(jc)/(r2*r1))
               end do
               block(ic, ic) = block(ic, ic) + dEdcnij*dcf/r1
            end do

            do ic = 1, 3
               ii = 3*(iat - 1) + ic
               do jc = 1, 3
                  jj = 3*(jat - 1) + jc
                  hessian(ii, 3*(iat - 1) + jc) = hessian(ii, 3*(iat - 1) + jc) + block(ic, jc)
                  hessian(3*(jat - 1) + ic, jj) = hessian(3*(jat - 1) + ic, jj) + block(ic, jc)
                  hessian(ii, jj) = hessian(ii, jj) - block(ic, jc)
                  hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) = &
                     & hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) - block(ic, jc)
               end do
            end do
         end do
      end do
   end do
end subroutine add_coordination_number_hessian

end module tblite_ncoord
