! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Small second-order hyper-dual scalar used for exact derivatives of the
! low-dimensional ALPB/GBSA radial kernels.  It is deliberately local to the
! solvation implementation; no geometry finite differences are involved.
module tblite_solvation_hyperdual
   use mctc_env, only : wp
   implicit none
   private

   type, public :: hyperdual
      real(wp) :: v = 0.0_wp
      real(wp) :: dx = 0.0_wp
      real(wp) :: dy = 0.0_wp
      real(wp) :: dxy = 0.0_wp
   end type hyperdual

   public :: hd, hd_exp, hd_log, hd_sqrt, hd_tanh
   public :: operator(+), operator(-), operator(*), operator(/)

   interface operator(+)
      module procedure add_hh, add_hr, add_rh
   end interface
   interface operator(-)
      module procedure sub_hh, sub_hr, sub_rh, neg_h
   end interface
   interface operator(*)
      module procedure mul_hh, mul_hr, mul_rh
   end interface
   interface operator(/)
      module procedure div_hh, div_hr, div_rh
   end interface

contains

pure elemental function hd(v, dx, dy, dxy) result(a)
   real(wp), intent(in) :: v
   real(wp), intent(in), optional :: dx, dy, dxy
   type(hyperdual) :: a
   a%v = v
   if (present(dx)) a%dx = dx
   if (present(dy)) a%dy = dy
   if (present(dxy)) a%dxy = dxy
end function hd

pure elemental function add_hh(a,b) result(c)
   type(hyperdual), intent(in) :: a,b
   type(hyperdual) :: c
   c%v=a%v+b%v; c%dx=a%dx+b%dx; c%dy=a%dy+b%dy; c%dxy=a%dxy+b%dxy
end function
pure elemental function add_hr(a,b) result(c)
   type(hyperdual), intent(in) :: a; real(wp), intent(in) :: b; type(hyperdual) :: c
   c=a; c%v=c%v+b
end function
pure elemental function add_rh(a,b) result(c)
   real(wp), intent(in) :: a; type(hyperdual), intent(in) :: b; type(hyperdual) :: c
   c=b; c%v=c%v+a
end function

pure elemental function sub_hh(a,b) result(c)
   type(hyperdual), intent(in) :: a,b; type(hyperdual) :: c
   c%v=a%v-b%v; c%dx=a%dx-b%dx; c%dy=a%dy-b%dy; c%dxy=a%dxy-b%dxy
end function
pure elemental function sub_hr(a,b) result(c)
   type(hyperdual), intent(in) :: a; real(wp), intent(in) :: b; type(hyperdual) :: c
   c=a; c%v=c%v-b
end function
pure elemental function sub_rh(a,b) result(c)
   real(wp), intent(in) :: a; type(hyperdual), intent(in) :: b; type(hyperdual) :: c
   c%v=a-b%v; c%dx=-b%dx; c%dy=-b%dy; c%dxy=-b%dxy
end function
pure elemental function neg_h(a) result(c)
   type(hyperdual), intent(in) :: a; type(hyperdual) :: c
   c%v=-a%v; c%dx=-a%dx; c%dy=-a%dy; c%dxy=-a%dxy
end function

pure elemental function mul_hh(a,b) result(c)
   type(hyperdual), intent(in) :: a,b; type(hyperdual) :: c
   c%v = a%v*b%v
   c%dx = a%dx*b%v + a%v*b%dx
   c%dy = a%dy*b%v + a%v*b%dy
   c%dxy = a%dxy*b%v + a%dx*b%dy + a%dy*b%dx + a%v*b%dxy
end function
pure elemental function mul_hr(a,b) result(c)
   type(hyperdual), intent(in) :: a; real(wp), intent(in) :: b; type(hyperdual) :: c
   c%v=a%v*b; c%dx=a%dx*b; c%dy=a%dy*b; c%dxy=a%dxy*b
end function
pure elemental function mul_rh(a,b) result(c)
   real(wp), intent(in) :: a; type(hyperdual), intent(in) :: b; type(hyperdual) :: c
   c=b*a
end function

pure elemental function inv_h(a) result(c)
   type(hyperdual), intent(in) :: a; type(hyperdual) :: c
   real(wp) :: f1, f2
   c%v = 1.0_wp/a%v
   f1 = -c%v*c%v
   f2 = 2.0_wp*c%v*c%v*c%v
   c%dx = f1*a%dx
   c%dy = f1*a%dy
   c%dxy = f2*a%dx*a%dy + f1*a%dxy
end function
pure elemental function div_hh(a,b) result(c)
   type(hyperdual), intent(in) :: a,b; type(hyperdual) :: c
   c=a*inv_h(b)
end function
pure elemental function div_hr(a,b) result(c)
   type(hyperdual), intent(in) :: a; real(wp), intent(in) :: b; type(hyperdual) :: c
   c=a*(1.0_wp/b)
end function
pure elemental function div_rh(a,b) result(c)
   real(wp), intent(in) :: a; type(hyperdual), intent(in) :: b; type(hyperdual) :: c
   c=a*inv_h(b)
end function

pure elemental function hd_exp(a) result(c)
   type(hyperdual), intent(in) :: a; type(hyperdual) :: c
   real(wp) :: f
   f=exp(a%v); c%v=f; c%dx=f*a%dx; c%dy=f*a%dy; c%dxy=f*(a%dx*a%dy+a%dxy)
end function

pure elemental function hd_log(a) result(c)
   type(hyperdual), intent(in) :: a; type(hyperdual) :: c
   real(wp) :: f1,f2
   c%v=log(a%v); f1=1.0_wp/a%v; f2=-f1*f1
   c%dx=f1*a%dx; c%dy=f1*a%dy; c%dxy=f2*a%dx*a%dy+f1*a%dxy
end function

pure elemental function hd_sqrt(a) result(c)
   type(hyperdual), intent(in) :: a; type(hyperdual) :: c
   real(wp) :: f1,f2
   c%v=sqrt(a%v); f1=0.5_wp/c%v; f2=-0.25_wp/(a%v*c%v)
   c%dx=f1*a%dx; c%dy=f1*a%dy; c%dxy=f2*a%dx*a%dy+f1*a%dxy
end function

pure elemental function hd_tanh(a) result(c)
   type(hyperdual), intent(in) :: a; type(hyperdual) :: c
   real(wp) :: t,f1,f2
   t=tanh(a%v); f1=1.0_wp-t*t; f2=-2.0_wp*t*f1
   c%v=t; c%dx=f1*a%dx; c%dy=f1*a%dy; c%dxy=f2*a%dx*a%dy+f1*a%dxy
end function

end module tblite_solvation_hyperdual
