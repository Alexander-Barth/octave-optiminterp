!
!  n-dimentional optimal interpolation
!  Copyright (C) 2005 Alexander Barth <abarth@marine.usf.edu>
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

!  Author: Alexander Barth <abarth@marine.usf.edu>

subroutine optiminterpwrapper(n,nf,gn,on,nparam,ox,of,ovar,    &
                               param,m,gx,gf,gvar)
 use optimal_interpolation
 implicit none
 integer :: m,n,nf,gn,on,nparam
 real(8) :: gx(n,gn),ox(n,on),of(nf,on),ovar(on),param(nparam)
 real(8) :: gf(nf,gn),gvar(gn)


 call optiminterp(ox,of,ovar,param,m,gx,gf,gvar)
!!$
!!$ integer :: i,j
!!$ real :: x=1
!!$ do i=1,1000
!!$   write(6,*) 'i',i
!!$   do j= 1,1000000
!!$    x = x + sin(x)
!!$   end do
!!$ end do

end subroutine optiminterpwrapper
