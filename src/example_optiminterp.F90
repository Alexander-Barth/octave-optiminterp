!  
! Example for using the n-dimentional optimal interpolation Fortran module
! Copyright (C) 2005 Alexander Barth <abarth@marine.usf.edu>
!  
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!  
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!  
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
! MA 02110-1301, USA.
!  

! Author: Alexander Barth <abarth@marine.usf.edu>

program test_optiminterp
 use optimal_interpolation
 implicit none

 ! n times n is the sized of the 2D background grid
 ! on: number of observations
 ! m: number of influential observations

 integer, parameter :: n=100, on=200, m=30

 ! xi(1,:) and xi(2,:) represent the x and y coordindate of the 
 ! grid of the interpolated field
 ! fi and vari are the interpolated field and its error variance resp.

 real(wp) :: xi(2,n*n), fi(1,n*n), vari(n*n)

 ! xi(1,:) and xi(2,:) represent the x and y coordindate of the 
 ! observations
 ! f and var are observations and their error variance resp.


 real(wp) :: x(2,on), var(on), f(1,on)

 ! param: inverse of the correlation length

 real(wp) :: param(2)

 integer :: i,j,l

 ! we use a simple random generator instead of Fortran 90's
 ! random_number in the hope that the results will be
 ! platform independent

 integer(8), parameter :: A = 1664525_8, B = 1013904223_8,
 &              Mo = 4294967296_8
 integer(8) :: next = 0_8


 ! create a regular 2D grid
 l = 1
 do j=1,n
   do i=1,n
     xi(1,l) = (i-1.)/(n-1.)
     xi(2,l) = (j-1.)/(n-1.)
     l = l+1
   end do
 end do

 ! param is the inverse of the correlation length
 param = 1./0.1

 ! the error variance of the observations
 var = 0.01

 ! location of observations

 do j=1,on
   do i=1,2
     ! simple random generator
     next = mod(A*next  + B,Mo)
     x(i,j) = real(next,8)/real(Mo,8)
   end do
 end do

 ! the underlying function to interpolate

 f(1,:) = sin( x(1,:)*6 ) * cos( x(2,:)*6);

 ! fi is the interpolated function and vari its error variance
 call optiminterp(x,f,var,param,m,xi,fi,vari)

 ! control value

 write(6,'(A,F10.6)') 'Expected:',2.2205936104348591E-002
 write(6,'(A,F10.6)') 'Obtained:',fi(1,1)

end program test_optiminterp

