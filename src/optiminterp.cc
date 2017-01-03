// Copyright (C) 2005 Alexander Barth <barth.alexander@gmail.com>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

// n-dimensional optimal interpolation

#include <octave/oct.h>
#include "f77-fcn.h"

#if defined (OCTAVE_HAVE_F77_INT_TYPE)
#  define TO_F77_INT(x) octave::to_f77_int (x)
#else
typedef octave_idx_type F77_INT;
#  define TO_F77_INT(x) (x)
#endif

extern "C"
{
  F77_RET_T 
  F77_FUNC (optiminterpwrapper, OPTIMINTERPWRAPPER) 
    (
     const F77_INT& n, const F77_INT& nf, const F77_INT& gn,
     const F77_INT& on, const F77_INT& nparam,  
     double* ox, double* of, double* ovar, double* param,
     const F77_INT& m, double* gx,  double* gf, double* gvar);
}

DEFUN_DLD (optiminterp, args, ,
"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{fi},@var{vari}] =} optiminterp(@var{x},@var{f},@var{var},@var{len},@var{m},@var{xi}) \n\
Performs an @var{n}-dimensional optimal interpolation. \n\
\n\
Every elements in @var{f} corresponds to a data point (observation)\n\
at location the @var{x(:,1)}, @var{x(:,2)},... @var{x(:,n)} in a @var{n}-dimensional space with the error variance @var{var}\n\
\n\
Each element of the @var{n} by 1 vector @var{len} is the correlation length in the corresponding direction.\n\
\n\
@var{m} represents the number of influential points.\n\
\n\
@var{xi(:,1)}, @var{xi(:,2)} ... @var{xi(:,n)} are the coordinates of the points where the field is\n\
interpolated. @var{fi} is the interpolated field and @var{vari} is \n\
its error variance.\n\
\n\
The background field of the optimal interpolation is zero.\n\
For a different background field, the background field\n\
must be subtracted from the observation, the difference \n\
is mapped by OI onto the background grid and finally the\n\
background is added back to the interpolated field.\n\
\n\
The error variance of the background field is assumed to \n\
have a error variance of one. The error variances of the observations need to be scaled accordingly. \n\
@end deftypefn\n\
@seealso{optiminterp1,optiminterp2,optiminterp3,optiminterp4}\n")
{
  octave_value_list retval;
  
  if (args.length() != 6) {
    print_usage ();
    return octave_value();
  }

  Matrix ox = args(0).matrix_value();
  Matrix of = args(1).matrix_value();
  Matrix ovar = args(2).matrix_value();
  Matrix param = args(3).matrix_value();
  F77_INT m = TO_F77_INT ((int)args(4).scalar_value());
  Matrix gx = args(5).matrix_value();

  // The Fortran routine expects the inverse of the correlation length
  param = 1/param;

  /*
   In Octave and Matlab convention (see gridatan), the positions of series 
   of points in n-D space is represented by an m-by-n matrix. However, it is more
   efficient to represent it as a n-ny-m matrix. The Fortran code uses
   the later convention. Therefore we need to transpose the following 
   matrices:
  */

  of = of.transpose();
  gx = gx.transpose();
  ox = ox.transpose();
  //octave_stdout << "shape(of) " << of.rows() << "x" << of.cols() << std::endl;

  F77_INT n = TO_F77_INT (gx.rows());
  F77_INT nf = TO_F77_INT (of.rows());
  F77_INT gn = TO_F77_INT (gx.cols());
  F77_INT on = TO_F77_INT (ox.cols());
  F77_INT nparam = TO_F77_INT (param.numel());
  Matrix gf = Matrix(nf,gn);
  Matrix gvar = Matrix(gn,1);

  /*
  octave_stdout << "nf " << nf << std::endl;
  octave_stdout << "n "  << n << std::endl;
  octave_stdout << "gn " << gn << std::endl;
  octave_stdout << "on " << on << std::endl;
  */

  OCTAVE_QUIT;
  if (gx.rows() != ox.rows()) {
    error ("optiminterp: number of rows in gx and ox must be the same");
    return octave_value();
  }

  if (ox.cols() != of.cols() || ox.cols() != ovar.numel()) {
    error ("optiminterp: number of columns in ox must be equal to the number of elements in of and ovar");
    return octave_value();
  }

  if (m > on) {
    warning("optiminterp: Number of influential points is larger than number of data points. Limiting the number of influential points to the number of data points");
    m = on;
  }


  F77_XFCN (optiminterpwrapper, OPTIMINTERPWRAPPER,
	    ( n,nf,gn,on,nparam,ox.fortran_vec(),of.fortran_vec(),ovar.fortran_vec(),
              param.fortran_vec(),m,gx.fortran_vec(),gf.fortran_vec(),gvar.fortran_vec()));

  if (f77_exception_encountered) {
      error ("unrecoverable error in optiminterp");
      return retval;
  }

  gf = gf.transpose();

  retval(0) = gf;
  retval(1) = gvar;

  return retval;
}
