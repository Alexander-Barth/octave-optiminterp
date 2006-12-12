/*
  n-dimentional optimal interpolation
  Copyright (C) 2005 Alexander Barth <abarth@marine.usf.edu>

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


  Author: Alexander Barth <abarth@marine.usf.edu>
*/									


#include <octave/oct.h>
#include "f77-fcn.h"

extern "C"
{
  int F77_FUNC (optiminterp_wrapper, OPTIMINTERP_WRAPPER) 
    (
     const int& n, const int& nf, const int& gn, const int& on, const int& nparam,  
     double* ox, double* of, double* ovar, double* param,
     const int& m, double* gx,  double* gf, double* gvar);
}

DEFUN_DLD (optiminterp, args, ,
	   "- Loadable Function: [gf,gvar] = optiminterp(ox,of,ovar,param,m,gx)\n\
     \n\
     Performs an N-dimensional optimal interpolation.")
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
  int m = (int)args(4).scalar_value();
  Matrix gx = args(5).matrix_value();

  /*
   In Octave and Matlab convention (see gridatan), the positions of series 
   of points in n-D space is represented by an m-by-n matrix. However, it is more
   efficient to represent it as a n-ny-m matrix. The Fortran code uses
   the later convension. Therefore we need to transpose the following 
   matrices:
  */

  of = of.transpose();
  gx = gx.transpose();
  ox = ox.transpose();
  //octave_stdout << "shape(of) " << of.rows() << "x" << of.cols() << std::endl;

  int n = gx.rows();
  int nf = of.rows();
  int gn = gx.cols();
  int on = ox.cols();
  int nparam = param.numel();
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
    error ("optiminterp: number of colums in ox must be equal to the number of elements in of and ovar");
    return octave_value();
  }

  if (m > on) {
    warning("optiminterp: Number of influential points is larger than number of data points. Limiting the number of influential points to the number of data points");
    m = on;
  }

  F77_XFCN (optiminterp_wrapper, OPTIMINTERP_WRAPPER,
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
