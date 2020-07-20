#ifndef FITS_GENERAL_NEW_H
#define FITS_GENERAL_NEW_H

//#include <cfitsio/fitsio.h>
#include <fitsio.h>

// General routines

//      FITS Image manipulation library
//      Copyright (C) 2000 Jeremy Sanders
//      Contact: jss@ast.cam.ac.uk
//               Institute of Astronomy, Madingley Road,
//               Cambridge, CB3 0HA, UK.

//      See the file COPYING for full licence details.

//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.

//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.

//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

#define USE_FITSIMAGE_CHECKS

//////////////////////////////////////////////////
// these need to be changed together to change the
// fitsio working float type

typedef double CFloatType;
const int CFITSFloatType = TDOUBLE;
//const CFloatType CNullValue = DOUBLENULLVALUE;
const CFloatType CNullValue = 0.0;   // does 0 work?
const CFloatType CSmallNumber = 1e-15;
const int CFITSDataFormat = FLOAT_IMG;

////////////
// Functions

extern int IsNull(CFloatType val);

#endif
