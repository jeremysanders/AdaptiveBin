#ifndef FITSPOSN_NEW_H
#define FITSPOSN_NEW_H

// Updated version
// Holds data on scales of images

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

class CFITSFile;

class CFITSPosn {

// public interface
public:
  CFITSPosn();                      // constructor: copy default values
  CFITSPosn(const CFITSPosn &copy); // copy constructor
  const CFITSPosn & operator=(const CFITSPosn &copy);

  void RescaleCoords(int factor);   // rescale coords by factor
  void WriteFITSHeader(CFITSFile &fileid);
  void ReadFITSHeader(CFITSFile &fileid);
  void Zero();

// public variables
public:
  char m_Telescop[9]; // telescope name
  char m_Instrume[9]; // instrument name
  char m_ObsMode[9];  // observing mode
  char m_RADECSys[9]; // Equatorial system reference

  int m_Equinox;      // equinox of observation

  double m_bZero;     // data scaling (zero posn)
  double m_bScale;    // data scale
  char m_bUnit[9];    // unit of data

  double m_crPix1;    // reference posn (axis 1)
  char m_cType1[9];   // projection (axis 1)
  double m_crVal1;    // value of reference (axis 1)
  double m_cDelt1;    // pixel size (axis 1)
  char m_cUnit1[9];   // unit (axis 1)

  double m_crPix2;    // reference posn (axis 2)
  char m_cType2[9];   // projection
  double m_crVal2;    // value of reference
  double m_cDelt2;    // pixel size
  char m_cUnit2[9];   // unit
};

#endif
