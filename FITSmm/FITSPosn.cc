// Routines for holding WCS Data (not very good stuff)

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

#include <string.h>
#include "FITSPosn.h"
#include "FITSFile.h"

CFITSPosn::CFITSPosn()
{
  Zero();
}

CFITSPosn::CFITSPosn(const CFITSPosn &copy)
{
  strcpy(m_Telescop, copy.m_Telescop);
  strcpy(m_Instrume, copy.m_Instrume);
  strcpy(m_ObsMode, copy.m_ObsMode);
  strcpy(m_RADECSys, copy.m_RADECSys);

  m_Equinox = copy.m_Equinox;
  m_bZero = copy.m_bZero;
  m_bScale = copy.m_bScale;
  strcpy(m_bUnit, copy.m_bUnit);

  m_crPix1 = copy.m_crPix1;
  strcpy(m_cType1, copy.m_cType1);
  m_crVal1 = copy.m_crVal1;
  m_cDelt1 = copy.m_cDelt1;
  strcpy(m_cUnit1, copy.m_cUnit1);

  m_crPix2 = copy.m_crPix2;
  strcpy(m_cType2, copy.m_cType2);
  m_crVal2 = copy.m_crVal2;
  m_cDelt2 = copy.m_cDelt2;
  strcpy(m_cUnit2, copy.m_cUnit2);
}

const CFITSPosn & CFITSPosn::operator=(const CFITSPosn &copy)
{
  strcpy(m_Telescop, copy.m_Telescop);
  strcpy(m_Instrume, copy.m_Instrume);
  strcpy(m_ObsMode, copy.m_ObsMode);
  strcpy(m_RADECSys, copy.m_RADECSys);

  m_Equinox = copy.m_Equinox;
  m_bZero = copy.m_bZero;
  m_bScale = copy.m_bScale;
  strcpy(m_bUnit, copy.m_bUnit);

  m_crPix1 = copy.m_crPix1;
  strcpy(m_cType1, copy.m_cType1);
  m_crVal1 = copy.m_crVal1;
  m_cDelt1 = copy.m_cDelt1;
  strcpy(m_cUnit1, copy.m_cUnit1);

  m_crPix2 = copy.m_crPix2;
  strcpy(m_cType2, copy.m_cType2);
  m_crVal2 = copy.m_crVal2;
  m_cDelt2 = copy.m_cDelt2;
  strcpy(m_cUnit2, copy.m_cUnit2);

  return *this;
}

void CFITSPosn::Zero()
{
  m_Telescop[0] = 0;
  m_Instrume[0] = 0;
  m_ObsMode[0] = 0;
  strcpy(m_RADECSys, "FK5");

  m_Equinox = 2000;

  m_bZero = 0.0;
  m_bScale = 1.0;
  m_bUnit[0] = 0;

  m_crPix1 = 0;
  strcpy(m_cType1, "RA---TAN");
  m_crVal1 = 0.0;
  m_cDelt1 = 1.0e10;
  strcpy(m_cUnit1, "deg");

  m_crPix2 = 0;
  strcpy(m_cType2, "DEC--TAN");
  m_crVal2 = 0.0;
  m_cDelt2 = 1.0;
  strcpy(m_cUnit2, "deg");
}

void CFITSPosn::RescaleCoords(int factor)
{
  double dFactor = factor;

  m_crPix1 /= dFactor;
  m_cDelt1 *= dFactor;
  m_crPix2 /= dFactor;
  m_cDelt2 *= dFactor;
}

void CFITSPosn::WriteFITSHeader(CFITSFile &file)
{
  if( m_cDelt1 == 1.0e10 ) return;

  file.UpdateKey("TELESCOP", CFITSFile::tstring, m_Telescop, "Telescope name");
  file.UpdateKey("INSTRUME", CFITSFile::tstring, m_Instrume, "Instrument name");
  file.UpdateKey("OBS_MODE", CFITSFile::tstring, m_ObsMode, "Observing mode");
  file.UpdateKey("RADECSYS", CFITSFile::tstring, m_RADECSys, "Equatorial system reference");

  file.UpdateKey("EQUINOX", CFITSFile::tint, &m_Equinox, "Equinox");
  file.UpdateKey("BZERO", CFITSFile::tint, &m_bZero, 0);
  file.UpdateKey("BSCALE", CFITSFile::tfloat, &m_bScale, 0);
  file.UpdateKey("BUNIT", CFITSFile::tstring, m_bUnit, "Units of data");

  file.UpdateKey("CRPIX1", CFITSFile::tfloat, &m_crPix1, "Reference pixel");
  file.UpdateKey("CTYPE1", CFITSFile::tstring, m_cType1, "Projection");
  file.UpdateKey("CRVAL1", CFITSFile::tfloat, &m_crVal1, "Right Ascension");
  file.UpdateKey("CDELT1", CFITSFile::tfloat, &m_cDelt1, "Pixel size");
  file.UpdateKey("CUNIT1", CFITSFile::tstring, m_cUnit1, "Units of coordinate");

  file.UpdateKey("CRPIX2", CFITSFile::tfloat, &m_crPix2, "Reference pixel");
  file.UpdateKey("CTYPE2", CFITSFile::tstring, m_cType2, "Projection");
  file.UpdateKey("CRVAL2", CFITSFile::tfloat, &m_crVal2, "Declination");
  file.UpdateKey("CDELT2", CFITSFile::tfloat, &m_cDelt2, "Pixel size");
  file.UpdateKey("CUNIT2", CFITSFile::tstring, m_cUnit2, "Units of coordinate");
}

void CFITSPosn::ReadFITSHeader(CFITSFile &file)
{
  Zero();

  if( file.ReadKey("TELESCOP", CFITSFile::tstring, m_Telescop, "") == 0 )
    return;

  const char *defObsMode = "POINTING";
  const char *defRADEC = "";
  file.ReadKey("INSTRUME", CFITSFile::tstring, m_Instrume, "UNKNOWN");
  file.ReadKey("OBS_MODE", CFITSFile::tstring, m_ObsMode, defObsMode);
  file.ReadKey("RADECSYS", CFITSFile::tstring, m_RADECSys, defRADEC);

  const double defZero = 0.0;
  const double defScale = 1.0;
  const char *defUnit = "";
  file.ReadKey("EQUINOX", CFITSFile::tint, &m_Equinox, "J2000");
  file.ReadKey("BZERO", CFITSFile::tint, &m_bZero, &defZero);
  file.ReadKey("BSCALE", CFITSFile::tfloat, &m_bScale, &defScale);
  file.ReadKey("BUNIT", CFITSFile::tstring, m_bUnit, defUnit);

  file.ReadKey("CRPIX1", CFITSFile::tfloat, &m_crPix1);
  file.ReadKey("CTYPE1", CFITSFile::tstring, m_cType1);
  file.ReadKey("CRVAL1", CFITSFile::tfloat, &m_crVal1);
  file.ReadKey("CDELT1", CFITSFile::tfloat, &m_cDelt1, &defZero);
  file.ReadKey("CUNIT1", CFITSFile::tstring, m_cUnit1, defUnit);

  file.ReadKey("CRPIX2", CFITSFile::tfloat, &m_crPix2);
  file.ReadKey("CTYPE2", CFITSFile::tstring, m_cType2);
  file.ReadKey("CRVAL2", CFITSFile::tfloat, &m_crVal2);
  file.ReadKey("CDELT2", CFITSFile::tfloat, &m_cDelt2, &defZero);
  file.ReadKey("CUNIT2", CFITSFile::tstring, m_cUnit2, defUnit);
}
