#ifndef FITSIMAGE_NEW_H
#define FITSIMAGE_NEW_H

// Definition of FITS Image

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

#include "FITSGeneral.h"

class CFITSImage {
public:                    // public methods
  CFITSImage();
  CFITSImage(int xw, int yw);
  CFITSImage(const CFITSImage &other);
  ~CFITSImage();

  void Resize(int xw, int yw, int clear = 1);
    // resize the image (removes all values).
    // if clear, init with null.
  void Crop(int x1, int y1, int x2, int y2);
    // crops the image, preserving information

  // pixel routines
  /////////////////
  CFloatType GetPixel(int x, int y) const;
    // get (x, y)
  CFloatType GetPixelInterpolated(double x,
				  double y) const;
    // interpolate for non-integer coords (bounded)
  void SetPixel(int x, int y, CFloatType value);
    // set (x, y) to value
  void SetPixelBounded(int x, int y, CFloatType value);
    // set (x, y) to value, ignoring out-of image errors
  CFloatType GetPixelBounded(int x, int y) const;
    // get (x, y). If outside image gets null
  CFloatType GetPixelTotal() const;
    // get the total of all the pixels
  void SetAll(CFloatType value);
    // set all the pixels to value

  void NullRegion(int x1, int y1, int x2, int y2);
    // set region to null value
  void NullIfLess(CFloatType lower);
    // null pixel if less than lower
  void NullIfGreater(CFloatType upper);
    // null pixel if greater than upper
  void Bin(int factor);
    // bin the image by the factor. Adds pixels together.

  // operators
  ////////////
  void operator +=(const CFITSImage &other);
    // adds images together, (null+x) = x
  void operator -=(const CFITSImage &other);
    // subtracts images
  void operator /=(const CFITSImage &other);
    // divides image, x/0 = x/null = null
  void operator *=(const CFITSImage &other);
    // multiplies images together
  void operator =(const CFITSImage &other);
    // copy image
  void Square();
    // square pixels, null^2 = null
  void Sqrt();
    // sqrt pixels, ignore nulls, sqrt(-ve) -> error
  void Gaussian(CFloatType xsigma, CFloatType ysigma);
    // apply a gaussian to the image

  // breakout box
  ///////////////
  CFloatType *GetImageBuffer();
  const CFloatType *GetConstImageBuffer() const;
  int GetXW() const;
  int GetYW() const;
  int GetNoPixels() const;

private:                   // private data
  CFloatType *m_buffer;
  int m_xw, m_yw;

};

inline CFloatType *CFITSImage::GetImageBuffer()
{
  return m_buffer;
}

inline const CFloatType *CFITSImage::GetConstImageBuffer() const
{
  return m_buffer;
}

inline int CFITSImage::GetXW() const
{
  return m_xw;
}

inline int CFITSImage::GetYW() const
{
  return m_yw;
}

#endif

