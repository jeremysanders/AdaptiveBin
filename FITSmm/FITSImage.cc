// Routines for the manipulation of images

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FITSImage.h"

CFITSImage::CFITSImage()
{
  m_buffer = 0;
  m_xw = m_yw = -1;
}

CFITSImage::CFITSImage(int xw, int yw)
{
  m_buffer = 0;

  Resize(xw, yw);
}

CFITSImage::CFITSImage(const CFITSImage &other)
{
  m_buffer = 0;
  //  Resize(other.m_xw, other.m_yw);

  operator=(other);
}

CFITSImage::~CFITSImage()
{
  if(m_buffer != 0)
    delete[] m_buffer;
}

////////////////////////////////////////////////////////////////////////////

int CFITSImage::GetNoPixels() const
{
#ifdef USE_FITSIMAGE_CHECKS
  if(m_xw <= 0 || m_yw <= 0) {
    fprintf(stderr, "*  CFITSImage::GetNoPixels(): Invalid image size\n");
    exit(-1);
  }
#endif

  return m_xw*m_yw;
}

void CFITSImage::Resize(int xw, int yw, int clear)
{
  if(m_buffer != 0) {
    delete[] m_buffer;
    m_buffer = 0;
  }

  m_xw = xw; m_yw = yw;

  m_buffer = new CFloatType[GetNoPixels()];

  if( clear ) {
    // blank image
    CFloatType *iter = m_buffer;
    for(int i = GetNoPixels() - 1; i >=0; i--) {
      *iter = CNullValue;
      iter++;
    }
  }
}

void CFITSImage::Crop(int x1, int y1, int x2, int y2)
{
  if( x1 > x2 ) { int t=x1; x1=x2; x2=t; } // swap if wrong way round
  if( y1 > y2 ) { int t=y1; y1=y2; y2=t; }

  if( x1 >= m_xw || y1 >= m_yw || x1 < 0 || y1 < 0
      || x2 >= m_xw || y2 >= m_yw || x2 < 0 || y2 < 0 ) {
    fprintf(stderr, "* CFITSImage::Crop - area exceeds xw or yw\n");
    exit(-1);
  }

  int new_xw = x2-x1+1;
  int new_yw = y2-y1+1;
  CFloatType *newBuffer = new CFloatType[new_xw*new_yw];

  for(int y = y1; y <= y2; y++)
    for(int x = x1; x <= x2; x++)
      newBuffer[ (y-y1)*new_xw + (x-x1) ] = m_buffer[y*m_xw + x];

  delete[] m_buffer;
  m_xw = new_xw;
  m_yw = new_yw;
  m_buffer = newBuffer;
}

CFloatType CFITSImage::GetPixel(int x, int y) const
{
#ifdef USE_FITSIMAGE_CHECKS
  if( x >= m_xw || y >= m_yw || x < 0 || y < 0 ) {
    fprintf(stderr, "*  CFITSImage::GetPixel: Position exceeds xw or yw"
	    " (x=%i,y=%i)\n", x, y);
    abort();
  }
#endif

  return m_buffer[x + y*m_xw];
}

CFloatType CFITSImage::GetPixelInterpolated(double x, double y) const
{
  // make origin 0.,0.
  x += 0.5; y += 0.5;

  // need to find nearest integer+1/2
  // int+0.5 <= x < int+1.5
  // equiv  int <= x-0.5 < int+1
  int lx = int(floor(x-0.5));
  int ly = int(floor(y-0.5));
  int ux = lx + 1;
  int uy = ly + 1;

  // check bounds
  if( ux >= m_xw ) ux = m_xw-1;
  if( uy >= m_yw ) uy = m_yw-1;

  if( lx < 0 ) lx = 0;
  if( ly < 0 ) ly = 0;

  // find nearest pixels
  const double pix[4] = { GetPixel(lx, ly), GetPixel(ux, ly),
			  GetPixel(lx, uy), GetPixel(ux, uy) };

  // find x weighted average
  const double delleft = x-(lx+0.5);
  const double deltop = y-(ly+0.5);

  const double avbot = pix[2]*(1.-delleft) + pix[3]*delleft;
  const double avtop = pix[0]*(1.-delleft) + pix[1]*delleft;

  // find y weighted average
  const double av = avtop*(1.-deltop) + avbot*deltop;

  return av;
}

CFloatType CFITSImage::GetPixelBounded(int x, int y) const
{
  if( x >= m_xw || y >= m_yw || x < 0 || y < 0 )
    return CNullValue;
  else
    return m_buffer[x + y*m_xw];
}

CFloatType CFITSImage::GetPixelTotal() const
{
  CFloatType total = 0.0;
  for(int i = m_xw*m_yw-1; i >= 0; i--)
    total += m_buffer[i];
  return total;
}

void CFITSImage::SetAll(CFloatType value)
{
  for(int i = m_xw*m_yw-1; i >= 0; i--)
    m_buffer[i] = value;
}

void CFITSImage::SetPixel(int x, int y, CFloatType value)
{
#ifdef USE_FITSIMAGE_CHECKS
  if( x >= m_xw || y >= m_yw || x < 0 || y < 0 ) {
    fprintf(stderr, "*  CFITSImage::SetPixel: Position exceeds xw or yw"
	    " (x=%i,y=%i)\n", x, y);
    abort();
  }
#endif

  m_buffer[x + y*m_xw] = value;
}

void CFITSImage::SetPixelBounded(int x, int y, CFloatType value)
{
  if( x >= m_xw || y >= m_yw || x < 0 || y < 0 )
    return;
  else
    m_buffer[x + y*m_xw] = value;
}

void CFITSImage::NullRegion(int x1, int y1, int x2, int y2)
{
#ifdef USE_FITSIMAGE_CHECKS
  if( x1 >= m_xw || y1 >= m_yw || x2 >= m_xw || y2 >= m_yw
     || x1 < 0 || y1 < 0 || x2 < 0 || y2 < 0 ) {
    fprintf(stderr, "*  CFITSImage::NullRegion: Region exceeds image");
    exit(-1);
  }
#endif

  if(x2 < x1) { int t=x1; x1=x2; x2=t; }
  if(y2 < y1) { int t=y1; y1=y2; y2=t; }

  // slow, but works
  for(int x = x1; x <= x2; x++)
    for(int y = y1; y <= y2; y++)
      SetPixel(x, y, CNullValue);
}

void CFITSImage::Bin(int factor)
{
  if(factor == 1)
    return;

  if(factor != 2 && factor != 4) {
    fprintf(stderr, "*  CFITSImage::Bin: Invalid binning factor: %i\n",
	    factor);
    exit(-1);
  }

  if(m_xw % factor != 0 || m_yw % factor != 0) {
    fprintf(stderr, "*  CFITSImage::Bin: Image size not divisible by factor\n");
    exit(-1);
  }

  int new_xw = m_xw / factor; 
  int new_yw = m_yw / factor;

  CFITSImage *newImage = new CFITSImage(new_xw, new_yw);

  for(int x = new_xw - 1; x >= 0; x--) {
    for(int y = new_yw - 1; y >= 0; y--) {
      CFloatType sum = 0.0;
      int added = 0;
      
      for(int delx = factor - 1; delx >= 0; delx--)
	for(int dely = factor - 1; dely >= 0; dely--) {
	  CFloatType pixel = GetPixel( x*factor + delx, y*factor + dely);
	  if( ! IsNull(pixel) ) {
	    sum += pixel;
	    added = 1;
	  }
	}

      if( added )
	newImage -> SetPixel(x, y, sum);
      else
	newImage -> SetPixel(x, y, CNullValue);
    }
  }

  // copy new image over
  Resize(new_xw, new_yw, 0);
  CFloatType *to = m_buffer;
  CFloatType *from = newImage->m_buffer;

  for(int i = GetNoPixels()-1; i >= 0; i--) {
    *to = *from;
    to++; from++;
  }

  delete newImage;
}

//////////////////////////////////////////////////////////////

void CFITSImage::operator +=(const CFITSImage &other)
{
  if( m_xw != other.m_xw || m_yw != other.m_yw ) {
    fprintf(stderr, "*  CFITSImage::operator +=: "
	    "images are of different size\n");
    exit(-1);
  }

  CFloatType *to = m_buffer;
  CFloatType *from = other.m_buffer;

  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( IsNull(*to) )
      *to = *from;
    else {
      if( ! IsNull(*from) )
	*to += *from;
    }
    to++; from++;
  }

}

void CFITSImage::operator -=(const CFITSImage &other)
{
  if( m_xw != other.m_xw || m_yw != other.m_yw ) {
    fprintf(stderr, "*  CFITSImage::operator +=: "
	    "images are of different size\n");
    exit(-1);
  }

  CFloatType *to = m_buffer;
  CFloatType *from = other.m_buffer;

  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( IsNull(*to) )
      *to = *from;
    else {
      if( ! IsNull(*from) )
	*to -= *from;
    }
    to++; from++;
  }

}

void CFITSImage::operator /=(const CFITSImage &other)
{
  if( m_xw != other.m_xw || m_yw != other.m_yw ) {
    fprintf(stderr, "*  CFITSImage::operator /=: "
	    "images are of different size\n");
    exit(-1);
  }

  CFloatType *to = m_buffer;
  CFloatType *from = other.m_buffer;

  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( IsNull(*from) )
      *to = CNullValue;
    else {
      if( fabs(*from) < CSmallNumber )
	*to = 0.0;
      else
	*to /= *from;
    }

    to++; from++;
  }
}

void CFITSImage::operator *=(const CFITSImage &other)
{
    if( m_xw != other.m_xw || m_yw != other.m_yw ) {
    fprintf(stderr, "*  CFITSImage::operator /=: "
	    "images are of different size\n");
    exit(-1);
  }

  CFloatType *to = m_buffer;
  CFloatType *from = other.m_buffer;

  for(int i = GetNoPixels()-1; i >= 0; i--) {
    *to *= *from;
    to++; from++;
  }
}

void CFITSImage::operator =(const CFITSImage &other)
{
  if(other.m_buffer == 0) {
    if(m_buffer != 0) delete[] m_buffer;
    m_buffer = 0;
    m_xw = -1;
    m_yw = -1;
    return;
  }

  Resize(other.m_xw, other.m_yw, 0);

  CFloatType *to = m_buffer;
  CFloatType *from = other.m_buffer;

  for(int i = GetNoPixels()-1; i >= 0; i--) {
    *to = *from;
    to++; from++;
  }
}

void CFITSImage::Square()
{
  CFloatType *iter = m_buffer;
  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( ! IsNull(*iter) )
      (*iter) = (*iter)*(*iter);
    iter++;
  }
}

void CFITSImage::Sqrt()
{
  CFloatType *iter = m_buffer;
  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( ! IsNull(*iter) ) {
      if( *iter < 0.0 ) {
	fprintf(stderr, "*  CFITSImage::Sqrt() - negative pixel\n");
	exit(-1);
      }

      *iter = sqrt(*iter);
    }

    iter++;
  }
}

void CFITSImage::NullIfLess(CFloatType lower)
{
  CFloatType *iter = m_buffer;
  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( ! IsNull(*iter) && *iter < lower )
      *iter = CNullValue;
    iter++;
  }
}

void CFITSImage::NullIfGreater(CFloatType upper)
{
  CFloatType *iter = m_buffer;
  for(int i = GetNoPixels()-1; i >= 0; i--) {
    if( ! IsNull(*iter) && *iter > upper )
      *iter = CNullValue;
    iter++;
  }
}

void CFITSImage::Gaussian(CFloatType xsigma, CFloatType ysigma)
{
  printf(":  Applying Gaussian... ");
  fflush(stdout);

  CFloatType xsigmasqd = xsigma*xsigma;
  CFloatType ysigmasqd = ysigma*ysigma;

  // add up 3 sigmas worth of pixels around
  int xc = (int)ceil(xsigma*3.0);
  int yc = (int)ceil(ysigma*3.0);

  int bxw = 2*xc;
  int byw = 2*yc;

  // calculate a box for adding around
  CFloatType *gBox = new CFloatType[m_xw*m_yw];

  for(int x=0; x<bxw; x++)
    for(int y=0; y<byw; y++) {
      CFloatType delx = (CFloatType)(x-xc);
      CFloatType dely = (CFloatType)(y-yc);
      CFloatType pixel = exp( -(delx*delx)/(2.0*xsigmasqd)
			      -(dely*dely)/(2.0*ysigmasqd) );

      gBox[x + y*bxw] = pixel;
    }

  CFloatType *newBuffer = new CFloatType[ GetNoPixels() ];

  for(int yim=0; yim<m_yw; yim++) {
    for(int xim=0; xim<m_xw; xim++) {          // iterate over image
      CFloatType sum = 0.0;
      int changed = 0;

      CFloatType pixnorm = 0.;
      for(int yb=0; yb<byw; yb++)
	for(int xb=0; xb<bxw; xb++) {          // iterate over box
	  CFloatType pixel = GetPixelBounded(xim + xb - xc,
					     yim + yb - yc);
	  if( ! IsNull(pixel) ) {
	    sum += pixel * gBox[xb + yb*bxw];
	    pixnorm += gBox[xb + yb*bxw];
	    changed = 1;
	  }
	}

      if( changed )
	newBuffer[ xim + yim*m_xw ] = sum/pixnorm;
      else
	newBuffer[ xim + yim*m_xw ] = CNullValue;
    }
    
    if(yim % 48 == 0) {
      printf("%i ", yim);
      fflush(stdout);
    }
  }
  
  delete[] m_buffer;
  m_buffer = newBuffer;

  delete[] gBox;

  printf(" ... done\n");
}

