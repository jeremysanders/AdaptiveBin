// SigCalc.cc

#include <list>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "../FITSImage.h"
#include "Coord.hh"
#include "SigCalc.hh"

/// (Almost) Abstract class ///////////////

using namespace AdaptiveBin;
using std::sqrt;
using std::fabs;

CSigCalc::CSigCalc()
{
}

CSigCalc::~CSigCalc()
{
}

/// Count significance ///////////

CCountSig::CCountSig(const CFITSImage &image) :
  m_image(image),
  m_maskImage(GetXW(), GetYW())
{
}

int CCountSig::GetNoValues()
{
  return 1;
}

int CCountSig::GetXW()
{
  return m_image.GetXW();
}

int CCountSig::GetYW()
{
  return m_image.GetYW();
}

double CCountSig::GetValue(const cCoordList &points,
			   int value)
{
  assert( value == 0 );
  assert( ! points.empty() );

  cCoordList::const_iterator p = points.begin();

  double tot = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tot += m_image.GetPixel( p->m_x, p->m_y );
    ++count;
    ++p;
  }

  return tot/count;
}

double CCountSig::GetError(const cCoordList &points,
			   int value)
{
  assert( ! points.empty() );
  assert( value == 0 );

  cCoordList::const_iterator p = points.begin();

  double tot = 0.0;
  while( p != points.end() ) {
    tot += m_image.GetPixel( p->m_x, p->m_y );
    ++p;
  }
  if( fabs(tot) < 1e-5 ) return 1e5;

  return 1.0/sqrt(tot);
}

void CCountSig::Mask( CCoord p )
{
  m_maskImage.SetPixel(p.m_x, p.m_y, 1.0);
}

bool CCountSig::IsMasked( CCoord p )
{
  return( fabs(m_maskImage.GetPixel(p.m_x, p.m_y)) > 1e-5 );
}

// Ratio significance //////////////////////////////

CRatioSig::CRatioSig(const CFITSImage &imagea,
		     const CFITSImage &imageb  ) :
  m_imagea(imagea),
  m_imageb(imageb),
  m_maskImage(GetXW(), GetYW())
{
}

int CRatioSig::GetNoValues()
{
  return 3;  // ratio, count1, count2
}

int CRatioSig::GetXW()
{
  return m_imagea.GetXW();
}

int CRatioSig::GetYW()
{
  return m_imagea.GetYW();
}

double CRatioSig::GetValue(const cCoordList &points,
			   int value)
{
  assert( ! points.empty() );
  assert( value < 3 );

  cCoordList::const_iterator p = points.begin();

  double tota = 0.0;
  double totb = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tota += m_imagea.GetPixel( p->m_x, p->m_y );
    totb += m_imageb.GetPixel( p->m_x, p->m_y );
    ++p; ++count;
  }

  if( fabs(totb) < 1e-5) return 1e5;

  switch(value) {
  case 1: return tota/count;
  case 2: return totb/count;
  default: return tota/totb;
  }
}

double CRatioSig::GetError(const cCoordList &points,
			   int value)
{
  assert( ! points.empty() );

  cCoordList::const_iterator p = points.begin();

  double tota = 0.0;
  double totb = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tota += m_imagea.GetPixel( p->m_x, p->m_y );
    totb += m_imageb.GetPixel( p->m_x, p->m_y );
    ++count;
    ++p;
  }

  if( fabs(totb) < 1e-5 || fabs(tota) < 1e-5 )
    return 1e5;

  switch(value) {
  case 1: return 1.0/sqrt(tota);
  case 2: return 1.0/sqrt(totb);
  default: return sqrt( 1.0/tota + 1.0/totb );
  }
}

void CRatioSig::Mask( CCoord p )
{
  m_maskImage.SetPixel(p.m_x, p.m_y, 1.0);
}

bool CRatioSig::IsMasked( CCoord p )
{
  return( fabs(m_maskImage.GetPixel(p.m_x, p.m_y)) > 1e-5 );
}

//// Ratio file which also takes account of third colour

CRatioSig3::CRatioSig3(const CFITSImage &imagea,
		       const CFITSImage &imageb,
		       const CFITSImage &imagec)
  : CRatioSig(imagea, imageb),
    m_imagec(imagec)
{
}

int CRatioSig3::GetNoValues()
{
  return 7;
}

double CRatioSig3::GetValue(const cCoordList &points,
			    int value)
{
  assert( value < 7 );
  assert( ! points.empty() );

  cCoordList::const_iterator p = points.begin();

  double tota = 0.0;
  double totb = 0.0;
  double totc = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tota += m_imagea.GetPixel( p->m_x, p->m_y );
    totb += m_imageb.GetPixel( p->m_x, p->m_y );
    totc += m_imagec.GetPixel( p->m_x, p->m_y );
    ++count;
    ++p;
  }

  switch(value) {
  case 1: return tota/totb;
  case 2: return totb/totc;
  case 3: return totc/tota;
  case 4: return tota/count;
  case 5: return totb/count;
  case 6: return totc/count;
  default: return tota/totb;
  }
}

double CRatioSig3::GetError(const cCoordList &points,
			    int value)
{
  assert( value < 7 );
  assert( ! points.empty() );

  cCoordList::const_iterator p = points.begin();

  double tota = 0.0;
  double totb = 0.0;
  double totc = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tota += m_imagea.GetPixel( p->m_x, p->m_y );
    totb += m_imageb.GetPixel( p->m_x, p->m_y );
    totc += m_imagec.GetPixel( p->m_x, p->m_y );
    ++count;
    ++p;
  }

  if( fabs(totb) < 1e-5 || fabs(tota) < 1e-5 ||
      fabs(totc) < 1e-5 )
    return 1e5;

  double errab = sqrt(1.0/tota + 1.0/totb);
  double errbc = sqrt(1.0/totb + 1.0/totc);
  double errca = sqrt(1.0/totc + 1.0/tota);

  switch(value) {
  case 1: return errab;
  case 2: return errbc;
  case 3: return errca;
  case 4: return 1.0/sqrt(tota);
  case 5: return 1.0/sqrt(totb);
  case 6: return 1.0/sqrt(totc);
    //  default: return (errab+errbc+errca)/3.0;
  default: return sqrt(1.0/tota+1.0/totb+1.0/totc);
  }
}

//// Ratio file which also takes account of fourth(!) colour

CRatioSig4::CRatioSig4(const CFITSImage &imagea,
		       const CFITSImage &imageb,
		       const CFITSImage &imagec,
		       const CFITSImage &imaged)
  : CRatioSig(imagea, imageb),
    m_imagec(imagec), m_imaged(imaged)
{
}

int CRatioSig4::GetNoValues()
{
  return 9;
}

double CRatioSig4::GetValue(const cCoordList &points,
			    int value)
{
  assert( value < 9 );
  assert( ! points.empty() );

  cCoordList::const_iterator p = points.begin();

  double tota = 0.0;
  double totb = 0.0;
  double totc = 0.0;
  double totd = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tota += m_imagea.GetPixel( p->m_x, p->m_y );
    totb += m_imageb.GetPixel( p->m_x, p->m_y );
    totc += m_imagec.GetPixel( p->m_x, p->m_y );
    totd += m_imaged.GetPixel( p->m_x, p->m_y );
    ++count;
    ++p;
  }

  switch(value) {
  case 1: return tota/totb;
  case 2: return totb/totc;
  case 3: return totc/totd;
  case 4: return totd/tota;
  case 5: return tota/count;
  case 6: return totb/count;
  case 7: return totc/count;
  case 8: return totd/count;
  default: return tota/totb;
  }
}

double CRatioSig4::GetError(const cCoordList &points,
			    int value)
{
  assert( value < 9 );
  assert( ! points.empty() );

  cCoordList::const_iterator p = points.begin();

  double tota = 0.0;
  double totb = 0.0;
  double totc = 0.0;
  double totd = 0.0;
  int count = 0;
  while( p != points.end() ) {
    tota += m_imagea.GetPixel( p->m_x, p->m_y );
    totb += m_imageb.GetPixel( p->m_x, p->m_y );
    totc += m_imagec.GetPixel( p->m_x, p->m_y );
    totd += m_imaged.GetPixel( p->m_x, p->m_y );
    ++count;
    ++p;
  }

  if( fabs(totb) < 1e-5 || fabs(tota) < 1e-5 ||
      fabs(totc) < 1e-5 || fabs(totd) < 1e-5 )
    return 1e5;

  double errab = sqrt(1.0/tota + 1.0/totb);
  double errbc = sqrt(1.0/totb + 1.0/totc);
  double errcd = sqrt(1.0/totc + 1.0/totd);
  double errda = sqrt(1.0/totd + 1.0/tota);

  switch(value) {
  case 1: return errab;
  case 2: return errbc;
  case 3: return errcd;
  case 4: return errda;
  case 5: return 1.0/sqrt(tota);
  case 6: return 1.0/sqrt(totb);
  case 7: return 1.0/sqrt(totc);
  case 8: return 1.0/sqrt(totd);
  default: return sqrt(1.0/tota+1.0/totb+1.0/totc+1.0/totd);
  }
}
