//      Smoothing program for adaptively binned files

//      Described in Sanders and Fabian (submitted)
//      Routines for adaptively binning data
//      Copyright (C) 2000, 2001 Jeremy Sanders
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

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <parammm/parammm.hh>
#include <FITSFile.h>
#include <cassert>
#include <cmath>
#include "version.hh"

using std::cout;
using std::string;

extern "C" {
  double *c_natgridd(int npts, double *x, double *y, double *z,
                     int numx_out, int numy_out,
		     double *xvals_out, double *yvals_out,
		     int *retval);
  void c_nnseti(char *param, int val);
}

typedef std::vector<double> dbl_vec;

class smoother
{
public:
  smoother(const CFITSImage &binmap,
	   const CFITSImage &inimage);
  void smooth();
  const CFITSImage &out_image() { return m_out_image; }

private:
  unsigned get_no_bins() const;
  void calc_largest_near_dist();
  void find_centres();
  void find_bin_centre(int bin, double *xc, double *yc,
		       double *xsd, double *ysd) const;
  void smooth_part(unsigned xp, unsigned yp, unsigned wx, unsigned wy);

private:
  const CFITSImage m_binmap_image;
  CFITSImage m_in_image;
  const unsigned m_xw, m_yw;

  CFITSImage m_out_image;
  const unsigned m_no_bins;

  dbl_vec m_xc, m_yc, m_vals;
  double m_furthest_near;
};

smoother::smoother(const CFITSImage &binmap,
		   const CFITSImage &image)
  : m_binmap_image(binmap),
    m_in_image(image),
    m_xw(binmap.GetXW()), m_yw(binmap.GetYW()),
    m_out_image(m_xw, m_yw),
    m_no_bins(get_no_bins())
{
}

unsigned smoother::get_no_bins() const
{
  int maxbin = -1;
  for(unsigned y=0; y<m_yw; ++y)
    for(unsigned x=0; x<m_xw; ++x) {
      const int thebin = int(m_binmap_image.GetPixel(x, y));
      if(thebin > maxbin) maxbin = thebin;
    }
  assert(maxbin >= 0);

  std::cout << "Found " << maxbin+1 << " bins\n";
  return unsigned(maxbin+1);
}

class bindata
{
public:
  double val;
  double totx, toty;
  unsigned count;
  int miny, maxy;
  int minx, maxx;

  bindata() 
  {
    val = -1;
    totx = toty = 0.;
    count = 0;
    minx = miny = 100000;
    maxx = maxy = -1;
  }
};

void smoother::find_centres()
{
  std::cout << "Finding centres...\n";

  std::vector<bindata> totals(m_no_bins);

  // loop, adding data for pixels
  for(unsigned y=0; y<m_yw; ++y)
    for(unsigned x=0; x<m_xw; ++x)
      {
	const int bin = int(m_binmap_image.GetPixel(x, y));
	if( bin >= 0 )
	  {
	    assert(bin < int(m_no_bins));

	    bindata& bd = totals[bin];
	    bd.val = m_in_image.GetPixel(x, y);
	    bd.totx += x;
	    bd.toty += y;
	    bd.count ++;
	    if( int(y) < bd.miny ) bd.miny = y;
	    if( int(y) > bd.maxy ) bd.maxy = y;
	    if( int(x) < bd.minx ) bd.minx = x;
	    if( int(x) > bd.maxx ) bd.maxx = x;
	  }
      }

  for(unsigned i=0; i<m_no_bins; i++)
    {
      const bindata& bd = totals[i];
      m_xc.push_back( bd.totx/bd.count );
      m_yc.push_back( bd.toty/bd.count );
      m_vals.push_back( bd.val );
    }

  std::cout << "Done\n";
}

void smoother::calc_largest_near_dist()
{
  // find maximum distance between nearest neighbours

  std::cout << "Finding furthest nearest neighbour\n";

  double fur_near_dist_2 = -1e99;

  for(unsigned bin=0; bin<m_no_bins; ++bin)
    {
      double min_dist_2 = 1e99;

      // find nearest neighbour
      for(unsigned b2=0; b2<m_no_bins; ++b2)
	{
	  if( bin == b2 ) continue;

	  const double dx = m_xc[b2]-m_xc[bin];
	  const double dy = m_yc[b2]-m_yc[bin];
	  const double d_2 = dx*dx + dy*dy;
	  if( d_2 < min_dist_2 )
	    min_dist_2 = d_2;
	}

      // is it further away that any of the others?
      if( min_dist_2 > fur_near_dist_2 )
	fur_near_dist_2 = min_dist_2;
    }

  std::cout << "Done\n";

  m_furthest_near = std::sqrt(fur_near_dist_2);
}

void smoother::smooth_part(unsigned xp, unsigned yp, unsigned wx, unsigned wy)
{
  double xc[m_no_bins], yc[m_no_bins], vals[m_no_bins];

  // set return coord section
  double retx[wx], rety[wy];
  for(unsigned x=0; x<wx; ++x)
    retx[x] = x+xp;
  for(unsigned y=0; y<wy; ++y)
    rety[y] = y+yp;

  const double del_dist = m_furthest_near*3;
  unsigned count_bins = 0;
  // find appropriate bins for this section
  for(unsigned bin=0; bin<m_no_bins; ++bin)
    {
      const double x = m_xc[bin];
      const double y = m_yc[bin];

      if( x >= (xp-del_dist) && x <= (xp+wx+del_dist) &&
	  y >= (yp-del_dist) && y <= (yp+wy+del_dist) )
	{
	  xc[count_bins] = x;
	  yc[count_bins] = y;
	  vals[count_bins] = m_vals[bin];
	  count_bins++;
	}
    }

  if( count_bins < 8 )
    return;

  // smooth section
  std::cout << "Smoothing " << xp << "," << yp << " to "
	    << xp+wx-1 << "," << yp+wy-1 << '\n';

  int ret = 0;
  const double *out = c_natgridd(count_bins, xc, yc, vals,
				 wx, wy, retx, rety, &ret);

  std::cout << "Done. Return value " << ret << '\n';


  for(unsigned y=0; y<wy; ++y)
    for(unsigned x=0; x<wx; ++x) {
      double val = out[x*wx + y];

      const unsigned rx = x+xp;
      const unsigned ry = y+yp;

      if( m_binmap_image.GetPixel(rx, ry) < 0. ) // ignore masked areas
	val = CNullValue;

      m_out_image.SetPixel(rx, ry, val);
    }

}


void smoother::smooth()
{
  m_xc.clear();
  m_yc.clear();
  m_vals.clear();

  // find centres of each bin
  find_centres();
  // calculate the largest nearest neighbour
  //  calc_largest_near_dist();
  m_furthest_near = 16;
  // set to 1 for non-linear interpolation
  c_nnseti("igr",0);

  //  c_nnseti("sdi",1);
  //  c_nnseti("rad",1);

  const unsigned segment = 32;
  for(unsigned y=0; y<m_yw; y += segment)
    for(unsigned x=0; x<m_xw; x += segment)
      {
	unsigned wx = segment;
	unsigned wy = segment;

	if( wx+x >= m_xw ) wx = m_xw - x - 1;
	if( wy+y >= m_yw ) wy = m_yw - y - 1;

	smooth_part(x, y, wx, wy);
      }

}

int main(int argc, char *argv[])
{
  string file_out = "adbin_smoothed.fits";

  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("out", 'o',
				      parammm::pstring_opt(&file_out),
				      "set out file (def adbin_smoothed.fits)",
				      "FILE"));


  params.set_autohelp("Usage: ABPostSmooth [OPTIONS] binmap.fits image.fits\n"
		      "Smooths an adaptively-binned image using interpolation\n"
		      "Written by Jeremy Sanders, 2001.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(c_adbin_version,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();

  params.interpret_and_catch();

  if(params.args().size() != 2)
    params.show_autohelp();

  CFITSImage binmap, image;
  CFITSPosn posn;
  {
    CFITSFile f(params.args()[0].c_str(), CFITSFile::existingro);
    binmap = f.GetImage();
    posn = f.GetPosn();
  }{
    CFITSFile f(params.args()[1].c_str(), CFITSFile::existingro);
    image = f.GetImage();
  }

  smoother s(binmap, image);
  s.smooth();

  {
    CFITSFile f(file_out.c_str(), CFITSFile::create);
    f.SetImage( s.out_image() );
    f.SetPosn( posn );
    f.WriteImageInclNull();
  }
}
