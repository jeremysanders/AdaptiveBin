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
typedef std::vector<unsigned> unsigned_vec;

class smoother
{
public:
  smoother(const CFITSImage &binmap,
	   const CFITSImage &inimage);
  void smooth();
  const CFITSImage &out_image() { return m_out_image; }

private:
  unsigned get_no_bins() const;
  void find_centres();
  void find_bin_centre(int bin, double *xc, double *yc,
		       double *xsd, double *ysd) const;

private:
  const CFITSImage m_binmap_image;
  CFITSImage m_in_image;
  const unsigned m_xw, m_yw;

  CFITSImage m_out_image;
  const unsigned m_no_bins;

  dbl_vec m_xc, m_yc, m_vals;
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

void smoother::smooth_subset(unsigned x1, unsigned y1,
			     unsigned x2, unsigned y2)
{
  int bin_x1 = int(x1) - int(x2-x1);
  int bin_x2 = int(x2) + int(x2-x1);
  int bin_y1 = int(y1) - int(y2-y1);
  int bin_y2 = int(y2) + int(y2-y1);

  if( bin_x1 < 0 ) bin_x1 = 0;
  if( bin_y1 < 0 ) bin_y1 = 0;

  if( bin_x2 >= m_xw ) bin_x2 = m_xw - 1;
  if( bin_x2 >= m_yw ) bin_y2 = m_yw - 1;

  // get bin subset
  unsigned_vec binsinreg;

  for(unsigned i=0; i<m_no_bins; i++)
    {
      if( m_xc[i] >= unsigned(bin_x1) && m_xc[i] <= unsinged(bin_x2) &&
	  m_yc[i] >= unsigned(bin_y1) && m_yc[i] <= unsinged(bin_y2) )
	{
	  binsinreg.push_back(i);
	}
    }

  // copy subset
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
  // set to 1 for non-linear interpolation
  c_nnseti("igr",0);
  c_nnseti("non", 1);

  //  c_nnseti("sdi",1);
  //  c_nnseti("rad",1);

  double xc[m_no_bins], yc[m_no_bins], vals[m_no_bins];
  for(unsigned bin=0; bin<m_no_bins; ++bin)
    {
      xc[bin] = m_xc[bin];
      yc[bin] = m_yc[bin];
      vals[bin] = m_vals[bin];
    }


  // set return coord section
  double retx[m_xw], rety[m_yw];
  for(unsigned x=0; x<m_xw; ++x)
    retx[x] = x;
  for(unsigned y=0; y<m_yw; ++y)
    rety[y] = y;

  int ret = 0;
  const double *out = c_natgridd(m_no_bins, xc, yc, vals,
				 m_xw, m_yw, retx, rety, &ret);

  std::cout << "Done. Return value " << ret << '\n';

  for(unsigned y=0; y<m_yw; ++y)
    for(unsigned x=0; x<m_xw; ++x)
      {
	double val = out[x*m_yw + y];
	
	if( m_binmap_image.GetPixel(x, y) < 0. ) // ignore masked areas
	  val = CNullValue;

	//	val = m_in_image.GetPixel(x, y);
	
	m_out_image.SetPixel(x, y, val);
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
