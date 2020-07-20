// Program takes an existing pixel file and bins image using it
// run with --help to get command line options

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

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <FITSFile.h>
#include <FITSImage.h>
#include <FITSPosn.h>
#include <parammm/parammm.hh>
#include "binmodule.hh"
#include "version.hh"

using std::string;
using std::vector;
using std::ostringstream;
using std::endl;
using std::cerr;
using std::cout;

class abpixelcopy_program
{
public:
  abpixelcopy_program(int argc, char *argv[]);
  ~abpixelcopy_program();

  void run();

private:
  void bin_image();
  void add_history_list(CFITSFile *file);

private:
  string m_in_fname, m_binmap_fname, m_out_fname, m_err_fname;

  CFITSImage m_bin_image;
  CFITSImage m_in_image;
  CFITSImage m_out_image, m_err_image;
  CFITSPosn m_in_posn;
  double m_backgrnd;
  bool m_verbose;

  vector<string> m_history_list;
};

abpixelcopy_program::abpixelcopy_program(int argc, char *argv[])
  : m_binmap_fname("adbin_binmap.fits"),
    m_backgrnd(0.), m_verbose(false)
{
  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("image", 'i',
				      parammm::pstring_opt(&m_in_fname),
				      "set input image file (req)",
				      "FILE"));
  params.add_switch( parammm::pswitch("out", 'o',
				      parammm::pstring_opt(&m_out_fname),
				      "set out file (req)",
				      "FILE"));
  params.add_switch( parammm::pswitch("err", 'e',
				      parammm::pstring_opt(&m_err_fname),
				      "set out error map (req)",
				      "FILE"));
  params.add_switch( parammm::pswitch("binmap", 'n',
				      parammm::pstring_opt(&m_binmap_fname),
				      "set binmap file (def adbin_binmap.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("pixel", 'p',
				      parammm::pstring_opt(&m_binmap_fname),
				      "set input pixel file (DISCOURAGED)",
				      "FILE"));
  params.add_switch( parammm::pswitch("background", 'b',
				      parammm::pdouble_opt(&m_backgrnd),
				      "set background counts/pixel",
				      "VAL"));
  params.add_switch( parammm::pswitch("verbose", 0,
				      parammm::pbool_noopt(&m_verbose),
				      "display more information",
				      ""));

  params.set_autohelp("Usage: ABPixelCopy [OPTIONS] --image=in.fits "
		      "--out=out.fits --err=err.fits\n"
		      "Bins image using pixel file, creating output and error "
		      "file.\nDoes background subtraction.\n"
		      "Written by Jeremy Sanders, 2000, 2001.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(c_adbin_version,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();
  params.interpret_and_catch();

  if( m_in_fname.empty() || m_out_fname.empty() || m_err_fname.empty() )
    params.show_autohelp();
}

abpixelcopy_program::~abpixelcopy_program()
{
}

void abpixelcopy_program::add_history_list(CFITSFile *file)
{
  const int no = m_history_list.size();

  for(int i=0; i<no; ++i) {
    const string line = "adbin: " + m_history_list[i];
    file -> WriteHistory(line.c_str());
  }
}

void abpixelcopy_program::run()
{
  // load images
  {
    CFITSFile infile(m_in_fname.c_str(), CFITSFile::existingro);
    m_in_image = infile.GetImage();
    m_in_posn = infile.GetPosn();
  }{
    CFITSFile infile(m_binmap_fname.c_str(), CFITSFile::existingro);
    m_bin_image = infile.GetImage();
  }

  bin_image();

  m_history_list.push_back("file created by ABPixelCopy v. " +
			   string(c_adbin_version));
  m_history_list.push_back("input image: " + m_in_fname);
  m_history_list.push_back("bin map (input): " + m_binmap_fname);
  m_history_list.push_back("output image: " + m_out_fname);
  m_history_list.push_back("error map: " + m_err_fname);
  {
    ostringstream o;
    o << "background: " << m_backgrnd << '\0';
    m_history_list.push_back(o.str());
  }

  if(m_verbose) {
    cout << "\nHeaders written to output files:\n";
    for(int i=0; i<int(m_history_list.size()); ++i)
      cout << m_history_list[i] << endl;
    cout << endl;
  }

  // save image
  {
    CFITSFile outfile(m_out_fname.c_str(), CFITSFile::create);
    outfile.SetImage(m_out_image);
    //outfile.SetPosn(m_in_posn);
    outfile.WriteImage();
    outfile.WriteHistory("adbin: file is output image");
    add_history_list(&outfile);
  }{
    CFITSFile outfile(m_err_fname.c_str(), CFITSFile::create);
    outfile.SetImage(m_err_image);
    //outfile.SetPosn(m_in_posn);
    outfile.WriteImage();
    outfile.WriteHistory("adbin: file is error map");
    add_history_list(&outfile);
  }
}

void abpixelcopy_program::bin_image()
{
  // rewritten to avoid going through the image too many times

  const int xw = m_in_image.GetXW();
  const int yw = m_in_image.GetYW();

  m_out_image.Resize(xw, yw);
  m_out_image.SetAll( 0. / 0. );
  m_err_image.Resize(xw, yw);

  int max_bin = -1;
  for(int y=0; y<yw; y++)
    for(int x=0; x<xw; x++) {
      const int bin = (int)m_bin_image.GetPixel(x, y);
      if( bin > max_bin ) max_bin = bin;
    }

  // collect lists of pixels in bins
  vector<AdaptiveBin::pixlist> pixels_in_bins(max_bin+1);
  for(int y=0; y<yw; y++)
    for(int x=0; x<xw; x++) {
      const int pix = (int)m_bin_image.GetPixel(x, y);
      if( pix >= 0 and not isnan(m_in_image.GetPixel(x, y)) )
	pixels_in_bins[pix].push_back( AdaptiveBin::pixel(x, y) );
    }

  // do the binning
  for(int bin_no=max_bin; bin_no>=0; --bin_no) {
    double tot = 0.;
    AdaptiveBin::pixlist &p = pixels_in_bins[bin_no];

    if( p.size() == 0 ) {
      cerr << "Warning - missing bin "
	   << bin_no
	   << " in bin map\n";
      continue;
    }

    for(int i=p.size()-1; i>=0; --i) {
      const double pix = m_in_image.GetPixel( p[i].x(), p[i].y() );
      tot += pix;
    }

    const double av = tot/p.size() - m_backgrnd;
    const double ferr = sqrt(tot + m_backgrnd*p.size()) /
      (tot - m_backgrnd*p.size());

    for(int i=p.size()-1; i>=0; --i) {
      m_out_image.SetPixel(p[i].x(), p[i].y(), av);
      m_err_image.SetPixel(p[i].x(), p[i].y(), ferr);
    }
  } // loop over bins

}

int main(int argc, char *argv[])
{
  abpixelcopy_program program(argc, argv);
  program.run();
  return 0;
}
