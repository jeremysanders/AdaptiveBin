//      AdaptiveAnnuli Program

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

// program to take an image, and adjust the size of annuli
// until they contain a certain number of counts

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <parammm/parammm.hh>
#include <FITSFile.h>

#include "version.hh"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;
using std::ifstream;

const double c_masked = -1.;
const double c_notdone = -2.;

class annuli_prog
{
public:
  annuli_prog(int argc, char **argv);
  ~annuli_prog();

  void run();

private:
  void load_images();
  void write_output_image();

  void mask_output_binmap(const CFITSImage &mask_image);
  void create_annuli();
  void create_fixed_annuli();

  bool are_remaining_pixels();

  void set_pixel_with_bin(int x, int y, int bin);

  void sum_pixels_less_radius(double rad, double *counts,
			      int *number);
  void assign_pixels_less_radius(double rad, int bin);

private:
  int m_xw, m_yw;
  double m_min_counts;
  int m_xc, m_yc;
  double m_bg_counts;
  string m_binmap_filename, m_input_filename;
  string m_mask_filename;
  string m_fixedann_filename;

  int m_no_sectors;
  double m_start_angle;

  CFITSImage m_binmap_out;
  CFITSImage m_input_image;
  CFITSPosn m_posn;
};

annuli_prog::annuli_prog(int argc, char **argv)
  : m_min_counts(5000.), m_xc(0), m_yc(0),
    m_bg_counts(0.),
    m_binmap_filename("adannuli_binmap.fits"),
    m_no_sectors(1),
    m_start_angle(0.)
{
  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("binmap", 'n',
				      parammm::pstring_opt(&m_binmap_filename),
				      "set binmap out file (def adannuli"
				      "_binmap.fits)", "FILE"));
  params.add_switch( parammm::pswitch("mask", 'm',
				      parammm::pstring_opt(&m_mask_filename),
				      "set input mask filename (optional)",
				      "FILE"));
  params.add_switch( parammm::pswitch("background", 'b',
				      parammm::pdouble_opt(&m_bg_counts),
				      "set bg counts/pixel (def 0)",
				      "VAL"));
  params.add_switch( parammm::pswitch("xc", 'x',
				      parammm::pint_opt(&m_xc),
				      "set x-centre pixel ",
				      "INT"));
  params.add_switch( parammm::pswitch("yc", 'y',
				      parammm::pint_opt(&m_yc),
				      "set y-centre pixel ",
				      "INT"));
  params.add_switch( parammm::pswitch("counts", 'c',
				      parammm::pdouble_opt(&m_min_counts),
				      "set minimum no counts (def 5000)",
				      "VAL"));
  params.add_switch( parammm::pswitch("sectors", 's',
				      parammm::pint_opt(&m_no_sectors),
				      "set number of sectors (def 1)",
				      "VAL"));
  params.add_switch( parammm::pswitch("angle", 'a',
				      parammm::pdouble_opt(&m_start_angle),
				      "set sector starting angle (def 0)",
				      "DEG"));
  params.add_switch( parammm::pswitch("fixed", 'f',
				      parammm::pstring_opt
				      (&m_fixedann_filename),
				      "set filename for list of fixed annuli",
				      "FILE"));

  params.set_autohelp("Usage: AdaptiveAnnuli [OPTIONS] file\n"
		      "Adaptively creates annuli on an image\n"
		      "Written by Jeremy Sanders, 2001.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(c_adbin_version,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();

  params.interpret_and_catch();

  if(params.args().size() != 1)
    params.show_autohelp();

  m_input_filename = params.args()[0];

  if(!m_fixedann_filename.empty()) {
    cout << "Reading annuli from " << m_fixedann_filename << endl;
  } else {
    cout << "Creating annuli with a minimum of " << m_min_counts
	 << " counts\n";
  }
}

annuli_prog::~annuli_prog()
{
}

void annuli_prog::run()
{
  load_images();
  if( !m_fixedann_filename.empty())
    create_fixed_annuli();
  else
    create_annuli();
  write_output_image();
}

void annuli_prog::load_images()
{
  // load input file
  {
    CFITSFile f(m_input_filename.c_str(), CFITSFile::existingro);
    m_input_image = f.GetImage();
    m_posn = f.GetPosn();
    m_xw = m_input_image.GetXW();
    m_yw = m_input_image.GetYW();
  }

  // set output image
  m_binmap_out.Resize(m_xw, m_yw);
  m_binmap_out.SetAll(c_notdone);

  // load mask
  CFITSImage mask_image;
  if( ! m_mask_filename.empty() ) {
    CFITSFile f(m_mask_filename.c_str(), CFITSFile::existingro);
    mask_image = f.GetImage();
  } else {
    mask_image.Resize(m_xw, m_yw);
    mask_image.SetAll(0.);
  }

  // mask masked pixels
  mask_output_binmap( mask_image );
}

void annuli_prog::write_output_image()
{
  CFITSFile f(m_binmap_filename.c_str(), CFITSFile::create);
  f.SetImage(m_binmap_out);
  f.SetPosn(m_posn);
  f.WriteImage();

  f.UpdateKey("RAD_XC",  CFITSFile::tint, &m_xc,
	      "X centre pixel for radial bins");
  f.UpdateKey("RAD_YC",  CFITSFile::tint, &m_yc,
	      "Y centre pixel for radial bins");
  f.UpdateKey("RAD_SEC", CFITSFile::tint, &m_no_sectors,
	      "Number of radial sectors");
  f.UpdateKey("RAD_MCTS", CFITSFile::tfloat, &m_min_counts,
	      "Minimum number of counts in radial bins");
  f.UpdateKey("RAD_BCTS", CFITSFile::tfloat, &m_bg_counts,
	      "Background number of counts/pixel for radial bins");

  // add input file as history
  ostringstream o;
  o << "AdaptiveAnnuli binmap generated from " << m_input_filename
    << '\0';
  f.WriteHistory(o.str().c_str());
}

void annuli_prog::mask_output_binmap(const CFITSImage &mask_image)
{
  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {
      if( mask_image.GetPixel(x, y) > 0. )
	m_binmap_out.SetPixel(x, y, c_masked);
    }
}

bool annuli_prog::are_remaining_pixels()
{
  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {
      if( fabs(m_binmap_out.GetPixel(x, y) - c_notdone) < 1e-8 )
	return true;
    }

  return false;
}

void annuli_prog::sum_pixels_less_radius(double rad, double *counts,
					 int *number)
{
  *counts = 0.;
  *number = 0;

  const double radius_sqd = rad*rad;

  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {
      const int dx = x-m_xc;
      const int dy = y-m_yc;
      const double curr_radius_sqd = dx*dx + dy*dy;
      if( curr_radius_sqd < radius_sqd &&
	  fabs(m_binmap_out.GetPixel(x, y) - c_notdone) < 1e-8 ) {
	(*counts) += (m_input_image.GetPixel(x, y) - m_bg_counts);
	++(*number);
      }
    }
}

void annuli_prog::set_pixel_with_bin(int x, int y, int bin)
{
  if( fabs(m_binmap_out.GetPixel(x, y) - c_notdone) < 1e-8  ) {
    const int dx = x-m_xc;
    const int dy = y-m_yc;
    const double rad_start = m_start_angle / 180. * M_PI;
    const double rad_per_sector = (M_PI * 2.) / m_no_sectors;

    // find sector
    const double angle = atan2(double(dy), double(dx)) -
      rad_start + M_PI;

    // move sector into range 0 -> m_no_sectors-1
    int sector = int(angle / rad_per_sector);
    while( sector >= m_no_sectors )
      sector -= m_no_sectors;
    while(sector < 0)
      sector += m_no_sectors;

    m_binmap_out.SetPixel(x, y, bin*m_no_sectors + sector);
  }
}

void annuli_prog::assign_pixels_less_radius(double rad, int bin)
{
  const double radius_sqd = rad*rad;

  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {
      const int dx = x-m_xc;
      const int dy = y-m_yc;
      const double curr_radius_sqd = dx*dx + dy*dy;
      if( curr_radius_sqd < radius_sqd )
	set_pixel_with_bin(x, y, bin);
    }
}

void annuli_prog::create_annuli()
{
  double last_radius = 0.;
  int bin_no = 0;

  // maximum radius
  const double maxdist = sqrt(m_xw*m_xw + m_yw*m_yw);

  while( are_remaining_pixels() ) {
    last_radius += 0.5;

    double count;
    int pixels;

    sum_pixels_less_radius(last_radius, &count, &pixels);
    if( count >= m_min_counts || last_radius > maxdist) {
      cout << "Created bin " << bin_no << " at radius "
	   << last_radius << endl;
      assign_pixels_less_radius(last_radius, bin_no);
      ++bin_no;
    }
  }

  // finished
}

void annuli_prog::create_fixed_annuli()
{
  ifstream in_list(m_fixedann_filename.c_str());
  if(!in_list) {
    cerr << "Cannot open " << m_fixedann_filename << endl;
    exit(1);
  }

  int bin = 0;
  while(in_list) {
    double rad;
    in_list >> rad;
    if(!in_list) continue;

    cout << "Annulus " << bin << " inside " << rad << " pixels\n";
    const double rad_sqd = rad*rad;

    for(int y=0; y<m_yw; ++y)
      for(int x=0; x<m_xw; ++x) {
	const int dx = x-m_xc;
	const int dy = y-m_yc;
	const double curr_radius_sqd = dx*dx + dy*dy;
	if(curr_radius_sqd < rad_sqd)
	  set_pixel_with_bin(x, y, bin);
    }

    ++bin;
  }

  cout << "Remaining pixels are annulus " << bin << endl;
  for(int y=0; y<m_yw; ++y)
    for(int x=0; x<m_xw; ++x) {
      // function only sets "unset" pixels
      set_pixel_with_bin(x, y, bin);
    }
}

////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  annuli_prog program(argc, argv);
  program.run();
  return 0;
}
