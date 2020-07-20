// Program for creating a binmap consisting of just annuli

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

// program to create a set of annuli coming out from a point
// in a binmap

#include <iostream>
#include <string>
#include <math.h>

#include <parammm/parammm.hh>
#include <FITSFile.h>
#include <FITSImage.h>
#include "version.hh"

using std::string;

class annuli_prog
{
public:
  annuli_prog(int argc, char **argv);

  void make_annuli();

  void process();

private:
  int m_xc, m_yc;        // x centre and y centre of annuli
  double m_radius;       // annuli radii

  CFITSImage m_outimage;
  CFITSPosn m_posn;

  string m_binmap_fname;
  string m_input_fname;
};


annuli_prog::annuli_prog(int argc, char **argv)
  : m_xc(-1), m_yc(-1),
    m_radius(20.),

    m_binmap_fname("annuli_binmap.fits")
{
  parammm::param params(argc, argv);

  params.add_switch( parammm::pswitch("binmap", 'n',
				      parammm::pstring_opt(&m_binmap_fname),
				      "set binmap out file "
				      "(def annuli_binmap.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("xc", 'x',
				      parammm::pint_opt(&m_xc),
				      "X centre",
				      "INT"));
  params.add_switch( parammm::pswitch("yc", 'y',
				      parammm::pint_opt(&m_yc),
				      "Y centre",
				      "INT"));
  params.add_switch( parammm::pswitch("radius", 'r',
				      parammm::pdouble_opt(&m_radius),
				      "Annuli widths",
				      "PIX"));

  params.set_autohelp("Usage: AnnuliMap [OPTIONS] file\n"
		      "Produces a annulus map on the same scale as input binmap\n"
		      "Written by Jeremy Sanders, 2001.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(c_adbin_version,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();

  params.interpret_and_catch();

  if(params.args().size() < 1)
    params.show_autohelp();

  m_input_fname = params.args()[0];
}


void annuli_prog::make_annuli()
{
  const int xw = m_outimage.GetXW(),
    yw = m_outimage.GetYW();

  for(int y=0; y<yw; ++y)
    for(int x=0; x<xw; ++x) {

      const int dx = x-m_xc, dy = y-m_yc;

      const double radius = sqrt(dx*dx + dy*dy);

      const int binno = int( radius / m_radius);

      m_outimage.SetPixel(x, y, binno);
    }
}

void annuli_prog::process()
{
  {
    CFITSFile f(m_input_fname.c_str(), CFITSFile::existingro);
    CFITSImage temp = f.GetImage();
    m_outimage.Resize(temp.GetXW(), temp.GetYW());
    m_posn = f.GetPosn();

    if(m_xc < 0 || m_yc < 0) {
      m_xc = m_outimage.GetXW() / 2;
      m_yc = m_outimage.GetYW() / 2;
    }
  }

  make_annuli();

  {
    CFITSFile f(m_binmap_fname.c_str(), CFITSFile::create);
    f.SetImage(m_outimage);
    f.SetPosn(m_posn);
    f.WriteImage();
  }
}


int main(int argc, char *argv[])
{
  annuli_prog prog(argc, argv);

  prog.process();

  return 0;
}
