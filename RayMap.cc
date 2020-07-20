// program to create a set of rays coming out from a point
// in a binmap

#include <iostream>
#include <string>
#include <math.h>

#include <parammm/parammm.hh>
#include <FITSFile.h>
#include <FITSImage.h>
#include "version.hh"

using std::string;

class ray_prog
{
public:
  ray_prog(int argc, char **argv);

  void make_rays();

  void process();

private:
  int m_xc, m_yc;        // x centre and y centre of rays
  int m_number_rays;     // number of rays
  double m_start_angle;  // angle for rays to start

  CFITSImage m_outimage;
  CFITSPosn m_posn;

  string m_binmap_fname;
  string m_input_fname;
};


ray_prog::ray_prog(int argc, char **argv)
  : m_xc(-1), m_yc(-1),
    m_number_rays(4),
    m_start_angle(0.),

    m_binmap_fname("ray_binmap.fits")
{
  parammm::param params(argc, argv);

  params.add_switch( parammm::pswitch("binmap", 'n',
				      parammm::pstring_opt(&m_binmap_fname),
				      "set binmap out file "
				      "(def ray_binmap.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("xc", 'x',
				      parammm::pint_opt(&m_xc),
				      "X centre",
				      "INT"));
  params.add_switch( parammm::pswitch("yc", 'y',
				      parammm::pint_opt(&m_yc),
				      "Y centre",
				      "INT"));
  params.add_switch( parammm::pswitch("norays", 'o',
				      parammm::pint_opt(&m_number_rays),
				      "Number of rays",
				      "INT"));
  params.add_switch( parammm::pswitch("angle", 'a',
				      parammm::pdouble_opt(&m_start_angle),
				      "Start angle",
				      "DEG"));

  params.set_autohelp("Usage: RayMap [OPTIONS] file\n"
		      "Produces a ray map on the same scale as input binmap\n"
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


void ray_prog::make_rays()
{
  const int xw = m_outimage.GetXW(),
    yw = m_outimage.GetYW();

  const double rayangle = 360./m_number_rays;

  for(int y=0; y<yw; ++y)
    for(int x=0; x<xw; ++x) {

      const int dx = x-m_xc, dy = y-m_yc;

      // get angle between -180 -> 180
      const double angle = atan2(dy, dx) * 180.0 / M_PI;

      double da = angle-m_start_angle;
      if(da < 0) da += 360.;

      const int binno = int( da / rayangle );

      m_outimage.SetPixel(x, y, binno);
    }
}

void ray_prog::process()
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

  make_rays();

  {
    CFITSFile f(m_binmap_fname.c_str(), CFITSFile::create);
    f.SetImage(m_outimage);
    f.SetPosn(m_posn);
    f.WriteImage();
  }
}


int main(int argc, char *argv[])
{
  ray_prog prog(argc, argv);

  prog.process();

  return 0;
}
