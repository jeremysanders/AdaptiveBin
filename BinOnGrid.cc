// simple program to bin by a factor, but preserves image size

#include <string>
#include <vector>
#include <map>

#include <parammm/parammm.hh>
#include "../FITSFile.h"

using std::string;
using std::vector;
using std::map;

class binongrid {
public:
  binongrid(int argc, char **argv);
  ~binongrid();
  void Run();

private:
  string m_outfname, m_infname, m_maskfname;
  int m_binsize;
};

binongrid::binongrid(int argc, char **argv)
  : m_outfname("pixel.fits"),
    m_binsize(1)
{
  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("out", 'o',
				      parammm::pstring_opt(&m_outfname),
				      "set out pixel file",
				      "FILE"));
  params.add_switch( parammm::pswitch("in", 'i',
				      parammm::pstring_opt(&m_infname),
				      "input file to get size from",
				      "FILE"));
  params.add_switch( parammm::pswitch("size", 's',
				      parammm::pint_opt(&m_binsize),
				      "set bin size "
				      "(def. 1)",
				      "INT"));
  params.add_switch( parammm::pswitch("mask", 'm',
				      parammm::pstring_opt(&m_maskfname),
				      "set mask file",
				      "FILE"));


  params.set_autohelp("Usage: BinOnGrid [options] --in=in.fits\n"
		      "Creates a pixel file with certain size bins\n"
		      "on an image with the same size as the input image\n"
		      "Written by Jeremy Sanders (2000-2005).",
		      "Report bugs to <jss@ast.cam.ac.uk>");

  params.enable_autohelp();
  params.interpret_and_catch();

  if(params.args().size() != 0)
    params.show_autohelp();

  if( m_infname.size() == 0 )
    params.show_autohelp();
}

binongrid::~binongrid()
{
}

void binongrid::Run()
{
  int xw, yw;
  {
    CFITSFile infile(m_infname.c_str(), CFITSFile::existingro);
    xw = infile.GetImage().GetXW();
    yw = infile.GetImage().GetYW();
  }

  CFITSImage maskimage(xw, yw, 1.);
  if( ! m_maskfname.empty() ) {
    CFITSFile maskfile(m_maskfname.c_str(), CFITSFile::existingro);
    maskimage = maskfile.GetImage();
  }

  int w = xw/m_binsize;
  if( xw % m_binsize != 0 ) w++;

  CFITSImage binpiximage(xw, yw);
  for(int y=0; y<yw; y++)
    for(int x=0; x<xw; x++) {
      const int pix = (x/m_binsize)+(y/m_binsize)*w;
      binpiximage.SetPixel(x, y, pix);
    }

  {
    CFITSFile outfile(m_outfname.c_str(), CFITSFile::create);
    outfile.SetImage(binpiximage);
    outfile.WriteImage();
  }
}

int main(int argc, char *argv[])
{
  binongrid prog(argc, argv);
  prog.Run();

  return 0;
}
